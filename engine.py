"""Pipeline profile calculations. All lengths are measured in the selected metric CRS."""
import io, json, math, zipfile, tempfile, re
from pathlib import Path
import numpy as np
from shapely.geometry import shape, Point, LineString
from shapely.ops import transform, unary_union
from pyproj import Transformer, CRS


def kp(m):
    return f'{int(m//1000)}+{m%1000:07.3f}'


def geo_features(data):
    if data.get('type') == 'FeatureCollection': return data['features']
    if data.get('type') == 'Feature': return [data]
    return [{'type':'Feature','geometry':data,'properties':{}}]


def read_vectors(path):
    """Read KML/KMZ, GeoJSON, or zipped shapefile using Fiona/GDAL."""
    import fiona
    ext = Path(path).suffix.lower()
    if ext == '.geojson' or ext == '.json':
        obj = json.loads(Path(path).read_text(encoding='utf-8'))
        return [(shape(f['geometry']),f.get('properties') or {}) for f in geo_features(obj) if f.get('geometry')], 'EPSG:4326'
    if ext == '.kmz':
        with zipfile.ZipFile(path) as z:
            names=[n for n in z.namelist() if n.lower().endswith('.kml')]
            if not names: raise ValueError('KMZ has no KML')
            with tempfile.TemporaryDirectory() as td:
                target=Path(td)/'route.kml'; target.write_bytes(z.read(names[0])); return read_vectors(target)
    if ext == '.kml':
        # Parse KML without requiring Fiona KML driver.
        from xml.etree import ElementTree as ET
        root=ET.parse(path).getroot(); out=[]
        for placemark in root.findall('.//{*}Placemark'):
            name=placemark.findtext('{*}name') or ''
            for tag in ('LineString','Point'):
                for el in placemark.findall('.//{*}'+tag):
                    raw=el.findtext('{*}coordinates') or ''
                    coords=[]
                    for token in raw.split():
                        parts=token.split(',')
                        if len(parts)>=2: coords.append((float(parts[0]),float(parts[1])))
                    if tag=='LineString' and len(coords)>=2: out.append((LineString(coords),{'name':name}))
                    if tag=='Point' and coords: out.append((Point(coords[0]),{'name':name}))
        return out,'EPSG:4326'
    if ext == '.zip':
        with zipfile.ZipFile(path) as z:
            names=z.namelist()
            if any(n.startswith('/') or '..' in Path(n).parts for n in names): raise ValueError('Unsafe archive')
            shp=[n for n in names if n.lower().endswith('.shp')]
            if len(shp)!=1: raise ValueError('Crossing ZIP must contain exactly one shapefile')
            with tempfile.TemporaryDirectory() as td:
                z.extractall(td)
                return read_vectors(Path(td)/shp[0])
    with fiona.open(str(path)) as ds:
        if not ds.crs: raise ValueError('Crossing shapefile CRS is missing (.prj required)')
        return [(shape(f['geometry']),dict(f['properties'])) for f in ds if f['geometry']],ds.crs_wkt or ds.crs


def project_features(features,source,target):
    t=Transformer.from_crs(source,target,always_xy=True)
    return [(transform(t.transform,g),a) for g,a in features]


def get_route(features):
    lines=[]
    for g,_ in features:
        if g.geom_type=='LineString': lines.append(g)
        elif g.geom_type=='MultiLineString': lines.extend(g.geoms)
    if not lines: raise ValueError('No pipeline line geometry found')
    # Multiple disconnected routes cannot be safely ordered automatically.
    if len(lines)!=1: raise ValueError('Select/export one continuous pipeline LineString per run')
    return lines[0]


def intersection_points(g):
    """Return every distinct point crossing; retain overlap segments for special handling.

    A linear overlap has no unique intersection chainage, so its entry and exit
    are recorded separately by crossings() rather than choosing an arbitrary
    midpoint. Geometry collections and multipart crossings are fully traversed.
    """
    if g.is_empty:
        return []
    if g.geom_type == 'Point':
        return [g]
    if g.geom_type == 'MultiPoint':
        return list(g.geoms)
    if g.geom_type in ('LineString', 'LinearRing'):
        coords = list(g.coords)
        return [Point(coords[0]), Point(coords[-1])] if len(coords) > 1 else []
    if hasattr(g, 'geoms'):
        return [p for part in g.geoms for p in intersection_points(part)]
    return []


def _attribute(attrs, field, fallback=''):
    """Shapefile field names are often truncated or differ in capitalization."""
    if not field:
        return fallback
    if field in attrs and attrs[field] not in (None, ''):
        return attrs[field]
    lookup = {str(k).casefold(): v for k, v in attrs.items()}
    value = lookup.get(str(field).casefold())
    return fallback if value in (None, '') else value


def _crossing_events(route, geom):
    """Return (chainage, event_kind) along the original route vertex order.

    Point intersections are exact. For an overlapping line, the beginning and
    end chainages are reported as separate events; a midpoint is not a crossing.
    Polygon crossings are reduced to their intersection with the pipeline and
    similarly represented by entry/exit events.
    """
    intersection = route.intersection(geom)
    if intersection.is_empty:
        return []
    events = []

    def visit(part):
        if part.is_empty:
            return
        kind = part.geom_type
        if kind == 'Point':
            events.append((float(route.project(part)), 'point'))
        elif kind == 'MultiPoint':
            for p in part.geoms:
                visit(p)
        elif kind in ('LineString', 'LinearRing'):
            coords = list(part.coords)
            if len(coords) >= 2:
                events.append((float(route.project(Point(coords[0]))), 'overlap_endpoint'))
                events.append((float(route.project(Point(coords[-1]))), 'overlap_endpoint'))
        elif hasattr(part, 'geoms'):
            for sub in part.geoms:
                visit(sub)
        else:
            # For unsupported geometry, fail explicitly instead of dropping it.
            raise ValueError(f'Unsupported intersection geometry: {kind}')

    visit(intersection)
    # Remove only duplicate points generated by the SAME crossing feature.
    # Distinct features at the same chainage must all be retained.
    events.sort(key=lambda item: item[0])
    unique = []
    for ch, kind in events:
        if unique and abs(ch - unique[-1][0]) < 0.001:
            if kind == 'overlap_endpoint':
                unique[-1] = (unique[-1][0], kind)
            continue
        unique.append((ch, kind))
    return unique


def crossings(route, features, id_field, description_field, reverse_crossing_chainage=False):
    """Compute crossing chainages with optional reverse crossing-only stationing.

    Each input crossing is intersected with the original ordered pipeline line.
    Multi-intersections are retained. A crossing coincident with the pipeline
    yields entry and exit chainages, not a misleading single midpoint.
    """
    if route.geom_type != 'LineString' or route.is_empty:
        raise ValueError('Crossing chainage requires one ordered pipeline LineString')
    found = []
    for i, (geom, attrs) in enumerate(features, 1):
        if geom is None or geom.is_empty:
            continue
        attrs = attrs or {}
        identifier = str(_attribute(attrs, id_field, _attribute(attrs, 'name', f'X-{i:03d}')))
        desc = str(_attribute(attrs, description_field, ''))
        events = _crossing_events(route, geom)
        for event_index, (ch, event_kind) in enumerate(events, 1):
            geometric_ch = max(0.0, min(float(route.length), ch))
            ch = (float(route.length) - geometric_ch) if reverse_crossing_chainage else geometric_ch
            ch = max(0.0, min(float(route.length), ch))
            label = f'CH: {kp(ch)} - {identifier} - {desc}'.rstrip(' -')
            if event_kind == 'overlap_endpoint' and len(events) > 1:
                label += ' - ' + (('EXIT' if event_index == 1 else 'ENTRY' if event_index == len(events) else 'BOUNDARY') if reverse_crossing_chainage else ('ENTRY' if event_index == 1 else 'EXIT' if event_index == len(events) else 'BOUNDARY'))
            found.append({
                'chainage_m': round(ch, 3), 'kp': kp(ch), 'id': identifier,
                'geometric_chainage_m': round(geometric_ch, 3),
                'reverse_crossing_chainage': bool(reverse_crossing_chainage),
                'description': desc, 'label': label,
                'overlap_midpoint': False,
                'event_type': event_kind,
                'source_feature_index': i,
            })
    return sorted(found, key=lambda item: (item['chainage_m'], item['source_feature_index']))


def sample_route(route,step):
    chainages=list(np.arange(0,route.length,step,dtype=float))+[route.length]
    return [(float(ch),route.interpolate(float(ch))) for ch in chainages]


def raster_elevations(samples,raster_path,crs):
    import rasterio
    from rasterio.warp import transform as convert
    with rasterio.open(raster_path) as ds:
        if not ds.crs: raise ValueError('DEM raster has no CRS')
        xs,ys=convert(crs,ds.crs,[p.x for _,p in samples],[p.y for _,p in samples])
        result=[]
        for vals in ds.sample(zip(xs,ys),masked=True):
            v=vals[0]
            result.append(None if np.ma.is_masked(v) or not np.isfinite(v) or (ds.nodata is not None and float(v)==ds.nodata) else round(float(v),3))
        return result


def online_elevations(samples,crs,source):
    """SRTM via OpenTopoData public endpoint. ASTER requires a configured compatible endpoint."""
    import requests,os
    if source=='aster':
        endpoint=os.getenv('ASTER_ELEVATION_URL','')
        if not endpoint: raise ValueError('ASTER online provider not configured. Set ASTER_ELEVATION_URL to a compatible OpenTopoData API endpoint, or upload ASTER GeoTIFF.')
    else: endpoint=os.getenv('SRTM_ELEVATION_URL','https://api.opentopodata.org/v1/srtm30m')
    t=Transformer.from_crs(crs,4326,always_xy=True)
    elevations=[]
    for start in range(0,len(samples),75):
        batch=samples[start:start+75]
        locations='|'.join(f'{t.transform(p.x,p.y)[1]:.7f},{t.transform(p.x,p.y)[0]:.7f}' for _,p in batch)
        r=requests.get(endpoint,params={'locations':locations},timeout=45);r.raise_for_status()
        payload=r.json()
        if len(payload.get('results',[]))!=len(batch):raise ValueError('Elevation provider returned incomplete results')
        elevations.extend([None if x.get('elevation') is None else float(x['elevation']) for x in payload['results']])
    return elevations


def make_dxf(samples,crossings,output,hscale,vscale):
    import ezdxf
    valid=[(s['chainage_m'],s['elevation_m']) for s in samples if s['elevation_m'] is not None]
    if not valid: raise ValueError('No valid elevations to export')
    datum=math.floor(min(y for _,y in valid)/5)*5-5
    exag=hscale/vscale
    doc=ezdxf.new('R2010');msp=doc.modelspace()
    for name,color in [('GROUND',3),('CROSSING_LINES',1),('CROSSING_LABELS',2),('PROFILE_GRID',8)]:doc.layers.new(name,dxfattribs={'color':color})
    segment=[]
    for s in samples:
        if s['elevation_m'] is None:
            if len(segment)>1:msp.add_lwpolyline(segment,dxfattribs={'layer':'GROUND'})
            segment=[]
        else:segment.append((s['chainage_m'],(s['elevation_m']-datum)*exag))
    if len(segment)>1:msp.add_lwpolyline(segment,dxfattribs={'layer':'GROUND'})
    top=max((y-datum)*exag for _,y in valid)+30
    for c in crossings:
        x=c['chainage_m'];msp.add_line((x,0),(x,top),dxfattribs={'layer':'CROSSING_LINES'})
        msp.add_text(c['label'],dxfattribs={'layer':'CROSSING_LABELS','height':2.5,'rotation':90,'insert':(x+1,top-4)})
    msp.add_line((0,0),(samples[-1]['chainage_m'],0),dxfattribs={'layer':'PROFILE_GRID'})
    doc.saveas(output)


def make_excel(samples,crossings,path):
    from openpyxl import Workbook
    wb=Workbook();ws=wb.active;ws.title='Ground Profile';ws.append(['Chainage (m)','KP','Elevation (m)','Easting','Northing'])
    for s in samples:ws.append([s['chainage_m'],s['kp'],s['elevation_m'],s['easting'],s['northing']])
    wc=wb.create_sheet('Crossings');wc.append(['Chainage (m)','KP','Crossing ID','Description','Label','Overlap midpoint'])
    for c in crossings:wc.append([c['chainage_m'],c['kp'],c['id'],c['description'],c['label'],c['overlap_midpoint']])
    for sheet in wb:
        sheet.freeze_panes='A2';sheet.auto_filter.ref=sheet.dimensions
        for col in sheet.columns:sheet.column_dimensions[col[0].column_letter].width=min(50,max(15,max(len(str(cell.value or '')) for cell in col)+2))
    wb.save(path)
