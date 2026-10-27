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
    if g.is_empty:return []
    if g.geom_type=='Point':return [g]
    if g.geom_type=='MultiPoint':return list(g.geoms)
    if g.geom_type in ('LineString','LinearRing'):
        return [g.interpolate(0.5,normalized=True)] # overlapping alignment: midpoint with warning
    if hasattr(g,'geoms'):
        return [p for x in g.geoms for p in intersection_points(x)]
    return []


def crossings(route,features,id_field,description_field,reverse_crossing_chainage=False):
    """Crossing events; optionally station from the opposite end of the route.

    Only crossing station values are changed; route geometry and ground samples
    are untouched. Set the option explicitly per project, not globally.
    """
    found=[]
    for i,(g,attrs) in enumerate(features,1):
        intersection=route.intersection(g)
        pts=intersection_points(intersection)
        for j,p in enumerate(pts,1):
            geometric_ch=route.project(p)
            ch=(route.length-geometric_ch) if reverse_crossing_chainage else geometric_ch
            ch=max(0.0,min(route.length,ch))
            identifier=str(attrs.get(id_field) or attrs.get('name') or f'X-{i:03d}')
            desc=str(attrs.get(description_field) or '')
            found.append({'chainage_m':round(ch,3),'kp':kp(ch),'id':identifier,'description':desc,
                          'label':f'CH: {kp(ch)} - {identifier} - {desc}'.rstrip(' -'),
                          'overlap_midpoint':intersection.geom_type in ('LineString','MultiLineString'),
                           'geometric_chainage_m':round(geometric_ch,3),
                           'reverse_crossing_chainage':bool(reverse_crossing_chainage)})
    return sorted(found,key=lambda x:x['chainage_m'])


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
