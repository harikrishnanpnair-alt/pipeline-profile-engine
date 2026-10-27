"""Standalone engineering profile API. Does not alter the existing GDB backend."""
import os,uuid,tempfile,shutil,zipfile
from pathlib import Path
from fastapi import FastAPI,UploadFile,File,Form,HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from engine import read_vectors,project_features,get_route,crossings,sample_route,raster_elevations,online_elevations,kp,make_dxf,make_excel

app=FastAPI(title='Pipeline Engineering Profile Engine',version='0.1.0')
origins=[x.strip() for x in os.getenv('ALLOWED_ORIGINS','http://localhost:8000,http://127.0.0.1:8000').split(',') if x.strip()]
app.add_middleware(CORSMiddleware,allow_origins=origins,allow_methods=['GET','POST'],allow_headers=['*'])
ROOT=Path(tempfile.gettempdir())/'pipeline_profile_jobs';ROOT.mkdir(exist_ok=True)
MAX_BYTES=int(os.getenv('MAX_UPLOAD_MB','100'))*1024*1024

@app.get('/health')
def health():return {'ready':True,'service':'pipeline-profile-engine','version':'0.1.0'}

async def save_upload(file,directory,allowed):
    if file is None:return None
    suffix=Path(file.filename or '').suffix.lower()
    if suffix not in allowed:raise HTTPException(400,f'Unsupported file extension: {suffix}')
    target=directory/('input_'+uuid.uuid4().hex+suffix)
    size=0
    with target.open('wb') as f:
        while chunk:=await file.read(1024*1024):
            size+=len(chunk)
            if size>MAX_BYTES:raise HTTPException(413,'Upload too large')
            f.write(chunk)
    return target

@app.post('/profile/generate')
async def generate(route:UploadFile=File(...),crossing:UploadFile|None=File(None),surface:UploadFile|None=File(None),
                   source:str=Form('srtm'),epsg:int=Form(32640),interval:float=Form(25),
                   id_field:str=Form('Crossing_ID'),description_field:str=Form('Description'),
                   hscale:int=Form(2500),vscale:int=Form(100),
                   reverse_crossing_chainage:bool=Form(True)):
    if source not in ('srtm','aster','dem'):raise HTTPException(400,'Source must be srtm, aster or dem; TIN conversion is not yet implemented')
    if epsg not in (32639,32640,27039,27040):raise HTTPException(400,'Unsupported projected CRS')
    if not 1<=interval<=1000:raise HTTPException(400,'Interval must be 1–1000 m')
    if not 100<=hscale<=20000 or not 10<=vscale<=2000:raise HTTPException(400,'Invalid drawing scale')
    job=uuid.uuid4().hex;directory=ROOT/job;directory.mkdir()
    try:
        r=await save_upload(route,directory,{'.kmz','.kml','.geojson','.json','.shp'})
        c=await save_upload(crossing,directory,{'.zip','.geojson','.json','.kml','.kmz'})
        d=await save_upload(surface,directory,{'.tif','.tiff','.img'})
        if source=='dem' and d is None:raise ValueError('Upload a georeferenced DEM')
        rfeatures,rcrs=read_vectors(r);line=get_route(project_features(rfeatures,rcrs,f'EPSG:{epsg}'))
        if line.length>500000:raise ValueError('Limit each profile run to 500 km')
        cs=[]
        if c:
            features,ccrs=read_vectors(c)
            cs=crossings(line,project_features(features,ccrs,f'EPSG:{epsg}'),id_field,description_field,
                         reverse_crossing_chainage=reverse_crossing_chainage)
        samples=sample_route(line,interval)
        if len(samples)>20001:raise ValueError('Too many samples; increase interval')
        elev=raster_elevations(samples,d,f'EPSG:{epsg}') if source=='dem' else online_elevations(samples,f'EPSG:{epsg}',source)
        rows=[{'chainage_m':round(ch,3),'kp':kp(ch),'elevation_m':z,'easting':round(p.x,3),'northing':round(p.y,3)} for (ch,p),z in zip(samples,elev)]
        if not any(s['elevation_m'] is not None for s in rows):raise ValueError('No elevations returned within surface coverage')
        make_dxf(rows,cs,directory/'Engineering_Profile.dxf',hscale,vscale)
        make_excel(rows,cs,directory/'Engineering_Profile.xlsx')
        import json
        result={'job_id':job,'route_length_m':round(line.length,3),'source':source,'epsg':epsg,'reverse_crossing_chainage':reverse_crossing_chainage,'samples':rows,'crossings':cs,
                'downloads':{'dxf':f'/profile/download/{job}/Engineering_Profile.dxf','excel':f'/profile/download/{job}/Engineering_Profile.xlsx'},
                'note':'Ground elevations from terrain sources are approximate; no pipe design level is inferred.'}
        (directory/'result.json').write_text(json.dumps(result),encoding='utf-8')
        return result
    except HTTPException:raise
    except Exception as e:raise HTTPException(400,str(e))

@app.get('/profile/download/{job}/{filename}')
def download(job:str,filename:str):
    if len(job)!=32 or any(c not in '0123456789abcdef' for c in job) or filename not in ('Engineering_Profile.dxf','Engineering_Profile.xlsx'):
        raise HTTPException(400,'Invalid download')
    p=ROOT/job/filename
    if not p.is_file():raise HTTPException(404,'File unavailable')
    return FileResponse(p,filename=filename)
