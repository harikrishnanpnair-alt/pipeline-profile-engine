# Pipeline Engineering Profile Engine — initial MVP

Independent FastAPI backend for the existing GIS Processing Portal. Does not replace `gdb_qc_cloud`.

## Deploy
Create a GitHub repository named `pipeline-profile-engine`, upload these files to its root, and create a **Render Docker Web Service** pointing to it. Configure `ALLOWED_ORIGINS=https://harikrishnanpnair-alt.github.io` and, optionally, `SRTM_ELEVATION_URL`. Start command is in Dockerfile. Visit `/health` and `/docs` to test.

## Current features
- Route KMZ/KML/GeoJSON upload (single continuous LineString).
- Crossing zipped shapefile (with `.shp`, `.shx`, `.dbf`, `.prj`) or GeoJSON.
- WGS84/Nahrwan UTM zones 39/40, chainage and intersection labels.
- SRTM online using public OpenTopoData SRTM30m endpoint (subject to external rate limits).
- ASTER online **only when** a compatible elevation API is configured in `ASTER_ELEVATION_URL`; otherwise upload an ASTER GeoTIFF using DEM mode.
- Uploaded georeferenced DEM `.tif`, `.tiff`, `.img`.
- AutoCAD DXF with green ground line, red crossing lines, yellow vertical labels; Excel schedule.

## Limitations / planned
- **TIN is not implemented in this MVP**. Convert LandXML/TIN to raster DEM first, or add a dedicated triangulation module.
- Crossing line overlap is represented by midpoint and flagged; polygon intersections may need entry/exit handling.
- DXF is a starting profile, not yet a fully templated engineering-data-band drawing. Web portal integration is pending.
- Uploaded files stored in ephemeral temp storage, without scheduled cleanup. Use restricted access for production.
- SRTM/ASTER terrain is not survey-grade; verify vertical datum and units.
- No built-in authentication; protect service before public production deployment.

## API
`POST /profile/generate` multipart: `route`, optional `crossing`, optional `surface`, `source=srtm|aster|dem`, `epsg=32640`, `interval=25`, `id_field=Crossing_ID`, `description_field=Description`, `hscale=2500`, `vscale=100`. Returns samples, crossings, and DXF/Excel download paths.

## Run locally
`pip install -r requirements.txt "fiona>=1.9,<2"` then `uvicorn main:app --reload`.
