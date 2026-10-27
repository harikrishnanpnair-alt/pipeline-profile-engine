# Sy-175 corrected pipeline profile backend

Replace **both** `main.py` and `engine.py` in `pipeline-profile-engine`, then redeploy Render.

## Important change

The original Sy-175 KMZ's last vertex is the engineering start near RX-01.
The engine now orients the **single route** before generating BOTH ground profile samples and crossing chainages.
No independent crossing-only mirroring is performed.

The API accepts `route_direction=reverse` (default for this package, Sy-175) or `route_direction=forward` (for other KMZs whose first vertex is CH 0+000).

POST `/profile/generate` multipart form fields include:
- `route`: KMZ/KML/GeoJSON
- `crossing`: ZIP shapefile/GeoJSON (optional)
- `route_direction`: `reverse` or `forward` (default `reverse`)
- `id_field`: `CrossingNa` (default)
- `description_field`: `CrossingTy` (default)
- existing `epsg`, `source`, `interval`, `hscale`, `vscale`, `surface` unchanged

The response now includes `route_direction`, `start_easting`, `start_northing`, `end_easting`, and `end_northing` for diagnostics.

## Verified with your supplied files
- Pipeline_Sy_175.kmz
- Crossings_Line.zip
- EPSG:27039
- Route length approximately 3,847.459 m
- RX-01 near CH 0+009.558 and CH 0+015.558 (overlap endpoints)

The original KMZ is not edited. Ground elevations require a live DEM/online source and were not fetched in this offline test.

**Other projects:** Set `route_direction=forward` if their KMZ already starts at engineering CH 0+000. Do not send any old `reverse_crossing_chainage` flag.
