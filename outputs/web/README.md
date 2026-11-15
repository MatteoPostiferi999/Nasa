# Norilsk v0 - data contract for the interface

All rasters: Cloud-Optimized GeoTIFF, EPSG:4326, float32, NaN = no data.
Units mm. Sign: vertical + = uplift, - = subsidence; east + = eastward.
Reference: mean of the plateau bedrock patches (zone code 1) = 0.

| file | content |
|---|---|
| s1_seasonal_2020.tif | S1 vertical displacement 20200615 -> 20200919 |
| s1_seasonal_2025.tif | S1 vertical displacement 20250613 -> 20250905 |
| nisar_early_vertical.tif | NISAR 2026 early window, vertical component |
| nisar_late_vertical.tif | NISAR 2026 late window, vertical component |
| nisar_early_east.tif | NISAR 2026 early window, east-west component |
| nisar_late_east.tif | NISAR 2026 late window, east-west component |
| zones.tif | uint8 zone code (see below), 0 = outside all zones |
| zones.geojson | one MultiPolygon per zone, properties `zone`, `code` |
| zone_timeseries.csv | zone, summer, date, disp_mm (zone mean), disp_median_mm, se_mm, n_px, tdd_cum |
| s2_spill_2020*.png | Sentinel-2 true colour of the 2020 spill area (no georeference; context image) |

Zone codes: 1 = reference, 2 = industrial, 3 = accident, 4 = urban, 5 = control_NW, 6 = control_SE

Notes for display
- S1 seasonal = displacement between the first and last summer acquisition,
  only on slopes < 5 deg (LOS converted to vertical assuming no horizontal motion).
- NISAR vertical/east come from one asc + one desc 12-day pair (dates in
  nisar_meta.json); they are a cross-check, not a time series.
- tdd_cum: cumulative thawing degree days (degC day) since 1 Jan, ERA5-Land,
  at the zone's climate point.
- Credit: contains modified Copernicus Sentinel data (2020, 2025); NISAR data
  courtesy NASA/ISRO via ASF DAAC; ERA5-Land (Copernicus C3S) via Open-Meteo;
  OpenStreetMap contributors (zones); ESA WorldCover 2021; Copernicus DEM.
