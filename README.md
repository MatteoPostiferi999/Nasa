# Norilsk permafrost — data download (v0)

Question: is the ground under Norilsk moving in a widespread way because of
thaw (summers 2017–2026)? Three hypotheses, each with an observable signature:

| Hypothesis | Expected ground behaviour |
|---|---|
| Climatic thaw | Subsidence also in natural tundra far from structures; seasonal cycle + multi-year trend tied to temperature |
| Anthropogenic thaw | Subsidence concentrated around structures / heat sources; nearby tundra stable |
| Structural defect / corrosion | Ground stable (consistent with, not proof of) |

Design: industrial area vs natural tundra controls with similar soil and
slope; InSAR referenced to bedrock on the plateau, not the lowland.

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Earthdata Login (free, <https://urs.earthdata.nasa.gov>), needed for HyP3 and
NISAR downloads. Put it in `~/.netrc`, never in the code:

```
machine urs.earthdata.nasa.gov login YOUR_USER password YOUR_PASS
```
`chmod 600 ~/.netrc`. Alternative: env vars `EARTHDATA_USERNAME` / `EARTHDATA_PASSWORD`.

## Steps

| # | Script | Data | Login | Size |
|---|---|---|---|---|
| 01 | `01_s1_search.py` | Sentinel-1 SLC scene list, picks the best stack | no | CSV |
| 02 | `02_s1_hyp3.py plan/submit/status/download` | S1 interferograms (HyP3 INSAR_GAMMA, 40 m) | yes | a few hundred MB/pair |
| 03 | `03_nisar.py search/download/subset` | NISAR L-band GUNW/GCOV, summer 2026 | download only | GB per frame; subset → MB |
| 04 | `04_climate_era5.py` | ERA5-Land air/soil temperature, snow, TDD | no | CSV |
| 05 | `05_dem_landcover.py` | Copernicus DEM 30 m, slope, WorldCover, reference/control masks | no | ~300 MB |
| 06 | `06_osm_infrastructure.py` | OSM buildings, tanks, pipelines, industrial land | no | GeoJSON |
| 07 | `07_cci_permafrost.py` | ESA CCI ground temperature + active layer (2015–2023) | no | cropped NetCDF |

Suggested order for a first version: 04, 05, 06 (fast, no login) → 01 →
02 `plan` → 02 `submit --years 2020 2025` → 03 → 07.

```bash
python 01_s1_search.py                     # prints suggested --path/--frame
python 02_s1_hyp3.py plan   --path P --frame F --years 2020 2025
python 02_s1_hyp3.py submit --path P --frame F --years 2020 2025
python 02_s1_hyp3.py status --path P       # wait until SUCCEEDED (~1 h)
python 02_s1_hyp3.py download --path P
```

## Before analysis

- `config.py`: the HPP-3 tank location differs between sources (two points
  kept). Check on Google Earth/OSM. Control and reference points are
  placeholders: choose them from `ref_candidate.tif` / `tundra_candidate.tif`.
- HyP3 free quota: 8,000 credits/month; 40 m InSAR = 15 credits/pair.
  All 10 summers ≈ 150–200 pairs ≈ 2,500–3,000 credits.
- Sentinel-1B failed in Dec 2021: 2022–2024 has 12-day revisit only.
- NISAR PROVISIONAL products start 2026-06-17 → one summer, no trend.
- Next step (time series): load `data/hyp3/` into MintPy
  (`load_data.py` with the `hyp3` processor), reference on `ref_candidate`.

## Evidence strength (state it in the submission)

1. Air warming (ERA5): strong. 2. Ground warming (CCI): medium, modelled.
3. Ice-rich ground under structures: weak, low-res maps only.
4. Surface deformation (InSAR): strong, summer only, area-wide.
5. Differential settlement of a single tank: not observable (a tank is one
   bright point in radar).
