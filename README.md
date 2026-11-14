# Norilsk permafrost — data download (v0)

Question: is the ground under Norilsk moving in a widespread way because of
thaw (summers 2017–2026)? Three hypotheses, each with an observable signature:

| Hypothesis | Expected ground behaviour |
|---|---|
| Climatic thaw | Subsidence also in natural tundra far from structures; seasonal cycle + multi-year trend tied to temperature |
| Anthropogenic thaw | Subsidence concentrated around structures / heat sources; nearby tundra stable |
| Structural defect / corrosion | Ground stable (consistent with, not proof of) |

Design: industrial area vs natural tundra controls with similar elevation
and slope (not floodplain, not wetland); InSAR referenced to bedrock on the
Kharaelakh plateau NE of Talnakh, away from mines.

### Hypothesis test (what is actually measured)

Snow decorrelates C-band, so each summer is a separate InSAR stack, and
Sentinel-1 has **no acquisitions at all over Norilsk in 2022-2024**. A direct
multi-year subsidence trend is therefore **not observable**. The test is:

> **seasonal subsidence vs sqrt(cumulative TDD), per summer, per zone**
> (industrial / urban terraces vs tundra control), referenced to plateau rock.

- Climatic thaw: tundra and urban both follow sqrt(TDD); the slope (a proxy
  for active-layer ice content / thaw depth) grows across summers.
- Anthropogenic thaw: urban/industrial slope clearly larger than tundra in
  the same summer; tundra stable or small.
- Structural defect: no significant seasonal subsidence in either zone.
Annual (late-summer to late-summer) pairs are an optional extra, coherent only
on built-up and rock.

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
| 01 | `01_s1_search.py` | Sentinel-1 scenes + bursts per path (frame-agnostic) | no | CSV |
| 02 | `02_s1_hyp3.py plan/submit/status/download` | S1 interferograms (HyP3 INSAR_ISCE_MULTI_BURST, 10x2 looks = 40 m) | yes | ~100 MB/pair |
| 03 | `03_nisar.py search/download/subset` | NISAR L-band GUNW/GCOV, summer 2026 | download only | GB per frame; subset → MB |
| 04 | `04_climate_era5.py` | ERA5-Land air temperature, TDD, snow depth (soil T: not usable) | no | CSV |
| 05 | `05_dem_landcover.py` | Copernicus DEM 30 m, slope, WorldCover, reference/control masks | no | ~300 MB |
| 06 | `06_osm_infrastructure.py [--tanks]` | OSM landuse=industrial/quarry/landfill -> mask (+100 m); tanks/plants | no | GeoJSON + tif |
| 07 | `07_cci_permafrost.py` | ESA CCI ground temperature + active layer (2015–2023): **only ground-temperature source** | no | cropped NetCDF |
| 08 | `08_masks.py` | v2 reference + tundra-control masks (needs 05, 06) | no | tif + PNG |
| – | `quicklook.py`, `tank_amplitude.py` | QA plots; tank-location check from S1 amplitude | no | PNG |

Suggested order: 04, 05, 06, 08, 07 (no login) → 01 → 02 `plan` →
02 `submit --years 2020 2025` → 03.

```bash
python 01_s1_search.py                     # path 136 DESC, 6 bursts
python 02_s1_hyp3.py plan   --path 136 --years 2020 2025   # 30 pairs, 150 credits
python 02_s1_hyp3.py submit --path 136 --years 2020 2025
python 02_s1_hyp3.py status --path P       # wait until SUCCEEDED (~1 h)
python 02_s1_hyp3.py download --path P
```

## v0 analysis (one command)

```bash
python 02_s1_hyp3.py download --path 136     # after jobs SUCCEEDED
python 03_nisar.py download --match _111_A_037_ --match SH_20260621T --max 1   # x4
./run_v0.sh                # 09 time series -> 10 NISAR -> 11 analysis -> 12 web
./run_v0.sh --synthetic    # same chain on synthetic data with a known signal
```

| # | Script | Output |
|---|---|---|
| 09 | `09_timeseries.py` | numpy SBAS per summer -> `data/ts/s1_<year>_{los,vert}.tif` |
| 10 | `10_nisar.py` | NISAR LOS per pair, asc/desc -> vertical + east |
| 11 | `11_analysis.py` | zone x date table, sqrt(TDD) fits, comparisons, figures |
| 12 | `12_export_web.py` | `outputs/web/`: COGs (EPSG:4326), zones.geojson, CSV, data contract |

MintPy was not used for v0: a small least-squares SBAS is enough for
12/24-day summer networks and is verified end-to-end on synthetic data.

## Before analysis

- `config.py`: the HPP-3 tank location differs between sources (two points
  kept, unresolved; see `tank_amplitude.py`). `tundra_control` is still a
  placeholder: choose it from `tundra_candidate_v2.tif`.
- `plateau_reference` sits at the east edge of BBOX: most of the plateau is
  outside the AOI. Extending BBOX east (~88.9 E) is likely free in bursts.
- Mask v1 (step 05) put the reference on the mining massif south of Norilsk
  and WorldCover maps most industry as "bare": use the v2 masks.
- ERA5-Land soil temperature (1-2.5 m) is positive in continuous permafrost:
  **not usable** as evidence of ground warming. Use CCI (step 07).
  Open-Meteo `era5_land` has no snowfall; snow depth is used instead.
- HyP3 free quota: 8,000 credits/month; multi-burst 10x2 with 4-9 bursts =
  5 credits/pair. All 7 summers with data: 99 pairs ≈ 495 credits.
- Sentinel-1 path 136 DESC is the only usable stack (no ascending). Frame
  numbers drift (360 vs 358+363), hence burst processing. **No S1 data at
  all in 2022-2024** (S1B failed Dec 2021; S1A did not cover Norilsk);
  2025 = S1C, 2026 = S1C/D.
- NISAR PROVISIONAL products start 2026-06-17 → one summer, no trend. Two
  tracks cover the AOI fully: 111 ASC (frame 037) and 6 DESC (frame 053),
  so 2026 allows an ascending/descending (vertical vs east-west) check.
- Next step (time series): load `data/hyp3/` into MintPy
  (`load_data.py` with the `hyp3` processor), reference on `ref_candidate`.

## Evidence strength (state it in the submission)

1. Air warming (ERA5): strong. 2. Ground warming (CCI): medium, modelled,
   ~1 km, ends 2023 (ERA5-Land soil temperature not usable).
3. Ice-rich ground under structures: weak, low-res maps only.
4. Surface deformation (InSAR): strong, summer only, area-wide.
5. Differential settlement of a single tank: not observable (a tank is one
   bright point in radar).
