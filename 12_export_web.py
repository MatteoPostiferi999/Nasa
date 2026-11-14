"""Step 12: web-ready files for the interface (outputs/web/).

  s1_seasonal_<year>.tif       COG EPSG:4326, float32 mm, + up, last date - first date
  nisar_<window>_vertical.tif  COG EPSG:4326, mm, + up   (12-day pair)
  nisar_<window>_east.tif      COG EPSG:4326, mm, + east
  zones.tif                    COG EPSG:4326, uint8 zone codes (0 = none)
  zones.geojson                zone polygons (EPSG:4326), property zone/code
  zone_timeseries.csv          written by step 11
  README.md                    data contract
Also outputs/fig_nisar.png (asc LOS, desc LOS, vertical, east per window).

    python 12_export_web.py
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio.features import shapes
from rasterio.warp import Resampling, calculate_default_transform, reproject
from shapely.geometry import mapping, shape
from shapely.ops import unary_union

from common import ZONE_NAMES
from config import DIRS, OUTPUTS

WEB = OUTPUTS / "web"


def to_cog_4326(src_path, dst_path, band=None, resampling=Resampling.bilinear):
    with rasterio.open(src_path) as s:
        idx = band or s.count
        tr, w, h = calculate_default_transform(s.crs, "EPSG:4326", s.width, s.height, *s.bounds)
        dtype = s.dtypes[idx - 1]
        nodata = s.nodata if s.nodata is not None else (np.nan if "float" in dtype else 0)
        out = np.full((h, w), nodata, dtype)
        reproject(rasterio.band(s, idx), out, dst_transform=tr, dst_crs="EPSG:4326",
                  dst_nodata=nodata, resampling=resampling)
        desc = s.descriptions[idx - 1]
    prof = dict(driver="COG", width=w, height=h, count=1, dtype=dtype, crs="EPSG:4326",
                transform=tr, nodata=nodata, compress="deflate")
    with rasterio.open(dst_path, "w", **prof) as d:
        d.write(out, 1)
        if desc:
            d.set_band_description(1, desc)
    return dst_path


def seasonal_tifs():
    out = []
    for f in sorted(DIRS["ts"].glob("s1_*_vert.tif")):
        year = f.stem.split("_")[1]
        with rasterio.open(f) as s:
            d0, d1 = s.descriptions[0], s.descriptions[-1]
        dst = to_cog_4326(f, WEB / f"s1_seasonal_{year}.tif")   # last band = cumulative
        out.append({"file": dst.name, "from": d0, "to": d1})
    return out


def zones_files():
    src = DIRS["zones"] / "zones_utm30.tif"
    to_cog_4326(src, WEB / "zones.tif", resampling=Resampling.nearest)
    with rasterio.open(WEB / "zones.tif") as s:
        z = s.read(1)
        feats = {}
        for geom, val in shapes(z, mask=z > 0, transform=s.transform):
            feats.setdefault(int(val), []).append(shape(geom))
    fc = {"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {"code": c, "zone": ZONE_NAMES[c]},
         "geometry": mapping(unary_union(g).buffer(0.0003).buffer(-0.0003).simplify(0.0004))}
        for c, g in sorted(feats.items())]}
    (WEB / "zones.geojson").write_text(json.dumps(fc))


def nisar_files():
    out = []
    for f in sorted(DIRS["ts"].glob("nisar_*_vertical.tif")) + sorted(DIRS["ts"].glob("nisar_*_east.tif")):
        to_cog_4326(f, WEB / f.name)
        out.append(f.name)
    los = sorted(DIRS["ts"].glob("nisar_*_los.tif"))
    wins = ["early", "late"]
    if not los:
        return out
    fig, axs = plt.subplots(len(wins), 4, figsize=(18, 4.5 * len(wins)), squeeze=False)
    for r, win in enumerate(wins):
        panels = []
        for tag, label in [("111A", "asc 111 LOS"), ("6D", "desc 6 LOS")]:
            cand = [p for p in los if f"nisar_{tag}_" in p.name]
            cand = cand[:1] if win == "early" else cand[1:2]
            panels.append((cand[0] if cand else None, label))
        panels += [(DIRS["ts"] / f"nisar_{win}_vertical.tif", "vertical (+ up)"),
                   (DIRS["ts"] / f"nisar_{win}_east.tif", "east (+ east)")]
        for ax, (p, label) in zip(axs[r], panels):
            ax.set_xticks([]); ax.set_yticks([])
            if p is None or not p.exists():
                ax.set_title(f"{win}: {label} (missing)", fontsize=9)
                continue
            with rasterio.open(p) as s:
                a = s.read(1)
            lim = np.nanpercentile(np.abs(a), 98) if np.isfinite(a).any() else 1
            im = ax.imshow(a, cmap="RdBu", vmin=-lim, vmax=lim)
            ax.set_title(f"{win}: {label} [mm] {p.stem.replace('nisar_', '')}", fontsize=8)
            fig.colorbar(im, ax=ax, shrink=0.7)
    fig.suptitle("NISAR 2026 L-band GUNW (12-day pairs), referenced to plateau patches", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUTPUTS / "fig_nisar.png", dpi=120)
    plt.close(fig)
    return out


README = """# Norilsk v0 - data contract for the interface

All rasters: Cloud-Optimized GeoTIFF, EPSG:4326, float32, NaN = no data.
Units mm. Sign: vertical + = uplift, - = subsidence; east + = eastward.
Reference: mean of the plateau bedrock patches (zone code 1) = 0.

| file | content |
|---|---|
{rows}
| zones.tif | uint8 zone code (see below), 0 = outside all zones |
| zones.geojson | one MultiPolygon per zone, properties `zone`, `code` |
| zone_timeseries.csv | zone, summer, date, disp_mm (zone mean), disp_median_mm, se_mm, n_px, tdd_cum |
| s2_spill_2020*.png | Sentinel-2 true colour of the 2020 spill area (no georeference; context image) |

Zone codes: {codes}

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
"""


def main():
    WEB.mkdir(parents=True, exist_ok=True)
    seas = seasonal_tifs()
    zones_files()
    nis = nisar_files()
    rows = [f"| {s['file']} | S1 vertical displacement {s['from']} -> {s['to']} |" for s in seas]
    rows += [f"| {n} | NISAR 2026 {n.split('_')[1]} window, "
             f"{'vertical' if 'vertical' in n else 'east-west'} component |" for n in nis]
    codes = ", ".join(f"{c} = {n}" for c, n in ZONE_NAMES.items())
    (WEB / "README.md").write_text(README.format(rows="\n".join(rows), codes=codes))
    print("web files:", sorted(p.name for p in WEB.iterdir()))


if __name__ == "__main__":
    main()
