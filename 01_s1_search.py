"""Step 01: list Sentinel-1 IW SLC scenes over the AOI, summers 2017-2026.

No login needed. Output: data/s1/s1_scenes.csv + a summary of stacks
(same relative orbit + frame), so you can pick one stack for InSAR.

    python 01_s1_search.py
"""
import asf_search as asf
import pandas as pd
from shapely import wkt
from shapely.geometry import shape

from config import BBOX, DIRS, SUMMER, YEARS, bbox_wkt

aoi = wkt.loads(bbox_wkt(BBOX))
rows = []
for y in YEARS:
    res = asf.geo_search(
        dataset=asf.DATASET.SENTINEL1,
        processingLevel=asf.PRODUCT_TYPE.SLC,
        beamMode=asf.BEAMMODE.IW,
        intersectsWith=bbox_wkt(BBOX),
        start=f"{y}-{SUMMER[0]}",
        end=f"{y}-{SUMMER[1]}T23:59:59Z",
    )
    for p in res:
        pr = p.properties
        geom = shape(p.geometry)
        rows.append({
            "granule": pr["sceneName"],
            "time": pr["startTime"],
            "platform": pr["platform"],
            "path": pr["pathNumber"],
            "frame": pr["frameNumber"],
            "direction": pr["flightDirection"],
            "aoi_coverage": round(geom.intersection(aoi).area / aoi.area, 3),
        })
    print(f"{y}: {len(res)} scenes")

df = pd.DataFrame(rows)
if df.empty:
    raise SystemExit("No scenes found: check BBOX / dates.")
df["time"] = pd.to_datetime(df["time"])
df["date"] = df["time"].dt.date
df["year"] = df["time"].dt.year
df = df.sort_values("time")
out = DIRS["s1"] / "s1_scenes.csv"
df.to_csv(out, index=False)

full = df[df["aoi_coverage"] >= 0.95]
summary = (full.groupby(["path", "frame", "direction", "year"])
           .size().unstack("year", fill_value=0))
summary["total"] = summary.sum(axis=1)
summary = summary.sort_values("total", ascending=False)
summary.to_csv(DIRS["s1"] / "s1_stacks_summary.csv")

print(f"\nSaved {len(df)} scenes -> {out}")
print("\nStacks covering >=95% of the AOI (scenes per summer):")
print(summary.head(10).to_string())
if not summary.empty:
    best = summary.index[0]
    print(f"\nSuggested stack: path={best[0]} frame={best[1]} ({best[2]})")
    print(f"Next: python 02_s1_hyp3.py submit --path {best[0]} --frame {best[1]}")
else:
    print("\nNo single frame covers the AOI: shrink BBOX or use bursts.")
