"""Step 08: InSAR reference and tundra-control masks, v2 (needs 05 + 06).

Writes *_v2.tif next to the step-05 masks (v1 kept for comparison) and
outputs/masks_v2.png.

Reference (plateau bedrock NE of Talnakh, away from mines):
  WorldCover bare (60), elevation > 300 m, slope < 10 deg,
  inside REF_SECTOR, > 3 km from OSM industrial/quarry/landfill,
  patches >= 0.1 km2.
Tundra control (natural ground comparable to the urban terraces):
  WorldCover shrub/grass/moss-lichen (20, 30, 100); no wetland (90),
  water (80) or trees (10); elevation and slope within the urban-terrace
  range (WorldCover built-up p5-p95 elevation, p90 slope); height above
  nearest large water body (>1 km2: rivers, big lakes) >= HAND_MIN m,
  i.e. not floodplain; > 2 km from built-up (WorldCover 50 or OSM land use);
  patches >= 0.5 km2.

    python 08_masks.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import ListedColormap
from pyproj import Transformer
from scipy import ndimage as nd

from config import DIRS, POINTS, UTM_EPSG

RES = 30.0
PX_KM2 = (RES / 1000) ** 2
BUILT, BARE, WATER, WETLAND, TREE = 50, 60, 80, 90, 10
TUNDRA_CLASSES = [20, 30, 100]
REF_SECTOR = (88.40, 69.45, 88.65, 69.58)  # lon/lat box: Kharaelakh plateau NE of Talnakh
REF_MINE_DIST_KM = 3.0
TUNDRA_BUILT_DIST_KM = 2.0
HAND_MIN = 5.0          # m above nearest large water body
BIG_WATER_KM2 = 1.0
OUT = Path(__file__).parent / "outputs"


def read(path):
    with rasterio.open(path) as s:
        return s.read(1), s.profile


def write(path, arr, prof):
    with rasterio.open(path, "w", **dict(prof, dtype="uint8", nodata=None)) as d:
        d.write(arr.astype("uint8"), 1)


def min_patch(mask, km2):
    lab, n = nd.label(mask)
    if n == 0:
        return mask
    size = nd.sum(mask, lab, range(1, n + 1)) * PX_KM2
    return np.isin(lab, 1 + np.where(size >= km2)[0])


def lonlat_grid(prof):
    h, w = prof["height"], prof["width"]
    cols, rows = np.meshgrid(np.arange(w) + 0.5, np.arange(h) + 0.5)
    x, y = prof["transform"] * (cols, rows)
    return Transformer.from_crs(UTM_EPSG, 4326, always_xy=True).transform(x, y)


def main():
    dem, prof = read(DIRS["dem"] / "dem_utm30.tif")
    slope, _ = read(DIRS["dem"] / "slope_deg_utm30.tif")
    lc, _ = read(DIRS["landcover"] / "worldcover_utm30.tif")
    osm, _ = read(DIRS["osm"] / "osm_landuse_mask.tif")
    valid = dem > 0
    osm = osm.astype(bool)

    # urban terraces: main WorldCover built-up patches (> 0.5 km2)
    urban = min_patch(nd.binary_closing(lc == BUILT, iterations=3), 0.5) & (lc == BUILT)
    z_lo, z_hi = np.percentile(dem[urban], [5, 95])
    s_hi = np.percentile(slope[urban], 90)

    # height above nearest large water body (HAND proxy)
    big_water = min_patch(lc == WATER, BIG_WATER_KM2)
    _, (iy, ix) = nd.distance_transform_edt(~big_water, return_indices=True)
    hand = dem - dem[iy, ix]

    km = lambda m: nd.distance_transform_edt(~m) * RES / 1000
    dist_built = km((lc == BUILT) | osm)
    dist_osm = km(osm)
    lon, lat = lonlat_grid(prof)

    in_sector = ((lon >= REF_SECTOR[0]) & (lat >= REF_SECTOR[1])
                 & (lon <= REF_SECTOR[2]) & (lat <= REF_SECTOR[3]))
    ref = (valid & (lc == BARE) & (dem > 300) & (slope < 10) & in_sector
           & (dist_osm > REF_MINE_DIST_KM))
    ref = min_patch(ref, 0.1)

    tundra = (valid & np.isin(lc, TUNDRA_CLASSES)
              & (dem >= z_lo) & (dem <= z_hi) & (slope <= s_hi)
              & (hand >= HAND_MIN) & (dist_built > TUNDRA_BUILT_DIST_KM) & ~osm)
    tundra = min_patch(tundra, 0.5)

    write(DIRS["landcover"] / "ref_candidate_v2.tif", ref, prof)
    write(DIRS["landcover"] / "tundra_candidate_v2.tif", tundra, prof)

    # suggested reference point: centre pixel of the largest ref patch
    lab, n = nd.label(ref)
    sizes = nd.sum(ref, lab, range(1, n + 1))
    big = 1 + int(np.argmax(sizes))
    rr, cc = np.where(lab == big)
    k = np.argmin((rr - rr.mean()) ** 2 + (cc - cc.mean()) ** 2)
    r0, c0 = rr[k], cc[k]

    print(f"urban terraces: elev {z_lo:.0f}-{z_hi:.0f} m (p5-p95), slope p90 {s_hi:.1f} deg")
    print(f"ref v2: {ref.sum() * PX_KM2:.1f} km2 in {n} patches; largest "
          f"{sizes.max() * PX_KM2:.2f} km2, elev {dem[lab == big].mean():.0f} m")
    print(f"  suggested plateau_reference = ({lat[r0, c0]:.4f}, {lon[r0, c0]:.4f}), "
          f"{dist_osm[r0, c0]:.1f} km from nearest OSM industrial/quarry")
    print(f"tundra v2: {tundra.sum() * PX_KM2:.1f} km2 "
          f"(v1 was {read(DIRS['landcover'] / 'tundra_candidate.tif')[0].sum() * PX_KM2:.1f})")
    lt, nt = nd.label(tundra)
    print(f"  {nt} patches; classes %: " + ", ".join(
        f"{c}:{np.mean(lc[tundra] == c) * 100:.0f}" for c in TUNDRA_CLASSES))

    # figure
    OUT.mkdir(exist_ok=True)
    b = rasterio.transform.array_bounds(prof["height"], prof["width"], prof["transform"])
    ext = (b[0], b[2], b[1], b[3])
    to_utm = Transformer.from_crs(4326, UTM_EPSG, always_xy=True)
    fig, ax = plt.subplots(figsize=(10, 9))
    ax.imshow(np.ma.masked_where(~valid, dem), extent=ext, cmap="Greys_r")
    for m, col in [(osm, "orange"), (lc == BUILT, "yellow"),
                   (tundra, "green"), (ref, "red")]:
        ax.imshow(np.ma.masked_where(~m, m), extent=ext, cmap=ListedColormap([col]),
                  alpha=0.75, interpolation="nearest")
    pts = dict(POINTS, plateau_reference_NEW=(lat[r0, c0], lon[r0, c0]))
    for name, (la, lo) in pts.items():
        x, y = to_utm.transform(lo, la)
        ax.plot(x, y, "o", ms=7, mfc="magenta", mec="black")
        ax.annotate(name, (x, y), xytext=(5, 5), textcoords="offset points", fontsize=8,
                    bbox=dict(fc="white", alpha=0.7, lw=0, pad=1))
    ax.set_title("v2 masks: ref (red), tundra control (green), "
                 "OSM industrial/quarry/landfill +100 m (orange), WorldCover built (yellow)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "masks_v2.png", dpi=130)


if __name__ == "__main__":
    main()
