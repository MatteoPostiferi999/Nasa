"""Step 10: NISAR GUNW (L-band) -> LOS displacement, then vertical + east.

Input: the 4 GUNW .h5 files from step 03 (111 ASC and 6 DESC, early + late
summer 2026). Only the AOI window is read from each file.

Per file:
  unwrapped phase (frequencyA, HH) -> LOS displacement (mm)
      d_los = SIGN * phase * wavelength / (4 pi),  wavelength = c / centerFrequency
      SIGN = -1: positive = motion toward the satellite (same as step 09).
      ASSUMPTION to verify against the NISAR GUNW product spec: if the city
      shows uplift in summer while S1 shows subsidence, flip SIGN.
  coherence < COH_MIN masked; referenced to the mean of the reference patches.
  LOS unit vector (target -> sensor, east/north components) from the
  metadata radarGrid cube, at the layer closest to the local DEM height.

Asc/desc decomposition: for each window (early, late) the asc and desc pairs
(12 days each, ~5 days offset) are combined per pixel, neglecting north
motion (LOS is nearly insensitive to it):
    d_asc  = e_a * dE + u_a * dU
    d_desc = e_d * dE + u_d * dU
Key check: horizontal (east) vs vertical on flat ground (slope < 5 deg).

Outputs (data/ts/): nisar_<track><dir>_<ref>_<sec>_los.tif,
nisar_<window>_vertical.tif, nisar_<window>_east.tif, nisar_meta.json.

    python 10_nisar.py
"""
from pathlib import Path

import h5py
import numpy as np
from pyproj import CRS, Transformer
from rasterio.transform import from_origin
from scipy.ndimage import map_coordinates

from common import grid, read, save_json, to_grid, write_stack, zones
from config import BBOX, DIRS

C = 299_792_458.0
SIGN = -1.0
COH_MIN = 0.3
SLOPE_MAX = 5.0
WINDOWS = {"early": ("SH_20260621T", "SH_20260626T"),
           "late": ("SH_20260820T", "SH_20260825T")}   # (asc 111, desc 6) ref dates


def find(f, suffix, must=""):
    hits = []
    f.visit(lambda n: hits.append(n) if n.endswith(suffix) and must in n else None)
    return sorted(hits, key=len)


def parse_name(name):
    f = name.split("_")
    return {"track": int(f[5]), "dir": f[6], "ref": f[11][:8], "sec": f[13][:8]}


def epsg_of(f, group):
    """EPSG of a geocoded grid group: 'projection' dataset/attr in it or above."""
    g = group
    while g:
        if f"{g}/projection" in f:
            p = f[f"{g}/projection"]
            for k in ("epsg_code", "epsg"):
                if k in p.attrs:
                    return int(p.attrs[k])
            return int(p[()])
        g = g.rsplit("/", 1)[0] if "/" in g else ""
    raise KeyError(f"no projection above {group}")


def coords_of(f, group):
    g = group
    while g:
        if f"{g}/xCoordinates" in f:
            return f[f"{g}/xCoordinates"][:], f[f"{g}/yCoordinates"][:], g
        g = g.rsplit("/", 1)[0] if "/" in g else ""
    raise KeyError(f"no coordinates above {group}")


def aoi_slices(x, y, epsg):
    tr = Transformer.from_crs(4326, epsg, always_xy=True)
    xs, ys = tr.transform([BBOX[0], BBOX[2], BBOX[0], BBOX[2]],
                          [BBOX[1], BBOX[1], BBOX[3], BBOX[3]])
    ix = np.where((x >= min(xs)) & (x <= max(xs)))[0]
    iy = np.where((y >= min(ys)) & (y <= max(ys)))[0]
    if ix.size == 0 or iy.size == 0:
        raise ValueError("AOI not inside this product")
    return slice(ix[0], ix[-1] + 1), slice(iy[0], iy[-1] + 1)


def los_vectors(f, xg, yg, height):
    """East/north/up LOS components (target -> sensor) on the product grid."""
    base = find(f, "radarGrid/losUnitVectorX")
    if not base:
        raise KeyError("losUnitVectorX not found")
    rg = base[0].rsplit("/", 1)[0]
    cx, cy = f[f"{rg}/xCoordinates"][:], f[f"{rg}/yCoordinates"][:]
    h = f[f"{rg}/heightAboveEllipsoid"][:]
    k = int(np.argmin(np.abs(h - np.nanmedian(height))))   # one layer: AOI relief << cube spacing effect
    ex = f[f"{rg}/losUnitVectorX"][k]
    ny = f[f"{rg}/losUnitVectorY"][k]
    # fractional indices of the product grid in the (coarse) cube grid
    fx = np.interp(xg, cx, np.arange(cx.size)) if cx[0] < cx[-1] else \
        np.interp(xg, cx[::-1], np.arange(cx.size)[::-1])
    fy = np.interp(yg, cy, np.arange(cy.size)) if cy[0] < cy[-1] else \
        np.interp(yg, cy[::-1], np.arange(cy.size)[::-1])
    FX, FY = np.meshgrid(fx, fy)
    e = map_coordinates(ex, [FY, FX], order=1)
    n = map_coordinates(ny, [FY, FX], order=1)
    u = np.sqrt(np.clip(1 - e ** 2 - n ** 2, 0, 1))
    return e, n, u


def read_gunw(path, prof, dem_grid):
    with h5py.File(path, "r") as f:
        unw_path = find(f, "/unwrappedPhase", "frequencyA")[0]
        coh_path = unw_path.rsplit("/", 1)[0] + "/coherenceMagnitude"
        x, y, cgroup = coords_of(f, unw_path.rsplit("/", 1)[0])
        epsg = epsg_of(f, cgroup)
        sx, sy = aoi_slices(x, y, epsg)
        unw = f[unw_path][sy, sx].astype("float32")
        coh = f[coh_path][sy, sx].astype("float32") if coh_path in f else np.ones_like(unw)
        fpath = find(f, "frequencyA/centerFrequency")
        wl = C / float(f[fpath[0]][()]) if fpath else 0.2384
        xs, ys = x[sx], y[sy]
        e, n, u = los_vectors(f, xs, ys, dem_grid)
    unw[(unw == 0) | ~np.isfinite(unw)] = np.nan
    dx, dy = xs[1] - xs[0], ys[1] - ys[0]
    tr = from_origin(xs[0] - dx / 2, ys[0] - dy / 2, dx, -dy) if dy < 0 else None
    if tr is None:  # y increasing: flip to north-up
        unw, coh, e, n, u = (a[::-1] for a in (unw, coh, e, n, u))
        tr = from_origin(xs[0] - dx / 2, ys[-1] + dy / 2, dx, dy)
    crs = CRS.from_epsg(epsg)
    d = SIGN * unw * wl / (4 * np.pi) * 1000
    d[coh < COH_MIN] = np.nan
    g = lambda a: to_grid(a, tr, crs, prof)
    return {"los": g(d), "coh": g(coh), "e": g(e), "n": g(n), "u": g(u),
            "wavelength_m": wl, "epsg": epsg}


def main():
    prof = grid()
    z, ref_patches, slope = zones()
    dem, _ = read(DIRS["dem"] / "dem_utm30.tif")
    files = sorted(DIRS["nisar"].glob("NISAR_L2_*GUNW*.h5"))
    if not files:
        raise SystemExit("No GUNW .h5 in data/nisar: run 03 download.")
    DIRS["ts"].mkdir(parents=True, exist_ok=True)
    prods, meta = {}, {"sign": "LOS + toward satellite; vertical + up; east + east",
                       "sign_assumed": True, "products": {}}
    ref = ref_patches > 0
    flat = slope < SLOPE_MAX
    for p in files:
        info = parse_name(p.name)
        r = read_gunw(p, prof, dem)
        off = np.nanmean(r["los"][ref])
        r["los"] -= off
        key = f"{info['track']}{info['dir']}_{info['ref']}_{info['sec']}"
        write_stack(DIRS["ts"] / f"nisar_{key}_los.tif", [r["los"]], prof, [key])
        prods[p.name] = (info, r)
        meta["products"][key] = {**info, "wavelength_m": r["wavelength_m"],
                                 "epsg": r["epsg"], "ref_offset_mm": round(float(off), 2),
                                 "valid_px": int(np.isfinite(r["los"]).sum()),
                                 "mean_incidence_deg": float(np.degrees(np.arccos(np.nanmean(r["u"]))))}
        print(f"{key}: lambda {r['wavelength_m'] * 100:.2f} cm, valid "
              f"{np.isfinite(r['los']).mean() * 100:.0f}%, ref offset {off:.1f} mm")

    for win, (asc_tag, desc_tag) in WINDOWS.items():
        a = next((v for k, v in prods.items() if asc_tag in k and v[0]["dir"] == "A"), None)
        d = next((v for k, v in prods.items() if desc_tag in k and v[0]["dir"] == "D"), None)
        if not (a and d):
            print(f"{win}: asc/desc pair missing, no decomposition")
            continue
        (ia, ra), (id_, rd) = a, d
        det = ra["e"] * rd["u"] - rd["e"] * ra["u"]
        dE = (ra["los"] * rd["u"] - rd["los"] * ra["u"]) / det
        dU = (ra["e"] * rd["los"] - rd["e"] * ra["los"]) / det
        write_stack(DIRS["ts"] / f"nisar_{win}_vertical.tif", [dU], prof, [f"{win} up mm"])
        write_stack(DIRS["ts"] / f"nisar_{win}_east.tif", [dE], prof, [f"{win} east mm"])
        m = flat & np.isfinite(dU) & np.isfinite(dE)
        stats = {"asc": f"{ia['ref']}-{ia['sec']}", "desc": f"{id_['ref']}-{id_['sec']}",
                 "flat_px": int(m.sum()),
                 "flat_abs_east_p50_mm": float(np.nanmedian(np.abs(dE[m]))) if m.any() else None,
                 "flat_abs_up_p50_mm": float(np.nanmedian(np.abs(dU[m]))) if m.any() else None,
                 "flat_east_over_up_p50": float(np.nanmedian(np.abs(dE[m]) / (np.abs(dU[m]) + 1)))
                 if m.any() else None,
                 "per_zone": {}}
        from common import ZONE_NAMES
        for c, name in ZONE_NAMES.items():
            mz = m & (z == c)
            if mz.sum() > 10:
                stats["per_zone"][name] = {"up_mean_mm": round(float(np.nanmean(dU[mz])), 2),
                                           "east_mean_mm": round(float(np.nanmean(dE[mz])), 2),
                                           "n_px": int(mz.sum())}
        meta[win] = stats
        print(f"{win}: flat ground median |east| {stats['flat_abs_east_p50_mm']:.1f} mm vs "
              f"|up| {stats['flat_abs_up_p50_mm']:.1f} mm")
    save_json(DIRS["ts"] / "nisar_meta.json", meta)


if __name__ == "__main__":
    main()
