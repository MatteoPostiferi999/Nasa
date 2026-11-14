"""Step 11: zone x date table and the sqrt(TDD) test, per zone and per summer.

For each summer stack from step 09 (vertical, mm, relative to first date and
referenced to the plateau patches) and each zone:
  disp(t) = mean vertical displacement of the zone's pixels at date t
  fit     disp(t) = a + b * sqrt(TDD_cum(t))      (Stefan-type thaw model)
TDD_cum is the cumulative thawing degree days (ERA5-Land, from 1 Jan) at the
zone's climate point (config.ZONE_CLIMATE). b < 0 means subsidence that
grows with thaw; b is in mm per sqrt(degC day).

Uncertainty of b: OLS standard error (dates are few, so this is a lower
bound: spatially correlated atmospheric noise is not in it). The per-date
zone uncertainty is the spatial std / sqrt(n_eff), n_eff = n_px / PX_PER_IND
(neighbouring pixels are not independent).

Comparisons: industrial vs each control, control_NW vs control_SE,
difference of b with combined standard error (z-score).

Gap model (PREDICTION, flagged): per zone, b averaged over summers,
predicted end-of-summer subsidence for 2022-2024 from ERA5 TDD;
leave-one-summer-out check of b.

Outputs: outputs/web/zone_timeseries.csv, outputs/zone_fits.csv,
outputs/zone_comparisons.csv, outputs/gap_prediction.csv,
outputs/fig_disp_vs_sqrt_tdd.png, outputs/fig_subsidence_<year>.png.

    python 11_analysis.py
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio

from common import ZONE_NAMES, zones
from config import DIRS, OUTPUTS, SUMMER, ZONE_CLIMATE

PX_PER_IND = 25          # ~150 m correlation length on a 30 m grid
GAP_YEARS = [2022, 2023, 2024]
COLORS = {"reference": "black", "industrial": "tab:orange", "accident": "tab:red",
          "urban": "gold", "control_NW": "tab:green", "control_SE": "tab:cyan"}


def tdd_series(point):
    df = pd.read_csv(DIRS["climate"] / f"era5land_daily_{point}.csv",
                     index_col="date", parse_dates=["date"])
    return df["tdd_cum"]


def fit(x, y):
    """OLS y = a + b x. Returns b, se_b, a, r2, n."""
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = len(x)
    if n < 3 or np.ptp(x) == 0:
        return dict(b=np.nan, se_b=np.nan, a=np.nan, r2=np.nan, n=n)
    X = np.c_[np.ones(n), x]
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ coef
    s2 = res @ res / (n - 2)
    cov = s2 * np.linalg.inv(X.T @ X)
    r2 = 1 - (res @ res) / np.sum((y - y.mean()) ** 2) if np.ptp(y) > 0 else np.nan
    return dict(b=coef[1], se_b=np.sqrt(cov[1, 1]), a=coef[0], r2=r2, n=n)


def zone_table(z):
    rows = []
    for f in sorted(DIRS["ts"].glob("s1_*_vert.tif")):
        year = int(f.stem.split("_")[1])
        with rasterio.open(f) as s:
            stack = s.read().astype("float32")
            dates = [s.descriptions[i] for i in range(s.count)]
        for c, name in ZONE_NAMES.items():
            m = z == c
            tdd = tdd_series(ZONE_CLIMATE[name])
            for i, d in enumerate(dates):
                v = stack[i][m]
                v = v[np.isfinite(v)]
                if v.size < 10:
                    continue
                t = pd.Timestamp(d)
                rows.append({"zone": name, "summer": year, "date": t.date(),
                             "disp_mm": float(np.mean(v)), "disp_median_mm": float(np.median(v)),
                             "se_mm": float(np.std(v) / np.sqrt(max(v.size / PX_PER_IND, 1))),
                             "n_px": int(v.size), "tdd_cum": float(tdd.get(t, np.nan))})
    return pd.DataFrame(rows)


def main():
    z, _, _ = zones()
    (OUTPUTS / "web").mkdir(parents=True, exist_ok=True)
    tab = zone_table(z)
    if tab.empty:
        raise SystemExit("No s1_*_vert.tif in data/ts: run 09 first.")
    tab["sqrt_tdd"] = np.sqrt(tab["tdd_cum"])
    tab[["zone", "summer", "date", "disp_mm", "disp_median_mm", "se_mm", "n_px",
         "tdd_cum"]].to_csv(OUTPUTS / "web" / "zone_timeseries.csv", index=False)

    fits = []
    for (zone, yr), g in tab.groupby(["zone", "summer"]):
        r = fit(g["sqrt_tdd"].values, g["disp_mm"].values)
        end = g.sort_values("date").iloc[-1]
        fits.append({"zone": zone, "summer": yr, **r,
                     "seasonal_disp_mm": float(end["disp_mm"]),
                     "tdd_end": float(end["tdd_cum"])})
    fits = pd.DataFrame(fits)
    fits.round(4).to_csv(OUTPUTS / "zone_fits.csv", index=False)

    comps = []
    for yr, g in fits.groupby("summer"):
        g = g.set_index("zone")
        for a, b in [("industrial", "control_NW"), ("industrial", "control_SE"),
                     ("accident", "control_NW"), ("accident", "control_SE"),
                     ("urban", "control_NW"), ("urban", "control_SE"),
                     ("control_NW", "control_SE")]:
            if a in g.index and b in g.index:
                d = g.at[a, "b"] - g.at[b, "b"]
                se = np.hypot(g.at[a, "se_b"], g.at[b, "se_b"])
                comps.append({"summer": yr, "zone_a": a, "zone_b": b, "diff_b": d,
                              "se": se, "z": d / se if se > 0 else np.nan})
    comps = pd.DataFrame(comps)
    comps.round(4).to_csv(OUTPUTS / "zone_comparisons.csv", index=False)

    # gap model (prediction): b averaged over summers; LOSO check
    gap = []
    for zone, g in fits.groupby("zone"):
        b_mean = g["b"].mean()
        tdd = tdd_series(ZONE_CLIMATE[zone])
        for y in GAP_YEARS:
            s0, s1 = (pd.Timestamp(f"{y}-{SUMMER[0]}"), pd.Timestamp(f"{y}-{SUMMER[1]}"))
            pred = b_mean * (np.sqrt(tdd.get(s1, np.nan)) - np.sqrt(tdd.get(s0, np.nan)))
            gap.append({"zone": zone, "year": y, "predicted_seasonal_mm": pred,
                        "b_used": b_mean, "type": "PREDICTION (not measured)"})
        for _, row in g.iterrows():
            others = g[g.summer != row.summer]["b"]
            if len(others):
                gap.append({"zone": zone, "year": int(row.summer),
                            "loso_b_pred": others.mean(), "b_obs": row.b,
                            "type": "leave-one-summer-out check"})
    pd.DataFrame(gap).round(3).to_csv(OUTPUTS / "gap_prediction.csv", index=False)

    # figure: disp vs sqrt(TDD)
    years = sorted(tab.summer.unique())
    fig, axs = plt.subplots(1, len(years), figsize=(6 * len(years), 4.5), squeeze=False)
    for ax, yr in zip(axs[0], years):
        for zone, g in tab[tab.summer == yr].groupby("zone"):
            ax.errorbar(g.sqrt_tdd, g.disp_mm, yerr=g.se_mm, fmt="o", ms=4,
                        color=COLORS.get(zone), label=zone, capsize=2)
            f = fits[(fits.zone == zone) & (fits.summer == yr)].iloc[0]
            if np.isfinite(f.b):
                xx = np.linspace(g.sqrt_tdd.min(), g.sqrt_tdd.max(), 10)
                ax.plot(xx, f.a + f.b * xx, "-", color=COLORS.get(zone), lw=1)
        ax.axhline(0, color="grey", lw=0.5)
        ax.set_title(f"Summer {yr}: S1 vertical vs sqrt(TDD)")
        ax.set_xlabel("sqrt(cumulative TDD) [sqrt(degC day)]")
        ax.set_ylabel("vertical displacement since first date [mm] (+ up)")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUTPUTS / "fig_disp_vs_sqrt_tdd.png", dpi=130)
    plt.close(fig)

    # figure: seasonal subsidence map per summer (last date)
    for f in sorted(DIRS["ts"].glob("s1_*_vert.tif")):
        with rasterio.open(f) as s:
            last = s.read(s.count)
            d0, d1 = s.descriptions[0], s.descriptions[-1]
        fig, ax = plt.subplots(figsize=(10, 7.5))
        lim = np.nanpercentile(np.abs(last), 98) if np.isfinite(last).any() else 1
        im = ax.imshow(last, cmap="RdBu", vmin=-lim, vmax=lim)
        ax.contour(z == 2, levels=[0.5], colors="orange", linewidths=0.5)
        ax.contour(z == 3, levels=[0.5], colors="magenta", linewidths=0.8)
        ax.contour(z == 1, levels=[0.5], colors="black", linewidths=0.8)
        ax.set_title(f"S1 vertical displacement {d0} -> {d1} [mm], + up "
                     "(slope < 5 deg only). Orange industrial, magenta TPP-3, black reference",
                     fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=ax, shrink=0.7)
        fig.tight_layout()
        fig.savefig(OUTPUTS / f"fig_subsidence_{f.stem.split('_')[1]}.png", dpi=130)
        plt.close(fig)

    pd.set_option("display.width", 140)
    print(fits[["zone", "summer", "b", "se_b", "r2", "n", "seasonal_disp_mm"]].round(3).to_string(index=False))
    print()
    print(comps.round(3).to_string(index=False))
    (OUTPUTS / "analysis_summary.json").write_text(json.dumps(
        {"fits": fits.round(4).to_dict("records"),
         "comparisons": comps.round(4).to_dict("records")}, default=str, indent=1))


if __name__ == "__main__":
    main()
