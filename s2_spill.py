"""Presentation asset: Sentinel-2 L2A true colour of the 2020 TPP-3 diesel spill.

Planetary Computer (anonymous). Mosaics the two MGRS tiles of each date over
TPP-3 -> Daldykan -> Ambarnaya. Credit: contains modified Copernicus
Sentinel data 2020.

    python s2_spill.py   -> outputs/web/s2_spill_<date>.png (+ titled copy)
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import planetary_computer as pc
import pystac_client
import rasterio
from pyproj import Transformer
from rasterio.windows import from_bounds

BOX = (87.75, 69.30, 88.10, 69.55)   # lon/lat
SCENES = {"20200531": ["S2B_MSIL2A_20200531T060639_R134_T45WWS_20201027T072154",
                       "S2B_MSIL2A_20200531T060639_R134_T45WWT_20200919T124641"],
          "20200608": ["S2A_MSIL2A_20200608T061641_R034_T45WWS_20200919T124640",
                       "S2A_MSIL2A_20200608T061641_R034_T45WWT_20200919T124645"]}
TPP3 = (69.3227, 87.9499)
OUT = Path(__file__).parent / "outputs"


def main():
    (OUT / "web").mkdir(parents=True, exist_ok=True)
    col = pystac_client.Client.open(
        "https://planetarycomputer.microsoft.com/api/stac/v1",
        modifier=pc.sign_inplace).get_collection("sentinel-2-l2a")
    for d, ids in SCENES.items():
        out = None
        for iid in ids:
            with rasterio.open(col.get_item(iid).assets["visual"].href) as s:
                t = Transformer.from_crs(4326, s.crs, always_xy=True)
                xs, ys = t.transform([BOX[0], BOX[2], BOX[0], BOX[2]],
                                     [BOX[1], BOX[1], BOX[3], BOX[3]])
                w = from_bounds(min(xs), min(ys), max(xs), max(ys), s.transform)
                a = s.read(window=w, boundless=True, fill_value=0)
                tr = s.window_transform(w)
            out = a if out is None else np.where(out.sum(0) > 0, out, a)
        rgb = np.clip(np.moveaxis(out, 0, -1) / 255 * 1.3, 0, 1)
        plt.imsave(OUT / "web" / f"s2_spill_{d}.png", rgb)
        c, r = ~tr * t.transform(TPP3[1], TPP3[0])
        fig, ax = plt.subplots(figsize=(9, 10))
        ax.imshow(rgb)
        ax.plot(c, r, "o", mfc="none", mec="magenta", ms=14, mew=2)
        ax.annotate("TPP-3", (c, r), xytext=(10, -10), textcoords="offset points",
                    color="magenta", fontsize=10, weight="bold")
        ax.axis("off")
        ax.set_title(f"Sentinel-2 L2A true colour {d[:4]}-{d[4:6]}-{d[6:]}, "
                     "TPP-3 -> Daldykan -> Ambarnaya\n"
                     "Contains modified Copernicus Sentinel data 2020", fontsize=9)
        fig.tight_layout()
        fig.savefig(OUT / f"s2_spill_{d}_titled.png", dpi=130)
        plt.close(fig)
        print(d, rgb.shape, f"valid {(out.sum(0) > 0).mean():.2f}")


if __name__ == "__main__":
    main()
