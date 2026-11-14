"""Shared settings: study area, points of interest, periods, output folders."""
from pathlib import Path

# ---------------------------------------------------------------- study area
# Bounding box (lon_min, lat_min, lon_max, lat_max), WGS84.
# Covers Norilsk city, Talnakh (N), Kayerkan + HPP-3 (W), Nadezhda plant,
# surrounding tundra and the edge of the plateau (bedrock reference).
BBOX = (87.55, 69.20, 88.65, 69.58)

# Metric CRS for derived rasters (slope etc.): UTM 45N covers 84-90 E.
UTM_EPSG = 32645

# Points of interest (lat, lon). APPROXIMATE: check every point on
# Google Earth / OSM before using it in analysis.
POINTS = {
    "norilsk_centre": (69.349, 88.201),
    # Sources disagree on the HPP-3 tank 5 location; both kept until checked.
    "hpp3_tank_a": (69.3278, 87.9644),
    "hpp3_tank_b": (69.3794, 87.7444),
    "talnakh": (69.49, 88.39),
    # Placeholders: pick real control sites from DEM + land cover (step 05).
    "tundra_control": (69.42, 87.70),
    "plateau_reference": (69.30, 88.35),
}

# ------------------------------------------------------------------ periods
YEARS = list(range(2017, 2027))       # summers 2017-2026
# Snow-free InSAR window (MM-DD). Snow decorrelates C-band; keep summer only.
SUMMER = ("06-10", "10-05")
NISAR_START = "2026-06-17"            # first PROVISIONAL (calibrated) products

# ------------------------------------------------------------------ outputs
DATA = Path(__file__).parent / "data"
DIRS = {k: DATA / k for k in
        ["s1", "hyp3", "nisar", "climate", "dem", "landcover", "osm", "cci"]}
for d in DIRS.values():
    d.mkdir(parents=True, exist_ok=True)


def bbox_wkt(b=BBOX) -> str:
    x0, y0, x1, y1 = b
    return f"POLYGON(({x0} {y0},{x1} {y0},{x1} {y1},{x0} {y1},{x0} {y0}))"
