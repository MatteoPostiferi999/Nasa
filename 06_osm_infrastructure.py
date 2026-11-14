"""Step 06: OpenStreetMap infrastructure in the AOI -> GeoJSON (Overpass API).

Buildings, industrial land, storage tanks, power plants, pipelines,
tailings. Used to separate built-up vs natural ground and to locate
heat sources (anthropogenic thaw hypothesis).

    python 06_osm_infrastructure.py
"""
import json

import requests

from config import BBOX, DIRS

URL = "https://overpass-api.de/api/interpreter"
s, w, n, e = BBOX[1], BBOX[0], BBOX[3], BBOX[2]
QUERY = f"""
[out:json][timeout:300];
(
  way["building"]({s},{w},{n},{e});
  way["landuse"~"industrial|construction|landfill"]({s},{w},{n},{e});
  way["man_made"~"storage_tank|pipeline|tailings_pond|wastewater_plant"]({s},{w},{n},{e});
  node["man_made"="storage_tank"]({s},{w},{n},{e});
  way["power"~"plant|generator|substation"]({s},{w},{n},{e});
  way["landuse"="quarry"]({s},{w},{n},{e});
);
out geom tags;
"""


def to_feature(el):
    tags = el.get("tags", {})
    props = {"osm_id": f"{el['type']}/{el['id']}", **tags}
    if el["type"] == "node":
        geom = {"type": "Point", "coordinates": [el["lon"], el["lat"]]}
    else:
        coords = [[p["lon"], p["lat"]] for p in el.get("geometry", [])]
        if len(coords) < 2:
            return None
        closed = coords[0] == coords[-1] and len(coords) >= 4
        linear = tags.get("man_made") == "pipeline"
        geom = ({"type": "Polygon", "coordinates": [coords]}
                if closed and not linear
                else {"type": "LineString", "coordinates": coords})
    return {"type": "Feature", "geometry": geom, "properties": props}


def layer(props):
    if props.get("man_made") == "storage_tank":
        return "tanks"
    if props.get("man_made") == "pipeline":
        return "pipelines"
    if "building" in props:
        return "buildings"
    return "industrial"


def main():
    r = requests.post(URL, data={"data": QUERY}, timeout=600)
    r.raise_for_status()
    feats = [f for f in map(to_feature, r.json()["elements"]) if f]
    groups = {}
    for f in feats:
        groups.setdefault(layer(f["properties"]), []).append(f)
    for name, fs in groups.items():
        path = DIRS["osm"] / f"osm_{name}.geojson"
        path.write_text(json.dumps({"type": "FeatureCollection", "features": fs}))
        print(f"{name}: {len(fs)} features -> {path.name}")


if __name__ == "__main__":
    main()
