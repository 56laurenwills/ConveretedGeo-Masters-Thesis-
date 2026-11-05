import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Tuple
import collections
from collections import Counter


def load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    else:
        with path.open("r", encoding="utf-8") as f:
            return json.load(f)

def extract_entries(data: Dict[str, Any]) -> List[Dict[str, Any]]:
#Format of json: json['Surfaces']['entries'], want the list of entries
    if isinstance(data.get("Surfaces"), dict) and "entries" in data["Surfaces"]:
        return data["Surfaces"]["entries"]

    raise KeyError(
        "Could not locate entries. Tried data['Surfaces']['entries'], "
        "data['surfaces']['entries'], and data['entries']."
    )

def is_list_of_5_numbers(x):
    return (
        isinstance(x, list)
        and len(x) == 5
        and all(isinstance(v, (int, float)) for v in x)
    )

def axis_kind(ax):
    """Map axis dict -> 'phi' | 'r' | 'z' (or None if unknown / missing)."""
    if not isinstance(ax, dict):
        return None
    tag = ax.get("value")

    if tag in ("binPhi", "AxisPhi"): return "phi"
    if tag in ("binR",   "AxisR"): return "r"
    if tag in ("binZ",   "AxisZ"): return "z"
    return None

def _is_missing(x):
    return x is None or x == "" or x == [] or x == {}
    
def extract_entry_meta(entry):
#extract the information from each entry
    v = entry.get("value")    or {}
    mat = v.get("material")     or {}
    bu  = mat.get("binUtility") or {}
    vol = entry.get("volume")
    lay = entry.get("layer")
    appr = entry.get("approach")
    bound = entry.get("boundary")
    binning = bu.get("binningdata")

    ax0 = binning[0] if isinstance(binning, list) and len(binning) > 0 else None
    ax1 = binning[1] if isinstance(binning, list) and len(binning) > 1 else None
    
    if not _is_missing(vol) and not _is_missing(lay) and not _is_missing(appr):
        id_type = "V/L/A"
        id_tuple = (int(vol), int(lay), int(appr))
    elif not _is_missing(vol) and not _is_missing(lay):
        id_type = "V/L"
        id_tuple = (int(vol), int(lay))
    elif not _is_missing(vol) and not _is_missing(bound):
        try:
            id_tuple = (int(vol), int(bound))
        except Exception:
            id_tuple = (int(vol), bound)
        id_type = "V/B"

    return {
        "volume": vol,
        "layer": lay, 
        "approach": appr,
        "boundary": bound,
        "binUtility" : bu,
        "binningdata": binning,
        "axes0": ax0,
        "axes1": ax1,
        "id_tuple": id_tuple,
        "id_type" : id_type
    }
def cell00_material(entry):
    v = entry.get("value") or {}
    m = v.get("material")
    if isinstance(m, list):
        m = next((mm for mm in m if isinstance(mm, dict)), None)
    if not isinstance(m, dict):
        return None

    try:
        cell00 = m["data"][0][0]
    except (KeyError, IndexError, TypeError):
        return None

    if isinstance(cell00, dict):
        return cell00.get("material")
    if isinstance(cell00, list):
        if cell00 and isinstance(cell00[0], dict):
            return cell00[0].get("material")
        return cell00
    return None

def read_axis(ax):
    pmin = float(ax.get("min"))
    pmax = float(ax.get("max"))
    nb = int(ax.get("bins", ax.get("nBins")))
    return pmin, pmax, nb

def extract_meta_data(entries):
    out = []
    uid = 1
    theoretical_total = 0 
    for entry_idx, entry in enumerate(entries):
        meta = extract_entry_meta(entry)
        vol= meta["volume"]
        lay = meta["layer"]
        appr = meta["approach"]
        bound = meta["boundary"]
        bu = meta["binUtility"]
        binning = meta['binningdata']
        ax0 = meta['axes0']
        ax1 = meta['axes1']
        k0, k1 = axis_kind(ax0), axis_kind(ax1)     
        kset = {k0, k1}
        id_tuple = meta['id_tuple']
        id_type = meta["id_type"]
        v = entry.get('value')
        mat = v.get('material')
        data = mat.get('data')
        mt = mat.get('mappingType')
        mm = mat.get('mapMaterial')
        t = mat.get('type')
        cell_mat = cell00_material(entry)
        has_cell00_material  = cell_mat is not None
        material_list_5 = is_list_of_5_numbers(cell_mat)
        material_list_len = len(cell_mat) if isinstance(cell_mat, list) else None

        phi_ax   = ax0 if k0 == "phi" else ax1
        other_ax = ax0 if k0 in ("z", "r") else ax1

        _, _, n_phi   = read_axis(phi_ax)
        _, _, n_other = read_axis(other_ax)
        theoretical_total += int(n_phi) * int(n_other)


        if kset == {"phi", "z"}:
            surface_type = "Cylindrical Surface"
        elif kset == {"phi", "r"}:
            surface_type = "Disc Surface"
        out.append({
            "id" : uid,
            "entry idx" : entry_idx,
            "id type" : id_type,
            "id_tuple" : id_tuple,
            "mapping type" : mt,
            "Does it have material?" : mm,
            "Type of material info" : t,
            "Surface Type" : surface_type
        })
        uid += 1
    return out, theoretical_total

def surface_summary(json_info):
    counts = Counter(row.get("Surface Type") for row in json_info)
    cylinder_counts = counts.get("Cylindrical Surface")
    disc_counts = counts.get("Disc Surface")
    total = len(json_info)
    return disc_counts, cylinder_counts, total

def volume_layer_summary(json_data):
    id_types = Counter(row.get("id type") for row in json_data)
    VLA_count = id_types.get("V/L/A")
    VL_count = id_types.get("V/L")
    VB_count = id_types.get("V/B")
    return VLA_count, VL_count, VB_count

def material_type_summary(json_data):
    mt = Counter(row.get("mapMaterial") for row in json_data)
    mt_false = mt.get(False)
    mt_true = mt.get(True)
    return mt_false, mt_true

def empty_material(json_data):
    material_empty = 0
    for row in json_data:
        if row.get('material_list_5') == False:
            material_empty += 1
    return material_empty
    

def main():
    parser = argparse.ArgumentParser(
        description="Extract and summarize 'entries' from an ACTS/ITk-style JSON."
    )
    parser.add_argument(
        "-f", "--file",
        required=True,
        help="Path to the JSON file"
    )
    args = parser.parse_args()
    json_path = Path(args.file)
    data = load_json(json_path)
    entries = extract_entries(data)
    json_data, total_bins = extract_meta_data(entries)
    disc, cylinder, total = surface_summary(json_data)
    VLA, VL, VB  = volume_layer_summary(json_data)
    mat_type_t, mat_type_f = material_type_summary(json_data)
    num_empty_mat = empty_material(json_data)
    print(f" Total Number of Surfaces: {total}, Type of Surface: N Cylindrical: {cylinder}, Number of Disc: {disc}")
    print(f" Number of Entries with V/L/A: {VLA}, with V/L: {VL}, with V/B: {VB}")
    print(f" Number of Entries with mapMaterial = False : {mat_type_f} vs. Entries with Null Material: {num_empty_mat}")
    print(f" Total number of bins: {total_bins}")


if __name__ == "__main__":
    main()
    