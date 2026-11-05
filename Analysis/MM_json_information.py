#A python script to quickly inspect material maps, you can 
import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Any, Dict, Iterable, Tuple, Set, List
from mpl_toolkits.mplot3d import Axes3D
#import pyg4ometry
import math
import itertools
import collections
import awkward as ak
import uproot
from collections import Counter
from collections import defaultdict


def load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    else:
        with path.open("r", encoding="utf-8") as f:
            return json.load(f)

def extract_entries(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    if isinstance(data.get("Surfaces"), dict) and "entries" in data["Surfaces"]:
        return data["Surfaces"]["entries"]

    raise KeyError(
        "Could not locate entries. Tried data['Surfaces']['entries'], "
        "data['surfaces']['entries'], and data['entries']."
    )
def _is_missing(x):
    return x is None or x == "" or x == [] or x == {}
    
def extract_entry_meta(entry):
    v   = entry.get("value") or {}
    mat = v.get("material") or {}
    bu  = mat.get("binUtility") or {}
    vol   = entry.get("volume") or {}
    lay   = entry.get("layer") or {}
    appr  = entry.get("approach") or {}
    bound = entry.get("boundary") or {}

    binning = bu.get("binningdata")

    ax0 = binning[0] if isinstance(binning, list) and len(binning) > 0 else None
    ax1 = binning[1] if isinstance(binning, list) and len(binning) > 1 else None

    if not _is_missing(vol) and not _is_missing(lay) and not _is_missing(appr):
        id_type  = "V/L/A"
        id_tuple = (int(vol), int(lay), int(appr))
    elif not _is_missing(vol) and not _is_missing(lay):
        id_type  = "V/L"
        id_tuple = (int(vol), int(lay))
    elif not _is_missing(vol) and not _is_missing(bound):
        try:
            id_tuple = (int(vol), int(bound))
        except Exception:
            id_tuple = (int(vol), bound)
        id_type = "V/B"

    return {
        "volume": entry.get("volume"),
        "layer": entry.get("layer"),
        "approach": entry.get("approach"),
        "boundary": entry.get("boundary"),
        "surface_type": entry.get("surface_type") or entry.get("type"),
        "binUtility" : bu,
        "binningData": binning,
        "axes0": ax0,
        "axes1": ax1,
        "id_tuple": id_tuple,
        "id_type" : id_type
    }


def axis_kind(ax):
    if not isinstance(ax, dict):
        return None
    tag = ax.get("value")

    if tag in ("binPhi", "AxisPhi"): return "phi"
    if tag in ("binR",   "AxisR"):   return "r"
    if tag in ("binZ",   "AxisZ"):   return "z"
    return None

def read_axis(ax):
    pmin = float(ax.get("min"))
    pmax = float(ax.get("max"))
    nb   = int(ax.get("bins", ax.get("nBins")))
    return pmin, pmax, nb


def get_material(mat_field):

    if mat_field is None:
        return(None, None, None, None, None)
    if isinstance(mat_field, str) and mat_field.strip().lower() == "Null":
        return 0.0, 0.0, 0.0, 0.0, 0.0
    if isinstance(mat_field, list) and mat_field:
        vals = []
        for idx in (0, 1, 2, 3, 4):
            vals.append(to_float_or_none(mat_field[idx]) if idx < len(mat_field) else None)
        return tuple(vals)
    return None, None, None, None, None


def to_float_or_none(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

        

def get_binning_info(mat_entries):
    out = []
    uid = 1
    
    for entry_idx, entry in enumerate(mat_entries):
        if not isinstance(entry, dict):
            continue

        meta = extract_entry_meta(entry)             
        bu = meta["binUtility"]
        id_tuple_mat = meta['id_tuple']
        id_type = meta["id_type"]
        ax0, ax1 = meta.get("axes0"), meta.get("axes1")
        k0, k1 = axis_kind(ax0), axis_kind(ax1)
        if k0 is None or k1 is None:
            continue
        kset = {k0, k1}
        if kset == {"phi", "z"}:
            surface_type = "Cylindrical Surface"

        elif kset == {"phi", "r"}: 
            surface_type = "Disc Surface"
            
        v   = (entry.get("value") or {})
        mat = (v.get("material") or {})
        grid = mat.get("data")
        if not isinstance(grid, list):
            continue
        amin, amax, anb = read_axis(ax0)
        bmin, bmax, bnb = read_axis(ax1)
        a_edges = np.linspace(amin, amax, anb + 1)
        b_edges = np.linspace(bmin, bmax, bnb + 1)
        n_a = len(a_edges) - 1
        n_b = len(b_edges) - 1
        for i_b in range(min(n_b, len(grid))):
            row = grid[i_b]
            if not isinstance(row, list):
                continue

            for i_a in range(min(n_a, len(row))):
                cell = row[i_a]

                if cell is None or (isinstance(cell, str) and cell.strip().lower() == "null"):
                    X0 = L0 = A = Z = rho = None
                    thickness = None
                elif not isinstance(cell, dict):
                    X0 = L0 = A = Z = rho = None
                    thickness = None
                else:
                    X0, L0, A, Z, rho = get_material(cell.get("material"))
                    thickness = to_float_or_none(cell.get("thickness"))

                out.append({
                    "id": uid, "entry_idx": entry_idx,
                    "surface_type": surface_type,
                    "id type" : id_type,
                    "id_tuple" : id_tuple_mat,
                    "X0 value" : X0, 
                    "L0 Value" : L0,
                    "A": A, 
                    "Z" : Z, 
                    "rho" : rho, 
                    "thickness": thickness
                })
                uid += 1
    return out

def surface_summary(json_info):
    counts = Counter(row.get("surface_type") for row in json_info)
    cylinder_counts = counts.get("Cylindrical Surface")
    disc_counts = counts.get("Disc Surface")
    total = len(json_info)
    return disc_counts, cylinder_counts, total

def volume_layer_summary(json_info):
    id_types = Counter(row.get("id type") for row in json_info)
    VLA_count = id_types.get("V/L/A")
    VL_count = id_types.get("V/L")
    VB_count = id_types.get("V/B")
    return VLA_count, VL_count, VB_count

def is_empty_thickness(v):
    if v is None:
        return True
    if isinstance(v, str):
        s = v.strip().lower()
        if s in ("none", "nan", ""):
            return True
        try:
            x = float(s)
        except ValueError:
            return True
        return (math.isnan(x) or x == 0.0)

    try:
        x = float(v)
    except (TypeError, ValueError):
        return True
    return (math.isnan(x) or x == 0.0)

def add_numeric(store, field, val):
    if val is None:
        return
    try:
        x = float(val)
    except (TypeError, ValueError):
        return
    if math.isnan(x):
        return
    store[field].append(x)

def split_id(id_key):
    vol = id_key[0] if len(id_key) > 0 else None
    lay = id_key[1] if len(id_key) > 1 else None
    app = id_key[2] if len(id_key) > 2 else None
    return vol, lay, app

def write_txt_file(entries, output):
    disc_count, cylinder_count, total_count = surface_summary(entries)
    VLA_count, VL_count, VB_count = volume_layer_summary(entries)
    material_by_id = defaultdict(lambda: {
        "thickness": [],
        "X0": [],
        "L0": [],
        "A": [],
        "Z": [],
        "rho": [],
        "n_cells_total": 0,
        "n_cells_empty": 0,
    })
    
    for e in entries:
        id_tuple = e.get('id_tuple')
        thick = e.get('thickness')
        x0_value = e.get('X0 value')
        l0_value = e.get('L0 Value')
        A_value = e.get('A')
        Z_value = e.get('Z')
        rho_value = e.get('rho') 

        id_key = tuple(id_tuple)

        d = material_by_id[id_key]
        d["n_cells_total"] += 1
        if is_empty_thickness(thick):
            d["n_cells_empty"] += 1
        
        add_numeric(d, "thickness", thick)
        add_numeric(d, "X0", x0_value)
        add_numeric(d, "L0", l0_value)
        add_numeric(d, "A", A_value)
        add_numeric(d, "Z", Z_value)
        add_numeric(d, "rho", rho_value)

    outfile = output

    def mean_std(vals):
        if len(vals) == 0:
            return 0, np.nan, np.nan
        arr = np.array(vals, dtype=float)
        return float(np.mean(arr)), float(np.std(arr))
    
    with open(outfile, "w") as f:
        f.write("Material Map Summary:\n")
        f.write(f"Total Surfaces: {total_count}\t Disc Surfaces: {disc_count}\t Cylinder Surfaces:{cylinder_count}\n")
        f.write(f"Formatting- Surfaces with: V/L/A:{VLA_count}\t V/L:{VL_count}\t V/B:{VB_count}\n")
        f.write(
            f"{'id_tuple':<10}"
            f"{'n_bins':^8}"
            f"{'thick':^10} {'std':^10} "
            f"{'X0':^10} {'std_X0':^10} "
            f"{'L0':^10} {'std_L0':^10} "
            f"{'A':^10} {'std_A':^10} "
            f"{'Z':^10}" #{'std_Z':^10} "
            f"{'rho':^10} {'std_rho':^10} "
            f"{'empty':^8}\n "
        )
    
        for id_key in sorted(material_by_id.keys()):
            d = material_by_id[id_key]
            n_cells = d["n_cells_total"]
            n_empty = d["n_cells_empty"]
            empty_frac = (n_empty / n_cells) if n_cells > 0 else np.nan
    
            mu_t, sd_t = mean_std(d["thickness"])
            mu_x, sd_x = mean_std(d["X0"])
            mu_l, sd_l = mean_std(d["L0"])
            mu_a, sd_a = mean_std(d["A"])
            mu_z, sd_z = mean_std(d["Z"])
            mu_r, sd_r = mean_std(d["rho"])
            id_str = str(id_key)
    
            f.write(
                f"{str(id_key):<10} "
                f"{n_cells:^8} "
                f"{mu_t:^10.2f} {sd_t:^10.2f} "
                f"{mu_x:^10.2e} {sd_x:^10.2e} "
                f"{mu_l:^10.2e} {sd_l:^10.2e} "
                f"{mu_a:^10.2f} {sd_a:^10.2f} "
                f"{mu_z:^10.2f}" # {sd_z:>10.2e} "
                f"{mu_r:^10.2e} {sd_r:^10.2e} "
                f"{empty_frac:^10.2f}\n"
        )
    
    print(f"Wrote stats for {len(material_by_id)} id_tuples to {outfile}")
        

def main():
    parser = argparse.ArgumentParser(
        description="Extract and summarize 'entries' from an ACTS/ITk-style JSON."
    )
    parser.add_argument(
        "-m", "--material_map",
        required=True,
        help="Path to material map JSON file"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to output txt"
    )
    args = parser.parse_args()
    mm_json = Path(args.material_map)
    out_txt = Path(args.output)

    mm_load = load_json(mm_json)
    mm_entry = extract_entries(mm_load)
    mm_binning = get_binning_info(mm_entry)
    write_txt_file(mm_binning, out_txt)

if __name__ == "__main__":
    main()



