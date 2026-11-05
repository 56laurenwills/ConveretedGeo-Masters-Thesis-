import json
from pathlib import Path

def _is_missing(x):
    return x is None or x == "" or x == [] or x == {}
    
def get_id_tuple(entry):
    v = entry.get("value") or {}
    vol = entry.get("volume") or {}
    lay = entry.get("layer") or {}
    appr = entry.get("approach") or {}
    bound = entry.get("boundary") or {}

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
        "volume": entry.get("volume"),
        "layer": entry.get("layer"),
        "approach": entry.get("approach"),
        "boundary": entry.get("boundary"),
        "id_tuple": id_tuple,
        #"id_type" : id_type
    }

def matches(pattern, eid):

    if eid is None or len(pattern) != len(eid):
        return False
    return all(p is None or p == v for p, v in zip(pattern, eid))

def should_remove(eid, remove_list):
    for pat in remove_list:
        if pat == eid:
            return True
        if isinstance(pat, tuple) and isinstance(eid, tuple) and matches(pat, eid):
            return True
    return False

def remove_entries(in_mm, out_mm, remove_list):
    in_json = Path(in_mm)
    out_json = Path(out_mm)

    with in_json.open("r", encoding="utf-8") as f:
            data = json.load(f)

    surfaces = data['Surfaces']
    entries = surfaces['entries']

    keep = []
    remove = []
    for e in entries:
        meta = get_id_tuple(e)
        id_tuple = meta["id_tuple"]
        if should_remove(id_tuple, remove_list):
            remove.append(id_tuple)
        else:
            keep.append(e)

    data["Surfaces"]["entries"] = keep
    with out_json.open("w") as f:
        json.dump(data, f, indent=2)

    # Summary
    print(f"Input entries: {len(entries)}")
    print(f"Kept entries: {len(keep)}")
    print(f"Removed entries: {len(remove)}")
    if remove:
        show = remove[:20]
        print("Example removed IDs:", show)
        
if __name__ == "__main__":
    remove_itk = [
        (2, None, 1),
        (8, None, 1),
        (9, None, 2),
        (10, None, 2),
        (13, None, 1),
        (14, None, 1),
        (15, None, 1),
        (16, None, 2),
        (18, None, 2),
        (19, None, 2),
        (20, None, 2),
        (22, None, 1),
        (23, None, 2),
        (24, None, 2),
        (25, None, 2)
    ]

    remove_entries(
        in_mm = "/path/to/material-map.json",
        out_mm = "filtered_material-map.json",
        remove_list = remove_ids_itk
    )


        
