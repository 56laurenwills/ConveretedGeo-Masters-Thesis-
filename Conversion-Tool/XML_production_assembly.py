from __future__ import annotations
 
import argparse
import math
import os
import warnings
from collections import OrderedDict
from typing import Dict, List, Optional, Tuple
from xml.etree import ElementTree as ET
from xml.dom import minidom
from collections import defaultdict
 
import numpy as np
import re

VolumeEntry = Tuple[object, object, np.ndarray, np.ndarray]

#COMMENT: Used Claude to help with the actual write out code, the formatting for gmx and such
#current format of this writer is to be used in conjunction with official ATLAS ITk gmx's, If you want to use it stand alone, need to add dtd path N2 back in

def parse_volume_id(name):
    m = re.search(r'vol(\d+)_bou(\d+)_lay(\d+)_app(\d+)', name)
    if m is None:
        raise ValueError(f"Could not parse volume id from {name!r}")
    vol, bou, lay, app = (int(g) for g in m.groups())
    return vol, bou, lay, app

def target_volume_for(vnum, bou, lay):
    return VOL_REDIRECT.get((vnum, bou, lay), vnum)

VOL_CONFIG = {
    2: ("OutermostBoundary_neg", 0, 1050, 50),
    3: ("OutermostBoundary_pos", 0, 1050, 50),
    4:  ("IS_EC_neg_converted", 26, 150, 2400),
    5:  ("IS_barrel_converted", 0, 115, 255),
    6:  ("IS_EC_pos_converted", 26, 150, 2400),
    9:  ("OEC_neg_converted", 130, 355, 2550),
    10: ("OP_barrel_converted", 160, 305, 1000),
    11: ("OEC_pos_converted", 130,  355, 2550),
    19: ("Strip_EC_neg_converted", 360, 1005, 2870),
    20: ("Strip_barrel_converted", 410, 790, 2330),
    21: ("Strip_EC_pos_converted", 360, 1005, 2870),
    94: ("OEC_neg_inner_converted", 0, 150, 3000),
    95: ("OEC_neg_mid_converted", 0, 150, 3000),
    96: ("OEC_neg_outer_converted", 0, 150, 30000),
    97: ("OEC_pos_inner_converted", 0, 150, 3000),
    98: ("OEC_pos_mid_converted", 0, 150, 3000),
    99: ("OEC_pos_outer_converted", 0, 150, 3000)
}

VOL_REDIRECT = { #these values are subject to change depending on which verion of material map
    #(1, 1, 1): 5,
    #(1, 2, 1): 6,
    #(1, 3, 1): 5,
    #(2, 2, 0): 19,
    #(3, 1, 0): 2,
    (7, 2, 0): 20, 
    (20, 1, 1): 19,
    (20, 2, 1): 21,
}

def write_split(all_volumes, output_dir, dtd_path, solid_caches=None, mode="assembly", reuse=True):
        import os
        os.makedirs(output_dir, exist_ok=True)
        buckets = defaultdict(list)
        skipped =[]
        temp_writer = GMXWriter(name="temp")
        for vol in all_volumes:
            if vol.material is not None:
                temp_writer.materials[vol.material.name] = vol.material

        mat_path = os.path.join(output_dir, "simplified_materials.gmx")
        temp_writer.write_material(mat_path)

        for vol in all_volumes:
            target_volume = getattr(vol ,"target_volume", None)
            if target_volume is None:
                try:
                    vnum, bou, lay, appr = parse_volume_id(vol.name)
                except ValueError:
                    skipped.append(vol.name)
                    continue
                target_volume = target_volume_for(vnum, bou, lay)
            if target_volume not in VOL_CONFIG:
                skipped.append(vol.name)
                continue
            buckets[target_volume].append(vol)
        if skipped:
            print(f"[GMX] Skipped AHHHH")
        for target_volume, vols in sorted(buckets.items()):
            stem, rmin, rmax, zhalflength = VOL_CONFIG[target_volume]
            out_path = os.path.join(output_dir, stem + ".gmx")
            print(f"[GMX] writing {stem}.gmx")
            writer = GMXWriter(name=stem)
            cache = solid_caches[target_volume] if (reuse and solid_caches is not None) else None
            writer.collect(vols, solid_cache=cache)
            use_per_layer = (target_volume == 20)
            writer.write_main(out_path, per_layer=use_per_layer, mode=mode, rmin=rmin, rmax=rmax, zhalflength=zhalflength)
            #print(f"[GMX] ahhhh")
            
def _eval(val):
    try:
        return float(val.eval())
    except AttributeError:
        return float(val)
    
def vol_number(vol_name):
    m = re.search(r'vol(\d+)', vol_name)
    if m is None:
        print ("ahhhh")
    return int(m.group(1))

class GMXWriter:
    def __init__(self, name):
        self.name = name
        self.materials: OrderedDict[str, object] = OrderedDict()
        self.shapes: OrderedDict[str, ET.Element] = OrderedDict()
        self.logvols: Dict[str, Tuple[str, str]] = {}

    @staticmethod
    def make_tubs(name, rmin, rmax, zhalflength, sphi, dphi):
        tubs = ET.Element("tubs")
        tubs.set("name", name)
        tubs.set("rmin", f"{rmin:.8g}")
        tubs.set("rmax", f"{rmax:.8g}")
        tubs.set("zhalflength", f"{zhalflength:.8g}")
        tubs.set("sphi", f"{sphi:.8g}")
        tubs.set("dphi", f"{dphi:.8g}")
        return tubs

    def collect(self, all_volumes, solid_cache=None):
        if solid_cache is not None:
            used = {vol.solid_name for vol in all_volumes}
            for s in solid_cache.unique_solids:
                if s["name"] in used:
                    self.shapes[s["name"]] = self.make_tubs(s["name"], s["rmin"], s["rmax"], s["zhalflength"], s["sphi"], s["dphi"])
            missing = used - set(self.shapes)
            if missing:
                    cache_names = [s["name"] for s in solid_cache.unique_solids]
                    print(f"[dbg] bins reference {len(used)} distinct solids; "
                      f"writer's cache holds {len(cache_names)}")
                    print(f"[dbg] sample referenced: {sorted(used)[:3]}")
                    print(f"[dbg] sample in cache:   {cache_names[:3]}")
                    print(f"[dbg] cache object id in writer: {id(solid_cache)}")
                    raise RuntimeError(f"[{self.name}] {len(missing)} solids referenced but not in this "
                                    f"file's cache, e.g. {next(iter(missing))}")

        for vol in all_volumes:
            if solid_cache is None:
                shape_name = f"{vol.name}_Tubs"
                self.shapes[shape_name] = self.make_tubs(shape_name, vol.rmin, vol.rmax, vol.zhalflength, vol.sphi, vol.sphi)
                phi_rot = 0.0
            else:
                shape_name = vol.solid_name
                phi_rot = vol.phi_rotation

            if vol.material is not None:
                mat_name = vol.material.name
                self.materials[mat_name] = vol.material

            self.logvols[vol.name] = {
                "shape": shape_name,
                "material": mat_name if vol.material is not None else "N2",
                "translation": vol.translation,
                "phi_rotation": phi_rot,
            }


    def write(self, output_path, dtd_path, rmin = 0, rmax = 1200, zhalflength = 3000):
        stem = output_path.replace(".gmx", "")
        self.dtd_path = dtd_path
        mat_path = stem + "material.gmx"
        self.write_material(mat_path)
        self.write_main(stem + ".gmx", mat_path, dtd_path, rmin = rmin, rmax = rmax, zhalflength=zhalflength)
    
    
    def write_material(self, path):
        ELEMENT_DEFS = {
            "hydrogen": ("H",  1,  1.00794),
            "helium": ("He", 2, 4.002602),
            "lithium": ("Li", 3, 6.941),
            "beryllium": ("Be", 4, 9.012182),
            "boron": ("B", 5, 10.811),
            "carbon": ("C",  6, 12.011),
            "nitrogen": ("N",  7, 14.00674),
            "oxygen": ("O",  8, 15.9994),
            "fluorine": ("F", 9, 18.9984032),
            "neon": ("Ne", 10, 20.1797),
            "sodium": ("Na", 11, 22.9897693),
            "magnesium": ("Mg", 12, 24.305),
            "aluminum": ("Al", 13, 26.981539),
            "silicon": ("Si", 14, 28.0855),
            "phosphorus": ("P", 15, 30.973762),
            "sulfur": ("S", 16, 32.07),
            "chlorine": ("Cl", 17, 35.453),
            "argon": ("Ar", 18, 39.948), 
            "titanium": ("Ti", 22, 47.867),
            "vanadium": ("V", 23, 50.9415),
            "chromium":    ("Cr", 24, 51.9961),
            "manganese":   ("Mn", 25, 54.938045),
            "iron": ("Fe", 26, 55.845),
            "cobalt": ("Co", 27, 58.93319),
            "nickel": ("Ni", 28, 58.6934),
            "copper": ("Cu", 29, 63.546),
            "zinc": ("Zn", 30, 65.38),
            "bromine": ("Br", 35, 79.90),
            "zirconium": ("Zr", 40, 91.224),
            "molybdenum": ("Mo", 42, 95.96),
            "silver": ("Ag", 47, 107.8682),
            "tin": ("Sn", 50, 118.710),
            "barium": ("Ba", 56, 137.327),
            "tungsten":    ("W",  74, 183.84),
            "gold": ("Au", 79, 196.966569),
            "lead":        ("Pb", 82, 207.2),
        }
        root = ET.Element("materials")

        used_elements = set()
        for mat_obj in self.materials.values():
            if hasattr(mat_obj, "components") and mat_obj.components:
                for comp, _, _ in mat_obj.components:
                    used_elements.add(comp.name.lower())
        for ename in used_elements:
            if ename not in ELEMENT_DEFS:
                print(f"[GMX] Warning: no definition for element, need to add: '{ename}' to ELEMENT_DFS list")
                continue
            symbol, Z, A = ELEMENT_DEFS[ename]
            e = ET.SubElement(root, "element")
            e.set("name", ename)
            e.set("shortname", symbol)
            e.set("Z", str(Z))
            e.set("A", f"{A:.6g}")
 
        for mat_name, mat_obj in self.materials.items():
            m = ET.SubElement(root, "material")
            m.set("name",    mat_name)
            m.set("density", f"{self.get_density(mat_obj):.6g}")
 
            if hasattr(mat_obj, "components") and mat_obj.components:
                for comp, frac, ftype in mat_obj.components:
                    mass_frac = _eval(frac)
                    tag = "elementref"
                    ref = ET.SubElement(m, tag)
                    ref.set("ref", comp.name.lower())
                    ref.set("fraction", f"{mass_frac:.6g}")
            else:
                ref = ET.SubElement(m, "elementref")
                ref.set("ref",      mat_name)
                ref.set("fraction", "1.0")
 
        if "nitrogen" not in used_elements:
            e = ET.SubElement(root, "element")
            e.set("name",      "nitrogen")
            e.set("shortname", "N")
            e.set("Z",         "7")
            e.set("A",         "14.00674")

        # always write N2, uncomment this section if you are not combining with ITkLayouts, also need to add the dtd path and such back in 
        #n2 = ET.SubElement(root, "material")
        #n2.set("name",    "N2")
        ##n2.set("density", "1.2504e-3 * 273.15/293.15 * 96265/100000")
        #ref = ET.SubElement(n2, "elementref")
        #ref.set("fraction", "1.00000")
        #ref.set("ref",      "nitrogen")
        
        self.write_pretty(root, path)
        print(f"[GMX] Materials written → {path}")

#COMMENT: This write_main was written with the help of Claude. 
    def write_main(self, path: str, mode: str = "envelope", per_layer: bool = False,
               rmin=0, rmax=1200, zhalflength=3000):

        assert mode in ("envelope", "assembly"), f"Unknown mode {mode!r}"

        shapes_elem = ET.Element("shapes")
        for shape_elem in self.shapes.values():
            shapes_elem.append(shape_elem)
        if mode == "envelope":
            env_shape_name = f"SimplifiedPassiveEnvelope_{self.name}"
            env = ET.SubElement(shapes_elem, "tube")
            env.set("name",        env_shape_name)
            env.set("rmin",        str(rmin))
            env.set("rmax",        str(rmax))
            env.set("zhalflength", str(zhalflength))

        if mode == "envelope":
            top_elem = ET.Element("logvol")
            top_elem.set("name",     self.name)
            top_elem.set("shape",    env_shape_name)
            top_elem.set("material", "N2")
            top_elem.set("envelope", "true")
        else:  # assembly
            top_elem = ET.Element("assembly")
            top_elem.set("name", self.name)

        def _emit_transform(parent, idx, vol_name, info):
            """Append a <transform>…</transform> placement for one bin."""
            tx, ty, tz = info["translation"]
            transform = ET.SubElement(parent, "transform")
            tf = ET.SubElement(transform, "transformation")
            tf.set("name", f"tf_{self.name}_{idx}")
            tr = ET.SubElement(tf, "translation")
            tr.set("x", f"{tx:.10g}")
            tr.set("y", f"{ty:.10g}")
            tr.set("z", f"{tz:.10g}")
            phi = info.get("phi_rotation", 0.0)
            if phi != 0.0:
                rot = ET.SubElement(tf, "rotation")
                rot.set("xcos", "0")
                rot.set("ycos", "0")
                rot.set("zcos", "1")
                rot.set("angle", f"{phi:.10g}")
            lv = ET.SubElement(transform, "logvol")
            lv.set("name",     vol_name)
            lv.set("shape",    info["shape"])
            lv.set("material", info["material"])

        sub_assemblies = [] 

        if per_layer and mode == "assembly":
            layer_groups = defaultdict(list)
            for idx, (vol_name, info) in enumerate(self.logvols.items()):
                m = re.search(r'lay(\d+)', vol_name)
                layer = int(m.group(1)) if m else -1
                layer_groups[layer].append((idx, vol_name, info))

            for layer, bins in sorted(layer_groups.items()):
                sub_name = f"{self.name}_layer{layer}"
                sub = ET.Element("assembly")
                sub.set("name", sub_name)
                for idx, vol_name, info in bins:
                    _emit_transform(sub, idx, vol_name, info)
                sub_assemblies.append(sub)
                ref = ET.SubElement(top_elem, "assemblyref")
                ref.set("ref", sub_name)
        else:
            for idx, (vol_name, info) in enumerate(self.logvols.items()):
                _emit_transform(top_elem, idx, vol_name, info)

        def serialise(elem):
            text = self.pretty_string(elem)
            if text.startswith("<?xml"):
                text = text[text.index("?>") + 2:].lstrip()
            return text

        parts = [serialise(shapes_elem)]
        parts.extend(serialise(sub) for sub in sub_assemblies)
        parts.append(serialise(top_elem))
        fragment = "\n".join(parts) + "\n"

        with open(path, "w", encoding="utf-8") as fh:
            fh.write(fragment)
        print(f"[GMX] {self.name} → {path}   (mode={mode}, per_layer={per_layer})")



    @staticmethod
    def get_density(mat) -> float:
        raw = getattr(mat, "density", None)
        if raw is None:
            return 1.0
        try:
            val = _eval(raw)
            return float(val)
        except Exception:
            return 1.0

    @staticmethod
    def pretty_string(root: ET.Element) -> str:
        raw = ET.tostring(root, encoding="unicode")
        dom = minidom.parseString(raw)
        return dom.toprettyxml(indent=" ")

    @staticmethod
    def write_pretty(root: ET.Element, path: str):
        raw = ET.tostring(root, encoding="unicode")
        dom = minidom.parseString(raw)
        pretty = dom.toprettyxml(indent="  ")
        # Strip the <?xml ...?> header – it will be written by the caller if needed
        lines = pretty.split("\n")
        if lines[0].startswith("<?xml"):
            lines = lines[1:]
        with open(path, "w", encoding="utf-8") as fh:
            fh.write('<?xml version="1.0" encoding="UTF-8" ?>\n')
            fh.write("\n".join(lines)) 
    


