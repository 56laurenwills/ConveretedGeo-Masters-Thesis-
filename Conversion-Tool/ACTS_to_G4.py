import pyg4ometry
import json
import itertools
import pyg4ometry.geant4 as _g4
from pyg4ometry import transformation as _tr
import numpy as np
from materialcreator import JsonMaterialCreator
from materialcreator import RootMaterialCreator
from solidcreator import SolidCreator
from XML_production_assembly import GMXWriter
from GMX_solid_creator import GMXBinVolume, GMX_solid_list
from Json_bin_collection import json_generate_layers
from Root_bin_collection import root_generate_layers, parse_root_material_map
import argparse
from dataclasses import dataclass, field

#This code was written with inspiration from Tong Qiu's original conversion code. 



json_Region_volumes = {
    "Pixel": [8, 9, 10, 13, 14, 15, 16, 18, 19, 20],
    "Strip": [22, 23, 24],
    "HGTD": [2, 25],
    "Barrel": [9, 16, 23],
    "Endcaps": [2, 8, 1, 13, 14, 15, 18, 19, 20, 22, 24, 25],
    "ITk" : [8, 9, 10, 13, 14, 15, 16, 18, 19, 20, 22, 23, 24],
    "All": [2, 6, 8, 9, 10, 13, 14, 15, 16, 18, 19, 20, 22, 23, 24, 25],
}

root_Region_volumes = {
    "Pixel": [4, 5, 6, 9, 10, 11],
    "Strip": [2, 3, 7, 19, 20, 21],
    "Barrel": [5, 10, 20],
    "Endcaps": [4, 6, 9, 11, 19, 21],
    "ITk" : [4, 5, 6, 9, 10, 11, 19, 20, 21],
}
#no sensitive name extraction for geometry regions yet, this is for GDML output
sensitive_name_extrations = {
    "Pixel" : ["InnerBarrelSingleMod_Sensor", "InnerBarrelQuadMod_Sensor", "InnerRingSingleMod_Sensor", 
               "InnerEndcapQuadMod_Sensor", "OuterBarrelQuadMod_Sensor","OuterEndcapQuadMod_Sensor",
                "InclinedQuadMod_Sensor"],
    "Strip" : ["ECSensor0", "ECSensor1", "ECSensor2", "ECSensor3", "ECSensor4", "ECSensor5", 
               "ECSensorBack0", "ECSensorBack1", "ECSensorBack2", "ECSensorBack3", "ECSensorBack4", 
               "ECSensorBack5", "BRLSensorSS", "BRLSensorMS"],
}
Detector_region = {"Pixel", "Strip", "HGTD"}
Geometry_region = {"Barrel", "Endcaps"}

def json_get_active_volumes(regions):
    volume = set()
    for r in regions:
        volume.update(json_Region_volumes[r])
    return volume

def root_get_active_volumes(regions):
    volume = set()
    for r in regions:
        volume.update(root_Region_volumes[r])
    return volume

def get_sensor_names(regions):
    names = set()
    for r in regions:
        names.update(sensitive_name_extrations[r])
    return names


class SingleBin:
    def __init__(self, bin_name, reg, world, solid_creator, gdml_status = False):
        self.bin_name = bin_name
        self.reg = reg
        self.world = world
        self.solid_creator = solid_creator
        self.gdml = gdml_status
        self.physical = None
        self.logical = None
        self.gmx_volume = None

    def set_volume(self, volume):
        self.volume = volume

    def set_material(self, material):
        self.material = material
    
    def set_PhiRange(self, phimin, phimax):
        self.phimin = phimin
        self.phimax = phimax
        self.dphi = abs(phimax - phimin)
    
    def set_ZRange(self, zmin, zmax):
        self.zmin = zmin
        self.zmax = zmax
        self.dz = abs(zmax - zmin)
    
    def set_Rrange(self, rmin, rmax):
        self.rmin = rmin
        self.rmax = rmax
        self.dr = abs(rmax - rmin)
    
    def set_translation(self, x, y, z):
        self.translation = [x, y, z]

    def set_rotation(self, x, y, z):
        self.rotation = [x, y, z]

    def disc_status(self, disc_status):
        self.disc_status = disc_status
    
    
    def create_volume(self):
        self.physical = None
        self.logical = None
        if self.dz == 0 or self.dr == 0 or self.dphi == 0 or self.material is None:
            return
        if self.rmin < 0.0001:
            self.rmin = 0.0001
        self.set_rotation(0, 0, self.phimin)
        if self.gdml:
            self.create_gdml_volume()
        else:
            self.create_gmx_volume()
            

    def create_gdml_volume(self):
        if self.disc_status:
            self.solid = self.solid_creator.create_tub(self.rmin+ 0.00001, self.rmax- 0.00001, self.dz, self.dphi- 0.001)
        else:
            self.solid = self.solid_creator.create_tub(self.rmin+ 0.00001, self.rmax- 0.00001, self.dz - 0.0001, self.dphi- 0.001)
        self.set_rotation(0, 0, self.phimin)
        self.logical = pyg4ometry.geant4.LogicalVolume(self.solid, self.material, f"{self.bin_name}_Logical", self.reg)
        self.physical = pyg4ometry.geant4.PhysicalVolume(self.rotation, self.translation, self.logical, f"{self.bin_name}_Physical", self.world, self.reg)

    def create_gmx_volume(self):
        self.gmx_volume = GMXBinVolume(
                        name = self.bin_name,
                        rmin = self.rmin,
                        rmax = self.rmax,
                        zhalflength = self.dz / 2.0,
                        sphi = self.phimin,
                        dphi = self.dphi,
                        translation = list(self.translation),
                        material = self.material,
                        rotation = self.rotation,
                        is_disc = self.disc_status,
            )



class Detector:
    def __init__(self, name, material_map_input, geometryR=None):
        self.name = name
        self.reg = pyg4ometry.geant4.Registry()
        self.solid_creator = SolidCreator(self.reg, tolerance=0)

        self.new_logicals = []
        self.new_physicals = []
        self.all_volumes = []
        self.solid_caches = {}

        if material_map_input.endswith('json'):
            with open(material_map_input) as mm:
                material_map = json.load(mm)
            material_map = material_map['Surfaces']['entries']
            self.materialmap = material_map
            self.input_type = 'json'
            self.material_creator = JsonMaterialCreator(self.reg)
            self.geometryR = geometryR
        elif material_map_input.endswith('root'):
            self.materialmap = parse_root_material_map(material_map_input)
            self.input_type = 'root'
            self.material_creator = RootMaterialCreator(self.reg)
            self.geometryR = None
    
    def create_world(self):
        self.world_solid = pyg4ometry.geant4.solid.Box("worldsolid", 50000, 50000, 100000, self.reg)
        self.create_air()
        self.world_logical = pyg4ometry.geant4.LogicalVolume(self.world_solid, self.air, "worldlogical", self.reg)
        self.reg.setWorld(self.world_logical.name)

    def use_world(self, logical):
        self.world_solid = logical.solid
        self.reg.transferSolid(logical.solid)
        self.reg.transferMaterial(logical.material)
        self.world_logical = pyg4ometry.geant4.LogicalVolume(self.world_solid, logical.material, "worldlogical", self.reg)
        self.reg.setWorld(self.world_logical.name)

    def add_volume_world(self, logical_sens, physical_sens, rot_sens, tra_sens):
        self.reg.transferSolid(logical_sens.solid)
        self.reg.transferMaterial(logical_sens.material)
        self.reg.transferLogicalVolume(logical_sens)

        tbxyz = _tr.matrix2tbxyz(np.linalg.inv(rot_sens))
        rot = [tbxyz[0], tbxyz[1], tbxyz[2], "rad"]
        pos = [float(tra_sens[0]), float(tra_sens[1]), float(tra_sens[2]), "mm"]

        new_physical = pyg4ometry.geant4.PhysicalVolume(rot, pos, logical_sens, physical_sens.name, self.reg.getWorldVolume(), self.reg)
        self.new_physicals.append(new_physical)
        self.all_volumes.append((new_physical, logical_sens, rot_sens, tra_sens))


    def add_sensitives_from_world(self, lv_world, reg, regions):
        self.materialnames = set()
        sensitive_hits = walk_sensitive(lv_world, reg, regions)
        print("Found sensor placements:", len(sensitive_hits))
        for pv, lv, rot_A, tra_A in sensitive_hits:
            self.add_volume_world(lv, pv, rot_A, tra_A)
            self.materialnames.add(lv.material.name)


    def create_air(self):
        self.c12 = _g4.Isotope("C12", 6, 12, 12)
        self.c13 = _g4.Isotope("C13", 6, 13, 13.0034)
        self.carbon = _g4.ElementIsotopeMixture("worldcarbon", "C", 2)
        self.carbon.add_isotope(self.c12, 0.9893)
        self.carbon.add_isotope(self.c13, 0.0107)

        self.n14 = _g4.Isotope("N14", 7, 14, 14.0031)
        self.n15 = _g4.Isotope("N15", 7, 15, 15.0001)
        self.nitrogen = _g4.ElementIsotopeMixture("worldnitrogen", "N", 2)
        self.nitrogen.add_isotope(self.n14, 0.99632)
        self.nitrogen.add_isotope(self.n15, 0.00368)

        self.o16 = _g4.Isotope("O16", 8, 16, 15.9949)
        self.o17 = _g4.Isotope("O17", 8, 17, 16.9991)
        self.o18 = _g4.Isotope("O18", 8, 18, 17.9992)
        self.oxygen = _g4.ElementIsotopeMixture("worldoxygen", "O", 3)
        self.oxygen.add_isotope(self.o16, 0.99757)
        self.oxygen.add_isotope(self.o17, 0.00038)
        self.oxygen.add_isotope(self.o18, 0.00205)

        self.ar36 = _g4.Isotope("Ar36", 18, 36, 35.9675)
        self.ar38 = _g4.Isotope("Ar38", 18, 38, 37.9627)
        self.ar40 = _g4.Isotope("Ar40", 18, 40, 39.9624)
        self.argon = _g4.ElementIsotopeMixture("worldargon", "Ar", 3)
        self.argon.add_isotope(self.ar36, 0.003365)
        self.argon.add_isotope(self.ar38, 0.000632)
        self.argon.add_isotope(self.ar40, 0.996003)

        self.air = _g4.MaterialCompound("worldair", 0.00120479, 4, self.reg)
        self.air.add_element_massfraction(self.carbon, 0.000124000124000124)
        self.air.add_element_massfraction(self.nitrogen, 0.755267755267755)
        self.air.add_element_massfraction(self.oxygen, 0.231781231781232)
        self.air.add_element_massfraction(self.argon, 0.0128270128270128)

        self.air.set_temperature(293.15)
    

    def generate_layers(self, regions, geo_tag, gdml_status):
        print("Generating detector")
        print(f"input_type: {self.input_type}")
        print(f"materialmap length: {len(self.materialmap)}") 
        if self.input_type == 'json':
            active_volumes = json_get_active_volumes(regions)
            json_generate_layers(self, active_volumes)
        elif self.input_type == 'root':
            active_volumes = root_get_active_volumes(regions)
            print(f"active_volumes: {active_volumes}")
            root_generate_layers(self, active_volumes, geo_tag, gdml_status)
        print("Detector generation complete")

    def save(self, filename, gdml = False): 
        if gdml:
            w = pyg4ometry.gdml.Writer()
            w.addDetector(self.reg)
            w.write(f"{filename}.gdml")
        else:
            dtd_path = ("/home/lwills/GeoModel/GeoModelTools/GeoModelXML/GeoModelXml/data/geomodel.dtd")
            from XML_production_assembly import write_split
            write_split(self.all_volumes, filename, dtd_path, self.solid_caches)
            

    def visualize(self):
        print("Visualizing detector")
        v = pyg4ometry.visualisation.VtkViewer()
        v.addLogicalVolume(self.world_logical)
        v.addAxes(20)
        v.view()

class CylinderR:
    def __init__(self, geometrymapfile):
        with open(geometrymapfile) as f:
            geometrymap = json.load(f)
        self.Rdic = {}
        geometrymap = geometrymap["Surfaces"]["entries"]
        for each in geometrymap:
            if "layer" not in each or "approach" not in each:
                continue
            if each["value"]["bounds"]["type"] != "CylinderBounds":
                continue
            key = self.getindex(int(each["volume"]), int(each["layer"]), int(each["approach"]))
            if key in self.Rdic:
                print(f"Error: CylinderR key {key} already exists in Rdic")
                exit(1)
            self.Rdic[key] = float(each["value"]["bounds"]["values"][0])

    def getR(self, volumeid, layerid, approachid):
        key = self.getindex(volumeid, layerid, approachid)
        if key not in self.Rdic:
            return -1
        return self.Rdic[key]

    def getindex(self, volumeid, layerid, approachid):
        return volumeid * 1000000 + layerid * 1000 + approachid
    
def build_sensitive_filter(regions):
    type_keywords = []
    if "All" in regions:
        return lambda pv: "Sensor" in pv.name
    if "Pixel" in regions:
        type_keywords.append("Pixel")
    if "Strip" in regions:
        type_keywords.append("Strip")
    if "HGTD" in regions:
        type_keywords.append("HGTD")

    geo_filter = None
    if "Barrel" in regions:
        geo_filter = lambda pv: "barrel_endcap_0" in pv.name
    elif "Endcap" in regions:
        geo_filter = lambda pv: "barrel_endcap_0" not in pv.name
    
    def pv_filter(pv):
        if "Sensor" not in pv.name:
            return False
        if type_keywords and not any (kw in pv.name for kw in type_keywords):
            return False
        if geo_filter is not None and not geo_filter(pv):
            return False
        return True
    return pv_filter


def pv_local_trans(pv):
        t = np.array(pv.position.eval(), dtype=float) if pv.position else np.zeros(3)
        R = np.linalg.inv(_tr.tbxyz2matrix(pv.rotation.eval())) if pv.rotation else np.eye(3)
        if getattr(pv, "scale", None):
            R = np.diag(pv.scale.eval()) @ R
        return R, t

def compose_trans(parent_Aa, child_Aa):
    A, a = parent_Aa
    B, b = child_Aa
    return (A @ B, A @ b + a)

def lv_has_sensdet(lv):
        if len(lv.auxiliary) > 0:
            if lv.auxiliary[0].auxtype == "SensDet":
                return True
        return False

def name_sensor_pv(pv, regions):
    name_list = get_sensor_names(regions)
    if any(name in pv.name for name in name_list):
        return True
    return False

def gdml_has_sensdet(reg):
    for lv in reg.logicalVolumeDict.values():
        if len(lv.auxiliary) > 0:
            if lv.auxiliary[0].auxtype == "SensDet":
                return True
    return False


def walk_sensitive(lv_world, reg, regions):
    use_aux = gdml_has_sensdet(reg)
    out = []
    
    stack = []
    for top_pv in lv_world.daughterVolumes:
        if "HGTD" in top_pv.name:
            continue
        stack.append((top_pv, (np.eye(3), np.zeros(3))))

    while stack:
        pv, Aa = stack.pop()
        lv = pv.logicalVolume
        if use_aux:
            match = lv_has_sensdet(lv)
        else:
            match = name_sensor_pv(pv, regions)
        if match:
            A, a = Aa
            out.append((pv, lv, A.copy(), a.copy()))
        for child in lv.daughterVolumes:
            Aa_child = compose_trans(Aa, pv_local_trans(child))
            stack.append((child, Aa_child))
    return out

def converting_ACTS_2_geant4(gdml, material_map_input, geometry_map, geo_tag, output_name, regions, include_sensitives, gdml_status):
    if gdml:
        print("GDML included- Reading")
        itk = pyg4ometry.gdml.Reader(gdml)
        reg = itk.getRegistry()
        world_volume = reg.worldVolume
    print("Importing material map")
    geometry_r = CylinderR(geometry_map) if geometry_map is not None else None

    converted_detector = Detector("output_name", material_map_input, geometry_r)
    converted_detector.create_world()
    print("Generating passive volumes from material map surfaces")
    converted_detector.generate_layers(regions, geo_tag, gdml_status)
    materialnames = set()
    if include_sensitives:
        print("Extracting sensitives turned on:")
        #filter_sensitives = build_sentive_filter(regions)
        converted_detector.add_sensitives_from_world(world_volume, reg, regions)
    print(f"Saving Detector to {output_name}")
    converted_detector.save(output_name, gdml_status)


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--gdml", type=str, default = "", help="The inputted GDML, needed to extract sensitives and the world volume"
    )
    p.add_argument(
        "--material", type=str, required = True, help="The material map to convert to Geant4- JSON or ROOT format accepted"
    )
    p.add_argument(
        "--geometry", type=str, required = False, default=None, help="The geometry used to create the material map (needed for the R value for cylinders if using json)"
    )
    p.add_argument(
        "--geotag", type=str, required = False, default="Run4-05-00", help="The Athena geometry tag used to create the material map"
    )
    p.add_argument(
        "--output", type=str, default="output_dir", help="The directory for outputted gmx's or gdml"
    )

    p.add_argument(
        "--regions", nargs="+", choices=list(json_Region_volumes.keys()), default=["ITk"], 
        help = "Which detector region that you want to convert, can select multiple. For Root no HGTD, no geometry descriptors either"
    )
    p.add_argument(
        "--gdml_output",action='store_true', default=False, help="Store the output as a gdml (need to provide an input), default if GeoModel XML"
        )

    p.add_argument(
        "--sensitives",action='store_true', default=False, help="Do you want to pull out the sensitives from the original GDML"
    )
    

    args = p.parse_args()

    chosen = set(args.regions)
    if chosen & Detector_region and chosen & Geometry_region:
        p.error("Cannot mix region types like ''Pixel'/ 'Strip'/ 'Endcap'")

    if args.sensitives and not args.gdml_output:
        p.error("Cannot pull out sensitives without gdml-output")

    if args.sensitives and not args.gdml:
        p.error("Cannot pull out sensitives without a provided gdml")
    converting_ACTS_2_geant4(args.gdml, args.material, args.geometry, args.geotag, args.output, args.regions, args.sensitives, args.gdml_output)
if __name__ == "__main__":
    main()
