
import json
import numpy as np
import pyg4ometry.geant4 as _g4
from pyg4ometry import transformation as _tr


class Json_DetectorBinCollection:
    def __init__(self, name, reg, world, materialmap, material_creator, solid_creator, r, approach):
        self.name = name
        self.reg = reg
        self.world = world
        self.materialmap = materialmap
        self.solid_creator = solid_creator
        self.material_creator = material_creator
        self.approach = int(approach)
        self.r = r 
        self.bins = []
        self.get_meta_physical()
        self.get_meta_binning()
        self.generate_allbins()
        
    
    def get_meta_physical(self):
        meta_map = self.materialmap["binUtility"]
        self.translation = None
        self.rotation = None
        if "transform" in meta_map:
            self.translation = meta_map["transform"]["translation"]
            self.rotation = meta_map["transform"]["rotation"]

        if self.translation is None:
            self.translation = [0, 0, 0]
        if self.rotation is None:
            self.rotation = [0, 0, 0]

    def get_meta_binning(self):
        meta_map = self.materialmap["binUtility"]
        self.binPhi = []
        self.binR = []
        self.binZ = []
        self.binorder = []
        self.bin_collection = []
        for each_binmetamap in meta_map["binningdata"]:
            each_binning_tem = np.linspace(each_binmetamap["min"], each_binmetamap["max"], each_binmetamap["bins"]+1)
            each_binning = [[each_binning_tem[i], each_binning_tem[i+1]] for i in range(len(each_binning_tem)-1)]
            if each_binmetamap["value"] == "binR" or each_binmetamap["value"] == "AxisR":
                self.binR = each_binning
                self.binorder.append(["binR", self.binR])
            elif each_binmetamap["value"] == "binPhi" or each_binmetamap["value"] == "AxisPhi":
                self.binPhi = each_binning
                self.binorder.append(["binPhi", self.binPhi])
            elif each_binmetamap["value"] == "binZ" or each_binmetamap["value"] == "AxisZ":
                self.binZ = each_binning
                self.binorder.append(["binZ", self.binZ])
            else:
                raise ValueError(f"Unknown binning type {each_binmetamap['value']}")


    def process_singlebin(self, name, binR, binPhi, binZ, material, disc_surface, which_approach, thickness, each_material_description_ZAF=None):
        from ACTS_to_G4 import SingleBin
        if len(each_material_description_ZAF) == 0:
            g4_material = None
        else:
            g4_material = self.material_creator.create_material(each_material_description_ZAF, material[4] * material[2] * 1000)
        singlebin = SingleBin(f"{self.name}_{name}", self.reg, self.world, self.solid_creator)
        singlebin.set_Rrange(binR[0], binR[1])
        singlebin.set_ZRange(binZ[0], binZ[1])
        singlebin.set_PhiRange(binPhi[0], binPhi[1])
        singlebin.set_material(g4_material)
        if disc_surface:
            if which_approach == 1:
                z_shift = -thickness/2
            elif which_approach == 2:
                z_shift = thickness /2
            else:
                z_shift = 0.0
        else:
            z_shift = 0.0
        if disc_surface == True:
            singlebin.set_translation(self.translation[0], self.translation[1], (self.translation[2]+ z_shift))
        else:
            z_bin_shift = (binZ[0] + binZ[1])/2
            if which_approach == 1:
                singlebin.set_translation(self.translation[0], self.translation[1], (self.translation[2] - z_bin_shift))
            elif which_approach == 2:
                singlebin.set_translation(self.translation[0], self.translation[1], (self.translation[2] + z_bin_shift))
            else:
                singlebin.set_translation(self.translation[0], self.translation[1], self.translation[2])
        singlebin.set_rotation(self.rotation[0], self.rotation[1], self.rotation[2])
        singlebin.create_volume()
        return singlebin


    def generate_allbins(self):
        material_map = self.materialmap["data"]
        approach = self.approach
        self.is_disc = False
        if len(material_map[0]) != len(self.binorder[0][1]):
            raise ValueError("Number of bins in binningdata and material data do not match")
        if len(material_map) != len(self.binorder[1][1]):
            raise ValueError("Number of bins in binningdata and material data do not match")

        for i in range(len(material_map)):
            each_bins = []
            for j in range(len(material_map[0])):
                binname = f"bin_{i}_{j}"
                each_material_map = material_map[i][j]
                each_material_discription = each_material_map["material"]
                each_material_discription_ZAF = each_material_map["ZAF"]
                #each_thickness = float(each_material_map["thickness"])
                each_thickness = 10.0
                binPhi = binZ = [-each_thickness/2, each_thickness/2]
                binR = [0, each_thickness]
                isrepresentative = False
                if self.binorder[0][0] == "binR":
                    binR = self.binorder[0][1][j]
                    self.is_disc = True
                if self.binorder[0][0] == "binPhi":
                    binPhi = self.binorder[0][1][j]
                if self.binorder[0][0] == "binZ":
                    binZ = self.binorder[0][1][j]
                if self.binorder[1][0] == "binR":
                    binR = self.binorder[1][1][i]
                    self.is_disc = True
                if self.binorder[1][0] == "binPhi":
                    binPhi = self.binorder[1][1][i]
                if self.binorder[1][0] == "binZ":
                    binZ = self.binorder[1][1][i]
                if self.binorder[0][0] in ("binPhi", "binZ") and self.binorder[1][0] in ("binPhi", "binZ"):
                    if approach == 1:
                        binR = [self.r - each_thickness,  self.r]
                    elif approach == 2:
                        binR = [self.r, self.r + each_thickness]
                    else:
                        binR = [self.r - each_thickness/2, self.r + each_thickness/2]
                if not isrepresentative:
                    eachbin = self.process_singlebin(binname, binR, binPhi, binZ, each_material_discription, self.is_disc, approach, each_thickness, each_material_discription_ZAF)
                    each_bins.append(eachbin)
            self.bin_collection.append(each_bins)

def json_check_layer(layer):
    json_bins = layer["value"]["material"]["data"]
    for each_bin in list(itertools.chain(*json_bins)):
        if float(each_bin["thickness"]) > 0:
            return True
    return False

def json_generate_layers(detector, active_volumes):
    detector.layers = {}
    i = 0
    count = 0
    for each_layer in detector.materialmap:
        name = f"dic_{i}_"
        if "volume" in each_layer:
            name += f"v{each_layer['volume']}_"

        if "layer" in each_layer:
            name += f"l{each_layer['layer']}_"

            if each_layer['volume'] not in active_volumes:
                continue
        if "approach" in each_layer:
            name += f"a{each_layer['approach']}_"
        else:
            continue

        if "boundary" in each_layer:
            name += f"{each_layer['boundary']}_"
            continue

        name = name[:-1]
        #if detector.json_check_layer(each_layer):
        if name in detector.layers:
            print(f"Layer {name} already exists")
        count += 1
        #if True:
        print(f"Creating layer {name}")
        r = detector.geometryR.getR(int(each_layer['volume']), int(each_layer['layer']), int(each_layer['approach']))
        approach = int(each_layer['approach'])
        each_bin_collection = Json_DetectorBinCollection(name, detector.reg, detector.world_logical, each_layer["value"]["material"], detector.material_creator, detector.solid_creator, r, approach)
        detector.layers[name] = each_bin_collection



