#Root_bin_collection.py

import uproot
import sys
import pyg4ometry.geant4 as _g4
from pyg4ometry import transformation as _tr
from Geometry_configs import Geometry_configs
import numpy as np
from collections import Counter
from GMX_solid_creator import GMX_solid_list
from XML_production_assembly import target_volume_for


class Surface:
    def __init__(self, name, n_bins, bin_names, option, mini, maxi, material):
        self.name = name
        self._parse_name()
        self.n_bins = n_bins
        self.bin_names = bin_names
        self.option = option
        self.min = mini
        self.max = maxi
        
        self.material = material

    def _parse_name(self):
        parts = self.name.replace('SurfaceMaterial_', '').split('_')
        self.volume = int(parts[0].replace('vol', ''))
        self.boundary = int(parts[1].replace('bou', ''))
        self.layer = int(parts[2].replace('lay', ''))
        self.approach = int(parts[3].replace('app', ''))

def parse_root_material_map(root_file_path):
    f = uproot.open(root_file_path)

    indexed_material = f['IndexedMaterial']
    all_elementZ = indexed_material['ElementZ'].array()
    all_elementFrac = indexed_material['ElementFrac'].array()
    all_thick = indexed_material['t'].array()
    eff_rho = indexed_material['rho'].array()
    
    if len(all_elementZ) == 0:
        print("Incompatible Root file, RecordedElementFrac from ACTS not included")
        sys.exit()
    
    grouped = {}
    for key in f.keys():
        key_clean = key.split(';')[0]
        if '/' in key_clean:
            surface_name, sub = key_clean.rsplit('/', 1)
        else:
            surface_name = key_clean
            sub = None
        if surface_name not in grouped:
            grouped[surface_name] = {}
        if sub:
            grouped[surface_name][sub] = f[key]
    
    surfaces = []
    surface_grouped = {name: hists for name, hists in grouped.items() if name.startswith('SurfaceMaterial')}
    for name, histograms in surface_grouped.items():
        n_bins = histograms['n'].to_numpy()
        bin_names = histograms['v'].to_numpy()
        option = histograms['o'].to_numpy()
        mini = histograms['min'].to_numpy()
        maxi = histograms['max'].to_numpy()
        i_hist = histograms['i'].to_numpy()
        i_array = i_hist[0]
        surface_material = {
            'elementZ' : all_elementZ[i_array],
            'elementFrac' : all_elementFrac[i_array],
            'thickness' : all_thick[i_array],
            'mass_density' : eff_rho[i_array]
        }

        surfaces.append(Surface(name, n_bins, bin_names, option, mini, maxi, surface_material))
    
    return surfaces

class Root_DetectorBinCollection:
    def __init__(self, name, reg, world_logical, surface, material_creator, solid_creator, gdml_status, geo_tag, solid_caches, all_volumes = None):
        self.name = name
        self.reg = reg
        self.world = world_logical
        self.materialmap = surface
        self.volume = surface.volume
        self.layer = surface.layer
        self.boundary = surface.boundary
        self.solid_creator = solid_creator
        self.r_value = []
        self.material_creator = material_creator
        self.gdml = gdml_status
        self.bins = []
        self.geo_tag = geo_tag
        self.all_volumes = all_volumes if all_volumes is not None else []
        self.cfg = Geometry_configs[geo_tag]
        self.z_value_list = self.cfg["z_lookup_endcaps"]
        self.reuse = True
        self.target_volume = target_volume_for(self.volume, self.boundary, self.layer)
        self.solid_caches = solid_caches
        self.get_binning_data()

    def get_z_translation(self):
        surface_name = (self.volume, self.boundary, self.layer)
        if surface_name not in self.z_value_list:
            print(f"name of surface not in z_value_list: {surface_name}")
            return 3000.0
        
        z = self.z_value_list[surface_name]
        if self.volume in (2, 3, 4, 6, 9, 11, 19, 21):
            return z if self.volume in (3, 6, 11, 21) else -z
        if self.volume in (1, 5, 10, 20):
            return z if self.boundary == 2 else -z 
        
    def cache_for (self, target):
        if target not in self.solid_caches:
            self.solid_caches[target] = GMX_solid_list(prefix=f"vol{target}")
        return self.solid_caches[target]

    def get_binning_data(self):
        n_vals = self.materialmap.n_bins[0]
        min_vals = self.materialmap.min[0]
        max_vals = self.materialmap.max[0]
        names = self.materialmap.bin_names[0]
        options = self.materialmap.option[0]

    
        self.bin_collection = []
        self.r_value = []
        self.is_disc = False
        truncate_r = True
        axis_0 = str(int(names[0]))
        axis_1 = str(int(names[1]))
        pair = (axis_0, axis_1)
        #axis names: 2 = Z, 3 = R, 4 = Phi, 5 = RPhi
        valid_pairs = {('3','4'), ('4','3'), ('5','2'), ('2','5')}
        for axis_name in (axis_0, axis_1):
            if axis_name in ('0', '1', '6', '7', '8'):
                print(f"Skipping {self.materialmap.name}: unsupported axis {axis_name}")
                return
        if pair not in valid_pairs:
            print(f"Skipping {self.materialmap.name}: unsupported axis pair {pair}")
            return
    
        #need to make sure that the phi values have a close boundary (-pi, pi) for the logic to work 
        for axis in range(2):
           if str(int(names[axis])) in ('5') and int(options[axis]) != 1:
               print(f"Skipping {self.materialmap.name}: axis {names[axis]} is not closed")
               return
    
        edges_0 = np.linspace(min_vals[0], max_vals[0], int(n_vals[0]) + 1)
        edges_1 = np.linspace(min_vals[1], max_vals[1], int(n_vals[1]) + 1)
        bins_0 = [[edges_0[i], edges_0[i+1]] for i in range(len(edges_0) - 1)]
        bins_1 = [[edges_1[i], edges_1[i+1]] for i in range(len(edges_1) - 1)]

        if pair in (('5', '2'), ('2', '5')):
            rphi_bins = bins_0 if axis_0 == '5' else bins_1
            z_bins    = bins_0 if axis_0 == '2' else bins_1

            rphi_min = min_vals[0] if axis_0 == '5' else min_vals[1]
            rphi_max = max_vals[0] if axis_0 == '5' else max_vals[1]
            R = (rphi_max - rphi_min) / (2 * np.pi)
            #add in an r shift to the barrel layers given the modulation of the sensors 
            if self.volume == 20: 
                R_shift = R + 15.0
            elif self.volume == 7:
                R_shift = R + 40.0
            else:
                R_shift = R + 5
            rphi_shift = -rphi_min 
    
        for i in range(len(self.materialmap.material['elementZ'])):
            each_bin = []
            for j in range(len(self.materialmap.material['elementZ'][i])):
                binname = f"bin_{i}_{j}"
                material = {key: self.materialmap.material[key][i][j] for key in self.materialmap.material}
                bin_thickness = self.materialmap.material['thickness'][i][j]
                target = self.target_volume
                if pair in (('5', '2'), ('2', '5')):
                    binR = [R_shift - (bin_thickness/2), R_shift + (bin_thickness/2)]
                    binZ = z_bins[j]
                    phi_min = (rphi_bins[i][0] + rphi_shift)/ R
                    phi_max = (rphi_bins[i][1] +rphi_shift) / R
                    binPhi = [phi_min, phi_max]
                    self.r_value.append(R)
    
                elif pair in (('3', '4'), ('4', '3')):
                    # disc (R, Phi)
                    r_bins = bins_0 if axis_0 == '3' else bins_1
                    phi_bins = bins_0 if axis_0 == '4' else bins_1
                    binR = r_bins[i]
                    binPhi = phi_bins[j]
                    binZ = [-bin_thickness / 2, bin_thickness / 2]
                    self.is_disc = True
                OuterMost_boundary_discs = [2, 3]
                target = self.target_volume
                if self.volume == 7: #adjusting where the volume 7 bins actually end up. This is just to make it work with the envelopes in ITkStrip
                    if j <= 24:
                        target = 2
                    if j >= 75:
                        target = 3 
                #this is for truncating the radial bins in the outer most boundary so that the inner part goes to ITkPixel and the outer part goes to ITkStrip
                if truncate_r:
                    if self.volume in OuterMost_boundary_discs and not (17 <= i):
                        continue
                    each_bin.append(self.create_binvolume(binname, binR, binZ, binPhi, material, bin_thickness, self.is_disc, target_volume = target))
                else:
                    each_bin.append(self.create_binvolume( binname, binR, binZ, binPhi, material, bin_thickness, self.is_disc, target_volume = target))
            self.bin_collection.append(each_bin)
        #debug print to figure out R value
        for r, count in sorted(Counter(self.r_value).items()):
            print(f"R = {r:.3f}, {count} counts")



    def create_binvolume(self, name, binR, binZ, binPhi, bin_material, thickness, disc_status, target_volume=None):
        from ACTS_to_G4 import SingleBin
        if target_volume is None:
            target_volume = self.target_volume
        if bin_material['thickness'] > 0.001:
            g4_material = self.material_creator.create_material(bin_material)
        else:
            g4_material = None
        singlebin = SingleBin(f"{self.name}_{name}", self.reg, self.world, self.solid_creator, self.gdml)
        singlebin.set_Rrange(binR[0], binR[1])
        singlebin.set_ZRange(binZ[0], binZ[1])
        singlebin.set_PhiRange(binPhi[0], binPhi[1])
        singlebin.set_material(g4_material)
        singlebin.disc_status(disc_status)

        #shifts added to minimize overlaps 
        if disc_status:
            z_translation = self.get_z_translation()
            if self.volume in (19, 21):
                z_shift = 25.0
            elif self.volume in (4, 6):
                 z_shift = 15.0
            elif self.volume in (9, 11):
                if self.boundary == 1:
                    z_shift = 10.0
                else:
                    z_shift = 0.0
            elif self.volume in (1, 3, 2, 5, 10, 20):
                 z_shift = 0.0
            if self.volume in (3, 5, 10, 6, 11, 21):
                 singlebin.set_translation(0, 0, z_translation + z_shift)
            elif self.volume in (2, 4, 9, 19, 20):
                 singlebin.set_translation(0, 0, z_translation - z_shift)
        else:
            singlebin.set_translation(0, 0, (binZ[0] + binZ[1]) / 2)
        singlebin.create_volume()
        vol = singlebin.gmx_volume
        if vol is not None and singlebin.material is not None:
            vol.target_volume = target_volume
            if self.reuse:
                cache = self.cache_for(target_volume)
                vol.solid_name, vol.phi_rotation = cache.get_params(vol)
            else:
                vol.solid_name = f"{vol.name}_Tubs"
                vol.phi_rotation = 0.0
            self.all_volumes.append(vol)
        return singlebin
    


def Root_check_layer(surface):
    return bool(np.any(surface.material['thickness'] > 0))

def root_generate_layers(detector, active_volumes, geo_tag, gdml_status):
    print(f"root_generate_layers called with {len(detector.materialmap)} surfaces")
    detector.layers = {}
    count = 0
    boundaries = False #to include the large encapsulating boundary surfaces or not
    cfg = Geometry_configs[geo_tag]
    boundary_vol = cfg["boundary_volumes"]
    
    for surface in detector.materialmap:
        #for specific layer debugging 
        #layer_list = [23, 24]
        sur_name = (surface.volume, surface.boundary, surface.layer)
        if surface.volume not in active_volumes:
            continue

        if not boundaries:
            if sur_name in boundary_vol:
                continue
            #if surface.layer not in layer_list:
                #continue

        name = f"dic_vol{surface.volume}_bou{surface.boundary}_lay{surface.layer}_app{surface.approach}"

        if Root_check_layer(surface):
            if name in detector.layers:
                print(f"Layer {name} already exists")
            count += 1
            print(f"Creating layer {name}")

            each_bin_collection = Root_DetectorBinCollection(
                name,
                detector.reg,
                detector.world_logical,
                surface,
                detector.material_creator,
                detector.solid_creator,
                gdml_status,
                geo_tag,
                detector.solid_caches,
                detector.all_volumes, 
            )
            detector.layers[name] = each_bin_collection

        else:
            print(f"Skipping layer {name}")
    for vol_id in sorted(detector.solid_caches):
        detector.solid_caches[vol_id].report()

