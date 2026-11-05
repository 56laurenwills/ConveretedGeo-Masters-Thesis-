
import numpy as np
import math
from dataclasses import dataclass

@dataclass
class GMXBinVolume:
    name: str
    rmin: float
    rmax: float
    zhalflength: float
    sphi: float
    dphi: float
    translation: list
    material: object
    rotation: object
    solid_name: str = None
    phi_rotation: float = 0.0
    is_disc: bool = False
    target_volume = None

class GMX_solid_list:

    def __init__(self, prefix, len_tol=1e-4, ang_tol=1e-6, canonical_phi=True):
        self.prefix = prefix
        self.len_tol = len_tol
        self.ang_tol = ang_tol
        self.canonical_phi = canonical_phi
        self.key_to_name = {}
        self.unique_solids = []
        self.n_requests = 0
        self.n_resused = 0

    @staticmethod
    def rounded(x, q, up):
        eps = 1e-6
        n = x/q #multiplying to get integer so can accurately round
        return math.ceil(n - eps) if up else math.floor(n+eps)


    @staticmethod
    def nearest(x, q):
        return round(x / q)
        

    def get_params(self, vol):
        L, A = self.len_tol, self.ang_tol

        r_rmin = self.nearest(vol.rmin, L) #rounding to prevent overlaps 
        r_rmax = self.nearest(vol.rmax, L)
        if vol.is_disc:
            r_dz = max(self.nearest(vol.zhalflength, L), 1) #avoids problems with rounding down
        else:
            r_dz = self.rounded(vol.zhalflength, L, up=False)
            if r_rmax <= r_rmin:
                r_rmax = r_rmin + 1
        if self.canonical_phi:
            r_sphi = 0
            r_dphi = self.rounded(vol.dphi, A, up=False)
            phi_rotation = vol.sphi
        else:
            r_sphi = self.rounded(vol.sphi, A, up=True)
            end = vol.sphi + vol.dphi
            r_dphi = self.rounded(end - r_sphi * A, A, up=False)
            phi_rotation = 0.0

        self.n_requests += 1
        key = (r_rmin, r_rmax, r_dz, r_sphi, r_dphi)
        name = self.key_to_name.get(key)
        if name is None:
            name = f"{self.prefix}_Tubs_{len(self.unique_solids)}"
            self.key_to_name[key] = name
            self.unique_solids.append(dict(
                name=name,
                rmin = r_rmin* L,
                rmax = r_rmax* L,
                zhalflength = r_dz * L,
                sphi = r_sphi * A,
                dphi = r_dphi * A,

            ))
        else:
            self.n_resused += 1
        return name, phi_rotation

    def report(self):
        n_unique = len(self.unique_solids)
        print(f"[{self.prefix}]: bins: {self.n_requests}, unique solids: {n_unique},"
              f"resused: {self.n_resused}")

