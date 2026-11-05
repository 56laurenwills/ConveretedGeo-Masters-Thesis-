import pyg4ometry.geant4 as _g4
import pyg4ometry

class SolidCreator:
    def __init__(self, reg, tolerance=0.001):
        self.tubs = {}
        self.tubs_lookup = {}
        self.tolerance = tolerance
        self.reg = reg
        self.reuse = True

    def compare_parameters(self, params1, params2):
        for key in params1:
            if params1[key] != 0:
                if abs((params1[key] - params2[key])/params1[key]) > self.tolerance:
                    return False
            else:
                if abs(params1[key] - params2[key]) != 0:
                    return False
        return True

    def create_tub(self, rmin, rmax, dz, dphi):
        parameters = {
            "rmin": rmin,
            "rmax": rmax,
            "dz": dz,
            "dphi": dphi
        }
        newname = f"Solid_Tubs_{len(self.tubs)}"
        #outputsolid = pyg4ometry.geant4.solid.Tubs(newname, parameters["rmin"], parameters["rmax"], parameters["dz"], 0, parameters["dphi"], self.reg)
        #self.tubs[newname] = outputsolid
        #return outputsolid

        ifcreate = True
        if self.reuse:
            for each_key in self.tubs:
                if self.compare_parameters(parameters, self.tubs[each_key][1]):
                    outputsolid = self.tubs[each_key[0]]
                    ifcreate = False
                    break
        


#my stuff?
        # if self.reuse:
        #     if self.tolerance == 0:
        #         param_key = (rmin, rmax, dz, dphi)
        #         if param_key in self.tubs_lookup:
        #             return self.tubs_lookup[param_key]
        #     else:
        #         for each_key in self.tubs:
        #             if self.compare_parameters(parameters, self.tubs[each_key][1]):
        #                 outputsolid = self.tubs[each_key][0]

        # ifcreate = True
        # if self.reuse:
        #     for each_key in self.tubs:
        #         if self.compare_parameters(parameters, self.tubs[each_key][1]):
        #             outputsolid = self.tubs[each_key][0]
        #             ifcreate = False
        #             break

        if ifcreate:
             outputsolid = pyg4ometry.geant4.solid.Tubs(newname, parameters["rmin"], parameters["rmax"], parameters["dz"], 0, parameters["dphi"], self.reg)
             self.tubs[newname] = [outputsolid, parameters]
        #     #print(f"Created new solid: {newname}")
        else:
            print(f"reuse existing solid: {each_key}")
            pass

        return outputsolid

if __name__ == "__main__":
    reg = pyg4ometry.geant4.Registry()
    solidobj = SolidCreator(reg)
    create_tub = solidobj.create_tub(10, 20, 0.3, 3)
    create_tub = solidobj.create_tub(10, 20, 0.3, 3)
    create_tub = solidobj.create_tub(10, 20, 0.3, 3.000000001)
    create_tub = solidobj.create_tub(10, 20, 0.3, 3.1)
    create_tub = solidobj.create_tub(10, 20, 0.3, 3.100000001)
