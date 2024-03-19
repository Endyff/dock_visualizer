import numpy as np
from rdkit import Chem
from spyrmsd.rmsd import symmrmsd

class Atom(dict):
    def __init__(self,):
        pass

    def dist(self, other):
        return ((self["x"] - other["x"])**2 + (self["y"] - other["y"])**2 + (self["z"] - other["z"])**2)**0.5

    def add_style(self, style:dict):
        self["style"] = {**self["style"], **style} if "style" in self else style

    def get_style(self):
        return self["style"] if "style" in self else {}
    
    def remove_style(self):
        self.pop("style", None)



class Molecule(list):
    def __init__(self):
        self.atoms = []
        self.bonds = []

    def __getitem__(self, idx):
        # return self.atoms[idx]
        return [x for x in self.atoms if int(x['idx']) == int(idx)][0]
    
    def __iter__(self):
        for elem in self.atoms:
            yield elem
    
    def __len__(self):
        return len(self.atoms)
    
    def select_distance(self, atom: Atom, distance):
        return [a for a in self.atoms if atom.dist(a) <= distance]
    
    def select_molecule_distance(self, other, distance):
        return [a for a in self.atoms if min([a.dist(b) for b in other.atoms]) <= distance]
    
    def get_mean(self):
        mean_x = sum([a["x"] for a in self.atoms])/len(self.atoms)
        mean_y = sum([a["y"] for a in self.atoms])/len(self.atoms)
        mean_z = sum([a["z"] for a in self.atoms])/len(self.atoms)
        return {"x": mean_x, "y": mean_y, "z": mean_z}

    def get_coords(self):
        return np.asarray([[a["x"], a["y"], a["z"]] for a in self.atoms])
    
    def get_elements(self):
        return np.asarray([a["element"] for a in self.atoms])
    
    def rmsd(self, other):
        coords1 = self.get_coords()
        coords2 = other.get_coords()
        elements1 = self.get_elements()
        elements2 = other.get_elements()
        try:
            suppl1 = Chem.SDMolSupplier()
            suppl1.SetData(str(self))
            am1 = Chem.rdmolops.GetAdjacencyMatrix(suppl1[0])
            suppl2 = Chem.SDMolSupplier()
            suppl2.SetData(str(other))
            am2 = Chem.rdmolops.GetAdjacencyMatrix(suppl2[0])

            cur_rmsd = symmrmsd(coords1, coords2, elements1, elements2, am1, am2)[0]
        except Exception as e:  
            return e
    
        return cur_rmsd
    
    @staticmethod
    def from_pdb(pdb):
        mol = Molecule()
        for line in pdb.split("\n"):
            if line.startswith("ATOM"):
                atom = Atom()
                atom["idx"] = int(line[6:11].strip())
    