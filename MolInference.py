"""
author:
    Parker Sornberger
email:
    pdso223@uky.edu
"""

from BestMol import BestMol
#from ocelot.schema.graph import MolGraph


class MolInference:
    __version__ = "1.0.0.1"
    def __init__(self, string, name = None, add_h=False):
        
        self.bm = BestMol(string, name)
        if add_h:
            self.bm.add_hydrogens_to_self()
        self.smiles = string
        self.fused_scaffold = self.bm.get_LGFR_BM()
        self.generic = self.bm.get_carbon_generic_LGFR()
        self.generic.name = name
        
        self.fragments_lgfr = self.bm.LGFR_fragments()
        
        
        self.connected_scaffold = self.bm.connected_ring_BM()
        self.ring_atoms = self.bm.ring_atoms
        
        self.mol = self.bm.mol
        self.lgfr_groups = self.bm.lgfr_atom_groups
        
        self.lgfr_atoms = self.bm.atoms_in_lgfr
        
        self.graph = self.bm.graph
        self.generic_graph = self.generic.graph
        self.name = name
        
        
    def __len__(self):
        return self.mol.GetNumAtoms()
    def __iter__(self):
        return self.mol.GetAtoms()
    def __str__(self):
        return F"Graphical anlysis for {self.smiles}"
    def __hash__(self):
        return hash(F"{self.smiles} Information")
    def __repr__(self):
        return F"Info for {self.smiles}"
    def __getitem__(self, name):
        if isinstance(name, int):
            return self.mol.GetAtomWithIdx(name)
        return getattr(self, name)
    @classmethod
    def from_bmol(cls,bm):
        return cls(bm.smiles, bm.name)




class ForCombo:
    __version__ = "1.0.0.1"
    def __init__(self, string, name = None, add_h=False):
        
        self.bm = BestMol(string, name)
        if add_h:
            self.bm.add_hydrogens_to_self()
        self.smiles = string

        
        self.mol = self.bm.mol

        
        self.graph = MolGraph.from_rdmol(self.mol).graph

        
    def __len__(self):
        return self.mol.GetNumAtoms()
    def __iter__(self):
        return self.mol.GetAtoms()
    def __str__(self):
        return F"Graphical anlysis for {self.smiles}"
    def __hash__(self):
        return hash(F"{self.smiles} Information")
    def __repr__(self):
        return F"Info for {self.smiles}"
    def __getitem__(self, name):
        if isinstance(name, int):
            return self.mol.GetAtomWithIdx(name)
        return getattr(self, name)
    @classmethod
    def from_bmol(cls,bm):
        return cls(bm.smiles, bm.name)








