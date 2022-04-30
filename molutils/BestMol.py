"""
author:
    Parker Sornberger
email:
    pdso223@uky.edu
"""

import rdkit.Chem as Chem
import rdkit
from collections import defaultdict
import selfies as sf
import rdkit.Chem.rdMolDescriptors as Rdesc
import rdkit.Chem.Descriptors as desc
import rdkit.Chem.Descriptors as Desc
import rdkit.Chem.Scaffolds.MurckoScaffold as MS
import warnings
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from tqdm import tqdm
from pprint import pprint
from copy import deepcopy
import math
from rdkit.Chem.Draw import IPythonConsole
from itertools import chain
from types import MappingProxyType
from PIL import Image
from collections.abc import Iterable
import networkx as nx

class BestMol(rdkit.Chem.rdchem.Mol):
    __version__="1.2.6.1"
    def __init__(self, smile, name = None):
        temp = Chem.MolFromSmiles(smile)
        better_str = Chem.MolToSmiles(temp)
        self.mol = Chem.MolFromSmiles(better_str)
        self.smiles = better_str
        self.name = name
        super(BestMol, self).__init__()

    def __len__(self):
        return self.mol.GetNumAtoms()
    
    def __iter__(self):
        return self.mol.GetAtoms()
    
    def __repr__(self):
        return F"BestMol {self.smiles}"
    def __copy__(self):
     
        return BestMol(self.smiles)
    
    def __getitem__(self, i):
        if isinstance(i, int):
            return self.mol.GetAtomWithIdx(i)
        elif isinstance(i, str):
            return getattr(self, i)
        elif isinstance(i, Iterable):
            return [self[thing] for thing in i]

        
    def __hash__(self):
        return hash(F"{self.smiles}+BestMolObject")
    
    def __eq__(self, other):
        if not isinstance(other, BestMol):
            print(F"You're trying to compare {type(other).__name__} to {type(self).__name__}!")
            return False
        self_smile = Chem.MolToSmiles(self.mol, kekuleSmiles=True)
        other_smile = Chem.MolToSmiles(other.mol, kekuleSmiles=True)
        return self_smile==other_smile

    
    def __str__(self):
        return F" BestMol object that has isomeric SMILES: {Chem.MolToSmiles(self.mol)}"
    
    def copy(self):
        return self.__copy__()
    
    

    def ring_tuple_as_symbols(self, method="LGFR"):
        def failure():
            print("Not a valid choice! \nWe have All Rings, LGFR, and LGCR \nGiving you all rings")
            return "All Rings"
        choices = defaultdict(failure)
        the_choices = {"All Rings": self.ring_tuple, "LGFR": self.lgfr_atom_groups}
        choices.update(the_choices)
        actions = choices[method]
        

    def make_own_core_generic(self, return_nonC=False, native_S=False):
        non_C_dict= {}
        valid = [index for ring in self.lgfr_atom_groups for index in ring]
        for i, atom in enumerate(self):
            if i not in valid:
                continue
            symbol = atom.GetSymbol()
            if symbol !="C":
                if symbol=="S" and native_S:
                    continue
                else:
                    non_C_dict[i]=symbol
            atom.SetAtomicNum(0)
            atom.SetIsAromatic(True)
            atom.SetFormalCharge(0)
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            atom.SetNoImplicit(0)
            atom.SetNumExplicitHs(0)

        for bond,b_tup in zip(self.bond_iterator, self.bond_tuples):
            beg, end = b_tup
            if beg in valid and end in valid:
                bond.SetBondType(Chem.BondType.AROMATIC)
                bond.SetIsAromatic(True)
        if return_nonC:
            return non_C_dict

    def get_generic_core_BM(self):
        gen = self.get_LGFR_BM()
        for atom in gen:
            atom.SetAtomicNum(0)
            atom.SetIsAromatic(True)
            atom.SetFormalCharge(0)
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            atom.SetNoImplicit(0)
            atom.SetNumExplicitHs(0)
        for bond,b_tup in zip(gen.bond_iterator, gen.bond_tuples):


            bond.SetBondType(Chem.BondType.AROMATIC)
            bond.SetIsAromatic(True)
        return gen

    def get_carbon_generic_LGFR(self):
        lgfr_mol=self.get_LGFR_BM().mol
        gen_lgfr =MS.MakeScaffoldGeneric(lgfr_mol)
        return BestMol(Chem.MolToSmiles(gen_lgfr))
        
        
    def return_generic_core_BM(self):
        mol = BestMol(self.smiles)
        valid = [index for ring in mol.lgfr_atom_groups for index in ring]
        for i, atom in enumerate(mol):
            if i not in valid:
                continue
            atom.SetAtomicNum(0)
            atom.SetIsAromatic(True)
            atom.SetFormalCharge(0)
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            atom.SetNoImplicit(0)
            atom.SetNumExplicitHs(0)

        for bond,b_tup in zip(mol.bond_iterator, mol.bond_tuples):
            beg, end = b_tup
            if beg in valid and end in valid:
                bond.SetBondType(Chem.BondType.AROMATIC)
                bond.SetIsAromatic(True)
        return mol

    def connected_ring_BM(self):
        
        rings = self.ring_atoms
        atoms = [k for k in range(len(self))]

        not_ring_atoms= [atom  for atom in atoms if atom not in rings]
        _m =deepcopy(self.mol)

        for i in not_ring_atoms:
            _m.GetAtomWithIdx(i).SetAtomicNum(0)
        mol3 = Chem.DeleteSubstructs(_m, Chem.MolFromSmarts('[#0]'))

        return BestMol(Chem.MolToSmiles(mol3)) 
    def CONNECTED_fragments(self, remove_dummies=True):
        bond_symbols={"=","%","#", "@"}
        if self.all_rings==0:
            #warnings.warn("This strucuture has no rings! -1 is returned")
            return 0
        chains = Chem.ReplaceCore(self.mol, self.connected_ring_BM().mol)
        
        as_frags = Chem.GetMolFrags(chains, asMols=True, sanitizeFrags = True)
        if remove_dummies:
            smiles = [Chem.MolToSmiles(frag) for frag in as_frags]
            cleansmiles = []
            for smile in smiles:
                
                if "*]" in smile:

                    where = smile.index("*]")
                    smile = smile[where+2:]
                    if smile[0] in bond_symbols:
                        smile = smile[1:]
                    cleansmiles.append(smile)
            
            cleaned = [Chem.MolFromSmiles(clean) for clean in cleansmiles]
            
            return cleaned 
        return as_frags
        
    def get_LGFR_BM(self):
        lgfr = self.lgfr_atom_groups
        rings = self.ring_tuple
        atoms = [k for k in range(len(self))]
        lg_list = [g for t in lgfr for g in t]
        tog = [atom  for atom in atoms if atom not in lg_list]
        _m =deepcopy(self.mol)

        for i in tog:
            _m.GetAtomWithIdx(i).SetAtomicNum(0)
        mol3 = Chem.DeleteSubstructs(_m, Chem.MolFromSmarts('[#0]'))
        bm = BestMol(Chem.MolToSmiles(mol3))
        bm.name = self.name
        return bm
    def LGFR_fragments(self, remove_dummies=True):
        bond_symbols={"=","%","#", "@"}
        if self.all_rings==0:
            #warnings.warn("This strucuture has no rings! -1 is returned")
            return 0
        chains = Chem.ReplaceCore(self.mol, self.get_LGFR_BM().mol)
        
        as_frags = Chem.GetMolFrags(chains, asMols=True, sanitizeFrags = True)
        if remove_dummies:
            smiles = [Chem.MolToSmiles(frag) for frag in as_frags]
            cleansmiles = []
            for smile in smiles:
                
                if "*]" in smile:

                    where = smile.index("*]")
                    smile = smile[where+2:]
                    if smile[0] in bond_symbols:
                        smile = smile[1:]
                    cleansmiles.append(smile)
            
            cleaned = [Chem.MolFromSmiles(clean) for clean in cleansmiles]
            
            return cleaned 
        return as_frags
    
    

    def get_LGFR_BM_aromatic(self):
        lgfr = self.lgfr_atom_groups
        rings = self.ring_tuple
        atoms = [k for k in range(len(self))]
        lg_list = [g for t in lgfr for g in t]
        tog = [atom  for atom in atoms if atom not in lg_list]
        _m =deepcopy(self.mol)
        for i in tog:
            _m.GetAtomWithIdx(i).SetAtomicNum(0)
        mol = Chem.DeleteSubstructs(_m, Chem.MolFromSmarts('[#0]'))
        mol = BestMol(Chem.MolToSmiles(mol))
        for bond,b_tup in zip(mol.bond_iterator, mol.bond_tuples):
            beg, end = b_tup
            if beg in lg_list and end in lg_list:
                bond.SetBondType(Chem.BondType.AROMATIC)
                bond.SetIsAromatic(True)
        return mol
    
    def add_C_to_get_para(self, inplace=False):
        """
        doesn't really work <3<3<3.

        """
        
        smiles = self.smiles
        colors = len(self.canonical_rank_atom_unique)
        last = first="C"
        after = "(C)"
        if not inplace:
            for i in range(len(self.smiles)):
                if i ==0:
                    try:
                        new=first+self.smiles
                        newbm = BestMol(new)
                        newcolors = len(newbm.canonical_rank_atom_unique)
                        if newcolors == (colors+1):
                            return newbm
                        else:
                            continue
                    except:
                        continue
                if i == (len(smiles)-1):
                    try:
                        new=self.smiles+last
                        newbm = BestMol(new)
                        newcolors = len(newbm.canonical_rank_atom_unique)
                        if newcolors == (colors+1):
                            return newbm
                        else:
                            continue
                    except:
                        continue
                else:
                    try:
                        new=smiles[:i]+after+smiles[i:]
                        newbm = BestMol(new)
                        newcolors = len(newbm.canonical_rank_atom_unique)
                        if newcolors == (colors+1):
                            return newbm
                        else:
                            continue
                    except:
                        continue
        else:
            print("implement later...")
            return None
        
    
    
    
    
    def make_self_bemis_mucko(self):
        self.mol = self.murcko_scaffold
    
    
    def add_hydrogens_to_self(self):
        self.mol = Chem.AddHs(self.mol)
    
    def has(self, other):
        try:
            does = self.mol.HasSubstructMatch(other)
            return does
        except:
            print("oops :(")
            return None
        
    def where_sub_at(self, other_as_mol):
        return self.mol.GetSubstructMatches(other_as_mol)
    
    @property
    def graph(self):
        graph = nx.Graph()
        for i, atom in self.element_dictionary.items():
            graph.add_node(i, symbol=atom)
        for (begin, end), bond in zip(self.bond_tuples, self.bond_list):
            graph.add_edge(begin, end, order = str(bond))
        return graph
            
    
    
    
    @property
    def atom_id_dict(self):
        return MappingProxyType({i:atom for (i, atom) in enumerate(self)})
    
    

    @property
    def selfies(self):
        return sf.encoder(self.smiles)

        
    @property
    def ideal_conf_num(self):
        conf_number = 0
        rotatable = int(self.rotatable_bond_count)
        if rotatable<3:
            conf_number = 50
        elif rotatable >6:
            conf_number = 300
        else:
            conf_number = rotatable**3
        return conf_number

    
    @property
    def weight(self):
        return desc.ExactMolWt(self.mol)

    @property
    def all_rings(self):
        return Rdesc.CalcNumRings(self.mol)


    @property
    def aromatic_rings(self):
        return Rdesc.CalcNumAromaticRings(self.mol)
    @property
    def murcko_scaffold(self):
        if self.all_rings==0:
            warnings.warn("This strucuture has no rings! NoneType is returned")
            return None
        return MS.GetScaffoldForMol(self.mol)
    @property
    def fragment_count(self):
        if self.all_rings==0:
            #warnings.warn("This strucuture has no rings! -1 is returned")
            return 0
        chains = Chem.ReplaceCore(self.mol, self.murcko_scaffold)
        
        as_frags = Chem.GetMolFrags(chains, asMols=True, sanitizeFrags = True)
        return len(as_frags)

    @property
    def number_of_heteroatoms(self):
        return Rdesc.CalcNumHeteroatoms(self.mol)
    @property
    def RingInfoObject(self):
        if self.all_rings == 0:
            print("No rings here! Returning NoneType")
            return None
        return self.mol.GetRingInfo()
    @property
    def ring_tuple(self):
        if self.all_rings == 0:
            print("No rings here! Returning empty list")
            return []
        return self.RingInfoObject.AtomRings()
    @property
    def ring_atoms(self):
        return frozenset(chain.from_iterable(self.ring_tuple))
    @property
    def ring_size_list(self):
        if self.all_rings == 0:
            print("No rings here! Returning empty list")
            return []
        return tuple([len(ring) for ring in self.ring_tuple])
    @property
    def unique_ring_size_list(self):
        if self.all_rings == 0:
            print("No rings here! Returning empty list")
            return []
        return frozenset(set(len(ring) for ring in self.ring_tuple))  
    
    @property
    def all_connected_rings(self):
        rings = self.ring_tuple
        size_list = self.ring_size_list
        tot_length = sum(size_list)
        unique = self.count_all_unique_ring_atoms
        if tot_length==unique:
            return True
        else:
            return False   
    @property
    def group_fused(self):
        rings = self.ring_tuple
        r_d = {i:set(ring) for (i,ring) in enumerate(rings)}
        fused_to = {i: [] for i in range(len(rings))}
        # if all_connected_rings(Bmol):
        #     return fused_to
        for index, ring in r_d.items():
            for sub_index, sub_ring in r_d.items():
                if index==sub_index:
                    continue
                temp = ring | sub_ring
                together = len(temp)
                total=len(sub_ring)+len(ring)
    
                if together != total:
    
                    fused_to[index].append(sub_index)
        return fused_to
    @property
    def chain_fused(self):
        visited = []
        chain = {}
        grouped = self.group_fused
        for index, neighbors in grouped.items():
            neighbors.append(index)
            neighbors = set(neighbors)
            copy = neighbors.copy()
            for second_place, n_list in grouped.items():
                no_visits = False
                
                for neighbor in copy:
                    if second_place==index:
                        continue
                    if neighbor in n_list:
                        neighbors |=set(n_list)
            if neighbors !=set():
                copy2=neighbors.copy()
                for n in copy2:
                    for _chain in chain.values():
                        if n in _chain:
                            neighbors |=_chain
                if neighbors not in chain.values():
                    chain[index]=neighbors
        return chain
    @property
    def all_fused(self):
        if self.single_bond_count !=0 or self.rotatable_bond_count !=0:
            return False
        ring_atoms = self.return_all_ring_atoms
    
        #unique = self.return_all_unique_ring_atoms
        #counts = set(ring_atoms.count(atom) for atom in unique)
        rings = self.ring_tuple
        counted = {i:sum(set(ring_atoms.count(atom) for atom in ring)) for (i, ring) in enumerate(rings) }

        if 1 in counted.values():
            return False
        else:
            return True
        # if sum(counts) !=1:
        #     return True
        # else:
        #     return False
    @property
    def lgfr(self):
        if self.all_fused:
            return frozenset(set(index for index in range(len(self.ring_tuple))))
        chained = self.chain_fused
        return frozenset(max(chained.values(), key=len))
    
    @property
    def lgfr_atom_groups(self):
        if self.all_fused:
            ring_tuple = self.ring_tuple
            lgfr = self.lgfr
            the_lgfr = [ring_tuple[group] for group in lgfr]
            return tuple(the_lgfr)
        
        ring_tuple = self.ring_tuple
        lgfr = self.lgfr
        the_lgfr = [ring_tuple[group] for group in lgfr]
        return tuple(the_lgfr)
    @property
    def atoms_in_lgfr(self):
        return frozenset(chain.from_iterable(self.lgfr_atom_groups))
        
    
    @property
    def ring_size_count_dict(self):
        if self.all_rings == 0:
            print("No rings here! Returning empty dict")
            return {}
        unique = self.unique_ring_size_list
        sizes = self.ring_size_list
        return MappingProxyType({uni:sizes.count(uni) for uni in unique})
    @property
    def count_all_ring_atoms(self):
        if self.all_rings == 0:
            print("No rings here! Returning 0")
            return 0
        return len([atom for ring in self.ring_tuple for atom in ring])
    @property
    def return_all_ring_atoms(self):
        if self.all_rings == 0:
            print("No rings here! Returning 0")
            return 0
        return tuple([atom for ring in self.ring_tuple for atom in ring])
    @property
    def count_all_unique_ring_atoms(self):
        if self.all_rings == 0:
            print("No rings here! Returning 0")
            return 0
        return len(set(atom for ring in self.ring_tuple for atom in ring))
    @property
    def return_all_unique_ring_atoms(self):
        if self.all_rings == 0:
            print("No rings here! Returning 0")
            return 0
        the_atoms = set(atom for ring in self.ring_tuple for atom in ring)
        return frozenset(the_atoms)
    @property
    def size_of_largest_ring(self):
        if self.all_rings == 0:
            print("No rings here! Returning 0")
            return 0
        try:
            info = self.mol.GetRingInfo()
            #max_ring_sizes = len(max(info.AtomRings(), key = len))
            
        
            curr_max = len(max(info.AtomRings(), key = len))
            
        
            return curr_max
        except:
            return 0
    @property
    def smallest_ring(self):
        if self.all_rings == 0:
            print("No rings here! Returning 0")
            return 0
        
        try:
            info = self.mol.GetRingInfo()
            #max_ring_sizes = len(min(info.AtomRings(), key = len))
            
        
            curr_min = len(min(info.AtomRings(), key = len))
            
        
            return curr_min
        except:
            return 0
    @property
    def radicals(self):
        return desc.NumRadicalElectrons(self.mol)
    
        
    @property
    def hetero_atoms_in_rings(self):
        if self.all_rings ==0:
            print("No rings in this structure! \nReturned -1 to show failure")
            return -1
        info = self.mol.GetRingInfo()
        atom_rings = info.AtomRings()
        
        for ring in atom_rings:
            diff = 0
            for atomid in ring:
                symbol = self.mol.GetAtomWithIdx(atomid).GetSymbol()
                if symbol.upper() != "C":
                    diff+=1

        return diff
    

    
    @property
    def rotatable_bond_count(self):
        return AllChem.CalcNumRotatableBonds(self.mol)
    
    @property
    def bond_iterator(self):
        return self.mol.GetBonds()
    @property
    def valence_electrons(self):
        return desc.NumValenceElectrons(self.mol)
    @property
    def radical_electrons(self):
        return desc.NumRadicalElectrons(self.mol)
    
    @property
    def has_symmetry(self):
        return len(set(rdmolfiles.CanonicalRankAtoms(self.mol, breakTies=False))) <=(self.__len__()//2)
    @property
    def canonical_rank_atom_literal(self):
        return tuple([canon for canon in rdmolfiles.CanonicalRankAtoms(self.mol)])
    @property
    def canonical_rank_atom_unique(self):
        return frozenset(set(rdmolfiles.CanonicalRankAtoms(self.mol,breakTies=False)))
    
    @property
    def atoms(self):
        return tuple([atom for atom in self])
    
    @property
    def elements(self):
        return frozenset(set(atom.GetSymbol() for atom in self))
    @property
    def element_dictionary(self):
        return MappingProxyType({i:atom.GetSymbol() for (i,atom) in enumerate(self)})
    
    @property
    def bonds_literal(self):
        return tuple([bond for bond in self.bond_iterator])
    @property
    def bond_list(self):
        return tuple([Chem.rdchem.Bond.GetBondType(bond) for bond in self.bonds_literal])
    @property
    def bond_tuples(self):
        literal = self.bonds_literal
        return tuple([(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in literal])

    @property
    def single_bond_count(self):

        single = rdkit.Chem.rdchem.BondType.SINGLE
        bonds = self.bond_list
        return bonds.count(single)

    @property
    def single_bonds_literal(self):

        single = rdkit.Chem.rdchem.BondType.SINGLE
        bond_type = Chem.rdchem.Bond.GetBondType


        return tuple([bond for bond in self.bonds_literal if bond_type(bond)==single])

    @property
    def single_bonds_tuples(self):
        literal = self.single_bonds_literal
        return tuple([(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in literal ])

    @property
    def double_bond_count(self):

        double = rdkit.Chem.rdchem.BondType.DOUBLE
        bonds = self.bond_list
        return bonds.count(double)
    @property
    def double_bonds_literal(self):
        double = rdkit.Chem.rdchem.BondType.DOUBLE
        bond_type = Chem.rdchem.Bond.GetBondType
        return tuple([bond for bond in self.bonds_literal if bond_type(bond)==double])
    @property
    def double_bonds_tuples(self):
        literal = self.double_bonds_literal
        return tuple([(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in literal ])
    @property
    def triple_bond_count(self):

        triple = rdkit.Chem.rdchem.BondType.TRIPLE
        bonds = self.bond_list
        return bonds.count(triple)
    @property
    def triple_bonds_literal(self):
        triple = rdkit.Chem.rdchem.BondType.TRIPLE
        bond_type = Chem.rdchem.Bond.GetBondType
        return tuple([bond for bond in self.bonds_literal if bond_type(bond)==triple])
    @property
    def triple_bonds_tuples(self):
        literal = self.triple_bonds_literal
        return tuple([(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in literal ])

    @property
    def aromatic_bond_count(self):

        single = rdkit.Chem.rdchem.BondType.AROMATIC
        bonds = self.bond_list
        return bonds.count(single)   
    @property
    def aromatic_bonds_literal(self):
        aromatic = rdkit.Chem.rdchem.BondType.AROMATIC
        bond_type = Chem.rdchem.Bond.GetBondType
        return tuple([bond for bond in self.bonds_literal if bond_type(bond)==aromatic])

    @property
    def aromatic_bonds_tuples(self):
        literal = self.aromatic_bonds_literal
        return tuple([(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in literal ])
    
    @property
    def symbols(self):
        return frozenset(set(self.smiles))
    
    @property
    def simple_lowest_energy_conf(self):
        xyz_func = rdkit.Chem.rdmolfiles.MolToXYZFile
        list_of_energies = []
        get_properties = AllChem.MMFFGetMoleculeProperties
        append_energy = list_of_energies.append
        obtain_conf_id = AllChem.EmbedMultipleConfs
        with_h = Chem.AddHs(self.mol)
        molecule_parameters = AllChem.ETKDG
        make_force = AllChem.MMFFGetMoleculeForceField
        _ids = obtain_conf_id(with_h, numConfs = self.ideal_conf_num, params = molecule_parameters() )
        for _id in tqdm(_ids,desc="Making conformers..."):
            properties = get_properties(with_h, mmffVariant = "MMFF94s")
            force_field = make_force(with_h, properties, confId=_id)
            force_field.Minimize()
            energy = float(force_field.CalcEnergy())
            append_energy((energy, with_h, _id))
        lowest = min(list_of_energies[:])[0]
        return lowest
    def summary(self, pretty:bool=True, file_dump:str=None, others:list=[]):
        selected_attributes = {"Name":"name", "Molecule Smiles": "smiles", "Number of Rings": "all_rings",
                               "LGFR Smiles": "get_LGFR_BM", "LGCR Smiles": "connected_ring_BM",
                               "LGFR Fragement Count" : "LGFR_fragments", "LGCR Fragment Count": "CONNECTED_fragments",
                               }
        
        
        if file_dump:
            
            with open(file_dump, "w") as file:
                def _write_factory(m, selected_attributes, dump=None):
                    
                        for name, attr in selected_attributes.items():
                            if "FRAG" in attr.upper():
                                frags = len(m[attr]())
                                file.write(F"{name} : {frags}\n")
                            elif "BM" in attr:
                                smi = Chem.MolToSmiles(m[attr]().mol)
                                file.write(F"{name} : {smi}\n")
                            else:
                                file.write(F"{name} : {m[attr]}\n")
                if not others:
                    
                    _write_factory(self, selected_attributes, file_dump)
                else:
                    
                    for other in [self]+others:
                        
                        _write_factory(other, selected_attributes, file_dump)
        if pretty:
            def _dict_factory():
                pass
        
            
        
    
    @property
    def flat_score(self):
        sym=False
        best=False
        # rings = self.all_rings
        # if rings >16 or rings <=1:
        #     dummy=1
            
        # else:
        #     bad_mol=True
        num_atoms = len(self)
        if self.all_fused:
            best=True

        if self.has_symmetry:
            sym = True
        if best and sym:
            return 3
        if best:
            return 2
        if sym:
            return 1
        else:
            return 0
                          
        
    
        
    
    @property
    def score(self):
        rings = self.all_rings
        
        if rings >5 or rings <=1:
            ring_score = -100
        else:
            ring_score = 1
        rotate = self.rotatable_bond_count
        if rotate>4:
            rot_score = 1/rotate
        else:
            rot_score = 1
        # if rings==len(self.lgfr):
        #     ring_score =ring_score+1
        
        if self.has_symmetry:
            sym_score = 2
        else:
            sym_score = -1
        if self.hetero_atoms_in_rings>1: 
            het_score =1
        else:
            het_score =0
        return rot_score + sym_score+ring_score+het_score
        
    def filter_self(self):
        bad = False
        if self.size_of_largest_ring>9 or self.smallest_ring<5:
            bad = True
        do_not_use = {"S", "Se", "Te"}
        intersection = self.elements & do_not_use
        if intersection:
            bad = True
        if self.radical_electrons:
            bad = True
        if bad:
            return 0
        else:
            return 1
        
        
        
    def to_sdf(self):
        pass
    def to_image(self, name="BM", use_png=True, dpi=300, size=(400,400), transparent =False):
        if use_png:
            Chem.Draw.MolToFile(self.mol, filename=F"{name}.png",dpi=dpi, size = size )
            if transparent:
                img = Image.open(F"{name}.png")
                rgba = img.convert("RGBA")
                datas = rgba.getdata()
                newData = []
                for item in datas:
                    if item[0] == 255 and item[1] == 255 and item[2] == 255:  # finding yellow colour
                # replacing it with a transparent value
                        newData.append((255, 255, 255, 0))
                    else:
                        newData.append(item)
                rgba.putdata(newData)
                rgba.save(F"{name}.png", "PNG")
        else:
            print("I am lazt af")
    
    def to_xyz(self, file_name="Your Awesome Molecule", clip=False,
               print_lowest_too=False, just_sdf=False, xyz_and_sdf=False):
        xyz_func = rdkit.Chem.rdmolfiles.MolToXYZFile
        list_of_energies = []
        get_properties = AllChem.MMFFGetMoleculeProperties
        append_energy = list_of_energies.append
        obtain_conf_id = AllChem.EmbedMultipleConfs
        with_h = Chem.AddHs(self.mol)
        molecule_parameters = AllChem.ETKDG
        make_force = AllChem.MMFFGetMoleculeForceField
        if not clip:
            confs = self.ideal_conf_num
        else:
            clip = 30
        _ids = obtain_conf_id(with_h, numConfs = confs, params = molecule_parameters() )
        for _id in tqdm(_ids, desc="Making conformers..."):
            properties = get_properties(with_h, mmffVariant = "MMFF94s")
            force_field = make_force(with_h, properties, confId=_id)
            force_field.Minimize()
            energy = float(force_field.CalcEnergy())
            append_energy((energy, with_h, _id))
        lowest = min(list_of_energies[:])
        if print_lowest_too:
            print(F"{file_name} has energy {lowest[0]}")
        its_name = F"{file_name}.xyz"
        the_good = lowest[1]
        its_id = lowest[2]
        xyz_func(the_good, its_name, confId = its_id)
    
     

class BestMolRL(BestMol):
    pass    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
