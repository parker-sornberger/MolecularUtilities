"""
author: Parker Sornberger
email: pdso223@uky.edu OR parkersornberger3@gmail.com 

pysmiles info:
    https://github.com/pckroon/pysmiles 
    It's a really awesome package

"""



from .Molecule import Molecule
#from pysmiles.pysmiles import read_smiles, write_smiles
import networkx as nx
import warnings
from types import MappingProxyType
from tqdm import tqdm

import matplotlib.pyplot as plt
import rdkit.Chem as Chem
import rdkit


class StaticMolGroupers:
    @staticmethod
    def group_by_lgfr(molecules):
        lg_dict = {}
        for mol in tqdm(molecules, desc = "getting LGFR groups..."):
            try:
                bm = Molecule(mol)
            except:
                continue
            else:
                lgfr = bm.get_LGFR_BM()
                if lgfr not in lg_dict:
                    lg_dict[lgfr] = [bm]
                else:
                    lg_dict[lgfr].append(bm)
        return lg_dict







class MolGrouper:
    def __init__(self, molecules, mol_names = None, method = "lgfr", substruct=None, 
                 group_singletons = True, generic = False, start_now=True):
        self.molecules = molecules
        
        self.allowed_methods = {"LGFR", "LGCR", "BEMIS"}
        if method.upper() not in self.allowed_methods and not substruct:
            warnings.warn(F"Howdy! Allowable methods are {self.allowed_methods} \nChoose on of these or search by substructure <3")
            print("Method not found... choosing LGFR")
            method = "lgfr"
        self.method = method
        if substruct:
            if isinstance(substruct, str):
                substruct = [substruct]
            self.method = "substruct"
        self.subs_to_search = substruct
        self.mol_names = mol_names
        self.failures = set()
        self.use_generic = generic
        self.group_singletons=group_singletons
        if start_now:
            self.started = True
            self._begin_groups()
        else:
            self._groups= {}
            self.started = False
    
    def __repr__(self):
        return F"{self.method} MolGrouper at {hex(id(self))} started is {self.started}"
    def __str__(self):
        return F"{self.method} MolGrouper with {len(self._groups)} groups and started is {self.started}"
    
    
    def __iter__(self):
        return iter(self._groups)
    def keys(self):
        return self._groups.keys()
    def values(self):
        return self._groups.values()
    def items(self):
        return self._groups.items()
    
    def __copy__(self):
        return self._groups.__copy__()
    def __len__(self):
        return len(self._groups)
    
    def copy(self):
        return self._groups.copy()
    
    
    def _by_lgfr(self, other=None):
        if other:
            molecules = other
        else:
            molecules = self.molecules
        lg_dict = {}
        for mol in tqdm(molecules, desc = "getting LGFR groups..."):
            try:
                bm = Molecule(mol)
            except:
                self.failures.add(mol)
                continue
            if self.use_generic:
                pass
            else:
                lgfr = bm.get_LGFR_BM()
                if lgfr not in lg_dict:
                    lg_dict[lgfr] = [bm]
                else:
                    lg_dict[lgfr].append(bm)
        return lg_dict
                    
        
        
        
    
    def _by_lgcr(self, other=None):
        lc_dict = {}
        if other:
            molecules = other
        else:
            molecules = self.molecules
        for mol in tqdm(molecules, desc = "getting LGFR groups..."):
            try:
                bm = Molecule(mol)
            except:
                self.failures.add(mol)
                continue
            if self.use_generic:
                pass
            else:
                lgcr = bm.connected_ring_BM()
                if lgcr not in lc_dict:
                    lc_dict[lgcr] = [bm]
                else:
                    lc_dict[lgcr].append(bm)
        return lc_dict
    
    def _by_bemis(self):
        pass
    
    
    def _by_substruct(self, other=None):
        
        subs = {b:Molecule(b) for b in self.subs_to_search}
        has_these = {"None": []}
        if other:
            molecules = other
        else:
            molecules = self.molecules 
        for mol in tqdm(molecules):
            try:
                bm = Molecule(mol)
            except:
                continue
            found_once = False
            for sub, submol in subs.items():
                if bm.has(submol.mol):
                    found_once = True
                    if submol not in has_these.keys():
                        has_these[submol] = [bm]
                    else:
                        has_these[submol].append(bm)
            if not found_once:
                has_these["None"].append(bm)
        if not len(has_these["None"]):
            has_these.pop("None")
        return has_these
        
    def _wrap_singletons(self):
        if self.subs_to_search:
            return 
        temp_groups = {"One Element Groups":[]}
        
        for key, group in self._groups.items():
            if len(group) <2:
                temp_groups["One Element Groups"].extend(group)
            else:
                temp_groups[key]=group
        self._groups=temp_groups
        

    def _begin_groups(self):
        if self.method.upper() == "LGFR":
            self._groups = self._by_lgfr()
        if self.method.upper() == "LGCR":
            self._groups = self._by_lgfr()
        if self.method =="substruct":
            self._groups = self._by_substruct()
        if self.group_singletons:
            self._wrap_singletons()
        

    @property
    def groups(self):
        if not self.started:
            print("starting...")
            self._begin_groups()
            self.started=True
        return MappingProxyType(self._groups)
    
    @property
    def largest_group(self):
        key, group= max(self.groups.items(), key=len)
        return MappingProxyType({key:group})
        #return MappingProxyType({key:val for (key,val) in max(self.groups.items(), key=len)})
    
    def update_groups(self, more_mols):
        if not self.started:
            self._begin_groups()
        if self.method.upper() == "LGFR":
            new = self._by_lgfr(other = more_mols)
        if self.method.upper() == "LGCR":
            new = self._by_lgfr(other = more_mols)
        if self.method =="substruct":
            new = self._by_substruct(other = more_mols)
        
        
        for key, group in new.items():
            if key not in self._groups:
                self._groups[key] = group
            else:
                
                existing_group = self._groups[key]
                for has in existing_group:
                    if has not in group:
                        group.append(has)
                self._groups[key] = group
                
        if self.group_singletons:
            self._wrap_singletons()
            
    def regroup(self, method, substruct = None):
        if method.upper() not in self.allowed_methods and not substruct:
            warnings.warn(F"Howdy! Allowable methods are {self.allowed_methods} \nChoose on of these or search by substructure <3")
            print("Method not found... choosing LGFR")
            method = "lgfr"
        self.method = method
        if substruct:
            if isinstance(substruct, str):
                substruct = [substruct]
            self.method = "substruct"
        self._begin_groups()
        
    def reverse(self):
        if not self.started:
            self._begin_groups()
            self.started = True
        self._groups = {val:key for (key, _val) in self._groups.items() for val in _val}
        
        

class NaivePolymerBuilder:
    def __init__(self, monomer, units, index = 0, return_as_smiles=True, get_XYZ=False, name = None):
        self.monomer = monomer
        self.units = units
        self.as_smiles = return_as_smiles
        self.write_xyz = get_XYZ        
        self.index = index 
        self.name =name
        self._build()
        if self.write_xyz:
            self._write()
    
    def plot(self, g):
        nx.draw(g)
        plt.show()
    def _conv_to_graph(self):
        self.mono_graph = Molecule(self.monomer).graph
    
    
    @staticmethod
    def greedy_LGFR_clockwise_traversal(lgfr, start, graph):
        paths = []
        left_neighbor = list(graph.neighbors(start))[0]
        if left_neighbor not in lgfr:
            left_neighbor = list(graph.neighbors(start))[-1]
        path_generator = nx.all_simple_paths(graph, source = start, target=left_neighbor)
        for path in path_generator:
            print(path)
            if len(path) != len(lgfr):
                continue
            temp_set = set(path)
            print(path)
            intersection = lgfr & temp_set
            if intersection == lgfr:
                print("asdasdasdas")
                return path
            
    
      
    def clockwise_traversal_about_LGFR(self,bm):
        graph = bm.graph
        first_lgfr_atom = bm.ring_tuple[0][0]
        
        lgfr_atoms = bm.ring_atoms
        starting_neighbors = list(graph.neighbors(first_lgfr_atom))
        right = starting_neighbors[1]
        
        little_path = [first_lgfr_atom]
    
        stop =len(bm)+5
        place = 0
        
        while right != first_lgfr_atom:
            little_path.append(right)
            try:
                temp_right = list(graph.neighbors(right))[1]
                if temp_right not in lgfr_atoms and len(list(graph.neighbors(right)))==3:
                    right = list(graph.neighbors(right))[-1]
                else:
                    right = temp_right
            except IndexError:
                right = list(graph.neighbors(right))[0]
            place +=1
            if place==stop:
                
                print("You fuckup")
                return self.greedy_LGFR_clockwise_traversal(lgfr_atoms, first_lgfr_atom, graph)
        if len(little_path) < len(lgfr_atoms):
            print("alternative branch")
            return self.greedy_LGFR_clockwise_traversal(lgfr_atoms, first_lgfr_atom, graph)
        return little_path


    @staticmethod
    def find_one_neighbor(graph):

        remove_bond = None
        for bond in graph.edges():
            if graph.edges[bond]["order"]=="AROMATIC":
                start, end = bond
                sneigh = list(graph.neighbors(start))
                s_neigh_l = [list(graph.neighbors(s)) for s in sneigh]
                if len(max(s_neigh_l, key=len))>2:
                    continue
                eneigh = list(graph.neighbors(end))
                e_neigh_l = [list(graph.neighbors(e)) for e in eneigh]
                if len(max(e_neigh_l, key=len))>2:
                    continue
                if len(sneigh)==2 and len(eneigh)==2:
                    remove_bond = bond
            if remove_bond:
                break
        return remove_bond
                

    def ensure_linear(self, b, f, where):
        
        temp = Molecule(self._to_build(b))
        clock = self.clockwise_traversal_about_LGFR(temp)
        print(clock)
        
        
    def _join_units(self):
        self._conv_to_graph()
        build_on = self.mono_graph.copy()
        for unit in range(self.units):
            b_len = len(build_on)
            span = range(b_len, b_len+1+len(self.mono_graph))
        
            newmap = {node:i for (node, i) in zip(self.mono_graph.nodes(), span)}
            old_map = {i:node for (node, i) in newmap.items()}
            #print(len(newmap)==len(frag.nodes()))
            
            end = newmap[self.index]
            adjusted = nx.relabel_nodes(self.mono_graph, newmap)
            find = self.find_one_neighbor(build_on)
            
            #self.ensure_linear(build_on, self.mono_graph, find)
            if not unit%2:
                find = find[0]
            else:
                find = find[0]
            comp = nx.compose(build_on, adjusted)
            #print(build_on.nodes[0])
            
            comp.add_edge(find,end)
           
            comp.edges[(find, end)]["order"] = "SINGLE"
            build_on = comp
            #self.plot(build_on)
            #print(build_on.nodes)
        return build_on
    @staticmethod
    def _to_build(dirty):
        bt_dict = {"SINGLE": rdkit.Chem.rdchem.BondType.SINGLE,
                   "DOUBLE" : rdkit.Chem.rdchem.BondType.DOUBLE,
                   "TRIPLE" : rdkit.Chem.rdchem.BondType.TRIPLE,
                   "AROMATIC" : rdkit.Chem.rdchem.BondType.AROMATIC}
        rwmol = Chem.RWMol()
        for node in dirty.nodes:
            curr_sym = dirty.nodes[node]["symbol"]
            atom = Chem.Atom(curr_sym)
            rwmol.AddAtom(atom)
            
        for bond in dirty.edges:
            
            bt = dirty.edges[bond]["order"]
            
            bt = bt_dict[bt]
            beg, end = bond
            rwmol.AddBond(beg, end, bt)
        return Chem.MolToSmiles(rwmol)


    def _validate_and_remove(self, dirty):
        self._poly = []

        try:
            dirty = self._to_build(dirty)
            bm = Molecule(dirty)
            if self.as_smiles:
                self._poly = bm.smiles
            else:
                self._poly = bm
        except:
            self._poly = None
    def _build(self):
        joined = self._join_units()
        
        self._validate_and_remove(joined)
    
    
    def _write(self):
        if not self.name:
            name = "polymer"
        else:
            name = self.name
        if self.as_smiles:
           mol = Molecule(self.polymer)
           mol.to_xyz(filename=name)
           
        else:

            self.polymer.to_xyz(file_name=name)    
    @property
    def polymer(self):
        return self._poly

        
class MolJoiner:
    def __init__(self, molecules, fragments, multi_frag, start_indices=None, eager = True, iterative=False,
                 frag_indices =None, infer_best=False, times_to_join = None, 
                 pick_indicies_manually = False, return_as_smiles = True):
        self.molecules = molecules
        self.mol_graphs = []
        
        self.fragments = fragments
        self.frag_graphs = []
        self.begin_at = start_indices
        self.end_at = frag_indices
        self.multi_frag = multi_frag
        
        if pick_indicies_manually:
            self._conv_to_graph()
            self._plot()
        
        
        if not self.begin_at and not self.end_at:
            print("No indices passed... Using iterative combination :)")
            self.iterative = True
        
        if infer_best:
            print("Attempting symmetric joining")
        self.infer_best = infer_best
        self.times_to_join = times_to_join
        self.eager_execution = eager
        self.as_smiles = return_as_smiles
        if eager:
            self._begin_join()
            self.started = True
        else:
            self.started = False
        
        
    def _plot(self):
        
        print("Input as 'start end' for scaffolds and 'start' for fragments")
        self.begin_at=[]
        self.end_at = []
        for g in self.mol_graphs:
            #labels = nx.get_node_attributes(g, "element")
            labels = {node:(node, g.nodes[node]["symbol"] ) for node in g.nodes()}
            nx.draw(g, with_labels=True, labels = labels)
            plt.show()
            start, end = input("Start (space) end: ").split()
            self.begin_at.append((start, end))
        for g in self.frag_graphs:
            #labels = nx.get_node_attributes(g, "element")
            labels = {node:(node, g.nodes[node]["symbol"] ) for node in g.nodes()}
            nx.draw(g, with_labels=True, labels = labels)
            plt.show()
            start= input("Start: ")
            self.end_at.append(start)
        
    def _infer_attempt(self):
        pass
    
    
    def _conv_to_graph(self):
        self.mol_graphs = [Molecule(smile).graph for smile in tqdm(self.molecules, desc = "Creating graphs...")]
        self.frag_graphs = [Molecule(smile).graph for smile in tqdm(self.fragments, desc = "Creating graphs...")]
        
   

    @staticmethod
    def _two_connect_graphs(g1, lg2, starts, ends):
        g1_len = len(g1)
        
        tot = g1_len
        
        comp = g1
        #print(len(g1), len(lg2), len(starts), len(ends))
        for start, end, frag in zip(starts, ends, lg2):
            #print(len(list(frag.neighbors(end))))
            span = range(tot, tot+1+len(frag))
            newmap = {node:i for (node, i) in zip(frag.nodes(), span)}
            #print(len(newmap)==len(frag.nodes()))
            
            end = newmap[end]
            adjusted = nx.relabel_nodes(frag, newmap)
            
            comp = nx.compose(comp, adjusted)
    
            # comp.nodes[start]["hcount"]=0
            # if comp.nodes[end]["aromatic"]:
            #     #print(444444)
            #     comp.nodes[end]["hcount"]=0
            # elif comp.nodes[end]["hcount"]>0:
            #     comp.nodes[end]["hcount"]=comp.nodes[end]["hcount"]-1
            comp.add_edge(start,end)
            comp.edges[(start, end)]["order"] = "SINGLE"
            tot +=(len(frag))
            
            #print(comp.nodes(data=True))
        return comp

    @staticmethod
    def _single_connect(g1, g2, start, end):
        label = list(range(len(g1), len(g2)+len(g1)+1))
        old_to_new = {old:lab for (old, lab) in zip(g2.nodes, label)}
        end = old_to_new[end]
        newg = nx.relabel_nodes(g2, old_to_new)
        comp = nx.compose(g1, newg)
        # comp.nodes[start]["hcount"]=0
        # if comp.nodes[end]["aromatic"]:
            
        #     comp.nodes[end]["hcount"]=0
        # elif comp.nodes[end]["hcount"]>0:
        #     comp.nodes[end]["hcount"]=comp.nodes[end]["hcount"]-1
        comp.add_edge(start,end)
        comp.edges[(start, end)]["order"] = "SINGLE"
        return comp



    def _begin_join(self):
        
        self._conv_to_graph()
        to_clean=[]
        print(self.multi_frag and not self.begin_at and not self.end_at)
        if not self.multi_frag:
            
            if not self.begin_at:
                self.begin_at = []
                for mol in self.mol_graphs:
                    cycles = nx.cycle_basis(mol, 0)
                    self.begin_at.append(cycles[0][0])
            if not self.end_at:
                self.end_at = []
                for mol in self.frag_graphs:
                    try:
                        cycles = nx.cycle_basis(mol, 0)
                        self.end_at.append(cycles[0][0])
                    except:
                        self.end_at.append(mol.nodes[0])
            for start, end, mol, frag in zip(self.begin_at, self.end_at, self.mol_graphs, self.frag_graphs):
                to_clean.append(self._single_connect(mol, frag, start, end))
        
        elif self.multi_frag and not self.begin_at and not self.end_at:
            
            if not self.begin_at:
                self.begin_at = [list(g.nodes) for g in self.mol_graphs]
            if not self.end_at:
                self.end_at = [list(f.nodes) for f in self.frag_graphs]
            for mol, nodes in zip(self.mol_graphs, self.begin_at):
                
                for node in nodes:
                    for _node in nodes:
                        if _node == node:
                            continue
                        starts = (node, _node)
                        
                        for end, frag in zip(self.end_at, self.frag_graphs):
                            
                            for _end, _frag in zip(self.end_at, self.frag_graphs):
                                e1 = list(frag.nodes)[0]
                                e2 = list(_frag.nodes)[0]
                                
                                ends = (e1, e2)
                                lg2 = (frag, _frag)
                                
                                to_clean.append(self._two_connect_graphs(mol, lg2, starts, ends))
        elif self.multi_frag and self.begin_at and self.begin_at:
            for mol, starts in zip(self.mol_graphs, self.begin_at):
                starts = list(map(int, starts))
                for frag, end in zip(self.frag_graphs, self.end_at):
                    for _frag, _end in zip(self.frag_graphs, self.end_at):
                        ends = (end, _end)
                        ends = list(map(int, starts))
                        lg2 = (frag, _frag)
                        to_clean.append(self._two_connect_graphs(mol, lg2, starts, ends))
                        
        elif not self.multi_frag and self.begin_at and self.begin_at:
            for start, end, mol, frag in zip(self.begin_at, self.end_at, self.mol_graphs, self.frag_graphs):
                to_clean.append(self._single_connect(mol, frag, start, end))
        
                                
                        
            
        
        
        self._validate_and_remove(to_clean)
    @staticmethod
    def _to_build(dirty):
        bt_dict = {"SINGLE": rdkit.Chem.rdchem.BondType.SINGLE,
                   "DOUBLE" : rdkit.Chem.rdchem.BondType.DOUBLE,
                   "TRIPLE" : rdkit.Chem.rdchem.BondType.TRIPLE,
                   "AROMATIC" : rdkit.Chem.rdchem.BondType.AROMATIC}
        rwmol = Chem.RWMol()
        for node in dirty.nodes:
            curr_sym = dirty.nodes[node]["symbol"]
            atom = Chem.Atom(curr_sym)
            rwmol.AddAtom(atom)
            
        for bond in dirty.edges:
            bt = dirty.edges[bond]["order"]
            bt = bt_dict[bt]
            beg, end = bond
            rwmol.AddBond(beg, end, bt)
        return Chem.MolToSmiles(rwmol)          


    
    def _validate_and_remove(self, to_clean):
        self._connect = []
        for dirty in to_clean:
            try:
                dirty = self._to_build(dirty)
                bm = Molecule(dirty)
                if self.as_smiles:
                    self._connect.append(bm.smiles)
                else:
                    self._connect.append(bm)
            except:
                continue
    
    

    
    
    @property
    def joined_molecules(self):
        if not self.started:
            self._begin_join()
            
        return self._connect
        
        
    
    
    def update_joining(self, molecules, fragments, start_indices=None,  iterative=False,frag_indices =None):
        pass
        
        
        


class MolFuser:
    def __init__(self, base_scaffolds, scaffs_to_attach, start_indices=None,  
                 frag_indices =None, eager=True, infer_best=False, times_to_fuse = 1, 
                 pick_indicies_manually=False, return_as_smiles = True, debug = True, keep_increasing = False):
        
        
        
        
        
        self.base = base_scaffolds
        self.mol_graphs = []
        
        self.to_attach = scaffs_to_attach
        self.frag_graphs = []
        self.begin_at = start_indices
        self.end_at = frag_indices
        
        if pick_indicies_manually:
            self._plot()
        
        
        if not self.begin_at and not self.end_at:
            print("No indices passed... Using iterative combination :)")
            self.iterative = True
        
        if infer_best:
            print("Attempting symmetric joining")
        self.infer_best = infer_best
        self.times_to_fuse = times_to_fuse
        self.joined_times = 1
        self.eager_execution = eager
        self.as_smiles = return_as_smiles
        self.debug = debug
        self.keep_increasing = keep_increasing
        self._old = []
        if eager:
            self._begin_to_fuse()
            self.started = True
        else:
            self.started = False
            
    
    def _plot():
        pass
    @staticmethod
    def pop_ring_edge(graph, aromatic_only = True, selected_indices=None,
                      first_atomatic_2_neighbor=True):
        if aromatic_only:
            if selected_indices:
                what_order = graph.edges[selected_indices]["order"]
                if what_order=="AROMATIC":
                    for selected in selected_indices:
                        graph.remove_node(selected)
            elif first_atomatic_2_neighbor:
                remove_bond = None
                for bond in graph.edges():
                    if graph.edges[bond]["order"]=="AROMATIC":
                        start, end = bond
                        sneigh = list(graph.neighbors(start))
                        s_neigh_l = [list(graph.neighbors(s)) for s in sneigh]
                        if len(max(s_neigh_l, key=len))>2:
                            continue
                        eneigh = list(graph.neighbors(end))
                        e_neigh_l = [list(graph.neighbors(e)) for e in eneigh]
                        if len(max(e_neigh_l, key=len))>2:
                            continue
                        if len(sneigh)==2 and len(eneigh)==2:
                            remove_bond = bond
                    if remove_bond:
                        break
                for rem in remove_bond:
                    graph.remove_node(rem)
                
        reorder = list(range(len(graph)))
        newmap = {old:new for (old, new) in zip(graph.nodes(), reorder)}
        graph = nx.relabel_nodes(graph, newmap)

        # def plot(g):
        #     nx.draw(g)
        #     plt.show()
        # plot(graph)
        return graph
    @staticmethod
    def find_aromatic_one_neighbor(shattered_graph):
        aromatic_one_neighbor = []
        for node in shattered_graph.nodes():
            if shattered_graph.nodes[node]["aromatic"] and len(list(shattered_graph.neighbors(node)))==1:
                aromatic_one_neighbor.append(node)
        return aromatic_one_neighbor
    @staticmethod
    def find_good_edge(graph):
        
        remove_bond = None
        for bond in graph.edges():
            if graph.edges[bond]["order"]=="AROMATIC":
                start, end = bond
                if len(list(graph.neighbors(start)))==2 and len(list(graph.neighbors(end)))==2:
                    remove_bond = bond
            if remove_bond:
                return remove_bond
    
    
    
    def _fuse_frag_to_scaffold(self,scaffold, frag, target_scaffold_edge = None):
        span = range(len(scaffold), len(scaffold)+len(frag)+1)
        newmap = {old:new for (old, new) in zip(frag.nodes(), span)}
        
        one_neighbor = self.find_aromatic_one_neighbor(frag)
        
        
        frag = nx.relabel_nodes(frag, newmap)
        
        first, second = one_neighbor
        first = newmap[first]
        second = newmap[second]
        
        
        compose = nx.compose(scaffold, frag)
        
        if target_scaffold_edge:
            
            beg, end = target_scaffold_edge
            # compose.nodes[beg]["hcount"]=compose.nodes[beg]["hcount"]-1
            # compose.nodes[end]["hcount"]=compose.nodes[end]["hcount"]-1
            compose.add_edge(beg, first)
            compose.add_edge(end, second)
        if not target_scaffold_edge:
            
            beg, end = self.find_good_edge(scaffold)
            # compose.nodes[beg]["hcount"]=compose.nodes[beg]["hcount"]-1
            # compose.nodes[end]["hcount"]=compose.nodes[end]["hcount"]-1
            compose.add_edge(beg, first)
            compose.edges[(beg, first)]["order"] = "SINGLE"
            compose.add_edge(end, second)
            compose.edges[(end, second)]["order"] = "SINGLE"
        return compose 
        

    
    
    def _conv_to_graph(self):
        self.mol_graphs = [Molecule(smile).graph for smile in tqdm(self.base, desc = "Creating graphs...")]
        self.frag_graphs = [Molecule(smile).graph for smile in tqdm(self.to_attach, desc = "Creating graphs...")]
        
        
        
    @staticmethod
    def _to_build(dirty):
        bt_dict = {"SINGLE": rdkit.Chem.rdchem.BondType.SINGLE,
                   "DOUBLE" : rdkit.Chem.rdchem.BondType.DOUBLE,
                   "TRIPLE" : rdkit.Chem.rdchem.BondType.TRIPLE,
                   "AROMATIC" : rdkit.Chem.rdchem.BondType.AROMATIC}
        rwmol = Chem.RWMol()
        for node in dirty.nodes:
            curr_sym = dirty.nodes[node]["symbol"]
            atom = Chem.Atom(curr_sym)
            rwmol.AddAtom(atom)
            
        for bond in dirty.edges:
            bt = dirty.edges[bond]["order"]
            bt = bt_dict[bt]
            beg, end = bond
            rwmol.AddBond(beg, end, bt)
        return Chem.MolToSmiles(rwmol)
    
    def _validate_and_remove(self, to_clean):
        if not self._old:
            self._joined = []
        else:
            self.joined=self._old
        
        for dirty in to_clean:
            try:
                dirty = self._to_build(dirty)
                bm = Molecule(dirty)
                
                if self.as_smiles:
                    self._joined.append(bm.smiles)
                else:
                    self._joined.append(bm)
            except:
                continue
        
    
    
    
    
    def _begin_to_fuse(self):
        self._conv_to_graph()
        if not self.end_at:
            self.end_at = [None] * len(self.to_attach)
        popped = []
        copy_frag = []
        copy_ind = []
        copy_fg = []
        for to_pop, ind, frag in tqdm(zip(self.frag_graphs, self.end_at, self.to_attach), desc = "Popping edges..."):
            try:
                popped.append(self.pop_ring_edge(to_pop, selected_indices=ind))
                copy_frag.append(frag)
                copy_ind.append(ind)
                copy_fg.append(to_pop)
            except:
                if self.debug:
                    print(F"Failed to pop edge for {frag} at indices {ind} \nRemoving from list")
                continue
        
        self.frag_graphs, self.end_at, self.to_attach = copy_fg, copy_ind, copy_frag
        if not self.begin_at:
            self.begin_at = [None] * len(self.mol_graphs)
        to_clean = []
        for scaff, begin, base in tqdm(zip(self.mol_graphs, self.begin_at, self.base), desc="fusing..."):
            for pop, end, att in zip(self.frag_graphs, self.end_at, self.to_attach):
                try:
                    to_clean.append(self._fuse_frag_to_scaffold(scaff, pop, target_scaffold_edge = begin))
                except:
                    if self.debug:
                        print(F"Failed to fuse for {base} at beginning indices {begin} and attached scaffold {att} at indices {end}")
                    continue
        
        self._validate_and_remove(to_clean)
        
        if self.joined_times< self.times_to_fuse:
            self.to_attach = self._joined +self._old
            if not self._old:
                self._old = self.to_attach
            else:
                
                self._old.extend(self.to_attach)
                
            if not self.as_smiles:
                self.to_attach = [bm.smiles for bm in self._joined]
            #self.b.append(self.to_attach)
            if self.keep_increasing:
                self.end_at=None
            self._conv_to_graph()   
            self.joined_times+=1
            
            #self._validate_and_remove()
            self._begin_to_fuse()
            
            
    
    
    
    
    @property
    def fused_molecules(self):
        if not self.started:
            self._begin_to_fuse()
        return tuple(frozenset(self._joined))
    
    
    












