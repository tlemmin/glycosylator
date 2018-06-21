#! /usr/bin/env python
'''
----------------------------------------------------------------------------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

2016 Thomas Lemmin
----------------------------------------------------------------------------
'''


import os
import sys
import re
import copy
import math
import numpy as np
import networkx as nx
from prody import *
from itertools import izip
from scipy.spatial import distance
from scipy.interpolate import interp1d
import sqlite3

from matplotlib.patches import Circle, Rectangle, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

SELF_BIN = os.path.dirname(os.path.realpath(sys.argv[0]))
#sys.path.insert(0, SELF_BIN + '/support')



#####################################################################################
#                                Support functions                                     #
#####################################################################################

def readLinesFromFile(fileName):
    """Reads all lines in a file
    Parameters:
        fileName: path to file
    Returns:
        lines: list with all the lines in a file
    """
    file = open(fileName,'r')                  # open the file
    lines = file.readlines()                  # read all the lines in the file to the list "lines"
    file.close()                                  # close the file
    return lines

def topological_sort(unsorted_graph):
    """Topological sorting of a graph 
    Parameters:
        unsorted_graph: dictionary representation of a graph
    Returns:
        sorted_graph: list of nodes and corresponding edges
    """
    sorted_graph = []
    #sort graph
    while unsorted_graph:
        acyclic = False
        for node, edges in unsorted_graph.items():
            for edge in edges:
                if edge in unsorted_graph:
                    break
            else:
                acyclic = True
                del unsorted_graph[node]
                sorted_graph.append((node, edges))

        if not acyclic:
            print "WARNING! Cyclique dependency occurred in ICs. Impossible to build residue"
            print unsorted_graph
            print sorted_graph
            return ''
            break
    return sorted_graph[::-1]    

def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return izip(a, a)

def rotation_matrix(axis, theta):
    '''Computes the rotation matrix about an arbitrary axis in 3D
    Code from: http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    Parameters:
        axis: axis
        theta: rotation angle
    Return: 
        rotation matrix

    '''

    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

    
def alphanum_sort(l):
    """Alphanumerical sort of a list from
    https://arcpy.wordpress.com/2012/05/11/sorting-alphanumeric-strings-in-python/
    Parameter:
        l: list
    Returns:
        alphanumerically sorted list
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)

aaa2a = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

#####################################################################################
#                                Topology functions                                     #
#####################################################################################

class CHARMMTopology:
    """Class for parsing and storing CHARMM topology files.
    Attributes:
        topology: dictionary storing topology (RES)
                        key: resname
                        value: dictionary with ATOM, BOND, CHARGE and IC
                        key: MASS contains all masses of atoms in topology
         patches: dictionary storing patches (PRES)
                        key: patchname
                        value: dictionary with dele, ATOM, BOND, CHARGE and IC
         atomnames_to_patch: dictionary storing patch name to connect two atoms
                        key: atom1-atom2
                        value: patchname
    """
    def __init__(self, fileIn):
         self.topology = {}
         self.patches = {}
         self.atomnames_to_patch = {}
         self.read_topology(fileIn)
    
    def reset(self):
        """Resets previsously read topology
        """
        self.topology = {}
        self.patches = {}
        self.atomnames_to_patch = {}
        
    def read_topology(self, fileIn):
        """Reads CHARMM topology file. 
        Parameters:
            fileIn: path to topology file
        Initialize:
            topology: dictionary storing topology (RES)
                        key: resname
                        value: dictionary with ATOM, BOND, CHARGE and IC
                        key: MASS contains all masses of atoms in topology
            patches: dictionary storing patches (PRES)
                        key: patchname
                        value: dictionary with dele, ATOM, BOND, CHARGE and IC
            atomnames_to_patch: dictionary storing patch name to connect two atoms
                        key: atom1-atom2
                        value: patchname
        """
        lines = readLinesFromFile(fileIn)
        topo_type=''
        residue={}
        if 'MASS' not in self.topology:
            self.topology['MASS'] = {}
        masses = self.topology['MASS']
        
        for line in lines:                                                             # Loop through each line 
            line = line.split('\n')[0].split('!')[0].split() #remove comments and endl
            if line:
                if line[0]=='RESI' or line[0]=='PRES':
                    if residue:
                        if topo_type == 'RESI':
                            self.topology[resname] = copy.copy(residue)
                        elif topo_type == 'PRES':
                            self.patches[resname] = copy.copy(residue)
                            key = '-'.join(sorted([residue['BOND'][0][1:], residue['BOND'][1][1:]]))
                            #######TODO adapt for multiple patches#####
                            #if atomnames_to_patch.has_key(key):
                            #    atomnames_to_patch[key].append(resname)
                            #else:
                            #    atomnames_to_patch[key]=[resname]
                            ###########################################
                            self.atomnames_to_patch[key] = resname
                    residue['dele'] = []
                    residue['ATOM'] = []
                    residue['BOND'] = []
                    residue['IC'] = []    
                    topo_type = line[0]
                    resname = line[1]
                    residue['CHARGE'] = line[2] 
                elif line[0] == 'ATOM':
                    self.read_atom(line, residue)
                elif line[0] == 'BOND':
                    self.read_bond(line, residue)
                elif line[0] == 'IC':
                    self.read_ICs(line, residue)
                elif line[0] == 'dele':
                    self.read_dele(line, residue)
                elif line[0] == 'MASS':
                    self.read_mass(line, masses)
                    
        if topo_type == 'RESI': 
            self.topology[resname] = copy.copy(residue)
        elif topo_type == 'PRES':
            self.patches[resname] = copy.copy(residue)
            key = '-'.join(sorted([residue['BOND'][0][1:], residue['BOND'][1][1:]]))
            #if atomnames_to_patch.has_key(key):
            #    atomnames_to_patch[key].append(resname)
            #else:
            #    atomnames_to_patch[key]=[resname]
            self.atomnames_to_patch[key] = resname

    def read_mass(self, mass, masses):
        mass[3]=float(mass[3])
        masses[mass[2]]=mass[3:]
         
    def read_atom(self, atom, residue):
        atom[3]=float(atom[3])
        residue['ATOM'].append(atom[1:])
        
    def read_dele(self, delatom, residue):
        residue['dele'] += (delatom[2:])    
        
    def read_bond(self, bond, residue):
        residue['BOND'] += (bond[1:])
        
    def read_ICs(self, ic, residue):
        ic[5:]=map(float, ic[5:])
        residue['IC'].append(ic[1:])

    def get_IC(self, ics, atom):
        atom_ic = ([ic for ic in ics if ic[3]==atom])
        return atom_ic

    def get_atom_name(self, ATOM):
        names=[]
        for a in ATOM:
            names.append(a[0])
        return names
        
    def get_atoms(self, resname):
        return self.topology[resname]['ATOM']
    
    def get_ICs(self, resname):
        return self.topology[resname]['IC']
        
    def get_bonds(self, resname):
        return self.topology[resname]['BOND']
        
class CHARMMParameters:
    """Class for parsing and storing CHARMM parameters files.
        Attributes:
            parameters: dictionary storing parameters
                keys: 'BONDS', 'ANGLES', 'DIHEDRALS', 'NONBONDED', 'IMPROPER', 'NBFIX', 'CMAP' and 'ATOMS'
                values: dictionary of parameters
                        BONDS: atom1-atom2                  ->  k0, d0
                        ANGLES: atom1-atom2-atom3           ->  k0, a0, kub, d0
                        DIHEDRALS: atom1-atom2-atom3-atom4  ->  k0, n, dela
                        NONBONDED: atom1                    ->  
                        IMPROPER: atom1-atom2-atom3-atom4   ->  
                        NBFIX: atom1-atom2                  ->
                        CMAP:
                        ATOM: atom1                         -> mass
    """
    def __init__(self, fileIn):
        self.parameters = {}
        self.read_parameters(fileIn)

    def read_parameters(self, fileIn):
        """Reads CHARMM parameter file. 
        Parameters:
            fileIn: path to parameter file
        Initializes:
            parameters: dictionary storing parameters
                keys: 'BONDS', 'ANGLES', 'DIHEDRALS', 'NONBONDED', 'IMPROPER', 'NBFIX', 'CMAP' and 'ATOMS'
                values: dictionary of parameters
                        BONDS: atom1-atom2                  ->  k0, d0
                        ANGLES: atom1-atom2-atom3           ->  k0, a0, kub, d0
                        DIHEDRALS: atom1-atom2-atom3-atom4  ->  k0, n, dela
                        NONBONDED: atom1                    ->  
                        IMPROPER: atom1-atom2-atom3-atom4   ->  
                        NBFIX: atom1-atom2                  ->
                        CMAP:
                        ATOM: atom1                         -> mass
        """
        lines = readLinesFromFile(fileIn)
        prm = {}
        prm_type = ''
        tags = ['BONDS', 'ANGLES', 'DIHEDRALS', 'NONBONDED', 'IMPROPER', 'NBFIX', 'CMAP', 'ATOMS']
        #initialize parameter dictionary
        for t in tags:
            if not t in self.parameters:
                self.parameters[t] = {}
        for line in lines:                                                             # Loop through each line 
            line = line.split('\n')[0].split('!')[0].split() #remove comments and endl
            if line:
                if line[0] in tags:
                    if prm:
                        if prm_type in self.parameters:
                            self.parameters[prm_type] = dict(self.parameters[prm_type].items() + prm.items()) 
                        else:
                         self.parameters[prm_type] = copy.copy(prm)
                    prm_type = line[0]
                    prm = {}
                    continue
                if prm_type:
                    eval('self.read_'+prm_type+'(line, prm)')
        self.parameters[prm_type] = copy.copy(prm)


    def read_BONDS(self, bond, prm):
        #    CC311D     NC2D1     320.00    1.430
        if len(bond)==4:
            prm['-'.join(bond[0:2])]=map(float, bond[2:])
        else:
            print "Invalid BOND: "+' '.join(bond)

    def read_ANGLES(self, angle, prm):
        #CT1            CC321        HCA2     33.430        110.10     !22.53     2.17900 
        if len(angle)==5 or len(angle)==7:
            prm['-'.join(angle[0:3])]=map(float, angle[3:])
        else:
            print "Invalid ANGLE: "+' '.join(angle)
            
    def read_DIHEDRALS(self, dihe, prm):
        #CC321C    OC3C61    CC311D    NC2D1        0.62    1        0.0
        if len(dihe)==7:
            key='-'.join(dihe[0:4])
            if key in prm:
                prm[key].append(map(float, dihe[4:]))
            else:
                prm[key] = [map(float, dihe[4:])]
        else:
            print "Invalid DIHEDRAL: "+' '.join(dihe)
            
    def read_IMPROPER(self, impr, prm):
        #NC2D1     CC2O1     CC311D    HCP1        20.00    0     0.00
        if len(impr)==7:
            prm['-'.join(impr[0:4])]=map(float, impr[4:])
        else:
            print "Invalid IMPROPER: "+' '.join(impr)

    def read_NONBONDED(self, vdw, prm):
        #CT3            0.0             -0.0780        2.040 ! 0.0 -0.01 1.9 ! alkane, 4/98, yin, adm jr.
        if len(vdw)==4 or len(vdw)==7:
            prm[vdw[0]]=map(float, vdw[1:])
        else:
            print "Invalid NONBONDED: "+' '.join(vdw)
            
    def read_NBFIX(self, nb, prm):
        #SOD        OCL            -0.07502        3.23
        if len(nb)==4:
            prm['-'.join(nb[0:2])]=map(float, nb[2:])
        else:
            print "Invalid NBFIX: "+' '.join(nb)
            
    def read_CMAP(self, cmap, prm):
        return -1
        
    def read_ATOMS(self, atom, prm):
        #MASS        31 H            1.00800
        if len(atom)==4:
            prm[atom[2]]=float(atom[3])
        else:
            print "Invalid ATOM/MASS: "+' '.join(atom)



#####################################################################################
#                                    Molecule                                        #
#####################################################################################
class Molecule:
    """Class for saving a molecule
    Attributes:
            name: structure name
            chain: chain id
            segn: segname
            id: glycosylator id
            key: string representation of connectivity
            atom_group: Prody AtomGroup
            bonds: list of all bonds
            angles: list of all angles
            dihedrals: list of all dihedrals
            connectivity: graph for connectivity of molecule (bonds)
            directed_connectivity: directed acyclique graph of molecule
            interresidue_connectivity: directed acyclique graph representing the interresidue bonds
    """
    def __init__(self, name, chain = 'X', segn = 'X'):
        """initialize AtomGroup used to build pdb from scratch
        Parameters:
            name: structure name (str)
            chain: chain id (str)
            segn: segname (str)
        Initializes:
            atom_group: AtomGroup 
            rootAtom: serial number of atom used as root for graphs
            bonds: list of all bonds
            angles: list of all angles
            dihedrals: list of all dihedrals
            connectivity: graph for connectivity of molecule (bonds)
            directed_connectivity: directed acyclique graph of molecule
            cycle_id: dictionary where keys are the serial number of atom in cycles and values the corresponding cycle in directed_connectivity
            torsionals: dihedral that can rotate (i.e. not in cycles)
            bond_length: dictionary of bond distance used to guess bonds. Keys are sorted by alphabetical order
        """
        self.name = name
        self.atom_group = AtomGroup(self.name)
        self.chain = chain
        self.segn = segn
        self.rootAtom = ','.join([segn, chain, 'O',]) 
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.connectivity = nx.Graph()
        self.directed_connectivity = nx.DiGraph()
        self.interresidue_connectivity = nx.DiGraph()
        self.cycle_id = {}
        self.torsionals = []
        self.bonded_uptodate = False

        self.prefix = ['segment', 'chain', 'resid', 'icode']
        #Defines distance for bond length between different element used in guess_bonds()
        self.elements = ['C', 'H', 'N', 'O']
        self.bond_length = {}
        self.bond_length['C-H'] = 1.20
        self.bond_length['H-N'] = 1.20
        self.bond_length['H-O'] = 1.20
        self.bond_length['C-C'] = 1.7
        self.bond_length['C-N'] = 1.7
        self.bond_length['C-O'] = 1.7

    def writePDB(self, filename, selection = 'all'):
        """Saves molecule to a PDB file
        Parameters:
            filename: path to PDB file
            selection: selection of a subset of the molecule (str)
        """
        writePDB(filename, self.atom_group.select(selection))

    def read_molecule_from_PDB(self, filename, rootAtom = 1, update_bonds = True, **kwargs):
        """Initialize molecule from a PDB file
        Parameters:
            filename: path to PDB file
            rootAtom: serial number of root atom
            update_bonds: guess bonds, angles, dihedrals and connectivity based on the distance between atoms in PDB
            **kwargs: any of the following which allows the selection of a subset of the PDB file
                subset: selection of a subset of the PDB file
                model: model number (int)
                chain: chain id (str)
        Initializes:
            atom_group
            chain: chain id
            segn: segment name
            connectivity: bonds, angles, dihedrals and graph 
        """
        PDBmolecule = parsePDB(filename, **kwargs)
        chain = set(PDBmolecule.getChids())
        segn = set(PDBmolecule.getSegnames())
        self.rootAtom = rootAtom
        if len(chain) == 1 and len(segn) == 1:
            self.atom_group = PDBmolecule
            self.chain = chain
            self.segn = segn
            a1 = self.atom_group.select('serial ' + str(self.rootAtom))
            self.rootRes = a1.getSegnames()[0] + ',' + a1.getChids()[0] + ',' + str(a1.getResnums()[0])+ ',' + a1.getIcodes()[0]
            if update_bonds:
                self.update_connectivity()
            else:
                self.build_connectivity_graph()
        else:
            print "several chains are present in PDB. Please select only one molecule"
            return -1
        return 0

    def set_id(self):
        segn = self.atom_group.getSegnames()
        chid = self.atom_group.getChids()
        res = self.atom_group.getResnums()
        ics = self.atom_group.getIcodes()
        at = self.atom_group.getNames()
        ser = self.atom_group.getSerials()
        ids = {}

        for i,s,c,r,ic,a in zip(ser,segn,chid,res,ics,at):
            ids[i] = {'id': ','.join([s,c,str(r),ic,a])}
        
        self.connectivity.add_nodes_from(ids.items())
    
    def get_ids(self, sel):
        segn = sel.getSegnames()
        chid = sel.getChids()
        res = sel.getResnums()
        ic = sel.getIcodes()
        rn = sel.getResnames()
        ids = []
        for s,c,r,i in zip(segn,chid,res,ic):
            ids.append(','.join([s, c, str(r), i]))
        return ids,rn

    def get_chain(self):
        return self.chain
    
    def get_segname(self):
        return self.segn
    
    def get_residue(self, res_id):
        """Returns an AtomGroup of given atom id; composed of 'segname,chain,resid,icode,atomName'
        """
        
        sel = []
        for p,s in zip(self.prefix, res_id.split(',')):
            if s:
                sel .append(p + ' ' + s)
        sel = ' and '.join(sel)
        return self.atom_group.select(sel)
        
    def get_atom(self, a_id, atom_name):
        """Returns an AtomGroup of given atom id; composed of 'segname,chain,resid,icode'
        """
        for p,s in zip(self.prefix, a_id.split(',')):
            if s:
                sel .append(p + ' ' + s)
        sel = ' and '.join(sel)
        sel += ' and name ' + atom_name
        return self.atom_group.select(sel)

    def set_atom_type(self, atom_type):
        """Assignes atom name, type and charge to each atom in the connectivity graph
        """
        self.connectivity.add_nodes_from(atom_type.items())

    def build_connectivity_graph(self):
        """Builds a connectivity graph for molecules (not protein) in AtomGroup
            Parameters:
                G: undirected graph of connected elements
            Returns:
                names: dictionary with residue id (get_id) as keys and residue name as value
        """
        kd = KDTree(self.atom_group.getCoords())
        kd.search(1.7)
        atoms = kd.getIndices()
        atom_names = self.atom_group.getNames()
        ids,rn = self.get_ids(self.atom_group)
        G = nx.Graph()
        for a1,a2 in atoms:
            id1 = ids[a1]
            rn1 = rn[a1]
            id2 = ids[a2]
            rn2 = rn[a2]
            
            if id1 != id2:
                e = (id1, id2)
                an1 = atom_names[a1]
                an2 = atom_names[a2]
                G.add_node(id1, resname=rn1)
                G.add_node(id2, resname=rn2)
                G.add_edge(id1, id2, patch = '', atoms = an1 + ':' + an2)
        #create directed graph and remove all unnecessary edges
        if G:
            self.interresidue_connectivity = G.to_directed()
            for edge in nx.dfs_edges(G,self.rootRes): 
                edge = list(edge)
                edge.reverse()
                self.interresidue_connectivity.remove_edge(*edge)

    def get_names(self):
        """returns a dictironary with the residue ids as keys and residue names as values
        """
        return nx.get_node_attributes(self.interresidue_connectivity, 'resname')


    def update_connectivity(self, update_bonds = True):
        """Updates all the connectivity (bond, angles, dihedrals and graphs)
        """
        if update_bonds:
            self.guess_bonds()
        #self.guess_angles()
        #self.guess_dihedrals()
        self.update_graphs()
        self.bonded_uptodate = True 

    def set_bonds(self, bonds, update_con = True):
        """ Define list of bonds in a molecule
        Parameters:
            bonds: list of bonds
            update_con: update the connectivity with these new bonds
        """
        inv_atom = {v: k for k, v in nx.get_node_attributes(self.connectivity, 'id').items()}
        newbonds = []
        for b1,b2 in bonds:
            if b1 in inv_atom and b2 in inv_atom:
                newbonds.append((inv_atom[b1], inv_atom[b2]))
            else:
                print 'Skipping bond', b1,b2
        self.connectivity = nx.Graph()
        self.connectivity.add_edges_from(newbonds)
        self.bonds = self.connectivity.edges()
        self.bonded_uptodate = False
        if update_con:
            self.update_connectivity(update_bonds =  False)

    def set_AtomGroup(self, AGmolecule, rootAtom = 1, bonds = None, update_bonds = False):
        """Creates a Molecule instance from AtomGroup. 
        Parameters:
            AGmolecule: prody AtomGroup object
            rootAtom: serial number of rootAtom. Default fist one
            bonds: list of bonds (e.g. generated with MoleculeBuilder)
            update_bonds: if bonds have to be guessed.
        """
        chain = set(AGmolecule.getChids())
        segn = set(AGmolecule.getSegnames())
        self.rootAtom = rootAtom
        if len(chain) == 1 and len(segn) == 1:
            self.atom_group = AGmolecule
            self.chain = chain
            self.segn = segn
            a1 = self.atom_group.select('serial ' + str(self.rootAtom))
            self.rootRes = a1.getSegnames()[0] + ',' + a1.getChids()[0] + ',' + str(a1.getResnums()[0]) + ',' + a1.getIcodes()[0]
            
            if bonds:
                self.set_id()
                self.set_bonds(bonds)

            if update_bonds:
                self.update_connectivity()
            else:
                self.build_connectivity_graph()

        else:
            print "Several chains are present in AtomGroup. Please select only one molecule"
            return -1
        return 0

    def add_residue(self, residue, newbonds, dele_atoms = []):
        """ Add a new residue to a molecule
        Parameters:
            residue: proDy AtomGroup
            newbonds: list of bonds to be added
        """
        if self.atom_group.select('resid ' + ri + 'and chain ' + chid):
            print 'WARNING! A residue with the same id (resid and chain) already exists. The new residue has not been added'
            return -1
            
        if dele_atoms:
            self.delete_atoms(dele_atoms)
            
        natoms = self.atom_group.numAtoms()
        self.atom_group += residue
        self.atom_group.setTitle(self.name)
        self.atom_group.setSerials(np.arange(natoms)+1)
        
        self.connectivity.add_edges_from(np.array(newbonds) + natoms)
        #self.connectivity.remove_edges_from(delete_bonds)
        self.bonds = self.connectivity.edges()
        self.update_connectivity(self, update_bonds = False)

    def delete_atoms(self, dele_atoms):
        """ removes atoms and bonds from molecule
        Parameter:
            del_atoms: serial number of atoms to be deleted
        """
        newbonds = []
        for a in sorted(dele_atoms, reverse = True):
            for b in self.bonds:
                if a in b:
                    continue
                elif a < b[0]:
                    b[0] -= 1
                elif a < b[1]:
                    b[1] -= 1
                newbonds.append(b)
        self.atom_group = self.atom_group.select('not serial ' + del_atoms.join(' ')).copy()
        #renumber atoms 
        self.atom_group.setSerial(np.arange(self.atom_group.numAtoms()))
        self.atom_group.setTitle(self.name)
        self.bonds = newbonds
        self.update_connectivity(self, update_bonds = False)    

    def guess_bonds(self, default_bond_length = 1.6):
        """Searches for all bonds in molecule
        Parameters:
            default_bond_length: maximum distance between two connected heavy atoms (Angstrom), if not present in bond_length dictionary
        """
        self.connectivity = nx.Graph()
        for a in self.atom_group:
            bonds = []
            sel = ''
            a_elem = a.getElement()
            if a_elem:
                # use predefined bond length for atom pairs 
                for e in self.elements:
                    key = '-'.join(sorted([a_elem, e]))
                    if key in self.bond_length:
                        if sel:
                            sel += ' or '
                        sel += '((element ' + e + ') and (within ' + str(self.bond_length[key]) + ' of serial ' + str(a.getSerial()) + '))'
            if not sel:
                sel = 'within ' + str(default_bond_length) + ' of serial ' + str(a.getSerial())
            sel = '(' + sel + ') and (not serial ' + str(a.getSerial()) + ')'
            
            # search for all neighboring atoms
            neighbors =  self.atom_group.select(sel)
            if neighbors:
                for aa in neighbors:
                        bonds.append((a.getSerial(), aa.getSerial()))

            self.connectivity.add_edges_from(bonds)
        self.bonds = self.connectivity.edges()
    

    def guess_angles(self):
        """Searches for all angles in a molecule based on the connectivity
        """
        self.angles = []
        for node in self.connectivity.nodes():
            self.angles.extend(self.find_paths(self.connectivity, node, 2))


    def guess_dihedrals(self):
        """Searches for all dihedrals in a molecule based on the connectivity
        """
        self.dihedrals = []
        for node in self.connectivity.nodes():
            self.dihedrals.extend(self.find_paths(self.connectivity, node, 3))
        
    def find_paths(self, G, node, length, excludeSet = None):
        """Finds all paths of a given length
        Parameters:
            G: graph (netwrokx)
            node: starting node
            length: length of path
            excludedSet: set
        Returns:
            paths: list of all paths of a length starting from node 
        """
        if excludeSet == None:
            excludeSet = set([node])
        else:
            excludeSet.add(node)
          
        if length == 0:
            return [[node]]
        paths = [[node]+path for neighbor in G.neighbors(node) if neighbor not in excludeSet for path in self.find_paths(G, neighbor, length-1, excludeSet)]
        excludeSet.remove(node)
        return paths

    def set_rootAtom(self, rootAtom):
        """Sets the rootAtom and updates all the directed graph
        """
        self.rootAtom = rootAtom
        a1 = self.atom_group.select('serial ' + str(self.rootAtom))
        self.rootRes = a1.getSegnames()[0] + ',' + a1.getChids()[0] + ',' + str(a1.getResnums()[0] + ',' + a1.getIcodes()[0]) 
        self.update_graphs()

    def update_graphs(self):
        """Updates connectivity and directed graphs.
            - seaches for all cycles in connectivity graph
            - rebuilts acyclique directed connectivity graph
                starts from rootAtom
                cycle are collapsed
                nodes have the parameters
                        isclycle: True/False
                        cycle: serial number of all atoms in cycle
                        name of node: 
                                        not cycle: serial number of atom
                                        cycle: string with all serial number joined by a '-'
        """
        cycles = nx.cycle_basis(self.connectivity, self.rootAtom)
        #flatten cycles
        self.cycle_id = {}
        for cycle in cycles:
            key = '-'.join(map(str, cycle))
            for a in cycle:
                self.cycle_id[a] = key
        self.directed_connectivity = nx.DiGraph()
        
        for edge in nx.dfs_edges(self.connectivity,self.rootAtom):
            directed_edge = []
            atoms = []
            for node in edge:
                atoms.append(self.atom_group.select('serial ' + str(node)))
                if node in self.cycle_id:
                    key = self.cycle_id[node]
                    if key not in self.directed_connectivity:
                        self.directed_connectivity.add_node(key, iscycle = True, cycle_id = map(int, self.cycle_id[node].split('-')))
                    directed_edge.append(key)
                else:
                    self.directed_connectivity.add_node(node, iscycle = False, cycle_id = [])
                    directed_edge.append(node)
            if directed_edge[0] == directed_edge[1]:
                continue
            self.directed_connectivity.add_edge(directed_edge[0], directed_edge[1])
            a1,a2 = atoms
            if a1.getResnums()[0] != a2.getResnums()[0]:
                r1 = a1.getSegnames()[0] + ',' + a1.getChids()[0] + ',' + str(a1.getResnums()[0]) + ',' + a1.getIcodes()[0]
                r2 = a2.getSegnames()[0] + ',' + a2.getChids()[0] + ',' + str(a2.getResnums()[0]) + ',' + a2.getIcodes()[0]
                self.interresidue_connectivity.add_node(r1, resname=a1.getResnames()[0])
                self.interresidue_connectivity.add_node(r2, resname=a2.getResnames()[0])
                self.interresidue_connectivity.add_edge(r1, r2, patch = '', atoms = a1.getNames()[0] + ':' + a2.getNames()[0])
                

    def define_torsionals(self, hydrogens=True):
        """Builds a list with all the torsional angles that can rotate
        Parameters:
            hydrogens: include torsional involving terminal hydrogens        
        Initializes:
            torsionals: a list of serial number of atom defining a torsional angle (quadruplet) 
        """
        self.torsionals = []
        #cycles = nx.cycle_basis(self.connectivity, self.rootAtom)
        if not hydrogens:
            elements = nx.get_node_attributes(self.connectivity, 'element')
            if not elements:
                print 'Elements have not been defined (use assign_atom_type). Hydrogens cannot be excluded.'
                hydrogens = True
        
        for dihe in self.dihedrals:
            #check that quadruplet is not in cycle
            if dihe[1] in self.cycle_id and dihe[2] in self.cycle_id:
                continue
            
            d_dihe = []
            for a in dihe[1:-1]:
                if a in self.cycle_id:
                    a = self.cycle_id[a]
                d_dihe.append(a)

            #check direction of dihedral
            if self.directed_connectivity.has_edge(d_dihe[0], d_dihe[1]):
                pass
            elif self.directed_connectivity.has_edge(d_dihe[1], d_dihe[0]):
                dihe.reverse() 
            else:
                continue

            #check if hydrogen
            if not hydrogens:
                if elements[dihe[0]] == 'H' or elements[dihe[-1]] == 'H':
                    continue
            #check if already in torsionals list
            exists = False
            if self.torsionals:
                for t in self.torsionals:
                    if dihe[1] == t[1] and dihe[2] == t[2]:
                        exists = True
                        break

            if exists:
                continue

            self.torsionals.append(dihe)

    def rotate_bond(self, torsional, theta, absolute =False):
        """Rotate the molecule around a torsional angle. Atom affected are in direct graph.
        Parameters:
            torsional: index of torsional angle (in torsionals) or list of serial number of atoms defining the torsional angle.
            theta: amount (degrees)
            absolute: theta 
        """
         
        if type(torsional) == int:
            torsional = self.torsionals[torsional]
        else:
            if torsional not in self.torsionals:
                print "Unknown torsional"
                return -1

        atoms = []
        a1 = torsional[-2]

        #check if last atom of torsional angle is in cycle
        if a1 in self.cycle_id:
            a1 = self.cycle_id[a1]
            atoms += a1.split('-')
        else:
            atoms.append(str(a1))

        for n in nx.descendants(self.directed_connectivity, a1):
            if type(n) == str:
                atoms += n.split('-')
            else:
                atoms.append(str(n))
        sel = self.atom_group.select('serial ' + ' '.join(atoms))
        axis_sel = self.atom_group.select('serial ' + ' '.join(map(str, torsional[1:-1])))
        v1,v2 = axis_sel.getCoords()
        axis = v2-v1
        c_angle = 0.
        if absolute:
            vec_sel = self.atom_group.select('serial ' + ' '.join(map(str, torsional)))
            c1,c2,c3,c4 = vec_sel.getCoords()
            c12 = c1-c2
            c12 = c12 / np.linalg.norm(c12)
            c43 = c4-c3
            c43 = c43 / np.linalg.norm(c43)
            c_angle = np.rad2deg(np.cos(np.dot(c12,c43)))
            theta = theta - c_angle
        coords = sel.getCoords() - v2
        M = rotation_matrix(axis, np.deg2rad(theta))
        coords = M.dot(coords.transpose())
        sel.setCoords(coords.transpose()+v2)
        return c_angle


#####################################################################################
#                                Builders                                            #
#####################################################################################
class MoleculeBuilder:
    """Class for building/modifying molecule
    """
    
    def __init__(self, topofile, paramfile, force_field = 'charmm'):
        """
        Parameters:
            topofile: path to topology file
            paramfile: path to parameters file
            force_field: force field name. Currently only CHARMM
        """
        if force_field == 'charmm':
            self.Topology = CHARMMTopology(topofile) 
            self.Parameters = CHARMMParameters(paramfile)
        else:
            print "unknown force field."


    def init_new_residue(self, resid, resname, chain, segname, i = 1):
        """Initializes a residue from scratch
        Parameters:
            resid: residue id (int)
            resname: name of residue (str)
            chain: chain id (char)
            segn: segname (str)
        Returns:
            residue: AtomGroup with all atoms initialized (from Topology)
            atoms_name: name of of all the atoms in residue
            bonds: list of bonds (segn,chid,resi,atom_name)
        """
        residue = AtomGroup(resname+str(resid))
        # add all the atoms
        atoms = self.Topology.get_atoms(resname)
        natoms = len(atoms)
        coords = np.zeros((natoms,3))
        residue.setCoords(coords)
        resn = []
        resi = []
        atoms_name = []
        chid = []
        segn = []
        occupancy = []
        beta = []
        serial = []
        element = []
        icode = []
        altloc = []
        for a in atoms:
            atoms_name.append(a[0])
            resn.append(resname)
            resi.append(resid)
            chid.append(chain)
            segn.append(segname)
            occupancy.append(1)
            beta.append(0)
            serial.append(i)
            i += 1
            icode.append('')
            element.append(a[0][0])
            altloc.append('')

        residue.setNames(atoms_name)
        residue.setResnums(resi)
        residue.setResnames(resn)
        residue.setChids(chid)
        residue.setSegnames(segn)
        residue.setOccupancies(occupancy)
        residue.setBetas(beta)
        residue.setSerials(serial)
        residue.setIcodes(icode)
        residue.setElements(element)
        residue.setAltlocs(altloc)
        bonds = []
        top_bonds = self.Topology.get_bonds(resname)

        id_r = '%s,%s,%d,,' % (segname,chain,resid)
        for a1,a2 in pairwise(top_bonds):
            bonds.append((id_r+a1, id_r+a2))

        return residue, atoms_name, bonds

    def copy_atom(self, src_atom, dst_atom):
        """copies all the attributes of one atom to another
        Parameters:
            src_atom: original atom
            dst_atom: copy atom
        """
        dst_atom.setCoords(src_atom.getCoords())
        dst_atom.setNames(src_atom.getNames())
        dst_atom.setResnums(src_atom.getResnums())
        dst_atom.setResnames(src_atom.getResnames())
        dst_atom.setChids(src_atom.getChids())
        dst_atom.setSegnames(src_atom.getSegnames())
        dst_atom.setOccupancies(src_atom.getOccupancies())
        dst_atom.setBetas(src_atom.getBetas())
        dst_atom.setSerials(src_atom.getSerials())
        dst_atom.setIcodes(src_atom.getIcodes())
        dst_atom.setElements(src_atom.getElements())
        dst_atom.setAltlocs(src_atom.getAltlocs())
        
    def add_missing_atoms(self, residue, resid = None):
        """Add all missing atoms to a ProDy residue from topology
        Parameters:
            residue: ProDy residue (AtomGroup)
        Retruns:
            complete_residue: Completed ProDy residue, with coordinates of missing atoms set to (0, 0, 0)
            missing_atoms: list of missing atom names
            bonds: list of new bonds
        """
        if not resid:
            resid = residue.getResnums()[0]
        complete_residue,atoms,bonds = self.init_new_residue(resid, residue.getResnames()[0], residue.getChids()[0], residue.getSegnames()[0])
        missing_atoms = []
        atoms_in_residue = residue.getNames()
        for a in atoms:
            if a in atoms_in_residue:
                atom = residue.select('name ' + a)
                catom = complete_residue.select('name ' + a)
                self.copy_atom(atom, catom)
            else:
                missing_atoms.append(a)
        complete_residue.setResnums([resid]*len(complete_residue)) 
        return complete_residue, missing_atoms, bonds

    def apply_patch(self, patch,  residue1, residue2):
        """
        Parameters:
            patch: name of patch (str) 
            residue1: first residue in patch (ProDy AtomGroup)
            residue2: second residue in patch (ProDy AtomGroup)
        Return:
            bonds: list of all new bonds
            dele_atoms: list of atoms which should be deleted
        """
        
        dele_atoms = self.dele_atoms(patch, residue1, residue2)
        bonds = self.patch_bonds(patch, residue1, residue2)
        return dele_atoms,bonds

    def build_IC_graph(self, atoms, ics):
        """Extracts ICs to build missing atoms
        Parameters:
            atoms: list of missing atoms
            ics: list of internal coordinates
        Returns:
            unsorted_graph: dictionay representation of dependency of ics
                            key: atom name
                            value: list of atoms connected to key
            required_atoms: list atoms required to build missing atoms
        """
        unsorted_graph = {}
        required_atoms = []
        atomsIC = [atom.replace('*', '') for ic in ics for atom in ic[0:4]]
        #Build graph
        for a in atoms:
            if a in atomsIC:
                required_atoms.append(a)
                if not a in unsorted_graph:
                    unsorted_graph[a] = []
                for ic in self.Topology.get_IC(ics, a):
                    for aic in ic[0:3]:
                        aic=aic.replace('*', '')
                        if aic in unsorted_graph:
                            unsorted_graph[aic].append(a)
                        else:
                            unsorted_graph[aic] = [a]
        return unsorted_graph,required_atoms


    def build_from_patch(self, link_residue, resid, resname, chain, segname, patch):
        """Build residue from a patch
        Parameters:
            link_residue: residue that the patch will use to build new residue
            resid: residue number
            resname: residue name
            chain: residue chain
            segname: residue segname
            patch: name of patch (str) 
        Returns:
            denovo_residue: complete residue (AtomGroup)
            dele_atoms: list of atoms which should be deleted
            bonds: list of all new bonds
        """
        denovo_residue, missing_atoms, bonds = self.init_new_residue(resid, resname, chain, segname)
        ics = self.Topology.patches[patch]['IC']
        #patch_atoms = sorted(set([atom.replace('*', '')[1:] for ic in ics for atom in ic[0:4] if atom.replace('*', '')[0]=='2']))
        patch_atoms = sorted(set([atom.replace('*', '') for ic in ics for atom in ic[0:4] if atom.replace('*', '')[0]=='2']))
        self.build_patch_missing_atom_coord(link_residue, denovo_residue, patch_atoms, ics)
        missing_atoms = [a for a in missing_atoms if '2' + a not in patch_atoms]
        ics = self.Topology.topology[resname]['IC']
        self.build_missing_atom_coord(denovo_residue, missing_atoms, ics)

        dele_atoms,b =  self.apply_patch(patch, link_residue, denovo_residue)
        bonds.extend(b)
        return denovo_residue, dele_atoms, bonds
    
    def patch_bonds(self, patch, residue1, residue2 = None):
        """
        Parameters:
            patch: name of patch (str)
            residue1: first residue in patch
            residue2: second residue in patch. None if not required
        Returns:
            bonds: list of bonds
        """
        bonds = []
        segn1 = residue1.getSegnames()[0] 
        chid1 = residue1.getChids()[0]
        resi1 = str(residue1.getResnums()[0])
        ic1 = str(residue1.getIcodes()[0])
        if residue2:
            segn2 = residue2.getSegnames()[0] 
            chid2 = residue2.getChids()[0]
            resi2 = str(residue2.getResnums()[0])
            ic2 = str(residue2.getIcodes()[0])
        for a1,a2 in pairwise(self.Topology.patches[patch]['BOND']):
            b = []
            for a in [a1,a2]:
                if a[0] == '1':
                    b.append('%s,%s,%s,%s,%s'%(segn1, chid1, resi1, ic1, a[1:]))
                if a[0] == '2':
                    if residue2:
                        b.append('%s,%s,%s,%s,%s'%(segn2, chid2, resi2, ic2, a[1:]))
                    else:
                        print "Warning BOND: missing residue2 required for patch " + patch
            bonds.append(b)
        return bonds


    def dele_atoms(self, patch, residue1, residue2 = None):
        """
        Parameters:
            patch: name of patch (str)
            residue1: first residue in patch
            residue2: second residue in patch. None if not required
        Returns:
            dele_atoms: list of atoms to delete. (segn, chid, resi, atom_name)
        """
        atoms = self.Topology.patches[patch]['dele']
        dele_atoms = []
        segn1 = residue1.getSegnames()[0] 
        chid1 = residue1.getChids()[0]
        resi1 = str(residue1.getResnums()[0])
        ic1 = str(residue1.getIcodes()[0])
        if residue2:
            segn2 = residue2.getSegnames()[0] 
            chid2 = residue2.getChids()[0]
            resi2 = str(residue2.getResnums()[0])
            ic2 = str(residue2.getIcodes()[0])

        for a in atoms:
            if a[0] == '1':
                dele_atoms.append('%s,%s,%s,%s,%s'%(segn1, chid1, resi1, ic1, a[1:]))
            if a[0] == '2':
                if residue2:
                    dele_atoms.append('%s,%s,%s,%s,%s'%(segn2, chid2, resi2, ic2, a[1:]))
                else:
                    print "Warning: missing residue2 required for patch " + patch
        return dele_atoms

    def delete_atoms(self, molecule, dele_atoms):
        """Removes all atoms that should be delete by patch
        Parameters:
            molecule: AtomGroup defining the molecule
            dele_atoms: list of atoms to be deleted. (segn, chid, resi, atom_name)
        """
        sel = []
        prefix = ['segment', 'chain', 'resid', 'icode', 'name']
        for a in dele_atoms:
            s1 = []
            for p,s in zip(prefix, a.split(',')):
                if s:
                    s1.append(p + ' ' + s)
           
            sel.append(' and '.join(s1))
        sel = 'not (' + ' or '.join(sel) + ')'
        
        return  molecule.select(sel).copy()

    def build_from_DUMMY(self, resid, resname, chain, segname, dummy_patch, dummy_coords = [[0, 0, 0], [0, 0, 1], [0, 1, 1]]):
        """Builds residue from DUMMY atoms
        Parameters:
            resid: residue id (int)
            chain: residue chain id (chr)
            segname: residue segname (str)
            dummy_patch: patch to build new residue 
            dummy_coords: coordinated of dummy atoms
        Returns:
            denovo_residue: complete residue (AtomGroup)
            bonds: list of all bonds
        """
        dummy_residue, dummy_atoms, bonds = self.init_new_residue(0, 'DUMMY', 'D', 'DUM')
        counter = 0
        for a in dummy_atoms:
            dummy_residue.select('name ' + a).setCoords([dummy_coords[counter]])
            counter += 1
        denovo_residue, dele_atoms, bonds = self.build_from_patch(dummy_residue, resid, resname, chain, segname, dummy_patch)
        del dummy_residue
        #bonds.extend(bonds_p)
        return denovo_residue, dele_atoms, bonds

    def build_patch_missing_atom_coord(self, link_residue, residue, missing_atoms, ICs):
        """Builds all missing atoms in residue from a patch linking it to link_residue
        Parameters:
            link_residue: first residue in patch (AtomGroup)
            residue: second residue in patch (AtomGroup)
            missing_atoms: list of missing atom in second residue
            ICs: list of internal coordinate to build missing atoms
        """
        unsorted_graph,required_atoms = self.build_IC_graph(missing_atoms, ICs)
        sorted_graph = topological_sort(unsorted_graph)
        atoms=[g[0] for g in sorted_graph if g[0] in required_atoms]
        for a in atoms:
            atomname = a[1:]
            atom = residue.select('name ' + atomname)
            ic = self.Topology.get_IC(ICs, a)
            if ic:
                ic = ic[0]
                c = 0
                for atom_ic in ic[0:4]:
                    if atom_ic.replace('*', '')[0]=='2':
                        exec('xa'+ str(c) +'= residue.select(\'name \' + atom_ic.replace(\'*\', \'\')[1:]).getCoords()[0]')        
                    else:
                        exec('xa'+ str(c) +'= link_residue.select(\'name \' + atom_ic.replace(\'*\', \'\')[1:]).getCoords()[0]')        
                    c += 1
                atom.setCoords(self.build_cartesian(xa0, xa1, xa2, ic[8], ic[7], ic[6]))

    def build_missing_atom_coord(self, residue, missing_atoms, ICs):
        """Builds all missing atoms based on the provided internal coordinates
            Parameters:
                residue: Prody residue (AtomGroup)
                missing_atoms: list with all missing atom name
                ICs: list of internal coordinates for building missing atoms
        """
        unsorted_graph,required_atoms = self.build_IC_graph(missing_atoms, ICs)
        sorted_graph = topological_sort(unsorted_graph)
        atoms=[g[0] for g in sorted_graph if g[0] in required_atoms]
        for a in atoms:
            atom = residue.select('name ' +a)
            ic = self.Topology.get_IC(ICs, a)
            if ic:
                ic = ic[0]
                xa1 = residue.select('name ' + ic[0]).getCoords()[0]
                xa2 = residue.select('name ' + ic[1]).getCoords()[0]
                xa3 = residue.select('name ' + ic[2].replace('*', '')).getCoords()[0]
                atom.setCoords(self.build_cartesian(xa1, xa2, xa3, ic[8], ic[7], ic[6]))    

    def build_cartesian(self, a1, a2, a3, r, theta, phi):
        """Builds missing atom from internal coordinates
            Parameters:
                a1: coordinates of atom1
                a2: coordinates of atom2
                a3: coordinates of atom3
                r: distance from atom3
                theta: angle between a2 a3 and missing atom
                phi: torsional angle formed by a1, a2, a3 and missing atom
            Returns:
                coordinates of missing atom
        """

        theta = np.radians(theta)
        phi  = np.radians(phi)
        cost = np.cos(theta)
        sint = np.sin(theta)
        cosp = np.cos(phi)
        sinp = np.sin(phi)
        rjk = a2 - a3
        rjk /= np.linalg.norm(rjk)
        rij = a1 - a2
        cross = np.cross(rij, rjk)
        cross /= np.linalg.norm(cross)
        cross2 = np.cross(rjk, cross)
        cross2 /= np.linalg.norm(cross2)
        wt = [r * cost, r * sint * cosp, r * sint * sinp]
        newc = rjk*wt[0] + cross2*wt[1] + cross*wt[2]
        return a3 + newc
    
    def get_bonds(self, residue):
        rn = residue.getRenames()[0]
        bonds = []
        for a1,a2 in pairwise(self.Topology[rn]['BOND']):
            i1 = residue.select('name ' + a1).getSerials[0]
            i2 = residue.select('name ' + a2).getSerials[0]
            bonds += (i1, i2)
        return bonds

    
#####################################################################################
#                                Glycosylator                                        #
#####################################################################################
class Glycosylator:
    def __init__(self, topofile, paramfile, force_field = 'charmm'):
        """
        Parameters:
            topofile: path to topology file
            paramfile: path to parameter file
            force_field: name of force field. Default charmm
        Initializes:
            builder: MoleculeBuilder
            connect_topology: dictionary describing the topology of known glycans
            glycan_keys: dictionary for identifying glycans (built from connect_topology) 
            glycoprotein: Atomgroup representing the glycoprotein
            protein: Atomgroup without the identified glycans
            sequences: dictionary with protein chain as keys and sequences as value
            sequons: dictionary with sequon residue id as keys and sequon triplet as value
            glycanMolecules: dictionary with protein residue as keys and Molecules instances as values
            glycans: dictionary with protein residue as keys and graph of connected glycan as values
            names: dictionary of residue names for linked glycans
        """
        self.builder = MoleculeBuilder(topofile, paramfile)
        self.connect_topology = {}
        self.glycan_keys = {}
        self.glycoprotein = None 
        self.protein = None
        self.sequences = {}
        self.sequons = {}
        self.glycanMolecules = {}
        self.glycans = {}
        self.names = {}
        self.prefix = ['segment', 'chain', 'resid', 'icode']

    def read_connectivity_topology(self, connectfile):
        """Parse file defining the topology of connectivity trees.
        This function will initialize connect_topology and glycan_keys

        Each connectivity is defined by
            RESI: name of the polyglycan
            UNIT: resname, [list of patches to residue]

        Parameter:
            fileName: path to connectivity file
        """
        lines = readLinesFromFile(connectfile)
        residue = {}
        self.connect_topology = {}
        nbr_units = 0
        for line in lines:                                                             # Loop through each line 
            line = line.split('\n')[0].split('!')[0].split() #remove comments and endl
            if line:
                if line[0] == 'RESI':    
                    if residue:
                            residue['#UNIT'] = nbr_units
                            self.connect_topology[resname] = copy.copy(residue)
                    residue['UNIT'] = []
                    resname = line[1]
                    nbr_units = 0
                elif line[0] == 'UNIT':
                    self.read_unit(line, residue)
                    nbr_units += 1
        residue['#UNIT'] = nbr_units
        self.connect_topology[resname] = copy.copy(residue)
        self.build_keys()
    
    def import_connectivity_topology(self, filename):
        """Import connectivity topology from sql database
        This function will initialize connect_topology
        Parameters:
            filename: path to database
        """
        try:
            conn = sqlite3.connect(filename)
        except:
            print "Error while connecting to the database " + filename
            return -1
        cursor =  conn.cursor()
        cursor.execute("SELECT * FROM glycans")
        glycans = cursor.fetchall()

        self.connect_topology = {}
        for glycan in glycans:
            name,tree =  glycan
            residue = {}
            residue['UNIT'] = []
            nbr_unit = 0
            for unit in tree.split('|'):
                unit = unit.split(' ')
                nbr_unit += 1
                if len(unit) > 2:
                    residue['UNIT'].append([unit[0],unit[1], unit[2:]])
                else:
                    residue['UNIT'].append([unit[0], ' ', []])
                
            residue['#UNIT'] = nbr_unit
            self.connect_topology[name] = residue
        self.build_keys()

    def add_glycan_to_connectivity_topology(self, name, connect_tree, overwrite = True):
        """Add new glycan to connect_topology dictionary
        Parameters:
            name: name of new glycan
            connect_tree: dictionary with connectivity tree
            overwrite: should an existing glycan be overwritten
        """
        if name in self.connect_topology and not overwrite:
            print "Glycan with same name " + name + "already exists. Please change name or allow overwritting"
            return -1
        self.connect_topology[name] = connect_tree
    
    def export_connectivity_topology(self, filename):
        """Export connectivity topology to sql database
        This function will 
        """
        try:
            conn = sqlite3.connect(filename)
        except:
            print "Error while connecting to the database " + filename
            return -1
        cursor =  conn.cursor()
        tn = 'glycans'
        gn = 'glycan_name'
        gt = 'glycan_tree'
        cursor.execute("DROP TABLE IF EXISTS {tn}".format(tn = tn))
        cursor.execute("CREATE TABLE {tn} ({gn} text, {gt} text)".format(tn =  tn, gn = gn, gt = gt))

        for key in self.connect_topology.keys():
            units = self.connect_topology[key]['UNIT']
            glycan = []
            for unit in units:
                v = []
                v.extend(unit[0:2])
                v.extend(unit[2])
                glycan.append(' '.join(v))
            glycan = '|'.join(glycan)
            
            cursor.execute("INSERT INTO {tn} VALUES (\'{gn}\', \'{gt}\')".format(tn = tn, gn = key, gt = glycan))

        conn.commit()
        conn.close()
    
    def build_keys(self):
        self.glycan_keys = {}
        for res in self.connect_topology:
            key = [r[0]+' '+' '.join(r[2]) for r in self.connect_topology[res]['UNIT']]
            self.glycan_keys[res] = key

    def read_unit(self, unit, residue):
        if len(unit)>2:
            residue['UNIT'].append([unit[1], unit[2], unit[3:]])
        else:
            residue['UNIT'].append([unit[1], ' ', []])
    
    def build_glycan_gaph(self, connect_tree):
        """Builds a graph representation of a connectivity tree
        Parameters:
            connect_tree: dictionnary of connectivity of a poly glycan
        Returns:
            unsorted_graph: dictionary repersenting the connectivity graph
        """
        inv_connect_tree = {v: k for k, v in connect_tree.items()} 
        unsorted_graph = {}
        for g in connect_tree:
            if not g in unsorted_graph:
                unsorted_graph[g]=[]
            key=connect_tree[g]
            if not key:
                continue
            gr=inv_connect_tree[' '.join(key.split()[:-1])]
            if gr in unsorted_graph:
                unsorted_graph[gr].append(g)
            else:
                unsorted_graph[gr]=[g]
        return unsorted_graph
    
    def load_glycoprotein(self, protein):
        """Load and detects glycans in a glycoprotein
        Parameters:
            protein: Atomgroup of 
        """
        if type(protein) == str:
            protein = parsePDB(protein)
        self.glycoprotein = protein
        sel  = self.glycoprotein.select('name CA')
        segn = sel.getSegnames()
        chids = sel.getChids()
       
        fragids = []
        for s,c in zip(segn,chids):
            fragids.append(s + ',' + c)
        fragids = set(fragids)
        for i in fragids:
            sel = ['name CA']
            for p,s in zip(['segment ', 'chain '], i.split(',')):
                if s:
                    sel.append(p+s)
            sel = ' and '.join(sel)
            sel  = self.glycoprotein.select(sel)
            #Old versions of Prody does have getSequence
            try:
                self.sequences[i] = sel.getSequence()
            except:
                self.sequences[i] = self.getSequence(sel)
                
            self.sequons.update(self.find_sequons(sel, self.sequences[i]))

        self.glycans,self.names = self.find_glycans('NGLB', 'ASN')
        self.extract_glycans()

    def save_glycoprotein(self, filename):
        """Saves glycoprotein to PDB file
        Parameter:
            filename: path 
        """
        glycosylated_protein = self.protein.copy()
        for g in self.glycanMolecules.values():
            glycosylated_protein += g.atom_group
            
        writePDB(filename, glycosylated_protein)

    def get_residue(self, res_id):
        """Returns an AtomGroup of given atom id; composed of 'segname,chain,resid,icode,atomName'
        """
        sel = []
        for p,s in zip(self.prefix, res_id.split(',')):
            if s:
                sel .append(p + ' ' + s)
        sel = ' and '.join(sel)
        return self.glycoprotein.select(sel)

    def getSequence(self, sel):
        res =  sel.getResnames()
        return "".join(aaa2a[r] for r in res)
    
    def get_start_resnum(self, chid):
        sel = ['name CA']
        for p,s in zip(['segment ', 'chain '], chid.split(',')):
            if s:
                sel.append(p+s)
        return self.glycoprotein.select(' and '.join(sel)).getResnums()[0]

    def get_glycans(self, sequons):
        """Returns a dictornary of glycans detected for sequons
        """
        glycans = {}
        for s in sequons:
            if s in self.glycans:
                glycans[s] = self.glycans[s]
        return glycans

    def get_sequons(self, chid):
        """Returns all sequons of a chain (chid)
        """
        return [k for k in self.sequons.keys() if chid in k[:len(chid)]]

    def find_sequons(self, sel, seq = None):
        """finds all sequon in AtomGroup
        """
        if not seq:
            seq = sel.getSequence()
        sequons = {}
        res = sel.getResnums()
        icodes = sel.getIcodes()
        segn =  sel.getSegnames()
        chid = sel.getChids()
        for m in re.finditer(r'(N[^P][S|T])', seq):
            idx = m.start()
            sequons[','.join([segn[idx], chid[idx], str(res[idx]), icodes[idx]])] = m.group()
        return sequons
    
    def find_glycans(self, patch, resname):
        """Looks from all molecules (not protein) which are linked (1.7A) to residue
            Parameters:
                patch: name of patch that connects molecule to residue
                resname: residue name to which the molecule is bound
            Returns:
                glycans: dictionary with residue id as key and a list of root unit and graph of glycan
                names: dictionary with glycan id as key and glycan resname as value
        """
        at = self.builder.Topology.patches[patch]['BOND'][0:2]
        for a in at:
            if a[0] == '1':
                a1 = a[1:]
            else:
                a2 = a[1:]
        #a2 is assumed to be the glycan atoms and a1 from protein
        sel = self.glycoprotein.select('(not protein and name %s) or (resname %s and name %s)' % (a2,resname,a1))
        kd = KDTree(sel.getCoords())
        kd.search(1.7)
        atoms = kd.getIndices()
        G = nx.Graph()
        ids,rn = self.get_ids(sel) 
        for a1,a2 in atoms:
            id1 = ids[a1]
            rn1 = rn[a1]
            id2 = ids[a2]
            rn2 = rn[a2]

            if id1 != id2:
                e = (id1, id2)
                G.add_edge(*e)
        names = self.connect_all_glycans(G)
        glycans = {}
        for graph in list(nx.connected_component_subgraphs(G)):
            for node in graph.nodes():
                if node not in names:
                    root_glycan = list(graph.neighbors(node))[0] 
                    graph.remove_node(node)
                    glycans[node] = [root_glycan, graph]
                    break
        return glycans,names
    
    def extract_glycans(self):
        sel_all = []    
        for k,v in self.glycans.items():
            r,g = v
            sel = []
            for node in g.nodes():
                selg = []
                for p,s in zip(self.prefix, node.split(',')):
                    if s:
                        selg.append(p + ' ' + s)
                if not sel:
                    # !!!hard coded!!! Should be patch dependent
                    rootAtom = self.glycoprotein.select(' and '.join(selg) + ' and name C1').getSerials()[0]
                sel.append(' and '.join(selg))
            sel = '(' + ') or ('.join(sel) +')'
            sel_all.append(sel)
            glycan = Molecule(k)
            glycan.set_AtomGroup(self.glycoprotein.select(sel).copy(), rootAtom = rootAtom)
            self.glycanMolecules[k] = glycan
        self.protein = self.glycoprotein.select('not ((' + ') or ('.join(sel_all) + '))').copy()
     
    def connect_all_glycans(self, G):
        """Builds a connectivity graph for molecules (not protein) in AtomGroup
            Parameters:
                G: undirected graph of connected elements
            Returns:
                names: dictionary with residue id (get_id) as keys and residue name as value
        """
        sel = self.glycoprotein.select('not protein')
        kd = KDTree(sel.getCoords())
        kd.search(1.7)
        atoms = kd.getIndices()
        atom_names = sel.getNames()
        names = {} 
        ids,rn = self.get_ids(sel)

        for a1,a2 in atoms:
            id1 = ids[a1]
            rn1 = rn[a1]
            id2 = ids[a2]
            rn2 = rn[a2]
            if id1 in G and id1 not in names:
                names[id1] = rn1

            if id2 in G and id2 not in names:
                names[id2] = rn2

            if id1 != id2:
                e = (id1, id2)
                an1 = atom_names[a1]
                an2 = atom_names[a2]
                patch = self.find_patch(an1, an2)
                G.add_edge(id1, id2, patch = patch, atoms = an1 + ':' + an2)
                #G.add_edge(*e)
                names[id1] = rn1
                names[id2] = rn2

        return names

    def get_ids(self, sel):
        segn = sel.getSegnames()
        chid = sel.getChids()
        res = sel.getResnums()
        ic = sel.getIcodes()
        rn = sel.getResnames()
        ids = []
        for s,c,r,i in zip(segn,chid,res,ic):
            ids.append(','.join([s, c, str(r), i]))
        return ids,rn

    def get_id(self, sel):
        segn = sel.getSegnames()[0]
        chid = sel.getChids()[0]
        res = sel.getResnums()[0]
        ic = sel.getIcodes()[0]
        rn = sel.getResnames()[0]
        return ','.join([segn, chid, str(res), ic]), rn


    def glycosylate(self, glycan_name, link_residue = None, link_patch = None, template_glycan_tree = {}, template_glycan = None, chain = 'X', segname = 'G1'):
        """Builds a polyglycan from a connectivity tree
        Parameters:
            glycan_name: name of poly glycan
            kwargs:
                segname: segname for new glyan (str)
                chain: chain for new
                link_residue: residue to link new glycan (AtomGroup)
                link_patch: name of patch (str) which should be used to link glycan
                template_glycan_tree: dictionary with the identified glycan that should be modified
                                key: identification of each unit in glycan
                                value: selection for residue; segname, chain, resid, icode
                template_glycan: AtomGroup containing the glycan that should be modified

        Returns:
            glycan: structure of the glycan (AtomGroup)
        """
        #if kwargs is not None:
        #    for key, value in kwargs.ashesteritems():

        resid = 1
        dummy_patch = 'DUMMY_MAN'
        prefix = ['segment', 'chain', 'resid', 'icode']
        glycan = None
        glycan_bonds = []
        built_glycan = {}
        inv_template_glycan_tree = {}
        dele_atoms = []
        if template_glycan_tree and template_glycan:
            inv_template_glycan_tree = {v: k for k, v in template_glycan_tree.items()}
            resid = template_glycan.getResnums()[-1] + 1
            chain = template_glycan.getChids()[0]
            segname = template_glycan.getSegnames()[0]
        
        if glycan_name in self.connect_topology:
            glycan_topo = self.get_connectivity_tree(glycan_name)
            sorted_units =  sorted(glycan_topo.keys(), key = len)
            for unit in sorted_units:
                new_residue = None
                if unit:
                    lunit = unit.split(' ')
                    previous = ' '.join(lunit[:-1])
                else:
                    lunit = []
                    previous = ''

                del_atom = []

                #check if residue exists
                if unit in inv_template_glycan_tree:
                    sel  = []
                    
                    for p,s in zip(prefix, inv_template_glycan_tree[unit].split(',')):
                        if s:
                            sel.append(p + ' ' + s)
                    sel = ' and '.join(sel)
                    sel_residue = template_glycan.select(sel)
                    if sel_residue.getResnames()[0] == glycan_topo[unit]:
                        built_glycan[unit] = ','.join([segname, chain, str(resid),])
                        #built_glycan[unit] = inv_template_glycan_tree[unit]
                        new_residue,missing_atoms,bonds = self.builder.add_missing_atoms(sel_residue, resid)
                        #new_residue.setResnums([resid]*len(new_residue))
                        ics = self.builder.Topology.topology[glycan_topo[unit]]['IC']
                        self.builder.build_missing_atom_coord(new_residue, missing_atoms, ics)
                    else:
                        print 'Error in connect tree!! Residue will be build de novo'

                #build first unit from DUMMY or from linked residue
                if not lunit and not new_residue:
                    if link_residue and link_patch:
                        new_residue, del_atom, bonds = self.builder.build_from_patch(link_residue, resid, glycan_topo[unit], chain, segname, link_patch)
                    else:
                        new_residue, del_atom, bonds = self.builder.build_from_DUMMY(resid, glycan_topo[unit], chain, segname, dummy_patch)
                elif previous in built_glycan and lunit:
                    patch = lunit[-1]
                    sel = []
                    for p,s in zip(prefix, built_glycan[previous].split(',')):
                        if s:
                            sel.append(p + ' ' + s)
                    sel = ' and '.join(sel)
                    previous_residue = glycan.select(sel)
                    if new_residue:
                        del_atom, b = self.builder.apply_patch(patch,previous_residue, new_residue)
                        bonds.extend(b)
                    else:
                        new_residue, del_atom, bonds = self.builder.build_from_patch(previous_residue, resid, glycan_topo[unit], chain, segname, patch)                           
                elif lunit:
                    print 'Error in connect tree!! Glycans will not be builts'
                    return [], []

                built_glycan[unit] = ','.join([segname, chain, str(resid),])
                dele_atoms += del_atom
                glycan_bonds.extend(bonds)
                
                if glycan:
                    glycan += new_residue
                else:
                    glycan = new_residue
                resid += 1
        else:
            print "Unkown Glycan"
            return [], []
            
        if dele_atoms:
            glycan = self.builder.delete_atoms(glycan, dele_atoms)
            # remove all non existing bonds
            tmp_bonds = []
#            print 'Bonds',glycan_bonds
#            print 'delele',dele_atoms
            for a1,a2 in glycan_bonds:
                if a1 in dele_atoms or a2 in dele_atoms:
                    continue
                else:
                    tmp_bonds.append((a1, a2))
            glycan_bonds = tmp_bonds

        #set serial number
        for i in range(len(glycan)):
            glycan.select('index ' + str(i)).setSerials([i+1])

        return glycan, glycan_bonds

    def get_connectivity_tree (self, glycan_name):
        """Builds connectivity tree of a glycan described in the connectivity topology
        Parameters:
            glycan_name: name of the glycan
        Returns:
            glycan_topo: dictionary with connectivy and resname of glycan
        """
        glycan_topo = {}
        units = self.connect_topology[glycan_name]['UNIT']
        for unit in units:
            rn,a,p = unit
            #need to add blank at the begining of each key except the root key, which is blank ('')
            if p:
                #key = ' ' + ' '.join(p)
                key = ' '.join(p)
            else:
                key = ''
            glycan_topo[key] = rn
        return glycan_topo

    def assign_patches(self, molecule):
        """Assignes patch names to each edge of interresidue digraph
        Parameters:
            molecule: Molecule object
        """
        G = molecule.interresidue_connectivity
        patches = {}
        atoms = nx.get_edge_attributes(G, 'atoms')
        for e in G.edges():
            a1,a2 = atoms[e].split(':')
            patch = self.find_patch(a1, a2)
            if not patch:
                print 'Unknown inter residue bond', a1, a2
            u,v=e
            G[u][v]['patch'] = patch
        
    def build_connectivity_tree(self, root_id, G):
        """Defines the connectivity within a glycan polymer
           Parameters:
            root_id: id of first residue of the glycan polymer ('segn,chid,resid')
            G: directed interresidue graph 
           Returns:
            connect_tree: dictionary of glycan connectivity
        """
        paths = nx.shortest_path(G, source = root_id)
        #put root residue in dict
        connect_tree = {}
        connect_tree[root_id] = ' '
        
        for n in paths:
            p = paths[n]
            edges = zip(p[1:],p[:-1])
            value = []
            for e2,e1 in edges:
                value += [G[e1][e2]['patch']]
            connect_tree[n] = ' ' .join(value)
        connect_tree.items().sort(key=lambda id:len(id[1]))
        return connect_tree

    def build_connect_topology(self, molecule):
        """Builds connectivity topology
        Parameters:
            molecule: Molecule object
        """
        G = molecule.interresidue_connectivity
        connect_tree = self.build_connectivity_tree(molecule.rootRes, G) 
        
        target=[]
        for r in connect_tree.keys():
            target.append(G.node[r]['resname'] + ' ' + connect_tree[r])
        len_target=len(target)
        residue = {}
        residue['UNIT'] =  target
        residue['#UNIT'] = len_target
        return residue


    def identify_glycan(self, molecule):
        """Identifies glycan name
        Parameters:
            molecule: Molecule object
        """
        G = molecule.interresidue_connectivity
        connect_tree = self.build_connectivity_tree(molecule.rootRes, G) 
        
        target=[]
        for r in connect_tree.keys():
            target.append(G.node[r]['resname'] + ' ' + connect_tree[r])
        len_target=len(target)
        
        for gk in self.glycan_keys:
            if len(set(target) & set(self.glycan_keys[gk])) == len_target and len(set(self.glycan_keys[gk])) == len_target:
                break
            else:
                gk = ''
        
        if gk:
            molecule.id = gk
        else:
            print 'Unknown glycan'
            gk = ''
            molecule.id = ''
        return gk
    
    def write_connectivity_topology(self, glycan_name, connect_tree, fileName):
        """Write connectivity tree of a glycan
            Parameters:
                glycan_name: name of glycan (str)
                fileName: path to output file (str)
        """
        file = open(fileName,'w')                  # open the file
        file.write('RESI ' + glycan_name + '\n')
        units = connect_tree.items()
        units.sort(key=lambda id:len(id[1]))
        for unit in units:
            file.write('UNIT ' + unit[0].get_resname() + ' C1 ' + unit[1] + '\n')
        file.close() 


    def find_patch(self, atom1, atom2):
        """Finds patch that connects two atoms
        Currently based only on atom names
        Parameters:
            atom1: name of first atom (str)
            atom2: name of second atom (str)
        Returns:
            patch name
        """
        key = '-'.join(sorted([atom1, atom2]))
        if key in self.builder.Topology.atomnames_to_patch:
            return self.builder.Topology.atomnames_to_patch[key]
        else:
             return ''
             
    def assign_atom_type(self, molecule, connect_tree=None):
        """ Returns a dictionary of atom types for a given molecule.
        Parameters:
           molecule: Molecule instance
           connect_tree: connectivity tree generated with "build_connectivity_tree". If not provided will be created in function
        Returns:
           atom_types: dictionary of atom types
                        key: atom serial number
                        value: atom name, atom type, charge, id, element
        """
        if not connect_tree:
           connect_tree = self.build_connectivity_tree(molecule.rootRes, molecule.interresidue_connectivity)
        inv_connect_tree = {v: k for k, v in connect_tree.items()}
        sorted_units =  sorted(inv_connect_tree.keys(), key = len)
        atom_type = {}
        masses = self.builder.Topology.topology['MASS']
        for unit in sorted_units:
            current = inv_connect_tree[unit]
            cur_res = molecule.get_residue(current)
            cur_rn = cur_res.getResnames()[0]
           
            cur_atoms = [a.strip() for a in cur_res.getNames()]
            cur_serial =  cur_res.getSerials()
            cur_atom_serial = {}
            for a,s in zip(cur_atoms, cur_serial):
                cur_atom_serial[a] = s
            atoms = self.builder.Topology.get_atoms(cur_rn)
            for a in atoms:
                an = a[0].strip()
                if an in cur_atom_serial:
                    atom_type[cur_atom_serial[an]] = {'name': an, 'type': a[1], 'charge': a[2], 'id': current+','+an, 'element': masses[a[1]][-1]}
            # correct atom type 
            lunit = unit.split(' ')
            patch = lunit[-1]
            if not patch:
                continue
            previous = inv_connect_tree[' '.join(lunit[:-1])]
            pre_res = molecule.get_residue(previous)
            
            pre_atoms = [a.strip(' ') for a in pre_res.getNames()]
            pre_serial =  pre_res.getSerials()
            pre_atom_serial = {}

            for a,s in zip(pre_atoms, pre_serial):
                pre_atom_serial[a] = s

            atoms = self.builder.Topology.patches[patch]['ATOM']
            for a in atoms:
                an = a[0].strip(' ')
                if an[0] == '1':
                    if an[1:] in pre_atom_serial: 
                        atom_type[pre_atom_serial[an[1:]]] = {'name': an[1:], 'type': a[1], 'charge': a[2], 'id':previous+','+an[1:], 'element': masses[a[1]][-1]}
                elif an[0] == '2':
                    if an[1:] in cur_atom_serial:
                        atom_type[cur_atom_serial[an[1:]]] = {'name': an[1:], 'type': a[1], 'charge': a[2], 'id':current+','+an[1:], 'element': masses[a[1]][-1]}
        return atom_type


#####################################################################################
#                                Sampler                                            #
#####################################################################################
    
class Sampler():
    """Class to sample conformations, based on a 
    """
    def __init__(self, molecules, envrionment, dihe_parameters, clash_dist = 1.5):
        """ 
        Parameters
            molecules: list of Molecules instances
            environment: AtomGroup that will not be samples (e.g. protein, membrane, etc.)
            dihe_parameters: dictionary of parameters for dihedrals from CHARMMParameters
            clash_dist = threshold for defining a clash (A)
        """
        self.molecules = molecules
        self.environment = envrionment
        self.clash_dist = clash_dist 
        self.energy = {}
        self.energy_lookup = []
        self.nbr_clashes = []
        mol_id =  0
        for molecule in self.molecules:
            types = nx.get_node_attributes(molecule.connectivity, 'type')
            lookup = []
            self.nbr_clashes.append(self.count_total_clashes(mol_id))
            for dihe in molecule.torsionals:
                atypes = []
                for d in dihe:
                    atypes.append(types[d])
                k1 = '-'.join(atypes)
                atypes.reverse()
                k2 = '-'.join(atypes)
                if k1 in self.energy:
                    lookup.append(k1)   
                    continue
                if k2 in self.energy:
                    lookup.append(k2)   
                    continue

                if k1 in dihe_parameters:
                    k = k1
                elif k2 in dihe_parameters:
                    k = k2
                else:
                    print 'Missing parameters for ' + k1
                    print 'This dihedral will be skipped'
                    lookup.append('skip')
                    continue
                par_list = dihe_parameters[k]
                self.energy[k] = self.compute_inv_cum_sum_dihedral(par_list)
                lookup.append(k)
            self.energy_lookup.append(lookup)
            mol_id += 1
        
    def compute_inv_cum_sum_dihedral(self, par_list, n_points = 100):
        """Computes the an interpolation of the inverse transform sampling of a CHARMM
        dihedral term: sum(k*[1-cos(n*phi -d)]).
        This allows to map a uniform distribution [0:1[ to the dihedral energy distribution
        Parameters:
            par_list: list of parameters [k, n, d] 
            n_points: number of points used for computing the energy function
        Returns:
            inv_cdf: interpolate object
        """
        phi = np.linspace(-180., 180., n_points)
        energy = np.zeros(phi.shape)
        for p in par_list:
            k,n,d = p
            energy += k*(1-np.cos((n*phi -d)*np.pi / 180.))
        energy = np.max(energy) - energy
        energy = energy / np.sum(energy[1:]*np.diff(phi))
        cum_values = np.zeros(energy.shape)
        cum_values[1:] = np.cumsum(energy[1:]*np.diff(phi))
        inv_cdf = interp1d(cum_values, phi)
        return inv_cdf
        
    
    def count_total_clashes(self, mol_id):
        """Counts the total number of clashes in the system. 
        Performed in two steps:
                - count the number of clashes within a molecules (KDTree)
                - count the number of clashes between the molecule and it's environment (Grid)
        Parameters:
            mol_id: id of the molecule
        """
        nbr_clashes,clashes = self.count_self_clashes(mol_id)
        nbr_clashes += self.count_environment_clashes(mol_id)
        return nbr_clashes

    def count_self_clashes(self, mol_id):
        """Counts the number of clashes for a molecule
        KDTree based
        Parameters:
            molecule: Molecule type object
        Returns
            nbr_clashes: the number of clashes
            clashes: list of clashing atom serial numbers
        """
        molecule = self.molecules[mol_id]
        kd = KDTree(molecule.atom_group.getCoords())
        kd.search(self.clash_dist)
        atoms = kd.getIndices()
        G = molecule.connectivity
        nbr_clashes = 0
        clashes = []
        for a1,a2 in atoms:
            if not G.has_edge(a1+1, a2+1):
                nbr_clashes += 1
                clashes += (a1+1,a2+1)
        return nbr_clashes,clashes
        
    def count_environment_clashes(self, mol_id):
        """Counts the number of a molecule and its environment
        Parameters:
            mol_id: id of molecule
        """
        XA = self.molecules[mol_id].atom_group.getCoords()
        nbr_clashes = 0
        if self.environment:
            XB = self.environment.getCoords()
            Y = distance.cdist(XA, XB, 'euclidean')
            nbr_clashes += np.sum(Y < self.clash_dist)

        for mol in np.arange(len(self.molecules)):
            if mol == mol_id:
                continue
            XB = self.molecules[mol].atom_group.getCoords()
            Y = distance.cdist(XA, XB, 'euclidean')
            nbr_clashes += np.sum(Y < self.clash_dist)
        return nbr_clashes 
    
    def remove_clashes(self, temp = 310, n = 1000, max_torsional = 0.3):
        """Monte Carlo sampling to remove clashes. The number of torsional anlges will decrease for max_torsional (fractions of total angles)
        to one single angle.
        Parameters:
            temp: temperature for 
            n: number of cycles
            max_torsional: maximal fractions of trosional angles that will be sampled.
        """
        R      = 1.9872156e-3    # Gas constant  kcal/mol/degree
        beta   = 1 / ( R * temp)
        
        i = 0
        n_mol = len(self.molecules)

        while i < n:
            mol_id = np.random.randint(n_mol)
            clashes =  self.nbr_clashes[mol_id]
            torsionals = self.molecules[mol_id].torsionals
            nbr_torsionals = int(max_torsional * len(torsionals))
            t_ids  = np.random.randint(len(torsionals), size = nbr_torsionals)
            r_s = np.random.random_sample((nbr_torsionals,))
            lenergy = []
            for t_id in t_ids:
                e = self.energy_lookup[mol_id][t_id]
                if e == 'skip':
                    break
                lenergy.append(e)
            if e == 'skip':
                print 'skipping'
                continue
            i += 1
            dtheta = []
            for e,r,t_id in zip(lenergy, r_s, t_ids):
                theta = self.energy[e](r)
                ctheta = self.molecules[mol_id].rotate_bond(int(t_id), theta, absolute = True)
                dtheta.append(ctheta-theta)
            clashes_new = self.count_total_clashes(mol_id)
            print clashes_new,clashes
            if clashes_new >= clashes:
                if np.exp( -beta *(clashes_new - clashes) ) <  np.random.random_sample():
                    # set back initial configuration
                    for t_id, t in zip(t_ids, dtheta):
                        self.molecules[mol_id].rotate_bond(int(t_id), t)
                    continue

            self.nbr_clashes[mol_id] = clashes_new
            max_torsional -= max_torsional/n


#####################################################################################
#                                Sampler                                            #
#####################################################################################
    
class Drawer():
    """Class to sample conformations, based on a 
    """
    def __init__(self):
        self.Symbols = {}
        self.Symbols['GLC'] = ['circle', 'blue']
        self.Symbols['MAN'] = ['circle', 'green']
        self.Symbols['BMA'] = ['circle', 'green']
        self.Symbols['GAL'] = ['circle', 'yellow']
        self.Symbols['NAG'] = ['square', 'blue']
        self.Symbols['NGA'] = ['square', 'yellow']
        self.Symbols['FUC'] = ['triangle', 'red']
        self.Symbols['SIA'] = ['rhombus', 'purple']

        self.Colors = {}
        self.Colors['blue'] = np.array([0, 0, 204])
        self.Colors['green'] = np.array([0, 204, 0])
        self.Colors['yellow'] = np.array([255, 255, 102])
        self.Colors['red'] = np.array([204, 0, 0])
        self.Colors['purple'] = np.array([167, 13, 206])
        self.side = .25
        self.scaling = 5. 

    def get_shape(self, name, pos, direction = 1, axis = 0):
        """Returns a Matplotlib patch 
            name: name of shape
            pos: vector defining where to draw object
            direction: direction for drawing glycan: 1 positive y; -1 negative y
            axis: axis along which to draw glycan. 0 for x and 1 for y.
        """
        pos = np.array(pos)
        if name in self.Symbols:
            shape,color = self.Symbols[name]
            color = self.Colors[color]/255.
            if shape == 'circle':
                patch = Circle(pos, self.side)
                return patch,color
            if shape == 'square':
                return Rectangle(pos-self.side, 2*self.side, 2*self.side),color
            if shape == 'triangle':
                return Polygon([(pos[0]-self.side, pos[1]), (pos[0]+self.side, pos[1]-self.side), (pos[0]+self.side, pos[1]+self.side)]),color
            if shape == 'rhombus':
                if axis:
                    angle = 135 if direction == -1 else -45.
                else:
                    angle = -135 if direction == -1 else 45.
                return Rectangle(pos, 2*self.side, 2*self.side, angle=angle),color

    def draw_protein_fragment(self, pos = [0., 0.], length = 10, protein_color = [.7, .7, .7], sequon_color = [.5, .5, .5], ax = None, axis = 0):
        """
        Parameters:
            pos: position of sequon
            length: length of fragment
            color: color of protein
            axis: axis along which to draw glycan. 0 for x and 1 for y.
        Returns:
            ax: axis of plot
        """
        if not ax:
            fig, ax = plt.subplots()
    	shape = [0, 0]
    	shape[axis] = length/self.scaling
    	shape[np.mod(axis+1, 2)] = 2*self.side
        p_fragment = np.array([0., 0.])
        p_fragment[axis] = -shape[axis]/2
        p_fragment[np.mod(axis+1, 2)] = -(self.side * 3.5) 
        protein = Rectangle(p_fragment - self.side, shape[0], shape[1])
        patches = []
        colors = []
        patches.append(protein)
        colors.append(protein_color)
        new_pos = np.array([0., 0.])
        new_pos[axis] = pos[axis] 
        new_pos[np.mod(axis+1, 2)] = pos[np.mod(axis+1, 2)] -(self.side * 3.5) 
        patches.append(Rectangle(new_pos - self.side, 2*self.side, 2*self.side)) 
        colors.append(sequon_color)

        p = PatchCollection(patches, zorder = 3) 
        p.set_edgecolors([.2, .2, .2])
        p.set_linewidth(2)
        p.set_facecolors(colors)
        ax.add_collection(p)
        
        return ax
        

    def draw_glycoprotein(self, length, start_resid, sequons, pos = [0., 0.], protein_color = [.7, .7, .7], sequon_color = {}, ax = None, axis = 0, trees = None, names =  None):
        """
        Parameters:
            length: length of the protein
            start_resid: Index of frist residue
            sequons: position of sequons
            pos: starting position of sequence
            protein_color: color of protein (list)
            sequon_color: dictionary of colors for sequons, with sequon id as key and color code as value. Default color [.5, .5, .5]
            axis: axis along which to draw glycan. 0 for x and 1 for y.
            trees: dictionary with sequon as key and list of root and glycan graph as value
            names: dictionary with glycan nodes as keys and resnames as values
        Returns:
            ax: axis of plot
        """
        if not ax:
            fig, ax = plt.subplots()
    	shape = [0, 0]
    	shape[axis] = length/self.scaling
    	shape[np.mod(axis+1, 2)] = 2*self.side
        protein = Rectangle(np.array(pos) - self.side, shape[0], shape[1])
        patches = []
        colors = []
        patches.append(protein)
        colors.append(protein_color)
        direction = -1
        sequons = alphanum_sort(sequons)
        for seq in sequons:
            s,c,rr,i = seq.split(',')
            if seq in sequon_color:
                color = sequon_color[seq]
            else:
                color = [.5, .5, .5]
            r = rr
            if i:
                r += ord(i.upper()) - 64 
            r = (float(r)-start_resid)/self.scaling
            new_pos = np.array([0., 0.])
            new_pos[axis] = pos[axis] + r
            new_pos[np.mod(axis+1, 2)] = pos[np.mod(axis+1, 2)] 
            patches.append(Rectangle(new_pos - self.side, 2*self.side, 2*self.side)) 
            colors.append(color)

            new_pos[np.mod(axis+1, 2)] = direction * (self.side * 3.5) 
            if axis:
                h_align = 'right' if direction == -1 else 'left'
                v_align = 'center' 
                rotation = 'horizontal'
            else:
                h_align = 'center' 
                v_align = 'top' if direction == -1 else 'bottom'
                rotation = 'vertical'
            ax.text(new_pos[0], new_pos[1], rr+i,
                    verticalalignment = v_align, horizontalalignment = h_align, 
                    rotation = rotation)
            direction *= -1

            if seq in trees:
                r_pos = [0, 0]
                r_pos[axis] = r 
                r_pos[np.mod(axis+1, 2)] = direction*(self.side * 3.5)
                root, tree =  trees[seq]
                self.draw_tree(tree, root, names, r_pos, direction = direction, ax= ax, axis=axis)


        p = PatchCollection(patches, zorder = 3) 
        p.set_edgecolors([.2, .2, .2])
        p.set_linewidth(2)
        p.set_facecolors(colors)
        ax.add_collection(p)
        
        return ax


    def draw_all_trees(self, trees, start_resid, names, ax = None, axis = 0):
        """Draws all glycans in a protein
        Parameters:
            trees: dictionary with sequon as key and list of root and glycan graph as value
            start_resid: first residue number of protein sequence
            names: dictionary with glycan nodes as keys and resnames as values
            ax: axis handle
            axis along which to draw glycan. 0 for x and 1 for y.
        Returns:
            ax: axis of plot
        """
        if not ax:
            fig, ax = plt.subplots()
        direction = 1
        keys = alphanum_sort(trees.keys())
        for k in keys:
            s,c,r,i = k.split(',')
            if i:
                r += ord(i.upper()) - 64
            r_pos = [0, 0]
            r_pos[axis] = (float(r)-start_resid)/self.scaling
            r_pos[np.mod(axis+1, 2)] = direction*(self.side * 3.5)
            root, tree =  trees[k]
            self.draw_tree(tree, root, names, r_pos, direction = direction, ax= ax, axis=axis)
            direction *= -1

        plt.axis('equal')
        plt.axis('off')
        fig = plt.gcf()
        fig.tight_layout()
        return ax
    

    def draw_tree(self, tree, root, names, root_pos = [0, 0], direction = 1, ax = None, axis = 0):
        """Draws a tree representation of a glycan
        Parameters:
            tree: netwrokX graph representation of glycan
            root: name of root node
            names: dictionary with nodes as keys and resnames as values
            root_pos: position of the root node (default [0,0])
            direction: direction for drawing glycan: 1 positive y; -1 negative y
            axis: axis along which to draw glycan. 0 for x and 1 for y.
            ax: axes handle
         Returns:
            ax: axis of plot
        """
        if not ax:
            fig, ax = plt.subplots()

        pos = self.compute_positions(tree, root, root_pos, direction, axis)
        #draw edges
        for edge in tree.edges():
            x1,y1 = pos[edge[0]]
            x2,y2 = pos[edge[1]]
            ax.add_line(mlines.Line2D([x1, x2],[y1, y2], lw=2., c=[.2, .2, .2]))
        #draw nodes
        patches = []
        colors = []
        for node in tree.nodes():
            patch,color = self.get_shape(names[node], pos[node], direction, axis)
            patches.append(patch)
            colors.append(color)
        p = PatchCollection(patches, zorder = 3) 
        p.set_edgecolors([.2, .2, .2])
        p.set_linewidth(2)
        p.set_facecolors(colors)
        ax.add_collection(p)

        return ax

    def tree_to_text(self, tree, root, names, glyPAC = '', visited = []):
        """Returns UIPAC like string of a glycan.
        []: branch
        (): connectivity
        Example: Man3 == NAG(14bb)NAG(14bb)BMA[(16ab)MAN](13ab)MAN
        Parameters:
            tree: directed graph of a rooted tree
            root: id of root vertex
        Returns:
            glyPAC: string representing the glycan
        """
        if not glyPAC:
            glyPAC = names[root] 
        while 1:           
            visited.append(root)
            successors = [n for n in tree.neighbors(root) if n not in visited]
            if not successors:
                break
            branches =  len(successors)
            for s in successors:
                patch = tree[root][s]['patch']
                if branches > 1:
                    glyPAC += '[(' + patch + ')' + names[s] 
                    glyPAC = self.tree_to_text(tree, s, names, glyPAC, visited)
                    glyPAC += ']'
                    branches -= 1
                else:
                    glyPAC += '(' + patch + ')' + names[s] 
            root = s
	return glyPAC
    
    def compute_positions(self, tree, root, root_pos = [0, 0], direction = 1, axis = 0 ):
        """Computes positions of vertices for drawing a tree
        1-2--3-8
           \-4--5
              \-6-7
        Parameters:
            tree: directed graph of a rooted tree
            root: id of root vertex
            root_pos: position of the root vertex (default: 0,0)
            direction: -1/1 direction to draw tree
            axis: axis along which to draw glycan. 0 for x and 1 for y.
        Returns:
            pos: dictionary with nodes as keys and positions as values
        """
        nodes = list(nx.dfs_preorder_nodes(tree, root))
        successors = nx.dfs_successors(tree, root)
        nodes.reverse()
        branches = {}
        for n in nodes:
            if n in successors:
                branches[n] = len(successors[n])
                for s in successors[n]:
                    if branches[s] > 1:
                        branches[n] += branches[s] - 1
            else:
                branches[n] = 1

        nodes.reverse()
        pos = {}

        for n in nodes:
            if not n in pos:
                pos[n] = root_pos

            if n in successors:
                nbr_branches = []
                for s in successors[n]:
                    nbr_branches.append(branches[s])
                c = 0
                for s in np.argsort(nbr_branches):
                    x = pos[n]
                    x_new = [0, 0]
                    x_new[axis] = x[axis] + c
                    x_new[np.mod(axis+1, 2)] = x[np.mod(axis+1, 2)] + direction
                    pos[successors[n][s]] = x_new
                    c += nbr_branches[s]

        return pos 



