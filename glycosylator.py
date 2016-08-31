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

SELF_BIN = os.path.dirname(os.path.realpath(sys.argv[0]))
#sys.path.insert(0, SELF_BIN + '/support')



#####################################################################################
#						        Support functions             						#
#####################################################################################

def readLinesFromFile(fileName):
    """Reads all lines in a file
    Parameters:
        fileName: path to file
    Returns:
        lines: list with all the lines in a file
    """
    file = open(fileName,'r') 				 # open the file
    lines = file.readlines() 				 # read all the lines in the file to the list "lines"
    file.close() 						         # close the file
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

#####################################################################################
#						        Topology functions             						#
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
        masses={}
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
        self.topology['MASS'] = copy.copy(masses)        

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
        self.read_parameters

    def read_parameters(fileIn):
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
                    eval('read_'+prm_type+'(line, prm)')
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
                prm[key]=map(float, dihe[4:])
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
#						            Molecule    									#
#####################################################################################
class Molecule:
    """Class for saving a molecule
    Attributes:
            molecule: Prody AtomGroup 
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
            molecule: AtomGroup 
            rootAtom: serial number of atom used as root for graphs
            bonds: list of all bonds
            angles: list of all angles
            dihedrals: list of all dihedrals
            connectivity: graph for connectivity of molecule (bonds)
            directed_connectivity: directed acyclique graph of molecule
            cycle_id: dictionary where keys are the serial number of atom in cacles and values the corresponding cycle in directed_connectivity
            torsionals: dihedral that can rotate (i.e. not in cycles)
            bond_length: dictionary of bond distance used to guess bonds. Keys are sorted by alphabetical order
        """
        self.molecule = AtomGroup(name)
        self.chain = chain
        self.segn = segn
        self.rootAtom = 1
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.connectivity = nx.Graph()
        self.directed_connectivity = nx.DiGraph()
        self.interresidue_connectivity = nx.DiGraph()
        self.cycle_id = {}
        self.torsionals = []

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
        writePDB(filename, self.molecule.select(selection))

    def read_molecule_from_PDB(self, filename, update_bonds = True, **kwargs):
        """Initialize molecule from a PDB file
        Parameters:
            filename: path to PDB file
            update_bonds: guess bonds, angles, dihedrals and connectivity based on the distance between atoms in PDB
            **kwargs: any of the following which allows the selection of a subset of the PDB file
                subset: selection of a subset of the PDB file
                model: model number (int)
                chain: chain id (str)
        Initializes:
            molecule
            chain: chain id
            segn: segment name
            connectivity: bonds, angles, dihedrals and graph 
        """
        PDBmolecule = parsePDB(filename, **kwargs)
        chain = set(PDBmolecule.getChids())
        segn= set(PDBmolecule.getSegnames())
        if len(chain) == 1 and len(segn) == 1:
            self.molecule = PDBmolecule
            self.chain = chain
            self.segn = segn
            if update_bonds:
                self.update_connectivity()
        else:
            print "several chains are present in PDB. Please select only one molecule"
    
    def get_chain(self):
        return self.chain
    
    def get_segname(self):
        return self.segn
    
    def get_residue(self, resid):
        return self.molecule[self.chain][resid]

    def update_connectivity(self):
        """Updates all the connectivity (bond, angles, dihedrals and graphs)
        """
        self.guess_bonds()
        self.guess_angles()
        self.guess_dihedrals()
        self.update_graphs()

    def guess_bonds(self, default_bond_length = 1.6):
        """Searches for all bonds in molecule
        Parameters:
            default_bond_length: maximum distance between two connected heavy atoms (Angstrom), if not present in bond_length dictionary
        """
        self.connectivity = nx.Graph()
        for a in self.molecule:
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
            neighbors =  self.molecule.select(sel)
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
                atoms.append(self.molecule.select('serial ' + str(node)))
                if node in self.cycle_id:
                    key = self.cycle_id[node]
                    if key not in self.directed_connectivity:
                        self.directed_connectivity.add_node(key, iscycle = True, cycle_id = map(int, self.cycle_id[node].split('-')))
                    directed_edge.append(key)
                else:
                    self.directed_connectivity.add_node(node, iscycle = False, cycle_id = [])
                    directed_edge.append(node)
                    
            self.directed_connectivity.add_edge(directed_edge[0], directed_edge[1])
            a1,a2 = atoms
            if a1.getResnums()[0] != a2.getResnums()[0]:
                r1 = a1.getSegnames()[0] + ',' + a1.getChids()[0] + ',' + str(a1.getResnums()[0])
                r2 = a2.getSegnames()[0] + ',' + a2.getChids()[0] + ',' + str(a2.getResnums()[0])
                self.interresidue_connectivity.add_node(r1)
                self.interresidue_connectivity.add_node(r2)
                self.interresidue_connectivity.add_edge(r1, r2, patch = '', atoms = a1.getNames()[0] + ':' + a2.getNames()[0])

    def define_torsionals(self, hydrogens=True):
        """Builds a list with all the torsional angles that can rotate
        Parameters:
            hydrogens: include torsional involving terminal hydrogens        
        Initializes:
            torsionals: a list of serial number of atom defining a trorsional angle (quadruplet) 
        """
        self.torsionals = []
        cycles = nx.cycle_basis(self.connectivity, self.rootAtom)
        for dihe in self.dihedrals:
            in_cycle = 0
            d_dihe = []
            for a in dihe:
                if a in self.cycle_id:
                    in_cycle += 1
                    a = self.cycle_id[a]
                d_dihe.append(a)

            #check that quadruplet is not in cycle
            if in_cycle > 1:
                continue
            #check direction of dihedral
            if self.directed_connectivity.has_edge(d_dihe[1], d_dihe[2]):
                pass
            elif self.directed_connectivity.has_edge(d_dihe[2], d_dihe[1]):
                dihe.reverse() 
            else:
                continue
            #check if already in torsionals list
            exists = False
            if self.torsionals:
                for t in self.torsionals:
                    if dihe[2] == t[2] and dihe[3] == t[3]:
                        exists = True
                        break

            if exists:
                continue

            self.torsionals.append(dihe)

#####################################################################################
#								Builders											#
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


    def new_residue(self, resid, resname, chain, segname):
        """Initializes a residue from scratch
        Parameters:
            resid: residue id (int)
            resname: name of residue (str)
            chain: chain id (char)
            segn: segname (str)
        Returns:
            residue: AtomGroup with all atoms initialized (from Topology)
            atoms_name: name of of all the atoms in residue
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
        i = 0
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
        return residue, atoms_name

    def copy_atom(self, src_atom, dst_atom):
        """copies all the attributes of one atom to another
        Parameters:
            src_atom: original atom
            dst_atom: copy atom
        """
        dst_atom.setCoords(src_atom.getCoords())
        dst_atom.setNames(src_atom.getName())
        dst_atom.setResnums(src_atom.getResnum())
        dst_atom.setResnames(src_atom.getResname())
        dst_atom.setChids(src_atom.getChid())
        dst_atom.setSegnames(src_atom.getSegname())
        dst_atom.setOccupancies(src_atom.getOccupancy())
        dst_atom.setBetas(src_atom.getBeta())
        dst_atom.setSerials(src_atom.getSerial())
        dst_atom.setIcodes(src_atom.getIcode())
        dst_atom.setElements(src_atom.getElement())
        dst_atom.setAltlocs(src_atom.getAltloc())
        
    def add_missing_atoms(self, residue):
        """Add all missing atoms to a ProDy residue from topology
        Parameters:
            residue: ProDy residue (AtomGroup)
        Retruns:
            missing_atoms: list of missing atom names
        """
        complete_residue,atoms = self.new_residue(residue.getResnum(), residue.getResname(), residue.getChid(), residue.getSegname())
        missing_atoms = []
        atoms_in_residue = residue.getNames()
        for a in atoms:
            if a in atoms_in_residue:
                atom = residue.getAtom(a)
                catom = complete_residue.select('name ' + a)
                self.copy_atom(atom, catom)
            else:
                missing_atoms.append(a)
        
        return complete_residue, missing_atoms


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
        """
        denovo_residue, missing_atoms = self.new_residue(resid, resname, chain, segname)
        ics = self.Topology.patches[patch]['IC']
        #patch_atoms = sorted(set([atom.replace('*', '')[1:] for ic in ics for atom in ic[0:4] if atom.replace('*', '')[0]=='2']))
        patch_atoms = sorted(set([atom.replace('*', '') for ic in ics for atom in ic[0:4] if atom.replace('*', '')[0]=='2']))
        self.build_patch_missing_atom_coord(link_residue, denovo_residue, patch_atoms, ics)
        missing_atoms = [a for a in missing_atoms if '2' + a not in patch_atoms]
        ics = self.Topology.topology[resname]['IC']
        self.build_missing_atom_coord(denovo_residue, missing_atoms, ics)
        
        dele_atoms = self.dele_atoms(patch, link_residue, denovo_residue)
        return denovo_residue, dele_atoms



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
        if residue2:
            segn2 = residue2.getSegnames()[0] 
            chid2 = residue2.getChids()[0]
            resi2 = str(residue2.getResnums()[0])

        for a in atoms:
            if a[0] == '1':
                dele_atoms.append((segn1, chid1, resi1, a[1:]))
            if a[0] == '2':
                if residue2:
                    dele_atoms.append((segn2, chid2, resi2, a[1:]))
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
        for a in dele_atoms:
            segn, chid, resi, atom_name = a
            sel.append('(segment %s and chain %s and resid %s and name %s)' % (segn, chid, resi, atom_name ))
        sel = 'not (' + ' or '.join(sel) + ')'
        aa = molecule.select(sel).copy()
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
        """
        dummy_residue, dummy_atoms = self.new_residue(0, 'DUMMY', 'D', 'DUM')
        counter = 0
        for a in dummy_atoms:
            dummy_residue.select('name ' + a).setCoords([dummy_coords[counter]])
            counter += 1
        denovo_residue, dele_atoms = self.build_from_patch(dummy_residue, resid, resname, chain, segname, dummy_patch)
        del dummy_residue
        return denovo_residue

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

#####################################################################################
#								Glycosylator										#
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
            glycan_keys: dictionary 
        """
        self.builder = MoleculeBuilder(topofile, paramfile)
        self.connect_topology = {}
        self.glycan_keys = {}

    def read_connectivity_topology(self, connectfile):
        """Parse file defining the topology of connectivity trees.
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

    def build_keys(self):
        self.glycan_keys = {}
        for res in self.connect_topology:
            key = [r[0]+' '+' '.join(r[2]) for r in connect_topology[res]['UNIT']]
            self.glycan_keys[res] = key
        return self.glycan_keys

    def read_unit(self, unit, residue):
        if len(unit)>2:
            residue['UNIT'].append([unit[1], unit[2], unit[3:]])
        else:
            residue['UNIT'].append([unit[1], [], []])
    
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


    def glycosylate(self, glycan_name, link_residue = None, patch = '', template_glycan_tree = {}, template_glycan = None, chain = 'X', segname = 'G1'):
        """Builds a polyglycan from a connectivity tree
        Parameters:
            glycan_name: name of poly glycan
            kwargs:
                segname: segname for new glyan (str)
                chain: chain for new
                link_residue: residue to link new glycan (AtomGroup)
                patch: name of patch (str)
                template_glycan_tree: dictionary with the identified glycan that should be modified
                                key: identification of each unit in glycan
                                value: selection for residue; segname, chain, resid
                template_glycan: AtomGroup containing the glycan that should be modified

        Returns:
            glycan: structure of the glycan (AtomGroup)
        """
        #if kwargs is not None:
        #    for key, value in kwargs.iteritems():
        


        resid = 1
        dummy_patch = 'DUMMY_MAN'
        glycan = None
        built_glycan = {}
        inv_template_glycan_tree = {}
        dele_atoms = []
        if template_glycan_tree and template_glycan:
            inv_template_glycan_tree = {v: k for k, v in template_glycan_tree.items()}
            resid = template_glycan.getResnums()[-1]

        if glycan_name in self.connect_topology:
            glycan_topo = self.get_connectivity_tree(glycan_name)
            sorted_units =  sorted(glycan_topo.keys(), key = len)
            for unit in sorted_units:
                new_residue = None
                if unit in inv_template_glycan_tree:
                    s,c,r = inv_template_glycan_tree[unit].split(',')
                    sel = '(segment %s) and (chain %s) and (resid %s)' % (s, c, r)
                    sel_residue = template_glycan.select(sel)
                    if sel_residue.getResnames()[0] == glycan_topo[unit]:
                        built_glycan[unit] = inv_template_glycan_tree[unit]
                        new_residue = sel_residue.copy()
                elif not new_residue:
                    if link_residue and patch:
                       print 'pass'   
                    else:
                        lunit = unit.split(' ')
                        previous = ' '.join(lunit[:-1])
                        #build first unit form dummy
                        if len(lunit) == 1:
                            new_residue = self.builder.build_from_DUMMY(resid, glycan_topo[unit], chain, segname, dummy_patch)
                            built_glycan[unit] = ','.join([segname, chain, str(resid)])
                        elif previous in built_glycan:
                            patch = lunit[-1]
                            s,c,r = built_glycan[previous].split(',')
                            sel = '(segment %s) and (chain %s) and (resid %s)' % (s, c, r)
                            previous_residue = glycan.select(sel)
                            new_residue, del_atom = self.builder.build_from_patch(previous_residue, resid, glycan_topo[unit], chain, segname, patch)
                            built_glycan[unit] = ','.join([segname, chain, str(resid)])
                            dele_atoms += del_atom

                if glycan:
                    glycan += new_residue
                else:
                    glycan = new_residue
                resid += 1
        if dele_atoms:
            glycan = self.builder.delete_atoms(glycan, dele_atoms)     
        return glycan

    def get_connectivity_tree (self, glycan_name):
        """Builds connectivity tree of a glycan described in the connectivity topology
        Parameters:
            glycan_name: name of the glycan
        Returns:
            glycan_topo: dictionary with connectivy and resname of glycon
        """
        glycan_topo = {}
        units = self.connect_topology[glycan_name]['UNIT']
        for unit in units:
            rn,a,p = unit
            #need to add blank at the begining of each key except the root key, which is blank ('')
            if p:
                key = ' ' + ' '.join(p)
            else:
                key = ''
            glycan_topo[key] = rn
        return glycan_topo

    def write_connectivity_topology(self, glycan_name, fileName):
        """Write connectivity tree of a glycan
            Parameters:
                glycan_name: name of glycan (str)
                fileName: path to output file (str)
        """
        file = open(fileName,'w') 				 # open the file
        file.write('RESI ' + glycan_name + '\n')
        units = self.connect_tree.items()
        units.sort(key=lambda id:len(id[1]))
        for unit in units:
            file.write('UNIT ' + unit[0].get_resname() + ' C1 ' + unit[1] + '\n')
        file.close() 
    
    def assign_unit_id(self, structure, rootAtom, bond_length, atomnames_to_patch):
        """Identify glycan based on connectivity
            Parameters:
                structure: glycan structure (AtomGroup)
                rootAtom: first atom of the glycan polymer (AtomGroup)
                bond_length: threshold for bond length (float)
                atomsnames_to_patch: dictionary of connecting atoms to patch (Topology)
            Returns:
                connect_tree: dictionary of connections
        """
        connectivity = self.find_interresidue_bond(structure, bond_length, atomnames_to_patch)
        connect_tree = {}
        connect_tree = build_connectivity_tree(rootAtom, connectivity, connect_tree)
        connect_tree[rootAtom.get_parent()] = ''
        connect_tree.items().sort(key=lambda id:len(id[1]))
        return connect_tree
    

    def build_connectivity_tree(self, root_atom, structure, bond_length = 1.55, connect_tree = {}):
        """Defines the connectivity within a glycan polymer
            Parameters:
                root_atom: first atom of the glycan polymer (AtomGoup)
                structure: AtomGroup containing structure 
                connect_tree: dictionary of glycan connectivity
            Returns:
                connect_tree: dictionary of glycan connectivity
        """
        #put root residue in connectivity tree
        if not connect_tree:
            root_resid = root_atom.getResnums()[0]
            root_chain = root_atom.getChids()[0]
            root_segname = root_atom.getSegnames()[0]
            connect_tree[','.join([root_segname, root_chain, str(root_resid)])] = ''

        while 1:
            connectivity = self.find_interresidue_bonds(root_atom, structure, bond_length)
            root_residue = ','.join([root_atom.getSegnames()[0], root_atom.getChids()[0], str(root_atom.getResnums()[0])])
            
            if not connectivity:
                break
            branch = 0
            for connect in connectivity:
                if connect[0] != root_residue:
                    connect_residue = connect[0]
                    connect_atom = structure.select('serial ' + connect[1])
                else:
                    connect_residue = connect[2]
                    connect_atom = structure.select('serial ' + connect[3])
                if root_residue in connect_tree:
                    connect_tree[connect_residue] = connect_tree[root_residue] + ' ' + connect[-1]
                else:
                    connect_tree[connect_residue] = connect[-1]
                if branch:
                    self.build_connectivity_tree(connect_atom, structure, bond_length=bond_length, connect_tree=connect_tree)
                else:
                    if connect[0] != root_residue:
                        new_root_atom = structure.select('serial ' + connect[1])
                    else:
                        new_root_atom = structure.select('serial ' + connect[3])
                    branch = 1
            root_atom = new_root_atom
        return connect_tree 

    
    def find_interresidue_bonds(self, root_atom, structure, bond_length = 1.55):
        """ Searches for all inter residue bonds and checks that a patch connecting two residues exists.
        Two bonds between the same pair of residue are not allowed
        Parameters:
            root_atom: AtomGroup for . bonds to this atom will not be considered
            structure: AtomGroup of the entire structure (e.g. PDB)
        Returns
            connectivity: list of bonds [atom1, atom2, patch]
        """
        idx = root_atom.getSerials()[0]
        root_name = root_atom.getNames()[0]
        root_resid = root_atom.getResnums()[0]
        root_chain = root_atom.getChids()[0]
        root_segname = root_atom.getSegnames()[0]
        sel = structure.select("(not protein) and (within " + str(bond_length) + " of same residue as serial " + str(idx) + ") and not (same residue as serial " + str(idx) + ")")
        names = sel.getNames()
        serials = sel.getSerials()
        resids = sel.getResnums()
        chains = sel.getChids()
        segnames = sel.getSegnames()
        resnames = sel.getResnames()
        connectivity = []
        for n,ss,r,c,s,rn in zip(names, serials, resids, chains, segnames, resnames):
            sel = structure.select("(same residue as serial " + str(idx) + ") and (within " + str(bond_length) + " of serial " + str(ss) + ")")
            root_a = sel.getNames()[0]
            root_serial = sel.getSerials()[0]
            #
            if root_a != root_name:
                patch = self.find_patch(root_a, n)
                if patch:
                    key1 = ','.join([root_segname, root_chain, str(root_resid)])
                    key2 = ','.join([s,c,str(r)])
                    connectivity.append([key1, str(idx), key2, str(ss), patch])
                else:
                    print 'Unknown inter residue bond'
    
        return connectivity



    def find_patch(self, atom1, atom2):
        key = '-'.join(sorted([atom1, atom2]))
        if key in self.builder.Topology.atomnames_to_patch:
            return self.builder.Topology.atomnames_to_patch[key]
        else:
             return ''

        
    def identify_glycan(self, connect_tree, glycan_keys):
        """Identify glycan polymer
            Parameters:
                connect_tree:
                glycan_keys:
            Returns:
                gk: glycan name (empty string if unknown)
        """
        target=[]
        for r in connect_tree.keys():
            target.append(r.get_resname() + ' ' + connect_tree[r])
        len_target=len(target)
        for gk in glycan_keys:
                if len(set(target) & set(glycan_keys[gk])) == len_target and len(set(glycan_keys[gk])) == len_target:
                    return gk
        print 'Unknown glycan'
        return ''
     
