#!/usr/bin/env python

from prody import *
import networkx as nx
import numpy as np
import os, sys, argparse

def readfile(fname):
    with open(fname) as f:
        content = [x.strip('\n') for x in f.readlines()]
    return content

def find_paths(G, node, length, excludeSet = None):
    """Finds all paths of a given length
    Parameters:
        G: graph (networkx)
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
    paths = [[node]+path for neighbor in G.neighbors(node) if neighbor not in excludeSet for path in find_paths(G, neighbor, length-1, excludeSet)]
    excludeSet.remove(node)
    return paths

def guess_ICs(structure, root_atom): 
    kd = KDTree(structure.getCoords())
    kd.search(1.7)
    atoms = kd.getIndices()
        
    atom_names = structure.getNames()
    root_idx = np.where(atom_names == root_atom)[0][0]
    G = nx.Graph()
    for id1,id2 in atoms:
        G.add_node(id1, name = atom_names[id1])
        G.add_node(id2, name = atom_names[id2])
        G.add_edge(id1, id2)

    t = nx.dfs_tree(G,root_idx)
    ics = []
    for node in t.nodes():
        ics_idx = find_paths(t, node, 3)
        for ic in ics_idx:
            ics.append(' '.join([atom_names[i] for i in ic]))
    return ics

def main():
    par = argparse.ArgumentParser()
    par.add_argument('--p', required = True, help = 'input PDB file')
    par.add_argument('--i', required = False, help = 'list of internal coordinates')
    par.add_argument('--r', required = False, help = 'root atom if no internal coordinates are provided')
    par.add_argument('--o', required = True, help = 'internal coordinate file')

    args = par.parse_args()
    structure = parsePDB(args.p)
    if args.i is None:
        ics = guess_ICs(structure, root_atom)
    else:
        ics = readfile(args.i)
    to_write = []
    for ic in ics:
        a1,a2,a3,a4 = ic.split()
        internal_coordinate = 'IC %s %s %s %s' % (a1,a2,a3,a4)
        a1 = structure.select('name ' + a1)
        a2 = structure.select('name ' + a2)
        a3 = structure.select('name ' + a3)
        a4 = structure.select('name ' + a4)
        bond1 = calcDistance(a1, a2)[0]
        angle1 = calcAngle(a1, a2, a3)[0]
        dihe = calcDihedral(a1, a2, a3, a4)[0]
        angle2 = calcAngle(a2, a3, a4)[0]
        bond2 = calcDistance(a3, a4)[0]
        internal_coordinate += '  %.2f %.2f %.2f %.2f %.2f' % (bond1, angle1, dihe, angle2, bond2)
        to_write.append(internal_coordinate)

    f = open(args.o, 'w')
    f.write('\n'.join(to_write))
    f.close


if __name__ == "__main__":
    main()
