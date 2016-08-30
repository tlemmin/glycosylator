#!/usr/bin/env python

import xml.etree.cElementTree as ET
from prody import *
import os, sys, argparse
import array
PDBx = '{http://deposit.pdb.org/pdbML/pdbx.xsd}'


def main():
    par = argparse.ArgumentParser()
    par.add_argument('--x', required = True, help = 'input XML file')
    par.add_argument('--o', required = True, help = 'output PDB')
    args = par.parse_args()

    tree = ET.parse(args.x)
    root = tree.getroot()
    residue = None
    for chem_atoms in root.iter(PDBx+'chem_comp_atomCategory'):
        for chem_atom in chem_atoms:
            resname = chem_atom.get('comp_id')
            name = chem_atom.get('atom_id')
            element = chem_atom.find(PDBx+'type_symbol')
            coords = []
            coords.append(float(chem_atom.find(PDBx+'pdbx_model_Cartn_x_ideal').text))
            coords.append(float(chem_atom.find(PDBx+'pdbx_model_Cartn_y_ideal').text))
            coords.append(float(chem_atom.find(PDBx+'pdbx_model_Cartn_z_ideal').text))
            
            atom = AtomGroup(name)
            atom.setCoords([coords])
            atom.setResnames([resname])
            atom.setResnums([1])
            atom.setNames([name])
            atom.setElements([element])
            if residue:
                residue += atom
            else:
                residue = atom
    residue.setTitle('PDBx ' + resname)
    writePDB(args.o, residue)


if __name__ == "__main__":
    main()

