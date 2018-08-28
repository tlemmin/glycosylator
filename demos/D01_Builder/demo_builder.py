#!usr/bin/env python 
 
  
################################################################################ 
# A GENERAL EXPLANATION 

""" 
demo_builder.py

The MoleculeBuilder class is used to build ab initio or add missing atoms to a molecule.
The topology and parameter files are based on CHARMM force field

"""

import glycosylator as gl
import prody as pd
import os


#Creates a builder instance
#This requires a topology and parameter in CHARMM format
myBuilder = gl.MoleculeBuilder(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.rtf'), os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.prm'))
myNAG = gl.Molecule('NAG')
myNAG.read_molecule_from_PDB(os.path.join(gl.GLYCOSYLATOR_PATH,'support/examples/NAG.pdb'))
NAG_complete,missing_atoms,bonds = myBuilder.add_missing_atoms(myNAG.atom_group)
#Add missing atoms for internal coordinates
myBuilder.build_missing_atom_coord(NAG_complete, missing_atoms, myBuilder.Topology.topology['NAG']['IC'])
myNAG.set_AtomGroup(NAG_complete, bonds = bonds, update_bonds = False) 
#Update angles and dihedrals based on bonds
myNAG.update_connectivity(update_bonds = False)
pd.writePDB('NAG_complete.pdb', NAG_complete)
#Connect a second NAG using a 14bb connectivity
NAG2, del_atoms, bonds2 = myBuilder.build_from_patch(myNAG.atom_group, 2, 'NAG', myNAG.get_chain(), myNAG.get_segname(), '14bb')
NAG_14bb_NAG = myBuilder.delete_atoms(NAG_complete, del_atoms)
NAG_14bb_NAG += myBuilder.delete_atoms(NAG2, del_atoms)
pd.writePDB('NAG_14bb_NAG.pdb', NAG_14bb_NAG)
#Build a MAN ab initio
#Load patch for building a molecule ab initio (DUMMY)
myBuilder.Topology.read_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/DUMMY.top'))
Man, del_atoms, bond_man = myBuilder.build_from_DUMMY(1, 'MAN', 'G', '1G', 'DUMMY_MAN')
pd.writePDB('MAN_DUMMY.pdb', Man)

