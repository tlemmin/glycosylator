#!usr/bin/env python 
 
  
""" 
demo_glycosylator.py

The Glycosylator class is used for modelling glycans.
The following code will illustrate how to manipulate single glycans
1. The initialization of a glycosylator requires a topology and parameter file
2. Additional topology files can be added
3. A connectivity tree file is used to identify and build glycans. A glycan is composed of "UNIT" (! are used for comments)
    A UNIT is defined by its residue name (3 letter code from PDB), the connecting atom, and the path from the root UNIT (first unit in glycan [connected to sequon])
    The path is defined by the patches that have applied to connect previous UNIT.

    Example for mannose 3: 
                                                -16ab-> MAN[Z4] -13ab-> MAN[Z7]
        NAG[Z1] -14bb-> NAG[Z2] -14bb-> BMA[Z3]-|
                                                -13ab-> MAN[Z7]    
    Corresponding connectivity topology:
    RESI MAN3                               !User defined, which has to be unique within a top file
    UNIT NAG                                !Z1 Root unit in glycan
    UNIT NAG C1 14bb                        !Z2 Unit connected to root with a 14bb patch
    UNIT BMA C1 14bb 14bb                   !Z3 Unit connected with 14bb to Z2 [path from root: 14bb]
    UNIT MAN C1 14bb 14bb 16ab              !Z4 Unit connected with 16ab to Z3 [path from root: 14bb 14bb]
    UNIT MAN C1 14bb 14bb 16ab 13ab         !Z7 Unit connected with 13ab to Z4 [path from root: 14bb 14bb 16ab]
    UNIT MAN C1 14bb 14bb 13ab              !Z9 Unit connected with 13ab to Z3 [path from root: 14bb 14bb]

4. The connectivity tree is atomatically guessed based on bonds and patches. This can be used to identify glycan in connectivity tree library
5. A template glycan can be modified to match the defined connectivity tree (e.g Mannose 6: MAN6_1;3,2)
6. A topology tree can also be directly defined as a dictionary
7. A glycan can be built ab initio
"""

import glycosylator as gl
import prody as pd
import os


################################################################################
# 1. Create a glycosylator
myGlycosylator = gl.Glycosylator(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.rtf'), os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.prm'))
# 2. Additional topology files (in this case for building glycans ab initio [DUMMY pathc])
myGlycosylator.builder.Topology.read_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/DUMMY.top'))
# 3. Load glycan connectivity tree library
myGlycosylator.read_connectivity_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/mannose.top'))


################################################################################
# 4. Identify glycan in mannose.top
myMan9 = gl.Molecule('mannose9')
myMan9.read_molecule_from_PDB(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/examples/man9.pdb'), update_bonds = True)
myGlycosylator.assign_patches(myMan9)
print "Identified glycan: ", myGlycosylator.identify_glycan(myMan9)

################################################################################
# 5. Trim man9 down to man6
connect_tree = myGlycosylator.build_connectivity_tree(myMan9.rootRes, myMan9.interresidue_connectivity)
man6, bonds6  = myGlycosylator.glycosylate('MAN6_1;3,2', template_glycan_tree = connect_tree, template_glycan = myMan9.atom_group)
myMan6 = gl.Molecule('mannose6')
myMan6.set_AtomGroup(man6, bonds = bonds6, update_bonds = False)
#Update angles and dihedrals based on bonds
myMan6.update_connectivity(update_bonds = False)
myGlycosylator.assign_patches(myMan6)
print "New glycan: ", myGlycosylator.identify_glycan(myMan6)
myMan6.writePDB('man6.pdb')


#Atom type and initialize from AtomGroup
atom_type = myGlycosylator.assign_atom_type(myMan6)
myMan6.set_atom_type(atom_type)

#Extend man8 from man6
connect_tree = myGlycosylator.build_connectivity_tree(myMan6.rootRes, myMan6.interresidue_connectivity)
man8,bonds8 = myGlycosylator.glycosylate('MAN8_2;4,2', template_glycan_tree = connect_tree, template_glycan = myMan6.atom_group)
pd.writePDB('man8_1.pdb', man8)

################################################################################
# 6. Manually define a new glycan and build it ab initio
# NAG-14bb->NAG-14bb->BMA-[13ab->MAN]16ab->MAN-[16ab->MAN]13ab->MAN-12aa->MAN
glycan_topo = {}
glycan_topo['#UNIT'] = 8
glycan_topo['UNIT'] = [['NAG', ' ', []],
                       ['NAG', 'C1', ['14bb']],
                       ['BMA', 'C1', ['14bb', '14bb']],
                       ['MAN', 'C1', ['14bb', '14bb', '16ab']],
                       ['MAN', 'C1', ['14bb', '14bb', '16ab', '16ab']],
                       ['MAN', 'C1', ['14bb', '14bb', '16ab', '13ab']],
                       ['MAN', 'C1', ['14bb', '14bb', '16ab', '13ab', '12aa']],
                       ['MAN', 'C1', ['14bb', '14bb', '13ab']]]
myglycan = {}
myglycan['NEW_Glycan'] = glycan_topo

################################################################################
# 7. Glycan built ab initio
myGlycosylator.connect_topology.update(myglycan)
myMan, bonds = myGlycosylator.glycosylate('NEW_Glycan')
pd.writePDB('myMan_ab_initio.pdb', myMan)
