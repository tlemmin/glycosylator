#!usr/bin/env python 
 
  
################################################################################ 
# A GENERAL EXPLANATION 

""" 
demo_glycosylator.py

The Glycosylator class is used for modelling glycans.
The following code will illustrate how to manipulate single glycans

"""

import glycosylator as gl
import prody as pd
import os

#Create a glycosylator
myGlycosylator = gl.Glycosylator(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.rtf'), os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.prm'))
myGlycosylator.builder.Topology.read_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/DUMMY.top'))
#Load topology information about glycans
myGlycosylator.read_connectivity_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/mannose.top'))


#Identify glycan in mannose.top
myMan9 = gl.Molecule('mannose9')
myMan9.read_molecule_from_PDB(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/examples/man9.pdb'), update_bonds = True)
myGlycosylator.assign_patches(myMan9)
print "Identified glycan: ", myGlycosylator.identify_glycan(myMan9)

#Trim man9 down to man6
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

#Manually define a new glycan and build it ab initio
#NAG-14bb->NAG-14bb->BMA-[13ab->MAN]16ab->MAN-[16ab->MAN]13ab->MAN-12aa->MAN
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

myGlycosylator.connect_topology.update(myglycan)
myMan, bonds = myGlycosylator.glycosylate('NEW_Glycan')
pd.writePDB('myMan_ab_initio.pdb', myMan)
