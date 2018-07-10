import glycosylator as GL
from prody import *
import matplotlib.pyplot as plt
import numpy as np
#Reads a NAG and adds all missing atoms

#Load a Man9 molecule
myGlycosylator = GL.Glycosylator('./support/toppar_charmm/carbohydrates.rtf', './support/toppar_charmm/carbohydrates.prm')
myGlycosylator.builder.Topology.read_topology('./support/topology/DUMMY.top')
myGlycosylator.read_connectivity_topology('./support/topology/mannose.top')

myMan9 = GL.Molecule('man9')
myMan9.read_molecule_from_PDB('./support/examples/man9.pdb',  update_bonds = False)
#build topology
man9, bonds9  = myGlycosylator.glycosylate('MAN3_1;1,1', glycan_molecule = myMan9)
myMan9.set_AtomGroup(man9, rootAtom=myMan9.rootAtom, bonds=bonds9)
myGlycosylator.assign_patches(myMan9)
atom_type = myGlycosylator.assign_atom_type(myMan9)

myMan9.set_atom_type(atom_type)
#Identify glycan
print "Identified glycan: ", myGlycosylator.identify_glycan(myMan9)
myMan9.define_torsionals(hydrogens =  False)
#print myMan9.get_all_torsional_angles()
#print myMan9.torsionals
#writePDB('Man3.pdb', myMan9.atom_group)
#exit()

#Detect clashes
dihe_parameters = myGlycosylator.builder.Parameters.parameters['DIHEDRALS']
vwd_parameters = myGlycosylator.builder.Parameters.parameters['NONBONDED']

mySampler = GL.Sampler([myMan9], None, dihe_parameters, vwd_parameters)
#torsionals = mySampler.get_all_torsional_angles()
mySampler.remove_clashes_GA()
writePDB('Man3.pdb', myMan9.atom_group)

#myGlycosylator.write_glycoprotein('HIV_test.pdb')

