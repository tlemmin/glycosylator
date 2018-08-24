import glycosylator as GL
from prody import *
import matplotlib.pyplot as plt
import numpy as np
#Reads a NAG and adds all missing atoms

#Load a Man9 molecule
myGlycosylator = GL.Glycosylator('./support/toppar_charmm/carbohydrates.rtf', './support/toppar_charmm/carbohydrates.prm')
myGlycosylator.builder.Topology.read_topology('./support/topology/DUMMY.top')
myGlycosylator.read_connectivity_topology('./support/topology/mannose.top')

myGlycosylator.load_glycoprotein('./test_sampler.pdb')
#myGlycosylator.load_glycoprotein('./env_4tvp.pdb')
myGlycosylator.build_glycan_topology(patch = 'NGLB')
#writePDB('88.pdb', myGlycosylator.glycanMolecules[',G,88,'].atom_group)
#Detect clashes
dihe_parameters = myGlycosylator.builder.Parameters.parameters['DIHEDRALS']
vwd_parameters = myGlycosylator.builder.Parameters.parameters['NONBONDED']
mySampler = GL.Sampler(myGlycosylator.glycanMolecules.values(), myGlycosylator.protein, dihe_parameters, vwd_parameters)
#torsionals = mySampler.get_all_torsional_angles()
#mySampler.compute_TotalEnergy(torsionals)
#mySampler.minimize_molecules()
mySampler.remove_clashes_GA(n_generation = 50, pop_size=30, mutation_rate=0.01)
#mySampler.remove_clashes()
myGlycosylator.write_glycoprotein('HIV_test.pdb')

