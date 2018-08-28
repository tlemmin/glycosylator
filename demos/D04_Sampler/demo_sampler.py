#!usr/bin/env python 
 
  
################################################################################ 
# A GENERAL EXPLANATION 

""" 
demo_sampler.py

The Sampler class is for sampling conformations of glycans

"""

import glycosylator as gl
import prody as pd
import os

#Create a glycosylator
myGlycosylator = gl.Glycosylator(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.rtf'), os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.prm'))
myGlycosylator.builder.Topology.read_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/DUMMY.top'))
#Load topology information about glycans
myGlycosylator.read_connectivity_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/mannose.top'))

#Build a mannose 9 ab initio
man9, bonds = myGlycosylator.glycosylate('MAN9_3;4,2')
myMan9 = gl.Molecule('mannose9')
myMan9.set_AtomGroup(man9, bonds = bonds, update_bonds = False)
myMan9.update_connectivity(update_bonds = False)
myGlycosylator.assign_patches(myMan9)
#Atom type and initialize from AtomGroup
atom_type = myGlycosylator.assign_atom_type(myMan9)
myMan9.set_atom_type(atom_type)
#Do not consider torsional angles invovling terminal hydrogen
myMan9.define_torsionals(hydrogens =  False)

#Remove clashes with genetic algorithm
dihe_parameters = myGlycosylator.builder.Parameters.parameters['DIHEDRALS']
vwd_parameters = myGlycosylator.builder.Parameters.parameters['NONBONDED']

mySampler = gl.Sampler([myMan9], None, dihe_parameters, vwd_parameters)
mySampler.remove_clashes_GA(n_generation = 15, pop_size=30, mutation_rate=0.01)
myMan9.writePDB('Man9_optimized.pdb')


