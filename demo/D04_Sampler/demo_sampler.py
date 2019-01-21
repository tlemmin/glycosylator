#!usr/bin/env python 
 
  
""" 
demo_sampler.py

The Sampler class is for sampling conformations of Molecules.
1. Initialize an instance of Glycosylator
2. Build a complete glycan (Mannose 9 here)
3. Assign atom types (CHARMM force field) to glycan
4. Torsional involving terminal hydrogens can be neglected for computation efficiency 
5. Initialization of a Sampler instance requires a list of molecules, eventually a AtomGoup which will not changed (e.g. protein), parameters for dihedral and VdW

"""

import glycosylator as gl
import prody as pd
import os
################################################################################
# 1. Create a glycosylator
myGlycosylator = gl.Glycosylator(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.rtf'), os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.prm'))
myGlycosylator.builder.Topology.read_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/DUMMY.top'))
#Load topology information about glycans
myGlycosylator.read_connectivity_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/mannose.top'))

################################################################################
# 2. Build a mannose 9 ab initio
man9, bonds = myGlycosylator.glycosylate('MAN9_3;4,2')
myMan9 = gl.Molecule('mannose9')
myMan9.set_AtomGroup(man9, bonds = bonds, update_bonds = False)
myMan9.update_connectivity(update_bonds = False)
myGlycosylator.assign_patches(myMan9)

################################################################################
# 3. Atom type and initialize from AtomGroup
atom_type = myGlycosylator.assign_atom_type(myMan9)
myMan9.set_atom_type(atom_type)
# 4. Do not consider torsional angles invovling terminal hydrogen
myMan9.define_torsionals(hydrogens =  False)

################################################################################
# 5. Remove clashes with genetic algorithm
dihe_parameters = myGlycosylator.builder.Parameters.parameters['DIHEDRALS']
vwd_parameters = myGlycosylator.builder.Parameters.parameters['NONBONDED']

mySampler = gl.Sampler([myMan9], None, dihe_parameters, vwd_parameters)
mySampler.remove_clashes_GA_iterative(n_iter = 1, n_generation = 30, pop_size=30, mutation_rate=0.01)
myMan9.writePDB('Man9_optimized.pdb')


