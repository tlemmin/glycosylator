#!usr/bin/env python 
 
  
################################################################################ 
# A GENERAL EXPLANATION 

""" 
demo_sampler.py

The Sampler class is for sampling conformations of glycans
The following code illustrates how to optimize all glycans in a glycoprotein

"""

import glycosylator as gl
import prody as pd
import os

#Create a glycosylator
myGlycosylator = gl.Glycosylator(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.rtf'), os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.prm'))
myGlycosylator.builder.Topology.read_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/DUMMY.top'))
#Load topology information about glycans
myGlycosylator.read_connectivity_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/mannose.top'))
#Load HIV-1 Env trimer
myGlycosylator.load_glycoprotein(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/examples/env_4tvp.pdb'))
myGlycosylator.build_glycan_topology(patch = 'NGLB')

#Remove clashes with genetic algorithm
dihe_parameters = myGlycosylator.builder.Parameters.parameters['DIHEDRALS']
vwd_parameters = myGlycosylator.builder.Parameters.parameters['NONBONDED']

mySampler = gl.Sampler(myGlycosylator.glycanMolecules.values(), myGlycosylator.protein, dihe_parameters, vwd_parameters)
mySampler.remove_clashes_GA(n_generation = 25, pop_size=30, mutation_rate=0.01)
myGlycosylator.write_glycoprotein('HIV_optimized.pdb')


