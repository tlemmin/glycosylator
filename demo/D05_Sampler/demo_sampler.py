#!usr/bin/env python 
 
  
""" 
demo_sampler.py

The Sampler class is for sampling conformations of glycans
The following code illustrates how to optimize all glycans in a glycoprotein
1. Initialization of a Glycosylator instance
2. Load and identify all glycans in a glycoprotein
3. Initialize a Sampler
4. List of all glycans, the protein will remain static, but will be considered for clashes


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
# 2. Load HIV-1 Env trimer
myGlycosylator.load_glycoprotein(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/examples/env_4tvp.pdb'))
myGlycosylator.build_glycan_topology(patch = 'NGLB')

################################################################################
# 3. Remove clashes with genetic algorithm
dihe_parameters = myGlycosylator.builder.Parameters.parameters['DIHEDRALS']
vwd_parameters = myGlycosylator.builder.Parameters.parameters['NONBONDED']
# 4. Consider all glycans in HIV-1 Env trimer
mySampler = gl.Sampler(myGlycosylator.glycanMolecules.values(), myGlycosylator.protein, dihe_parameters, vwd_parameters)
#mySampler.remove_clashes_GA(n_generation = 20, pop_size=30, mutation_rate=0.01)
mySampler.remove_clashes_GA_iterative(n_iter = 5, n_individues =  3, n_generation = 30, pop_size=30, mutation_rate=0.01)
myGlycosylator.write_glycoprotein('HIV_optimized.pdb')


