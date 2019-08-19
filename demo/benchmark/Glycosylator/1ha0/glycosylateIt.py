#!usr/bin/env python 
import glycosylator as gl
import prody as pd
import os
import glob

################################################################################
# 1. Create a glycosylator
myGlycosylator = gl.Glycosylator(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.rtf'), os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.prm'))
# 2. Additional topology files (in this case for building glycans ab initio [DUMMY pathc])
myGlycosylator.builder.Topology.read_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/DUMMY.top'))
# 3. Load glycan connectivity tree library
myGlycosylator.read_connectivity_topology('../glycans.dat')


################################################################################
# 4. Automatically detect all sequons and N-linked glycans
pdb_name = '1ha0'
print '#'*50
print 'Loading glycoprotein'
myGlycosylator.load_glycoprotein('../../' + pdb_name + '.pdb')
myGlycosylator.build_glycan_topology(patch = 'NGLB')


linked_glycans = {}
linked_glycanMolecules = {}
print 'Building glycans'
glycan_name = 'M9' 
for sequon_id in myGlycosylator.sequons:
    sequon = myGlycosylator.get_residue(sequon_id)
    seg,chain,resid,i =  sequon_id.split(',')
    
    segname = chain + resid
    #extend exiting glycans    
    if sequon_id in myGlycosylator.glycanMolecules: 
        template = myGlycosylator.glycanMolecules[sequon_id]
        template_tree =  myGlycosylator.build_connectivity_tree(template.rootRes, template.interresidue_connectivity)
        glycan, bonds = myGlycosylator.glycosylate(glycan_name, template_glycan_tree = template_tree, template_glycan = template.atom_group, link_residue = sequon, link_patch = 'NGLB', chain =  'G', segname=segname)
    #build glycan ab initio
    else:
        glycan, bonds = myGlycosylator.glycosylate(glycan_name, template_glycan_tree = None, template_glycan = None, link_residue = sequon, link_patch = 'NGLB', chain =  'G', segname=segname)

    new_glycan = gl.Molecule(sequon_id)
    new_glycan.set_AtomGroup(glycan, bonds = bonds)
    myGlycosylator.assign_patches(new_glycan)
    at = myGlycosylator.assign_atom_type(new_glycan) 
    new_glycan.set_atom_type(at)
    new_glycan.define_torsionals(hydrogens=False)
    linked_glycanMolecules[sequon_id] = new_glycan
    linked_glycans[sequon_id] = [new_glycan.rootRes, new_glycan.interresidue_connectivity]

################################################################################
# 5. Create a dictionary for storing new glycans
#Update glycan and save pdb
myGlycosylator.glycanMolecules.update(linked_glycanMolecules)
myGlycosylator.glycans.update(linked_glycans)

for g in myGlycosylator.glycanMolecules:
    myGlycosylator.names.update(myGlycosylator.glycanMolecules[g].get_names())


myGlycosylator.write_glycoprotein(pdb_name + '_step_0.pdb')
dihe_parameters = myGlycosylator.builder.Parameters.parameters['DIHEDRALS']
vwd_parameters = myGlycosylator.builder.Parameters.parameters['NONBONDED']

mySampler = gl.Sampler(myGlycosylator.glycanMolecules.values(), myGlycosylator.protein, dihe_parameters, vwd_parameters, clash_dist = 1.6)
mySampler.remove_clashes_GA_iterative(n_iter = 5, n_individues = 6, n_generation = 10, pop_size=30, mutation_rate=0.01)

myGlycosylator.write_glycoprotein(pdb_name + '_step_1.pdb')
mySampler.remove_clashes_GA_iterative(n_iter = 5, n_individues =  4, n_generation = 10, pop_size=30, mutation_rate=0.01)


myGlycosylator.write_glycoprotein(pdb_name + '_step_2.pdb')
mySampler.remove_clashes_GA_iterative(n_iter = 5, n_individues =  2, n_generation = 20, pop_size=30, mutation_rate=0.01)

myGlycosylator.write_glycoprotein(pdb_name + '_glycosylated_optimized.pdb')
