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
myGlycosylator.read_connectivity_topology('./glycans.dat')


################################################################################
# 4. Automatically detect all sequons and N-linked glycans
print '#'*50
print 'Loading glycoprotein'
myGlycosylator.load_glycoprotein('5fyl.pdb')
myGlycosylator.build_glycan_topology(patch = 'NGLB')

glycosylation = {}
glycosylation['88'] = 'M5'
glycosylation['133'] = 'FA3'
glycosylation['137'] = 'FA2'
glycosylation['156'] = 'M9'
glycosylation['160'] = 'M9'
glycosylation['197'] = 'M9'
glycosylation['234'] = 'M9'
glycosylation['262'] = 'M9'
glycosylation['276'] = 'M5'
glycosylation['295'] = 'M9'
glycosylation['301'] = 'M9'
glycosylation['332'] = 'M9'
glycosylation['339'] = 'M9'
glycosylation['355'] = 'M5'
glycosylation['363'] = 'M9'
glycosylation['386'] = 'M9'
glycosylation['392'] = 'M9'
glycosylation['396'] = 'M5'
glycosylation['411'] = 'M5'
glycosylation['448'] = 'M9'
glycosylation['462'] = 'FA3'
glycosylation['611'] = 'FA3'
glycosylation['618'] = 'FA2'
glycosylation['625'] = 'FA2'
glycosylation['637'] = 'M5'

linked_glycans = {}
linked_glycanMolecules = {}
print 'Building glycans'
for sequon_id in myGlycosylator.sequons:
    sequon = myGlycosylator.get_residue(sequon_id)
    seg,chain,resid,i =  sequon_id.split(',')
    
    glycan_name = glycosylation[resid]
    segname = chain + resid
    #extend exiting glycans    
    if sequon_id in myGlycosylator.glycanMolecules: 
        template = myGlycosylator.glycanMolecules[sequon_id]
        template_tree =  myGlycosylator.build_connectivity_tree(template.rootRes, template.interresidue_connectivity)
        glycan, bonds = myGlycosylator.glycosylate(glycan_name, template_glycan_tree = template_tree, template_glycan = template.atom_group, link_residue = sequon, link_patch = 'NGLB', chain =  'G', segname=segname)
    #build glycan ab initio
    else:
        glycan, bonds = myGlycosylator.glycosylate(glycan_name, link_residue = sequon, link_patch = 'NGLB', chain =  'G', segname=segname)

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


myGlycosylator.write_glycoprotein('5fyl_step_0.pdb')
dihe_parameters = myGlycosylator.builder.Parameters.parameters['DIHEDRALS']
vwd_parameters = myGlycosylator.builder.Parameters.parameters['NONBONDED']

mySampler = gl.Sampler(myGlycosylator.glycanMolecules.values(), myGlycosylator.protein, dihe_parameters, vwd_parameters, clash_dist = 1.6)
mySampler.remove_clashes_GA_iterative(n_iter = 10, n_individues = 6, n_generation = 10, pop_size=30, mutation_rate=0.01)

myGlycosylator.write_glycoprotein('5fyl_step_1.pdb')
mySampler.remove_clashes_GA_iterative(n_iter = 10, n_individues =  4, n_generation = 10, pop_size=30, mutation_rate=0.01)


myGlycosylator.write_glycoprotein('5fyl_step_2.pdb')
mySampler.remove_clashes_GA_iterative(n_iter = 10, n_individues =  2, n_generation = 20, pop_size=30, mutation_rate=0.01)

myGlycosylator.write_glycoprotein('5fyl_glycosylated_optimized.pdb')
