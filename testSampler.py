import glycosylator as GL
from prody import *
import matplotlib.pyplot as plt
import numpy as np
#Reads a NAG and adds all missing atoms

#Load a Man9 molecule
myGlycosylator = GL.Glycosylator('./support/toppar_charmm/carbohydrates.rtf', './support/toppar_charmm/carbohydrates.prm')
myGlycosylator.builder.Topology.read_topology('./support/topology/DUMMY.top')
myGlycosylator.read_connectivity_topology('./support/topology/mannose.top')

myGlycosylator

    myGlycosylator.assign_patches(myMan9)
    connect_tree = myGlycosylator.build_connectivity_tree(myMan9.rootRes, myMan9.interresidue_connectivity)
    man9, bonds9  = myGlycosylator.glycosylate('MAN9_3;4,2', template_glycan_tree = connect_tree, template_glycan = myMan9.atom_group)

    myMan9.set_AtomGroup(man9, rootAtom=myMan9.rootAtom, bonds=bonds9)
    myGlycosylator.assign_patches(myMan9)

    myMan9.define_torsionals()

    atom_type = myGlycosylator.assign_atom_type(myMan9)
    myMan9.set_atom_type(atom_type)
    myMan9.define_torsionals(hydrogens=False)


Nman8, Nbonds8 = myGlycosylator.glycosylate('MAN8_3;3,2', link_residue=ASN, link_patch = 'NGLB')


#Detect clashes
dihe_parameters = myGlycosylator.builder.Parameters.parameters['DIHEDRALS']
vwd_parameters = myGlycosylator.builder.Parameters.parameters['NONBONDED']
mySampler = GL.Sampler([myNman8], HIV, dihe_parameters, vwd_parameters)
print mySampler.compute_self_non_bonded_energy(0)
mySampler.remove_clashes(n = 50)
print mySampler.compute_self_non_bonded_energy(0)
writePDB('man8_min.pdb', myNman8.atom_group)
#nbr_count,clashes = mySampler.count_clashes(myMan9)
#print nbr_count, 'clashes were detected: ', clashes

