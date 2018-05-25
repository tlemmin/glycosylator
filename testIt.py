import glycosylator as GL
from prody import *
import matplotlib.pyplot as plt
import numpy as np
#Reads a NAG and adds all missing atoms
myMol = GL.Molecule('NAG')
myMol.read_molecule_from_PDB('support/examples/NAG.pdb')
myBuilder = GL.MoleculeBuilder('./support/toppar_charmm/carbohydrates.rtf', './support/toppar_charmm/carbohydrates.prm')
myBuilder.Topology.read_topology('./support/topology/DUMMY.top')
r,a,b = myBuilder.add_missing_atoms(myMol.atom_group['G'][1])
myBuilder.build_missing_atom_coord(r, a, myBuilder.Topology.topology['NAG']['IC'])
writePDB('NAG_complete.pdb', r)
#Connect a second NAG using a 14 connectivity
r2, d2, b2 = myBuilder.build_from_patch(r, 2, 'NAG', r.getChids()[0], r.getSegnames()[0], '14bb')
writePDB('NAG2_complete.pdb', r2)
#Build a MAN ab initio
rd, da, bd = myBuilder.build_from_DUMMY(1, 'MAN', 'G', '1G', 'DUMMY_MAN')
writePDB('MAN_DUMMY.pdb', rd)

#Load a Man9 molecule
myMan9 = GL.Molecule('man9')
myMan9.read_molecule_from_PDB('./support/examples/man9.pdb')
myGlycosylator = GL.Glycosylator('./support/toppar_charmm/carbohydrates.rtf', './support/toppar_charmm/carbohydrates.prm')
myGlycosylator.builder.Topology.read_topology('./support/topology/DUMMY.top')
myGlycosylator.read_connectivity_topology('./support/topology/mannose.top')

#Identify glycan
myGlycosylator.assign_patches(myMan9)
print "Identified glycan: ", myGlycosylator.identify_glycan(myMan9)

#Trim man9 down to man6
connect_tree = myGlycosylator.build_connectivity_tree(myMan9.rootRes, myMan9.interresidue_connectivity)
man6, bonds6  = myGlycosylator.glycosylate('MAN6_1;3,2', template_glycan_tree = connect_tree, template_glycan = myMan9.atom_group)
writePDB('man6.pdb', man6)


#Atom type and initialize from AtomGroup
myMan6 = GL.Molecule('man6')
myMan6.set_AtomGroup(man6, bonds = bonds6)
myGlycosylator.assign_patches(myMan6)
atom_type = myGlycosylator.assign_atom_type(myMan6)
myMan6.set_atom_type(atom_type)
#myMan6.set_bonds(bonds6)
#myMan6.read_molecule_from_PDB('man6.pdb')


connect_tree = myGlycosylator.build_connectivity_tree(myMan6.rootRes, myMan6.interresidue_connectivity)
man8,bonds8 = myGlycosylator.glycosylate('MAN8_2;4,2', template_glycan_tree = connect_tree, template_glycan = myMan6.atom_group)
writePDB('man8_1.pdb', man8)


#Build a man8 ab initio
man8, bonds8 = myGlycosylator.glycosylate('MAN8_3;3,2')
writePDB('man8.pdb', man8)

#Rotate man9
myMan9.define_torsionals()
myMan9.rotate_bond(8, 60)
myMan9.rotate_bond(59, 60)
writePDB('man9_rot.pdb', myMan9.atom_group)

#N-glycosylate
HIV = parsePDB('support/examples/4tvp.pdb')
ASN = HIV.select('resid 625 and chain B')
ASN.setSegnames(['B']*len(ASN))

Nman8, Nbonds8 = myGlycosylator.glycosylate('MAN8_3;3,2', link_residue=ASN, link_patch = 'NGLB')

myNman8 = GL.Molecule('N-man8')
myNman8.set_AtomGroup(Nman8, bonds=Nbonds8)
myGlycosylator.assign_patches(myNman8)
atom_type = myGlycosylator.assign_atom_type(myNman8)
myNman8.set_atom_type(atom_type)
myNman8.define_torsionals(hydrogens=False)


#Detect clashes
dihe_parameters = myGlycosylator.builder.Parameters.parameters['DIHEDRALS']
mySampler = GL.Sampler([myNman8], HIV, dihe_parameters)
#nbr_count,clashes = mySampler.count_clashes(myMan9)
#print nbr_count, 'clashes were detected: ', clashes

#Draw
protein = parsePDB('./support/examples/4tvp.pdb')
myGlycosylator.load_glycoprotein(protein)

myDrawer = GL.Drawer()
i = 1
for chid in  myGlycosylator.sequences.keys():
    l = len(myGlycosylator.sequences[chid])
    sequons = [k for k in myGlycosylator.sequons.keys() if chid in k[:len(chid)]]
    ax = myDrawer.draw_glycoprotein(l, myGlycosylator.get_start_resnum(chid), sequons, axis = i, trees=myGlycosylator.glycans, names = myGlycosylator.names)
    #glycans = {}
    #for s in sequons:
    #    if s in myGlycosylator.glycans:
    #        glycans[s] = myGlycosylator.glycans[s]
    #ax = myDrawer.draw_all_trees(glycans, myGlycosylator.get_start_resnum(chid), myGlycosylator.names, ax=ax, axis = i)
    i = np.mod(i+1, 2)
plt.show(block=False)

