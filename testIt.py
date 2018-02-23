import glycosylator as GL
from prody import *

#Reads a NAG and adds all missing atoms
myMol = GL.Molecule('NAG')
myMol.read_molecule_from_PDB('support/examples/NAG.pdb')
myBuilder = GL.MoleculeBuilder('./support/toppar_charmm/mannose.rtf', './support/toppar_charmm/mannose.prm')
myBuilder.Topology.read_topology('./support/topology/DUMMY.top')
r,a,b = myBuilder.add_missing_atoms(myMol.molecule['G'][1])
myBuilder.build_missing_atom_coord(r, a, myBuilder.Topology.topology['NAG']['IC'])
writePDB('NAG_complete.pdb', r)
#Connect a second NAG using a 14 connectivity
r2, da, b2 = myBuilder.build_from_patch(r, 2, 'NAG', r.getChids()[0], r.getSegnames()[0], '14bb')
writePDB('NAG2_complete.pdb', r2)
#Build a MAN ab initio
rd, bd = myBuilder.build_from_DUMMY(1, 'MAN', 'G', '1G', 'DUMMY_MAN')
writePDB('MAN_DUMMY.pdb', rd)

#Load a Man9 molecule
myMan9 = GL.Molecule('man9')
myMan9.read_molecule_from_PDB('./support/examples/man9.pdb')
myGlycosylator = GL.Glycosylator('./support/toppar_charmm/mannose.rtf', './support/toppar_charmm/mannose.prm')
myGlycosylator.builder.Topology.read_topology('./support/topology/DUMMY.top')
myGlycosylator.read_connectivity_topology('./support/topology/mannose.top')

#identify glycan
myGlycosylator.assign_patches(myMan9)
print myGlycosylator.identify_glycan(myMan9)

#trim it down to man6
connect_tree = myGlycosylator.build_connectivity_tree(myMan9.rootRes, myMan9.interresidue_connectivity)
man6, bonds6  = myGlycosylator.glycosylate('MAN6_1;3,2', template_glycan_tree = connect_tree, template_glycan = myMan9.molecule)
writePDB('man6.pdb', man6)


#extend man6 to man8
myMan6 = GL.Molecule('man6')
myMan6.read_molecule_from_PDB('man6.pdb')
myGlycosylator.assign_patches(myMan6)

connect_tree = myGlycosylator.build_connectivity_tree(myMan6.rootRes, myMan6.interresidue_connectivity)
man8,bonds8 = myGlycosylator.glycosylate('MAN8_2;4,2', template_glycan_tree = connect_tree, template_glycan = myMan6.molecule)
writePDB('man8_1.pdb', man8)


#Build a man8 ab initio
man8, bonds8 = myGlycosylator.glycosylate('MAN8_3;3,2')
writePDB('man8.pdb', man8)

#Rotate man9
myMan9.define_torsionals()
myMan9.rotate_bond(8, 60)
myMan9.rotate_bond(59, 60)
writePDB('man9_rot.pdb', myMan9.molecule)
