import glycosylator as GL

#Load a Man9 molecule
myMan9 = GL.Molecule('man9')
myMan9.read_molecule_from_PDB('./support/examples/man9.pdb')
myGlycosylator = GL.Glycosylator('./support/toppar_charmm/carbohydrates.rtf', './support/toppar_charmm/carbohydrates.prm')
myGlycosylator.builder.Topology.read_topology('./support/topology/DUMMY.top')
myGlycosylator.read_connectivity_topology('./support/topology/mannose.top')

#Identify glycan
myGlycosylator.assign_patches(myMan9)
print "Identified glycan: ", myGlycosylator.identify_glycan(myMan9)


myDrawer = GL.Drawer()



print myDrawer.tree_to_text(tree, root, self.names)
