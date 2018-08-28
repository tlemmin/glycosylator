#!usr/bin/env python 
 
  
################################################################################ 
# A GENERAL EXPLANATION 

""" 
demo_molecule.py

The Molecule class stores the coordiantes and connectivity of a molecule. 

"""

import glycosylator as gl
import prody as pd
import os

#Create a Molecule 
myMan9 = gl.Molecule('mannose9')
#Load PDB file and guess all bonds, angles and dihedral angles
myMan9.read_molecule_from_PDB(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/examples/man9.pdb'), update_bonds = True)
#Define all torsionals angles
#Torsionals are defined by serial number of atoms
myMan9.define_torsionals()
print '#'*50
print 'Number of torsional angles: {}'.format(len(myMan9.torsionals))
print '#'*50

#Changing torisonal angle
idx =  10 
torsional = myMan9.torsionals[idx]
print 'Current dihedral for the torsional angle #{} ({}): {}'.format(idx, torsional, myMan9.measure_dihedral_angle(idx))
print 'Rotating the angle by 55 degrees...'
#Set absolute angle
myMan9.rotate_bond(idx, 45., absolute = True)
#The torsional angle can also be defined as the corresponding list of atom serial numbers
myMan9.rotate_bond(torsional, 10., absolute = False)
print 'New dihedral angle: {}'.format(myMan9.measure_dihedral_angle(idx))
print '#'*50

############################################################################
#Creating molecule directly from AtomGroup
myGlycan =  gl.Molecule('glycan')
HIV_env = pd.parsePDB(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/examples/env_4tvp.pdb'))
myGlycan.set_AtomGroup(HIV_env.select('resid 1088 to 1094'), update_bonds = True)
myGlycan.writePDB('myglycan.pdb')
