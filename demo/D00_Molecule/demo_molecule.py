#!usr/bin/env python 
 

""" 
demo_molecule.py

The Molecule class stores the coordinates and connectivity of a molecule. It provides functions for detecting cycles and rotate bonds.
1. A Molecule can be initialize with a pdb file. The file should containt only a single molecule (i.e. single chain and segname)
Torsionals angles are defined by the serial number of the atoms. The torsional angle can be rotate/measured either by:
    2. by the index.
    3. explicitly defining it (list of four the serial number of atoms)
The rotation angle is given in degrees and can be used to set 
    4. a new diherdal angle 
    5. to add an increment to the current angle.

6. A Molecule can also be initilized with an AtomGroup (from Prody)

"""

import glycosylator as gl
import prody as pd
import os

############################################################################
#Create a Molecule 
myMan9 = gl.Molecule('mannose9')
# 1. Initialization with PDB file
#Load PDB file and guess all bonds, angles and dihedral angles
myMan9.read_molecule_from_PDB(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/examples/man9.pdb'), update_bonds = True)
#Define all torsionals angles
#Torsionals are defined by serial number of atoms
myMan9.define_torsionals()
print '#'*50
print 'Number of torsional angles: {}'.format(len(myMan9.torsionals))
print '#'*50

############################################################################
#Changing torsional angle

# 2. index of torsional angle
idx =  10 
# 3. list of 4 atoms
torsional = myMan9.torsionals[idx]
print 'Current dihedral for the torsional angle #{} ({}): {}'.format(idx, torsional, myMan9.measure_dihedral_angle(idx))
print 'Rotating the angle by 55 degrees...'
# 4. Set absolute angle (angle defined with index)
myMan9.rotate_bond(idx, 45., absolute = True)
# 5. Rotate by increment (angle explicitly defined with list of 4 atoms)
myMan9.rotate_bond(torsional, 10., absolute = False)
print 'New dihedral angle: {}'.format(myMan9.measure_dihedral_angle(idx))
print '#'*50

############################################################################
# 6. Creating molecule directly from AtomGroup

myGlycan =  gl.Molecule('glycan')
HIV_env = pd.parsePDB(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/examples/env_4tvp.pdb'))
myGlycan.set_AtomGroup(HIV_env.select('resid 1088 to 1094'), update_bonds = True)
myGlycan.writePDB('myglycan.pdb')
