#!usr/bin/env python 
 
  
""" 
demo_polymer.py

The following code will illustrate how to use glycosylator to build a polymer. In this case,
we will build a small peptide.
1. Load glycosylator with topology files for amino acids
"""

import glycosylator as gl
import prody as pd
import os


################################################################################
# 1. Create a glycosylator
myGlycosylator = gl.Glycosylator('peptide.str', 'peptide.prm')

################################################################################

################################################################################
# 2. Define sequence of peptide: 
# GLYCSYLATR
peptide_topo = {}
peptide_topo['#UNIT'] = 10
peptide_topo['UNIT'] = [['GLY', ' ', []],
                       ['LEU', 'C1', ['LINK']],
                       ['TYR', 'C1', ['LINK', 'LINK']],
                       ['CYS', 'C1', ['LINK', 'LINK', 'LINK']],
                       ['SER', 'C1', ['LINK', 'LINK', 'LINK', 'LINK']],
                       ['TYR', 'C1', ['LINK', 'LINK', 'LINK', 'LINK', 'LINK']],
                       ['LEU', 'C1', ['LINK', 'LINK', 'LINK', 'LINK', 'LINK', 'LINK']],
                       ['ALA', 'C1', ['LINK', 'LINK', 'LINK', 'LINK', 'LINK', 'LINK', 'LINK']],
                       ['THR', 'C1', ['LINK', 'LINK', 'LINK', 'LINK', 'LINK', 'LINK', 'LINK', 'LINK']],
                       ['ARG', 'C1', ['LINK', 'LINK', 'LINK', 'LINK', 'LINK', 'LINK', 'LINK', 'LINK', 'LINK']]]
                       
mypeptide = {}
mypeptide['NEW_peptide'] = peptide_topo

################################################################################
# 7. Glycan built ab initio
myGlycosylator.connect_topology.update(mypeptide)
peptide, bonds = myGlycosylator.glycosylate('NEW_peptide')
pd.writePDB('peptide.pdb', peptide)
