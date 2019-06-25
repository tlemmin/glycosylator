#!usr/bin/env python 
  
""" 
demo_glycolipid.py

The following code will illustrate how to build a glycolipid.
1. The initialization of a Glysosylator instance requires  topology and parameter files
2. Load lipid, in this case a CERAMIDE 16
3. Modify the CERAMIDE with two GALACTOSE (Galactolipid)
"""

import glycosylator as gl
import prody as pd
import os
################################################################################
# 1. Create a glycosylator
myGlycosylator = gl.Glycosylator(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.rtf'), os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.prm'))
myGlycosylator.builder.Topology.read_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/DUMMY.top'))

################################################################################
# 2. Load ceramide
cer160 = pd.parsePDB(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/examples/cer160.pdb'))
# Define connectivity of di-galactose
glycan_topo = {}
glycan_topo['#UNIT'] = 2
glycan_topo['UNIT'] = [['GAL', ' ', []],
                       ['GAL', 'C1', ['14bb']]]
myglycan = {}
myglycan['GAL2'] = glycan_topo

################################################################################
# 3. Glycan built ab initio
myGlycosylator.connect_topology.update(myglycan)

gal2,bonds2 = myGlycosylator.glycosylate('GAL2', link_residue=cer160, link_patch = 'CERB')
galactolipid = cer160 + gal2
pd.writePDB('galactolipid.pdb', galactolipid)