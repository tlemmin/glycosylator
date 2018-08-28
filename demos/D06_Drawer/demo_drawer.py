#!usr/bin/env python 
 
  
################################################################################ 
# A GENERAL EXPLANATION 

""" 
demo_sampler.py

The Sampler class is for sampling conformations of glycans
The following code illustrates how to optimize all glycans in a glycoprotein

"""

import glycosylator as gl
import prody as pd
import os
import matplotlib.pyplot as plt
import numpy as np
#Create a glycosylator
myGlycosylator = gl.Glycosylator(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.rtf'), os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.prm'))
myGlycosylator.builder.Topology.read_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/DUMMY.top'))
#Load topology information about glycans
myGlycosylator.read_connectivity_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/mannose.top'))
#Load HIV-1 Env trimer
myGlycosylator.load_glycoprotein(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/examples/env_4tvp.pdb'))
myGlycosylator.build_glycan_topology(patch = 'NGLB')

myDrawer = gl.Drawer()
for chid in  myGlycosylator.sequences.keys():
    l = len(myGlycosylator.sequences[chid])
    sequons = [k for k in myGlycosylator.sequons.keys() if chid in k[:len(chid)]]
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111)
    myDrawer.draw_glycoprotein(l, myGlycosylator.get_start_resnum(chid), sequons, ax = ax , trees=myGlycosylator.glycans, names = myGlycosylator.names)
    ax.axis('equal')
    ax.axis('off')
plt.show()
