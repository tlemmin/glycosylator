#!usr/bin/env python 
 
  
################################################################################ 
# A GENERAL EXPLANATION 

""" 
demo_glycosylator.py

The Glycosylator class is used for modelling glycans
The following code will illustrate how to manipulate entire glycoprotein

"""

import glycosylator as gl
import prody as pd
import os

#Create a glycosylator
myGlycosylator = gl.Glycosylator(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.rtf'), os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.prm'))
myGlycosylator.builder.Topology.read_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/DUMMY.top'))
#Load topology information about glycans
myGlycosylator.read_connectivity_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/mannose.top'))


#Manually N-glycosylate one sequon
HIV_Env = pd.parsePDB(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/examples/env_4tvp.pdb'))
sequon = HIV_Env.select('resid 625 and chain B')

man8,bonds8 = myGlycosylator.glycosylate('MAN8_3;3,2', link_residue=sequon, link_patch = 'NGLB')
myMan8 = gl.Molecule('N-man8')
myMan8.set_AtomGroup(man8, bonds=bonds8)
myGlycosylator.assign_patches(myMan8)
atom_type = myGlycosylator.assign_atom_type(myMan8)
myMan8.set_atom_type(atom_type)
myMan8.define_torsionals(hydrogens=False)

#Automatically detect all sequons and N-linked glycans
print '#'*50
print 'Loading glycoprotein'
myGlycosylator.load_glycoprotein(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/examples/env_4tvp.pdb'))
myGlycosylator.build_glycan_topology(patch = 'NGLB')
#Identify all glycans
print 'Detected glycans'
print '-'*50
for sequon,glycan in myGlycosylator.glycanMolecules.items():
    print sequon, '->', myGlycosylator.identify_glycan(glycan) 

print '#'*50
print 'Modifying sequon N262'
#Modify N-linked glycan at sequon ASN 262
#key: segn,chid,resid,icode
key = ',G,262,'
#ASN AtomGroup
sequon = myGlycosylator.get_residue(key)
#Get AtomGroup
original_glycan = myGlycosylator.glycanMolecules[key].atom_group
#Get Topology graph
r,t = myGlycosylator.glycans[key]
#Build topology tree from Topology graph
connect_tree = myGlycosylator.build_connectivity_tree(r, t)
segname = 'G'+key.split(',')[2]
glycan,bonds = myGlycosylator.glycosylate('MAN9_3;4,2', 
                                                template_glycan_tree = connect_tree, 
                                                template_glycan = original_glycan,
                                                link_residue=sequon, link_patch = 'NGLB', chain =  'G', segname=segname)
new_glycan = gl.Molecule(key)
new_glycan.set_AtomGroup(glycan, bonds = bonds)
myGlycosylator.assign_patches(new_glycan)
#Create a dictionary for storing new glycans
linked_glycanMolecules = {}
linked_glycanMolecules[key] = new_glycan
#Update glycan and save pdb
print 'Sequon new glycoprotein as HIV_N262_man9.pdb'
myGlycosylator.glycanMolecules.update(linked_glycanMolecules)
myGlycosylator.write_glycoprotein('HIV_N262_man9.pdb')
