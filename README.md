
# New implementation in python 3.0 available :)
https://github.com/ibmm-unibe-ch/glycosylator

Installation

```
pip install glycosylator
```


![alt text](./icons/glycosylator_logo.png)
# A Python framework for the rapid modeling of glycans
Glycosylator is a Python framework for the identification, modeling and
modification of glycans in protein structure. It can be used directly in a Python script
through its API or through its Graphical User Interface (GUI). The GUI provides a
straightforward 2D rendering of a glycoprotein that allows for a quick visual inspection of
the gylcosylation state of all the sequons on a protein structure. Modeled glycans can
be further refined by a genetic algorithm for removing clashes and sampling alternative
conformations. Glycosylator can also identify specific 3D glycans on a protein structure
using a library of predefined templates.
Glycosylator has been implemented in a generic way allowing the user to expand the library to incorporate other polymers.

Please cite:
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3097-6

## Dependencies
Glycosylator was developed with the following environment:
```
conda env create -f environment.yml
```
The dependencies for python2.7 are getting very old. networkx installed a newer version of `decoretor`, which is not compatible with python2.7. The older version can be found in the `support` folder.
## Usage 
To open the GUI:
```python glycosylator_GUI.py```
The GUI was developed on Linux, it seems to be a bit buggy on Mac OsX.

In Python script:
```
import glycosylator as gl
myGlycosylator = gl.Glycosylator(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.rtf'), os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.prm'))
myGlycosylator.builder.Topology.read_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/DUMMY.top'))
#Load topology information about glycans
myGlycosylator.read_connectivity_topology(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/topology/mannose.top'))
# 3. Automatically detect all sequons and N-linked glycans
print 'Loading glycoprotein'
myGlycosylator.load_glycoprotein(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/examples/env_4tvp.pdb'))
myGlycosylator.build_glycan_topology(patch = 'NGLB')
```

## Demo
The demo folder contains several examples, showing how to use the different classes provided by Glycosylator.

- The Molecule class is used for storing the coordinates (Prody’s AtomGroup) and
connectivity (NetworkX graph) of a molecule
- The MoleculeBuilder class is employed for building and editing molecules
- Glycosylator class was designed to deal specifically with glycans/glycoprotein and model glycans
- Sampler class implements a genetic algorithm for removing clashes between Molecules and their environment (e.g. protein)
- Drawer class is used for generating 2D symbolic representations of glycans according to the UIPAC standard

## Glycan text representation
Glycosylator uses a text representation to identify or build glycans.
A UNIT is defined by its residue name (3 letter code from PDB), the connecting atom, and the path from the root UNIT. The root UNIT is the first glycan that is connected to the sequon.
The path is defined by the patches (PRES in CHARMM ff) that are applied to connect to the previous UNIT.
Example for mannose 3:
```
                                         -16ab-> MAN[Z4] -13ab-> MAN[Z7]
NAG[Z1] -14bb-> NAG[Z2] -14bb-> BMA[Z3]-|
                                         -13ab-> MAN[Z7]    
```    
The corresponding text representation (connectivity topology):
```
    RESI MAN3                               !User defined, which has to be unique within a top file
    UNIT NAG                                !Z1 Root unit in glycan
    UNIT NAG C1 14bb                        !Z2 Unit connected to root with a 14bb patch
    UNIT BMA C1 14bb 14bb                   !Z3 Unit connected with 14bb to Z2 [path from root: 14bb]
    UNIT MAN C1 14bb 14bb 16ab              !Z4 Unit connected with 16ab to Z3 [path from root: 14bb 14bb]
    UNIT MAN C1 14bb 14bb 16ab 13ab         !Z7 Unit connected with 13ab to Z4 [path from root: 14bb 14bb 16ab]
    UNIT MAN C1 14bb 14bb 13ab              !Z9 Unit connected with 13ab to Z3 [path from root: 14bb 14bb]
```
## Adding new monosaccaride
Glycosylator uses the internal coordinates (IC) to build molecules. It requires the absence of circular dependency.
The easiest is to extract the coordinates from a PDB file of the optimized monomer. 
Two scripts are available for generating the coordinates: ./support/scripts/
 - buildICs.py: will compute all values for a list of ICs.
 - XML2PDB.py: extracts the structure from RCSB XML molecule file. 

## Setting up topology for Molecular Dynamics simulations
Glycosylator does not build directly the topology file for MD simulations. However, the output is fully compatible with several very good tools.
### VMD
Available from Theoretical and Computational Biophysics Group at UIUC: https://www.ks.uiuc.edu/Research/vmd/
1. autopsf: In most cases, glycans are dectected correctly and the topology file is generated automatically.
2. psfgen: Manually build the topology for glycans and combine the topology of the protein . Glycosylator provides a list of patches for connecting glycans. (see demo D10_psfgen)

### CHARMM-GUI
CHARMM-GUI requires the `CONECT` field for the glycans. @Vikasdubey0551 suggested two approaches (see issue 2):
#### Pymol
1. load glycan pdb to pymol
2. in pymol command line, write the following :
    > set pdb_retain_ids, on
    > set pdb_conect_all, off
    > save yourfile.pdb
3. This will write `CONECT` records in the end for only `HETATMs`.

#### Chimera
1. Upload glycan pdb to chimera by **clicking -> file -> open**.
2.  Go to **Select -> Chain -> choose chain containing glycans yourself**.
3. Go to **Tools-> Structure Editing -> Renumber Residues**.
4. Click on selected residues and renumber the glycan residues by typing a desired number.
5. Go to **File-> Save PDB -> save the pdb**.
6. This will write `CONECT` records in the end for only `HETATMs`.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)
