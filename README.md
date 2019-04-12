
![alt text](./icons/glycosylator_logo.png)
# Glycosylator: a Python framework for the rapid modeling of glycans
Glycosylator is a Python framework for the identification, modeling and
modification of glycans in protein structure. It can be used directly in a Python script
through its API or through its Graphical User Interface (GUI). The GUI provides a
straightforward 2D rendering of a glycoprotein that allows for a quick visual inspection of
the gylcosylation state of all the sequons on a protein structure. Modeled glycans can
be further refined by a genetic algorithm for removing clashes and sampling alternative
conformations. Glycosylator can also identify specific 3D glycans on a protein structure
using a library of predefined templates.
Glycosylator has been implemented in a generic way allowing the user to expand the library to incorporate other polymers.

## Computing Environment
Glycosylator was developed with the following packages:
```
prody == 1.9.4
numpy == 1.15.1
matplotlib == 2.1.2
networkx == 2.1
```
## Running glycosylator
To open the GUI:
`python glycosylator_GUI.py`

## Demo
The demo folder contains a 7 examples, showing how to use the different classes provided by Glycosylator.

- The Molecule class is used for storing the coordinates (Prody’s AtomGroup) and
connectivity (NetworkX graph) of a molecule.
- The MoleculeBuilder class is employed for building and editing molecules.
- Glycosylator class was designed to deal specifically with glycans/glycoprotein.
Sampler class implements a genetic algorithm for removing clashes between
- Molecules and their environment (e.g. protein).
- Drawer class is used for generating 2D symbolic representations of glycans
according to the UIPAC standard.

## Glycan text representation
A UNIT is defined by its residue name (3 letter code from PDB), the connecting atom, and the path from the root UNIT. The root UNIT is the first glycan is connected to the sequon.
The path is defined by the patches that have applied to connect previous UNIT.
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

