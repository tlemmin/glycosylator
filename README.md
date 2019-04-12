# Glycosylator: a Python framework for the rapid modeling of glycans



## Installation

## Computing Environment
The following packages are required:
```
prody
numpy
matplotlib
networkx
```

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

