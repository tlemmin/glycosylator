#!/usr/bin/env python

from prody import *
import os, sys, argparse

def readfile(fname):
    with open(fname) as f:
        content = [x.strip('\n') for x in f.readlines()]
    return content

def main():
    par = argparse.ArgumentParser()
    par.add_argument('--p', required = True, help = 'input PDB file')
    par.add_argument('--i', required = True, help = 'list of internal coordinates')
    par.add_argument('--o', required = True, help = 'internal coordinate file')
    args = par.parse_args()
    
    structure = parsePDB(args.p)
    ics = readfile(args.i)
    to_write = []
    for ic in ics:
        a1,a2,a3,a4 = ic.split()
        internal_coordinate = 'IC %s %s %s %s' % (a1,a2,a3,a4)
        a1 = structure.select('name ' + a1)
        a2 = structure.select('name ' + a2)
        a3 = structure.select('name ' + a3)
        a4 = structure.select('name ' + a4)
        bond1 = calcDistance(a1, a2)[0]
        angle1 = calcAngle(a1, a2, a3)[0]
        dihe = calcDihedral(a1, a2, a3, a4)[0]
        angle2 = calcAngle(a2, a3, a4)[0]
        bond2 = calcDistance(a3, a4)[0]
        internal_coordinate += '  %.2f %.2f %.2f %.2f %.2f' % (bond1, angle1, dihe, angle2, bond2)
        to_write.append(internal_coordinate)

    f = open(args.o, 'w')
    f.write('\n'.join(to_write))
    f.close


if __name__ == "__main__":
    main()
