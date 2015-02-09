#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
script to use(demonstrate) vaspy.doscar
"""

import argparse
import os
import sys
mypath = os.readlink(__file__) if os.path.islink(__file__) else __file__
sys.path.append(os.path.dirname(os.path.abspath(mypath)))
from vaspy.outcar import OUTCAR
from vaspy.doscar import DOSCAR, TDOS, PDOS
from vaspy import tools

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('doscar', metavar='DOSCAR_file')
parser.add_argument('-o', '--outcar', metavar='OUTCAR_file',
                    nargs='?', default=NotImplemented,
                    help="""Use OUTCAR file.
If file path is not given,
try to open OUTCAR file
in the same directory to DOSCAR.""")
# TIPS
# -o data/OUTCAR => args.outcar == "data/OUTCAR"
# -o             => args.outcar == None
# (not given)    => args.outcar == NotImplemented (default)
parser.add_argument('-f', '--fermi', metavar='value', type=float,
                    default=0.0,
                    help="""Fermi level correction.
Energy shifts by this value.
if --outcar is set, this option is ignored.""")
parser.add_argument('-s', '--site', metavar='atoms', action='append',
                    dest='atomset', type=tools.parse_AtomselectionNum,
                    help="""atom # specified with range.
Use "-" or ","
(ex.) --site 1,2,7-9""")
parser.add_argument('-a', '--as', metavar='name', action='append',
                    dest='atomsetname',
                    help="""the name of the range identified by --site.
(ex.) --as layer1
the name is used in the output filename.""")

args = parser.parse_args()

doscar = DOSCAR(args.doscar)
atomlist = list()
if args.outcar is not NotImplemented:
    if args.outcar is None:
        args.outcar = os.path.join(os.path.dirname(args.doscar), "OUTCAR")
    outcar = OUTCAR(args.outcar)
    atomlist = outcar.atom_names
    args.fermi = outcar.fermi
#
if atomlist == []:
    atomlist.extend("atom" + str(i) for i in range(1, doscar.nAtom + 1))
#
# construct TDOS & PDOS objects
#
tmp = doscar.dos_container
d = [TDOS(tmp.pop(0))]
#
d.extend(PDOS(*each) for each in zip(tmp, atomlist))  # tmp[1:] ?
#
if args.atomset is not None:
    if len(args.atomset) == len(args.atomsetname):
        sumPDOSs = list()
        for site, name in zip(args.atomset, args.atomsetname):
            each = PDOS()
            for atomNo in site:
                each += d[atomNo]
            each.site = name
            sumPDOSs.append(each)
        for summedPDOS in sumPDOSs:
            filename = summedPDOS.site + ".dat"
            try:  # Version safety
                file = open(filename, mode='w', newline='')
            except TypeError:
                file = open(filename, mode='wb')
            with file:
                summedPDOS.export_csv(file, delimiter='\t')
#
try:  # Version safety
    file = open("total.dat", mode='w', newline='')
except TypeError:
    file = open("total.dat", mode='wb')
with file:
    d[0].fermilevel_correction(args.fermi)
    d[0].export_csv(file, delimiter='\t')
for i, n in zip(d[1:], atomlist):
    i.fermilevel_correction(args.fermi)
    if isinstance(i, PDOS) and i.site == "":
        i.site = n
    filename = n + "_dos.dat"
    try:  # Version safety
        file = open(filename, mode='w', newline='')
    except TypeError:
        file = open(filename, mode='wb')
    with file:
        i.export_csv(file, delimiter='\t')
