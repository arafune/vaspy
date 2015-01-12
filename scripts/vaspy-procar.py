#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
script to use (demonstrate) vaspy.procar (Ver.3) module.
"""

import re
import functools as ft
import argparse
from vaspy import tools
from vaspy.outcar import OUTCAR
import vaspy.procar as procar

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--outcar', metavar='outcar_file',
                    help='''Use "OUTCAR" for the Fermi level correction.
outcar_file must be specified.
NOTE:  E-fermi of OUTCAR generated in
Band-calculation may NOT be reliable.''')
parser.add_argument('--fermi', metavar='value', type=float,
                    help='''Fermi level correction
Energy shifts by this value
if --outcar is set, this option is ignored''')
parser.add_argument('--site', metavar='atom_indices', dest='atomindex',
                    action='append',
                    type=tools.parse_AtomselectionNum,
                    help='''atom index specifed with range.
Use "-" or ","
 (ex.) --site 1,2,7-9''')
parser.add_argument('--as', metavar='name', nargs='+', dest='atomsetname',
                    action='append',
                    help='''the name of the sites identified
by --site option
the name is used in the title of the column''')
parser.add_argument('--orbital', metavar='orbitals', action='append',
                    type=ft.partial(re.split, r'[,:]'),
                    help='''orbital name deliminated by ":" or ",".
orbital names are:
s, p, pxpy, pz, d, dxy, dyz, dz2, dxz, dx2
 (ex.) --orbital s:pxpy:d''')
parser.add_argument('procar', metavar='PROCAR_file',
                    help='''PROCAR file''')

args = parser.parse_args()

# ---
if not (len(args.atomindex) == len(args.orbital) == len(args.atomsetname)):
    raise parser.error("--atom, --as and --orbital are mismatched.")
# ---

if args.outcar is not None:
    outcar = OUTCAR(args.outcar)
    fermi = outcar.fermi
elif args.fermi is not None:
    fermi = args.fermi
else:
    fermi = 0.0

sitenames = tuple(set([ e for inner in args.atomsetname for e in inner]))
flat_orbitals = tuple(set([ e for inner in args.orbital for e in inner]))

# As atomindex used here begins with "1", but siteindex used in procar.py
# internaly begins with "0".  (This is because VASP is fortran program !)
siteindex = [[ i-1 for i in internal] for internal in args.atomindex]

procar = procar.PROCAR(args.procar)
band = procar.band()
del procar  # for memory saving
if fermi != 0.0:
    band.energies -= fermi
for sites in siteindex:
    band.compose_sites(sites)
band.compose_orbital(flat_orbitals)
print (band.get_sitecomposed_data(sitenames, args.orbital))
