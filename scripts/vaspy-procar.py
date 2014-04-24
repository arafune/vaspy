#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
script to use (demonstrate) vaspy.procar module.
"""

import re
import functools as ft
import argparse
from vaspy import tools
from vaspy.outcar import OUTCAR
from vaspy.procar import PROCAR, Band

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
parser.add_argument('--site', metavar='atom_indices', dest='atomname',
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
if not (len(args.atomname) == len(args.orbital) == len(args.atomsetname)):
    raise parser.error("--atom, --as and --orbital are mismatched.")
# ---

if args.outcar is not None:
    outcar = OUTCAR(args.outcar)
    fermi = outcar.fermi
elif args.fermi is not None:
    fermi = args.fermi
else:
    fermi = 0.0

procar = PROCAR(args.procar)
band = procar.to_band()
    
if fermi != 0.0:
    band.fermilevel_correction(fermi)
    
output_band = Band()
for site, name, orbital in zip(args.atomname,
                               args.atomsetname,
                               args.orbital):
    tmpBand = band.site_integrate(site)
    tmpBand.rename_site(*name)
    tmpBand.extract_orbitals_in_place(orbital)
    output_band += tmpBand
output_band.sort()
print(output_band)
