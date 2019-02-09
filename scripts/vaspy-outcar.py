#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script file for use (demonstrate) vaspy.outcar module
"""

import argparse
from vaspy.outcar import OUTCAR
from vaspy import tools

parser = argparse.ArgumentParser(
    description='Show position/force evolution in VASP calculation',
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("-x", "--posx", action="store_true",
                    help="X-axis position", default=False)
parser.add_argument("-y", "--posy", action="store_true",
                    help="Y-axis position", default=False)
parser.add_argument("-z", "--posz", action="store_true",
                    help="Z-axis position", default=False)
parser.add_argument("-X", "--forcex", action="store_true",
                    help="Force along X-axis", default=False)
parser.add_argument("-Y", "--forcey", action="store_true",
                    help="Force along Y-axis", default=False)
parser.add_argument("-Z", "--forcez", action="store_true",
                    help="Force along Z-axis", default=False)
parser.add_argument("--site", type=tools.parse_AtomselectionNum,
                    help="""numbers deliminated by comma and hyphen
for SITE you want to see.  If not specified, all atoms are selected""")
parser.add_argument("outcarfiles", type=OUTCAR, nargs="+", metavar="OUTCAR",
                    help="OUTCAR_file(s)")
parser.add_argument("--plot", action="store_true", default=False,
                    help="""Show plot. Require matplotlib.
If this option is not selected,
the results are written into standard output.""")

args = parser.parse_args()

pff = [args.posx, args.posy, args.posz,
       args.forcex, args.forcey, args.forcez]  # pff is "posforce_flag"
if not any(pff):
    pff = [True] * 6

outcar = args.outcarfiles.pop(0)
if args.site == None:
    args.site = list(range(1, outcar.nions + 1))

headers = outcar.select_posforce_header(pff, args.site)
posforce = outcar.select_posforce(pff, args.site)
for outcar in args.outcarfiles:
    posforce.extend(outcar.select_posforce(pff, args.site))

if args.plot:
    import matplotlib.pyplot as plt
    for name, data in zip(headers, zip(*posforce)):
        plt.plot(data, label=name)
    plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left',
               ncol=len(headers) // 18 + 1)
    plt.show()
else:
    print("\t".join(headers))
    for item in posforce:
        print("\t".join([str(i) for i in item]))
