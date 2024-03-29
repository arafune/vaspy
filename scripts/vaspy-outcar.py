#!/usr/bin/env python
"""Script file for use (demonstrate) vaspy.outcar module."""

import argparse

import numpy as np
from vaspy import tools
from vaspy.outcar import OUTCAR

parser = argparse.ArgumentParser(
    description="Show position/force evolution in VASP calculation",
    formatter_class=argparse.RawTextHelpFormatter,
)

parser.add_argument(
    "-x",
    "--posx",
    action="store_true",
    help="X-axis position",
    default=False,
)
parser.add_argument(
    "-y",
    "--posy",
    action="store_true",
    help="Y-axis position",
    default=False,
)
parser.add_argument(
    "-z",
    "--posz",
    action="store_true",
    help="Z-axis position",
    default=False,
)
parser.add_argument(
    "-X",
    "--forcex",
    action="store_true",
    help="Force along X-axis",
    default=False,
)
parser.add_argument(
    "-Y",
    "--forcey",
    action="store_true",
    help="Force along Y-axis",
    default=False,
)
parser.add_argument(
    "-Z",
    "--forcez",
    action="store_true",
    help="Force along Z-axis",
    default=False,
)

parser.add_argument(
    "-r",
    "--relative",
    action="store_true",
    default=False,
    help="""In plot, the value is relative to the initial one.""",
)
parser.add_argument(
    "--site",
    type=tools.atom_selection_to_list,
    help="""numbers deliminated by comma and hyphen
for SITE you want to see.  If not specified, all atoms are selected""",
)
parser.add_argument(
    "outcarfiles",
    type=OUTCAR,
    nargs="+",
    metavar="OUTCAR",
    help="OUTCAR_file(s)",
)
parser.add_argument(
    "--plot",
    action="store_true",
    default=False,
    help="""Show plot. Require matplotlib.
If this option is not selected,
the results are written into standard output.""",
)

args = parser.parse_args()

pff = [
    args.posx,
    args.posy,
    args.posz,
    args.forcex,
    args.forcey,
    args.forcez,
]  # pff is "posforce_flag"
if not any(pff):
    pff = [True] * 6

outcar = args.outcarfiles.pop(0)
if args.site is None:
    args.site = list(range(1, outcar.n_atom + 1))

headers = outcar.select_posforce_header(pff, args.site)
posforce = outcar.select_posforce(pff, args.site)
for outcar in args.outcarfiles:
    posforce.extend(outcar.select_posforce(pff, args.site))

if args.plot:
    import matplotlib.pyplot as plt

    for name, data in zip(headers, zip(*posforce, strict=True), strict=True):
        if args.relative:
            data = np.array(data) - data[0]
        plt.plot(data, label=name)
    plt.legend(
        bbox_to_anchor=(-0.1, -0.1),
        loc="upper left",
        borderaxespad=0,
        ncol=len(headers) // 10 + 1,
    )
    plt.tight_layout()
    plt.show()
else:
    print("\t".join(headers))
    for item in posforce:
        print("\t".join([str(i) for i in item]))
