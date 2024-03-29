#! /usr/bin/env python
"""Script to use  vaspy.wavecar module."""

import argparse

from vaspy import poscar, wavecar

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--poscar",
    "-p",
    metavar="poscar_file",
    help="""POSCAR file  (default is "POSCAR")""",
    default="POSCAR",
)
parser.add_argument(
    "--spin",
    "-s",
    metavar="spin_index",
    type=int,
    default=0,
    help="""Spin index (0 or 1: default is 0)""",
)
parser.add_argument(
    "--k",
    "-k",
    metavar="k_index",
    type=int,
    default=1,
    help="""k index (>=1)  (default is 1)""",
)
parser.add_argument(
    "--band",
    "-b",
    metavar="band_index",
    type=int,
    default=1,
    help="""band index (>=1)""",
)
parser.add_argument(
    "--output",
    metavar="output_filename",
    type=str,
    help="""filename for output.
filename is 'output_s_(spin_index)_k_(k_index)_b_(bandindex).vasp'
If notspecified, standard output is used.
""",
)
parser.add_argument(
    "wavecar",
    metavar="WAVECAR_file",
    default="WAVECAR",
    help='WAVECAR file   (default is "WAVECAR")',
)

args = parser.parse_args()

try:
    pos: poscar.POSCAR = poscar.POSCAR(args.poscar)
except NameError:
    print("'POSCAR' file not exist")

wav: wavecar.WAVECAR = wavecar.WAVECAR(args.wavecar)
grid = wav.realspace_wfc(
    spin_i=args.spin,
    k_i=args.k - 1,
    band_i=args.band - 1,
    poscar=pos,
)
#
#
if args.output:
    filename = "{}_s_{:01d}_k_{:003d}_b_{:003d}.vasp".format(
        args.output,
        args.spin,
        args.k - 1,
        args.band - 1,
    )
    grid.save(filename)
else:
    print(grid)
