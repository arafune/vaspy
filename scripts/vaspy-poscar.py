#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
script to use(demonstrate) vaspy.poscar functions.
"""

import argparse
import functools as ft
from logging import DEBUG, INFO, Formatter, StreamHandler, getLogger
from typing import List
from vaspy import tools
from vaspy.poscar import POSCAR

LOGLEVEL = INFO
logger = getLogger(__name__)
fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
formatter = Formatter(fmt)
handler = StreamHandler()
handler.setLevel(LOGLEVEL)
logger.setLevel(LOGLEVEL)
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False


def split_to_float(string: str, n: int, name: str) -> List[float]:
    lis = string.split(",")
    if len(lis) != n:
        message = "--{0} option requires {1} numbers".format(name, n)
        raise argparse.ArgumentTypeError(message)
    return [float(i) for i in lis]


parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    epilog="""
NOTE: When you run this script on Windows Power Shell,
commas are regarded as delimiter of values.
So you must enclose values which contains commas with quotations.
(ex.) --site 1,2,3 -> failure / --site "1,2,3" -> OK""",
)
parser.add_argument(
    "--site",
    metavar="atoms",
    action="append",
    type=tools.atom_selection_to_list,
    help="""atoms specified with range using "-"
or comma-delimnated numbers.
(ex.) --site 1,2,7-9""",
)

group = parser.add_mutually_exclusive_group(required=False)
group.add_argument(
    "--translate",
    metavar="x,y,z",
    action="append",
    type=ft.partial(split_to_float, n=3, name="translate"),
    help="""displacement (AA unit) by three numbers
separated by comma.""",
)
group.add_argument(
    "--rotateX",
    metavar="theta,x,y,z",
    type=ft.partial(split_to_float, n=4, name="rotateX"),
    help="""Rotation around X-axis by "theta" at (x,y,z)
NOTE: this option is not taken into account
the periodic boundary.""",
)
group.add_argument(
    "--rotateY",
    metavar="theta,x,y,z",
    type=ft.partial(split_to_float, n=4, name="rotateY"),
    help="""Rotation around Y-axis by "theta" at (x,y,z)
NOTE: this option is not taken into account
the periodic boundary.""",
)
group.add_argument(
    "--rotateZ",
    metavar="theta,x,y,z",
    type=ft.partial(split_to_float, n=4, name="rotateZ"),
    help="""Rotation around Z-axis by "theta" at (x,y,z)
NOTE: this option is not taken into account
the periodic boundary.""",
)

parser.add_argument(
    "--output",
    metavar="file_name",
    help="""output file name
if not specified, use standard output""",
)
parser.add_argument("poscar", metavar="POSCAR_file (or CONTCAR_file)", type=POSCAR)
parser.add_argument(
    "--to_direct", action="store_true", help="""Change direct coordinates"""
)
parser.add_argument(
    "--to_cartesian", action="store_true", help="""Change cartesian coordinates"""
)
parser.add_argument(
    "--split",
    nargs=2,
    help="""Split into two POSCAR files.
The first is the file name for one POSCAR (molecule part).
The second is the file name for other POSCAR (substrate part).""",
)
#
args = parser.parse_args()
logger.debug("args: {}".format(args))

# translate option and rotate option are not set simulaneously.
if args.translate and any([args.rotateX, args.rotateY, args.rotateZ]):
    parser.error("Cannot set --translate and rotate option simultanaously.")
# rotate options are not set multiply.
if (args.rotateX, args.rotateY, args.rotateZ).count(None) < 2:
    parser.error("Cannot set multiple rotate options simultanaously.")
#

############
args.poscar.to_cartesian()

#
#  if "atom" option is not set, all atoms are concerned.
#
if not args.site:
    args.site = [
        tools.atom_selection_to_list("1-{0}".format(sum(args.poscar.atomnums)))
    ]
args.site = [i - 1 for i in args.site[0]]
#
#  Translation
#
if args.translate:
    if len(args.site) != len(args.translate):
        parser.error()
    for v, a in zip(args.translate, args.site):
        args.poscar.translate(v, a)
#
#  Rotation
#
if any([args.rotateX, args.rotateY, args.rotateZ]):
    if len(args.site) != 1:
        parser.error("--site option set once!")
    if args.rotateX:
        axis_name = "X"
        theta = args.rotateX[0]
        center = args.rotateX[1:]
    elif args.rotateY:
        axis_name = "Y"
        theta = args.rotateY[0]
        center = args.rotateY[1:]
    elif args.rotateZ:
        axis_name = "Z"
        theta = args.rotateZ[0]
        center = args.rotateZ[1:]
    args.poscar.atoms_rotate(args.site, axis_name, theta, center)

if args.to_direct and args.to_cartesian:
    raise ValueError("Error!!  Set either of --to_direct or --to_cartesian")

if args.to_direct:
    args.poscar.to_direct()

if args.to_cartesian:
    args.poscar.to_cartesian()
#
#  Output result
#
if args.output is not None:
    args.poscar.save(args.output)
elif args.split:
    logger.debug("args.site: {}".format(args.site))
    one, other = args.poscar.split(args.site)
    one.save(args.split[0])
    other.save(args.split[1])
else:
    print(args.poscar)
