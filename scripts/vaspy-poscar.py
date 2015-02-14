#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
script to use(demonstrate) vaspy.poscar functions.
"""

import argparse
import functools as ft
from vaspy import tools
from vaspy.poscar import POSCAR


def split_to_float(string, n, name):
    lis = string.split(',')
    if len(lis) != n:
        message = '--{0} option requires {1} numbers'.format(name, n)
        raise argparse.ArgumentTypeError(message)
    return [float(i) for i in lis]


parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    epilog="""
NOTE: When you run this script on Windows Power Shell,
commas are regarded as delimiter of values.
So you must enclose values which contains commas with quotations.
(ex.) --site 1,2,3 -> failure / --site "1,2,3" -> OK""")
parser.add_argument('--site', metavar='atoms', action='append',
                    type=tools.parse_AtomselectionNum,
                    help='''atoms specified with range using "-"
or comma-delimnated numbers.
 (ex.) --site 1,2,7-9''')

group = parser.add_mutually_exclusive_group(required=False)
group.add_argument('--translate', metavar='x,y,z', action='append',
                   type=ft.partial(split_to_float, n=3, name='translate'),
                   help='''displacement (AA unit) by three numbers
separated by comma.''')
group.add_argument('--rotateX', metavar='theta,x,y,z',
                   type=ft.partial(split_to_float, n=4, name='rotateX'),
                   help='''Rotation around X-axis by "theta" at (x,y,z)
NOTE: this option is not taken into account
the periodic boundary.''')
group.add_argument('--rotateY', metavar='theta,x,y,z',
                   type=ft.partial(split_to_float, n=4, name='rotateY'),
                   help='''Rotation around Y-axis by "theta" at (x,y,z)
NOTE: this option is not taken into account
the periodic boundary.''')
group.add_argument('--rotateZ', metavar='theta,x,y,z',
                   type=ft.partial(split_to_float, n=4, name='rotateZ'),
                   help='''Rotation around Z-axis by "theta" at (x,y,z)
NOTE: this option is not taken into account
the periodic boundary.''')

parser.add_argument('--output', metavar='file_name',
                    help='''output file name
if not specified, use standard output''')
parser.add_argument('poscar', metavar='POSCAR_file (or CONTCAR_file)',
                    type=POSCAR)

args = parser.parse_args()

# translate option and rotate option are not set simulaneously.
if args.translate and any([args.rotateX,
                           args.rotateY,
                           args.rotateZ]):
    parser.error("Cannot set --translate and rotate option simultanaously.")
# rotate options are not set multiply.
if (args.rotateX, args.rotateY, args.rotateZ).count(None) < 2:
    parser.error("Cannot set multiple rotate options simultanaously.")
#

############

# print(args.poscar) #DEBUG
args.poscar.to_Cartesian()
# print(args.poscar) #DEBUG

#
#  if "atom" option is not set, all atoms are concerned.
#
if not args.site:
    nAtoms = sum(args.poscar.ionnums)
    args.site = [tools.parse_AtomselectionNum('1-{0}'.format(nAtoms))]
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
        axis_name = 'X'
        theta = args.rotateX[0]
        center = args.rotateX[1:]
    elif args.rotateY:
        axis_name = 'Y'
        theta = args.rotateY[0]
        center = args.rotateY[1:]
    elif args.rotateZ:
        axis_name = 'Z'
        theta = args.rotateZ[0]
        center = args.rotateZ[1:]
    args.poscar.atoms_rotate(args.site, axis_name, theta, center)

#
#  Output result
#
if args.output is not None:
    args.poscar.save(args.output)
else:
    print(args.poscar)
