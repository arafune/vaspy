#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
script to use(demonstrate) vaspy.chgcar module
"""

import argparse
from vaspy.chgcar import CHGCAR

# HINT list methods for --spin option below.
# they are called with just one argument, CHGCAR instance.
spinmethods = {
        'mag':CHGCAR.magnetization,
        'magX':(lambda chg: chg.magnetiztion('x')),
        'magY':(lambda chg: chg.magnetiztion('y')),
        'magZ':(lambda chg: chg.magnetiztion('z')),
        'majority':CHGCAR.majorityspin,
        'minority':CHGCAR.minorityspin,
               }


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--add', action='store_true', default=False,
                 help="Add two CHGCAR files")
group.add_argument('--diff', action='store_true', default=False,
                 help="Get difference of two CHGCAR files")
group.add_argument('--spin', metavar='spin_operation',
                 help="""spin-relatated operation.
when this option is set --add, -diff are ignored,
and CHGCAR_file_2 must not be set.
spin operation is one of the followings:
mag : show the magnetisation
      density (for spin resolved calculations)
magX : show the magnetisation density in
       the X direction (for non collinear calc.)
magY : show the magnetisation density in
       the Y direction (for non collinear calc.)
magZ : show the magnetisation density in
       the Z direction (for non collinear calc.)
majority : extract the part for the
           majority spin (for spin resolved calc.)
minority : extract the part for the
           inority spin (for spin resolved calc.)""")
parser.add_argument('--output', metavar='file_name',
                 help="""output file name
if not specified, use standard output""")
parser.add_argument('CHGCAR_file_1', type=CHGCAR)
parser.add_argument('CHGCAR_file_2', type=CHGCAR, nargs='?')
#TIP: if CHGCAR_file_2 is not specified,
#*None* is stored in arguments.CHGCAR_file_2, not CHGCAR(None)


args = parser.parse_args()

#
if args.spin is not None:
    if args.CHGCAR_file_2 is not None:
        parser.error("Only one CHGCAR file for --spin operations")
    if args.spin in spinmethods:
        c = spinmethods[args.spin](args.CHGCAR_file_1)
    else:
        parser.error("Such spin operation parameter is not defined.")
#
if args.add or args.diff:
    if args.CHGCAR_file_2 is None:
        raise parser.error('Two CHGCAR files are required.')
    if args.add:
        c = args.CHGCAR_file_1 + args.CHGCAR_file_2
    else:
        c = args.CHGCAR_file_1 - args.CHGCAR_file_2
#
if args.output is not None:
    c.save(args.output)
else:
    print(c)