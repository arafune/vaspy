#! /usr/bin/env python
# -*- coding: utf-8 -*-
'''Script for manipulate XDATCAR files.'''

import argparse
import copy

from vaspy.poscar import POSCAR
from vaspy.xdatcar import XDATCAR

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--merge',
                    action='store_true',
                    default=False,
                    help='''Merge some XDATCAR files those''')
parser.add_argument('--split',
                    metavar='basenamePOSCAR',
                    help='''Split into several POSCAR files
The basename of POSCAR is provided by the argument''')
parser.add_argument('--poscar',
                    action='store_true',
                    default=False,
                    help='export into XDATCAR file from POSCAR files')
parser.add_argument("files",
                    nargs="+",
                    metavar="XDATCAR or POSCAR files",
                    help="XDATCAR file(s)")
#
args = parser.parse_args()
#
# args.merge, args.split, args.poscar はどれか一個。
#
assert (([args.merge, args.poscar].count(True) == 1 and args.split is None)
        or ([args.merge, args.poscar].count(True) == 0
            and args.split is not None)), "Option error"
#
if args.merge:
    xdatcars = []
    for xdatcar_file in args.files:
        xdatcars.append(XDATCAR(xdatcar_file))
    output_xdatcar = copy.deepcopy(xdatcars[0])
    configurations = []
    for xdatcar in xdatcars:
        configurations += xdatcar.configurations
    output_xdatcar.configurations = configurations
    print(output_xdatcar)
if args.split:
    assert len(args.files) == 1, '--split option takes the single XDATCAR file'
    xdatcar = XDATCAR(args.files[0])
    scaling_factor = xdatcar.scaling_factor
    cell_vecs = xdatcar.cell_vecs
    iontypes = xdatcar.iontypes
    ionnums = xdatcar.ionnums
    for frame, configuration in enumerate(xdatcar.configurations):
        poscar = POSCAR()
        poscar.system_name = xdatcar.system_name + '_frame_' + str(frame + 1)
        poscar.scaling_factor = scaling_factor
        poscar.cell_vecs = cell_vecs
        poscar.iontypes = iontypes
        poscar.ionnums = ionnums
        poscar.coordinate_type = 'Direct'
        poscar.positions = configuration
        poscar.save(args.split + '_' + str(frame + 1) + '.vasp')
if args.poscar:
    poscars = []
    for poscar_file in args.files:
        poscars.append(POSCAR(poscar_file))
    output_xdatcar = XDATCAR()
    output_xdatcar.system_name = poscars[0].system_name
    output_xdatcar.scaling_factor = poscars[0].scaling_factor
    output_xdatcar.cell_vecs = poscars[0].cell_vecs
    output_xdatcar.iontypes = poscars[0].iontypes
    output_xdatcar.ionnums = poscars[0].ionnums
    positions = []
    for poscar in poscars:
        poscar.to_direct
        positions.append(poscar.positions)
    output_xdatcar.configurations = positions
    print(output_xdatcar)
