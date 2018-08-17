#! /usr/bin/env python
# -*- coding: utf-8 -*-
'''
Script for manipulate XDATCAR files
'''


import argparse
import copy
from vaspy.poscar import POSCAR
from vaspy.xdatcar import XDATCAR

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--merge',action='store_true', default=False,
                    help='''Merge some XDATCAR files those''')
parser.add_argument('--split', metavar='basenamePOSCAR',
                   help='''Split into several POSCAR files
The basename of POSCAR is provided by the argument''')
parser.add_argument("xdatcarfiles", nargs="+", metavar="XDATCAR files",
                    help="XDATCAR file(s)")
#                    
args = parser.parse_args()
#
if args.merge:
    xdatcars = []
    for xdatcar_file in args.xdatcarfiles:
        xdatcars.append(XDATCAR(xdatcar_file))
    output_xdatcar = copy.deepcopy(xdatcars[0])
    configurations = []
    for xdatcar in xdatcars:
        configurations += xdatcar.configurations
    output_xdatcar.configurations = configurations
    print(output_xdatcar)
