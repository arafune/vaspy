#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
script to use  vaspy.wavecar module.
"""

import argparse
import vaspy.wavecar as wavecar
import vaspy.poscar as poscar

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--poscar', '-p', metavar='poscar_file',
                    help='''POSCAR file  (default is "POSCAR")''',
                    default='POSCAR')
parser.add_argument('--spin', '-s', metavar='spin_index',
                    type=int, default=0,
                    help='''Spin index (0 or 1: default is 0)''')
parser.add_argument('--k', '-k', metavar='k_index',
                    type=int, default=1,
                    help='''k index (>1)  (default is 1)''')
parser.add_argument('--band', '-b', metavar='band_index',
                    type=int, default=1,
                    help='''band index (>1)''')
parser.add_argument('--output', metavar='output_filename',
                    type=str,
                    help='''filename for output.
filename is 'output_s_(spin_index)_k_(k_index)_b_(bandindex).vasp'
If notspecified, standard output is used.
''')
parser.add_argument('wavecar', metavar='WAVECAR_file',
                    default='WAVECAR',
                    help='WAVECAR file   (default is "WAVECAR")')

args = parser.parse_args()

try:
    pos = poscar(args.poscar)
except NameError:
    print("'POSCAR' file not exist")

wav = wavecar(args.WAVECAR)
grid = wav.realspace_wfc(spin_i=args.spin,
                         k_i=args.k-1,
                         band_i=args.band-1,
                         poscar=pos)
#
#
if args.output:
    filename = '{0}_s_{1:01d}_k_{2:003d}_b_{3:003d}.vasp'.format(
        args.output, args.spin, args.k-1, args.band-1)
    grid.save(filename)
else:
    print(grid)
