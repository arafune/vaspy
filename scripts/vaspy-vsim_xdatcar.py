#! /usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Script to convert V_Sim ascii file to XDATCAR
'''
import numpy as np
import vaspy.vsim_asc as vsim_asc
import argparse
#
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--dim', type=str,
                    help='dimension of supercell, such as 2,2,1')
parser.add_argument('--frames', type=int, default=30,
                    help='Total number of frames, default is 30')
parser.add_argument('--mode', type=int, default=0,
                    help='phonon mode to calculate, default is 0')
parser.add_argument('--magnitude', type=float, default=1,
                    help='magnitude of the amplitude of the motion')
parser.add_argument('vsim', metavar='VSIM_file',
                    help='''VSIM ascii file''')
#
args = parser.parse_args()
#
dim = [int(x) for x in  args.dim.split(',')]
vsim_data = vsim_asc.VSIM_ASC(args.vsim)
n_modes = len(vsim_data.freqs)
positions = vsim_data.build_phono_motion(mode=args.mode, supercell=dim,
                                        n_frames=args.frames, magnitude=args.magnitude)
ionnums, iontypes = vsim_asc.ions_to_iontypes_ionnums(vsim_data.ions)
ionnums = [i * dim[0] * dim[1] * dim[2] for i in ionnums]
motion_frames = [[p[frame] for p in positions]
                 for frame in range(args.frames)]
#
supercell = vsim_asc.supercell_lattice_vectors(vsim_data.lattice_vectors, dim)
#
xdatcar_str = vsim_data.system_name + '_mode_'+str(args.mode)+'\n'
xdatcar_str += '       1.00\n'
for i in range(3):
    xdatcar_str += '      {:#.6f}   {:#.6f}    {:6f}\n'.format(
        supercell[i][0], supercell[i][1], supercell[i][2])
for element in iontypes:
    xdatcar_str += '    {}'.format(element)
xdatcar_str += '\n'
for ionnum in ionnums:
    xdatcar_str += '    {}'.format(ionnum)
xdatcar_str += '\n'
mat = np.linalg.inv(np.transpose(supercell))
for frame in range(1, args.frames+1):
    xdatcar_str += 'Direct configuration=    {}\n'.format(frame)
    positions = [mat.dot(v) for v in motion_frames[frame-1]]
    for position in positions:
        xdatcar_str += '    {:#.6f}    {:#.6f}    {:6f}\n'.format(
            position[0], position[1], position[2])

print(xdatcar_str)
