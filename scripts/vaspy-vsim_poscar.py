#! /usr/binenv python3
# -*- coding: utf-8 -*-
'''
Script to convert V_Sim ascii file to POSCAR
'''

import logging
logging.basicConfig(level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s -%(message)s')
logging.disable(logging.DEBUG)
#
import vaspy.poscar as poscar
import vaspy.vsim_asc as vsim_asc
import argparse
#
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--dim', type=str,
                    help='dimension of supercell, such as 2,2,1')
parser.add_argument('--frames', type=int,
                    help='Total number of frames, default is 30')
parser.add_argument('--mode', type=int,
                    help='phonon mode to calculate, default is 0')
parser.add_argument('--magnitude', type=float,
                    help='magnitude of the amplitude of the motion')
parser.add_argument('vsim', metavar='VSIM_file',
                    help='''VSIM ascii file''')
#
args = parser.parse_args()
#
dim = [int(x) for x in  args.dim.split(',')]
if not args.frames:
    args.frames = 30
if not args.magnitude:
    args.magnitude=1
vsim_asc = vsim_asc.VSIM_ASC(args.vsim)
n_modes = len(vsim_asc.freqs)
positions = vsim_asc.build_phono_motion(mode=args.mode, supercell=dim,
                                        n_frames=args.frames, magnitude=args.mag)
ionnums, iontypes = vsim_asc.ions_to_iontypes_ionnums(vsim_asc.ions)
ionnums = [i * dim[0] * dim[1] * dim[2] for i in ionnums]
motion_frames = [[p[i] for p in positions] for in range(args.frames)]
#
supercell = vsim_asc.supercell_lattice_vectors(vsim_asc.lattice_vectors, dim)
#
for frame in range(args.frames):
    vasp_poscar = poscar.POSCAR()
    vasp_poscar.system_name = vsim_asc.system_name + '_mode_'+str(args.mode)+'_frame_'+str(frame)
    vasp_poscar.scaling_factor = 1.0
    vasp_poscar.cell_vecs = supercell    
    vasp_poscar.iontypes = iontypes
    vasp_poscar.ionnums = ionnums
    vasp_poscar.coordinate_type = 'CARTESIAN'
    vasp_poscar.positions = motion_frames[frame]
    vaspy_poscar.save('POSCAR_mode_'+str(args.mode)+'_frame_'+str(frame)+'.vasp')
