# -*- codinng: utf-8 -*-
'''This module provides vsim_asc class

This module might not be a member of vaspy. But the input file .*ascii
is generated usually by phonopy, and the input files needed to phonopy
calculation is provided from the VASP calculations.

Purpose of this module
This module generates

1) POSCAR files for phonon displacement
2) POVRAY scene file for animation

The first is absolutely required.
 '''

import itertools
import os.path
import bz2
import numpy as np


class VSIM_ASC(object):
    '''VSIM_ASC

    Collection of phonon mode data from v_scim ascii file


    Attributes
    -----------

    system_name: str
        System name
    ions: list
        Atoms used
    positions: list
        List of ion positionn in static
    qpts: list
        List of qvectors
    freqs: list
        List of phonon frequencies
    d_vectors: np.array
        d_vectors[mode#][atom#] returns the displacement (complex) vector
    lattice.vectors: np.array
'''
    def __init__(self, filename=None):
        self.system_name = ""
        self.ions = []
        #
        self.qpts = []
        self.freqs = []
        #
        if filename:
            if os.path.splitext(filename)[1] == '.bz2':
                try:
                    thefile = bz2.open(filename, mode='rt')
                except AttributeError:
                    thefile = bz2.BZ2File(filename, mode='r')
            else:
                thefile = open(filename)
            self.load_file(thefile)

    def load_file(self, thefile):
        '''vsim.ascii parser

        Parameters
        ----------

        thefile: StringIO
            "VSIM.ascii" file
'''
        phonon_lines = []
        # the first line is system name
        self.system_name = next(thefile)[1:].strip()
        # the 2nd line represents dxx, dyx, dyy
        dxx, dyx, dyy = [float(x) for x in next(thefile).split()]
        # the 3rd line represents dzx, dzy, dzz
        dzx, dzy, dzz = [float(x) for x in next(thefile).split()]
        self.lattice_vector = np.array([[dxx, 0, 0],
                                        [dyx, dyy, 0],
                                        [dzx, dzy, dzz]])
        self.ions = []
        self.positions = []
        self.d_vectors = []
        for line in thefile:
            line = line.strip()
            if line[0] == "#" or line[0] == "!":
                phonon_lines.append(line[1:].strip())
            else:
                x, y, z, ion = line.split()
                self.ions.append(ion)
                self.positions.append(np.array(
                    [float(x), float(y), float(z)]))
        # self.ionnums, self.iontypes = ions_to_iontypes_ionnums(self.ions)
        #
        for line in phonon_lines:
            if 'metaData' in line:
                aline = line[15:]
                if aline[-1] == '\\':
                    aline = aline[:-1]
                modedata = aline.split(';')
                qpt = np.array([float(x) for x in modedata[0:3]])
                freq = float(modedata[3])
                self.qpts.append(qpt)
                self.freqs.append(freq)
            elif ']' in line:
                pass
            else:  # displacement vector
                vectors = [float(x) for x in line[1:-1].split(';')]
                self.d_vectors.append([vectors[0] + vectors[3]*1j,
                                       vectors[1] + vectors[4]*1j,
                                       vectors[2] + vectors[5]*1j])
        n_phonons = len(self.freqs)
        self.d_vectors = np.array(self.d_vectors).reshape(n_phonons,
                                                          len(self.ions),
                                                          3)
        self.freqs = np.array(self.freqs)

    def build_animation_frame(self, mode=0, supercell=(2, 2, 1), n_frames=30):
        '''Build data for creating POSCAR etc.,

        Parameters
        ----------

        mode: int
           mode number
        supercell: tuple
           supercell dimensions
        n_frames: int
           total number of animation frmaes
'''
        # カーテシアンのk をつくる
        #
        # 以下各原子についてsupercell の位置を求めそれぞれ
        # phonon mode に対応した原子位置変化を計算する
        for atom_index, position in enumerate(self.positions):
            #  変位ベクトルをもとめる
            #  mass
            #
            for cell_id in itertools.product(range(supercell[0]),
                                             range(supercell[1]),
                                             range(supercell[2])):
                # 原子位置を求める。
                # animate_atom_phononの実行
                pass


def animate_atom_phonon(position, qpt_cart, d_vector, mass=1.0,
                        n_frames=30, s_frame=0, e_frame=None,
                        magnitude=1):
    '''Return atom position series determined by d_vector and q

    Parameters
    ------------

    position: list, tuple, np.array
       position of atom in cartesian coordinate
    qpt_cart: np.array
       wavevector in cartesian coordinate
    d_vector: np.array
       displacement (complex) vecror
    mass: float
       mass of atom
    n_frames: int
       total number of animationn frames
    s_frame: int
       start number of frame
    e_frame: int or None
       end number of frame
    magnitude: float
        Scale factor for atom moving

    Returns
    ---------

    positions: list
        List of atom position representing animation
'''

    position = np.array(position)  # for safe
    positions = []
    if not e_frame:
        e_frame = s_frame + n_frames - 1
    for frame in range(s_frame, e_frame+1):
        exponent = np.exp(1.0j * np.dot(position, qpt_cart) -
                          2 * np.pi * frame/n_frames)
        normal_displ = np.array(list(map((lambda y: (y.real)),
                                         [x * exponent for x in d_vector])))
        pos = position + magnitude * normal_displ / np.sqrt(mass)
        positions.append(pos)
    return positions


def ions_to_iontypes_ionnums(ions):
    '''Return ionnums and iontypes list

    Returns
    --------

    ionnums
        list of number of ions

    iontypes
        list of ionnames


    Example
    -----------

    >>> ions_to_iontypes_ionnums(['Si', 'Si', 'Ag', 'Ag', 'Ag', \
                                  'Ag', 'H', 'H', 'Si'])
    ([2, 4, 2, 1], ['Si', 'Ag', 'H', 'Si'])
'''
    thelast = ''
    ionnums = []
    iontypes = []
    while ions:
        ion = ions.pop(0)
        if thelast == ion:
            ionnums[-1] = ionnums[-1] + 1
        else:
            ionnums.append(1)
            iontypes.append(ion)
        thelast = ion
    return ionnums, iontypes


if __name__ == '__main__':
    import doctest
    doctest.testmod()
