# -*- codinng: utf-8 -*-
"""This module provides vsim_asc class.

This module might not be a member of vaspy. But the input file .*ascii
is generated usually by phonopy, and the input files needed to phonopy
calculation is provided from the VASP calculations.

Purpose of this module
This module generates

1) POSCAR files for phonon displacement
2) POVRAY scene file for animation

The first is absolutely required.
 """
import bz2
import itertools
import logging
import os.path
from logging import Formatter, StreamHandler, getLogger

import numpy as np

import vaspy.const as const

logger = getLogger("LogTest")
logger.setLevel(logging.DEBUG)
stream_handler = StreamHandler()
handler_format = Formatter(' %(asctime)s - %(levelname)s - %(message)s ')
stream_handler.setFormatter(handler_format)


class VSIM_ASC(object):
    """Class for VSIM_ASC.

    Collection of phonon mode data from v_scim ascii file

    Attributes
    -----------
    system_name: str
        System name
    ions: list
        Atoms used
    positions: list
        List of ion position in static
    qpts: list
        List of qvectors
    freqs: list
        List of phonon frequencies
    d_vectors: np.array
        d_vectors[mode#][atom#] returns the displacement (complex) vector
    lattice.vectors: np.array

    """

    def __init__(self, filename=None):
        """Initialize."""
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
        """Parse vsim.ascii.

        Parameters
        ----------
        thefile: StringIO
            "VSIM.ascii" file

        """
        phonon_lines = []
        # the first line is system name
        self.system_name = next(thefile)[1:].strip()
        # the 2nd line represents dxx, dyx, dyy
        dxx, dyx, dyy = [float(x) for x in next(thefile).split()]
        # the 3rd line represents dzx, dzy, dzz
        dzx, dzy, dzz = [float(x) for x in next(thefile).split()]
        self.lattice_vectors = np.array([[dxx, 0, 0], [dyx, dyy, 0],
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
                self.positions.append(np.array([float(x), float(y), float(z)]))
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
                self.d_vectors.append([
                    vectors[0] + vectors[3] * 1j, vectors[1] + vectors[4] * 1j,
                    vectors[2] + vectors[5] * 1j
                ])
        n_phonons = len(self.freqs)
        self.d_vectors = np.array(self.d_vectors).reshape(
            n_phonons, len(self.ions), 3)
        self.freqs = np.array(self.freqs)

    def build_phono_motion(self,
                           mode=0,
                           supercell=(2, 2, 1),
                           n_frames=30,
                           magnitude=1):
        """Build data for creating POSCAR etc.

        Parameters
        ----------
        mode: int
           mode number
        supercell: tuple
           supercell dimensions
        n_frames: int
           total number of animation frmaes

        """
        qpt = self.qpts[mode]
        bmatrix = 2 * np.pi * np.linalg.inv(self.lattice_vectors).transpose()
        qpt_cart = qpt.dot(bmatrix)
        logger.debug(
            'qpt_cart[x] = {}, qpt_cart[y] = {}, qpt_cart[z] ={}'.format(
                qpt_cart[0], qpt_cart[1], qpt_cart[2]))
        #
        animation_positions = []
        for atom_i, position in enumerate(self.positions):

            for cell_id in itertools.product(
                    range(supercell[0]), range(supercell[1]),
                    range(supercell[2])):
                logger.debug(' cell_id:{}'.format(cell_id))
                abs_pos = position + (self.lattice_vectors[0] * cell_id[0]
                                      + self.lattice_vectors[1] * cell_id[1]
                                      + self.lattice_vectors[2] * cell_id[2])
                positions = animate_atom_phonon(
                    abs_pos,
                    qpt_cart,
                    self.d_vectors[mode][atom_i],
                    n_frames=n_frames,
                    magnitude=magnitude)
                animation_positions.append(positions)
        return animation_positions


def supercell_lattice_vectors(lattice_vectors, cell_id):
    """Return lattice vectors of supercell.

    Parameters
    ----------
    lattice_vectors: np.array
        3x3 matrix for original lattice vectors
    cell_id: tuple, list

    Returns
    -------
    np.array

    """
    supercell_vectors = []
    for x, x_i in zip(lattice_vectors, cell_id):
        supercell_vectors.append(x * x_i)
    return np.array(supercell_vectors)


def animate_atom_phonon(position,
                        qpt_cart,
                        d_vector,
                        n_frames=30,
                        s_frame=0,
                        e_frame=None,
                        magnitude=1.0):
    """Return atom position series determined by d_vector and q.

    Parameters
    ------------
    position: list, tuple, np.array
       position of atom in cartesian coordinate
    qpt_cart: np.array
       wavevector in cartesian coordinate
    d_vector: np.array
       displacement (complex) vecror
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

    """
    position0 = np.array(position)  # for safe
    positions = []
    if not e_frame:
        e_frame = s_frame + n_frames - 1
    for frame in range(s_frame, e_frame + 1):
        exponent = np.exp(
            1.0j * (np.dot(position0, qpt_cart)
                    - 2 * np.pi * frame / n_frames))
        logger.debug('r:{}, qpt_cart;{}, frame:{}, n_frames:{}'.format(
            position0, qpt_cart, frame, n_frames))
        logger.debug('arg_exponent:{}'.format(
            1.0j * (np.dot(position0, qpt_cart)
                    - 2 * np.pi * frame / n_frames)))
        logger.debug('exponent:{}'.format(exponent))
        normal_displ = np.array(
            list(map((lambda y: (y.real)), [x * exponent for x in d_vector])))
        logger.debug('normal_displ:{}'.format(normal_displ))
        # The displacement vector calculated by (at least) phonopy is
        # taken into account the mass of the atom.
        # If the calculated displacement vector
        # does not contain the mass effect,
        # the normal_displ should be devided by sqrt(mass)
        positions.append(position0 + magnitude * normal_displ)
        logger.debug('position.after_move:{}'.format(positions[-1]))
    return positions


def ions_to_iontypes_ionnums(ions):
    r"""Return ionnums and iontypes list.

    Returns
    --------
    ionnums
        list of number of ions
    iontypes
        list of ionnames


    Examples
    --------
    >>> ions_to_iontypes_ionnums(['Si', 'Si', 'Ag', 'Ag', 'Ag', \
                                  'Ag', 'H', 'H', 'Si'])
    ([2, 4, 2, 1], ['Si', 'Ag', 'H', 'Si'])

    """
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
