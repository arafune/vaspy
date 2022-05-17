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
from __future__ import annotations

import itertools
import logging
from logging import Formatter, StreamHandler, getLogger
from pathlib import Path
from typing import IO, Sequence

import numpy as np
from numpy.typing import NDArray

# import vaspy.const as const
from vaspy.tools import open_by_suffix

logger = getLogger("LogTest")
logger.setLevel(logging.DEBUG)
stream_handler = StreamHandler()
handler_format = Formatter(" %(asctime)s - %(levelname)s - %(message)s ")
stream_handler.setFormatter(handler_format)


class VSIM_ASC(object):
    """Class for VSIM_ASC.

    Collection of phonon mode data from v_scim ascii file

    Attributes
    -----------
    system_name: str
        System name
    atoms: list
        Atoms used
    positions: list
        list of atom position in static
    qpts: list
        list of qvectors
    freqs: list
        list of phonon frequencies
    d_vectors: np.array
        d_vectors[mode#][atom#] returns the displacement (complex) vector
    lattice.vectors: np.array

    """

    def __init__(self, filename: str | Path = "") -> None:
        """Initialize."""
        self.system_name: str = ""
        self.atoms: list[str] = []
        #
        self.qpts: list[NDArray[np.float64]] = []
        self.freqs: NDArray[np.float64]
        #
        if filename:
            self.load_file(open_by_suffix(str(filename)))

    def load_file(self, thefile: IO[str]) -> None:
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
        self.lattice_vectors: NDArray[np.float64] = np.array(
            [[dxx, 0, 0], [dyx, dyy, 0], [dzx, dzy, dzz]]
        )
        self.atoms = []
        self.positions = []
        d_vectors: list[list[complex]] = []
        freqs: list[float] = []
        for line in thefile:
            line = line.strip()
            if line[0] == "#" or line[0] == "!":
                phonon_lines.append(line[1:].strip())
            else:
                x, y, z, atom = line.split()
                self.atoms.append(atom)
                self.positions.append(np.array([float(x), float(y), float(z)]))
        # self.ionnums, self.iontypes = ions_to_iontypes_ionnums(self.ions)
        #
        for line in phonon_lines:
            if "metaData" in line:
                aline = line[15:]
                if aline[-1] == "\\":
                    aline = aline[:-1]
                modedata = aline.split(";")
                qpt: NDArray[np.float64] = np.array([float(x) for x in modedata[0:3]])
                freq = float(modedata[3])
                self.qpts.append(qpt)
                freqs.append(freq)
            elif "]" in line:
                pass
            else:  # displacement vector
                vectors = [float(x) for x in line[1:-1].split(";")]
                d_vectors.append(
                    [
                        vectors[0] + vectors[3] * 1j,
                        vectors[1] + vectors[4] * 1j,
                        vectors[2] + vectors[5] * 1j,
                    ]
                )
        n_phonons = len(freqs)
        self.d_vectors: NDArray[np.complex128] = np.array(d_vectors).reshape(
            n_phonons, len(self.atoms), 3
        )
        self.freqs = np.array(freqs)
        thefile.close()

    def build_phono_motion(
        self,
        mode: int = 0,
        supercell: tuple[int, int, int] = (2, 2, 1),
        n_frames: int = 30,
        magnitude: float = 1,
    ) -> list[list[NDArray[np.float64]]]:
        """Build data for creating POSCAR etc.

        Parameters
        ----------
        mode: int
            mode number
        supercell: tuple
            supercell dimensions
        n_frames: int
            total number of animation frames

        """
        qpt = self.qpts[mode]
        bmatrix: NDArray[np.float64] = (
            2 * np.pi * np.linalg.inv(self.lattice_vectors).transpose()
        )
        qpt_cart: NDArray[np.float64] = qpt.dot(bmatrix)
        logger.debug(
            "qpt_cart[x] = {}, qpt_cart[y] = {}, qpt_cart[z] ={}".format(
                qpt_cart[0], qpt_cart[1], qpt_cart[2]
            )
        )
        #
        animation_positions: list[list[NDArray[np.float64]]] = []
        for atom_i, position in enumerate(self.positions):
            for cell_id in itertools.product(
                range(supercell[0]), range(supercell[1]), range(supercell[2])
            ):
                logger.debug(" cell_id:{}".format(cell_id))
                abs_pos = position + (
                    self.lattice_vectors[0] * cell_id[0]
                    + self.lattice_vectors[1] * cell_id[1]
                    + self.lattice_vectors[2] * cell_id[2]
                )
                positions = animate_atom_phonon(
                    abs_pos,
                    qpt_cart,
                    self.d_vectors[mode][atom_i],
                    n_frames=n_frames,
                    magnitude=magnitude,
                )
                animation_positions.append(positions)
        return animation_positions


def supercell_lattice_vectors(
    lattice_vectors: NDArray[np.float64], cell_id: Sequence[int]
) -> NDArray[np.float64]:
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
    supercell_vectors: list[float] = []
    for x, x_i in zip(lattice_vectors, cell_id):
        supercell_vectors.append(x * x_i)
    return np.array(supercell_vectors)


def animate_atom_phonon(
    position: Sequence[float],
    qpt_cart: NDArray[np.float64],
    d_vector: NDArray[np.float64],
    n_frames: int = 30,
    s_frame: int = 0,
    e_frame: int | None = None,
    magnitude: float = 1.0,
) -> list[NDArray[np.float64]]:
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
        list of atom position representing animation

    """
    position0: NDArray[np.float64] = np.array(position)  # for safe
    positions: list[NDArray[np.float64]] = []
    if not e_frame:
        e_frame = s_frame + n_frames - 1
    for frame in range(s_frame, e_frame + 1):
        exponent = np.exp(
            1.0j * (np.dot(position0, qpt_cart) - 2 * np.pi * frame / n_frames)
        )
        logger.debug(
            "r:{}, qpt_cart;{}, frame:{}, n_frames:{}".format(
                position0, qpt_cart, frame, n_frames
            )
        )
        logger.debug(
            "arg_exponent:{}".format(
                1.0j * (np.dot(position0, qpt_cart) - 2 * np.pi * frame / n_frames)
            )
        )
        logger.debug("exponent:{}".format(exponent))
        normal_displ: NDArray[np.float64] = np.array(
            list(map((lambda y: (y.real)), [x * exponent for x in d_vector]))
        )
        logger.debug("normal_displ:{}".format(normal_displ))
        # The displacement vector calculated by (at least) phonopy is
        # taken into account the mass of the atom.
        # If the calculated displacement vector
        # does not contain the mass effect,
        # the normal_displ should be devided by sqrt(mass)
        positions.append(position0 + magnitude * normal_displ)
        logger.debug("position.after_move:{}".format(positions[-1]))
    return positions


if __name__ == "__main__":
    import doctest

    doctest.testmod()
