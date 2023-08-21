#! /usr/bin/env python
# -*- conding: utf-8 -*-
"""POSCAR class.

translate from poscar.rb of 2014/2/26, master branch

The below is an example of POSCAR file::

fcc (111) surface
    3.5200000
    0.707106780  0.000000000   0.000000000
    -0.353553390  0.612372400   0.000000000
    0.000000000  0.000000000  11.547000006
Cu
9
Selective dynamics
Direct
    0.000000000  0.000000000  0.000000000   F   F   F
    0.333333333  0.666666666  0.050000000   F   F   F
    0.666666666  0.333333333  0.100000000   F   F   F
    0.000000000  0.000000000  0.150000000   T   T   T
    0.333333333  0.666666666  0.200000000   T   T   T
    0.000000000  0.000000000  0.250000000   T   T   T
    0.333333333  0.666666666  0.300000000   T   T   T
    0.666666666  0.333333333  0.350000000   T   T   T
    0.000000000  0.000000000  0.400000000   T   T   T

    0.00000000E+00  0.00000000E+00  0.00000000E+00
    0.00000000E+00  0.00000000E+00  0.00000000E+00
    0.00000000E+00  0.00000000E+00  0.00000000E+00
    0.00000000E+00  0.00000000E+00  0.00000000E+00
    0.00000000E+00  0.00000000E+00  0.00000000E+00
    0.00000000E+00  0.00000000E+00  0.00000000E+00
    0.00000000E+00  0.00000000E+00  0.00000000E+00
    0.00000000E+00  0.00000000E+00  0.00000000E+00
    0.00000000E+00  0.00000000E+00  0.00000000E+00
"""

from __future__ import annotations

import copy
import itertools as it
import re
from logging import INFO, Formatter, StreamHandler, getLogger
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any

import numpy as np

from vaspy import tools
from vaspy.tools import open_by_suffix

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Sequence

    from numpy.typing import ArrayLike, NDArray

# logger
LOGLEVEL = INFO
logger = getLogger(__name__)
fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
formatter = Formatter(fmt)
handler = StreamHandler()
handler.setLevel(LOGLEVEL)
logger.setLevel(LOGLEVEL)
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False


class PosCarHead:
    """One of the parent classes of POSCAR class.

    Attributes
    ----------
    system_name: str
        system name
    scaling_factor: float
        scaling factor
    atom_types : list
        list of ion name
    atomnums : list
        list of number of atoms. Corresponding to `atom_types`

    """

    def __init__(self) -> None:
        """Initialize."""
        self.__cell_vecs: NDArray[np.float64] = np.array(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        )
        self.system_name: str = ""
        self.scaling_factor: float = 0.0
        self.atom_types: list[str] = []
        self.atomnums: list[int] = []
        self.__site_label: list[str] = []

    @property
    def cell_vecs(self) -> NDArray[np.float64]:
        """Return the matrix of the unit cell."""
        return self.__cell_vecs

    @cell_vecs.setter
    def cell_vecs(self, vec: Sequence[float]) -> None:
        """Setter of cell matrix.

        Parameters
        ----------
        vec: numpy.array or list or tuple
            3x3 matrix

        """
        if three_by_three(vec):
            self.__cell_vecs = np.array(vec)
        else:
            raise TypeError

    @property
    def realcell(self) -> NDArray[np.float64]:
        """Alias of cell_vecs to keep consistency with wavecar.py."""
        return self.__cell_vecs

    @realcell.setter
    def realcell(self, vec: Sequence[float]) -> None:
        """Alias of cell_vecs to keep consistency with wavecar.py.

        Parameters
        ----------
        vec: numpy.array or list or tuple
            3x3 matrix

        """
        if three_by_three(vec):
            self.__cell_vecs = np.array(vec)
        else:
            raise TypeError

    @property
    def site_label(self) -> list[str]:
        """Return list style of "site_label" (e.g.  "#0:Ag1")."""
        # for elm, n in zip(self.atom_types, self.atomnums):
        #     self.__site_label.extend(
        #         '#{0}:{1}{2}'.format(ii + m, elm, m + 1) for m in range(n))
        self.__site_label = []
        atomnames: list[str] = []
        for elm, atomnums in zip(self.atom_types, self.atomnums, strict=True):
            for j in range(1, atomnums + 1):
                elem_num = elm + str(j)
                if elem_num not in atomnames:
                    atomnames.append(elem_num)
                else:
                    while elem_num in atomnames:
                        j = j + 1
                        elem_num = elm + str(j)
                    else:
                        atomnames.append(elem_num)
        self.__site_label = ["#" + str(s) + ":" + a for s, a in enumerate(atomnames)]
        return self.__site_label

    @site_label.setter
    def site_label(self, value: list[str]) -> None:
        self.__site_label = value


class PosCarPos:
    """POSCAR_DOS Class.

    Attributes
    ----------
    coordinate_type: str
        "Direct" or "Cartesian"
    positions: list
        list of positions (np.array)

    """

    def __init__(self) -> None:
        """Initialize."""
        self.coordinate_type = ""
        self.positions: list[NDArray[np.float64]] = []
        self.coordinate_changeflags: list[str] = []
        self.selective: bool = False

    def is_cartesian(self) -> bool:
        """Return True if Cartesian coordinate is set.

        Returns
        -------
        boolean
            True if coordinate is cartesian

        """
        return bool(re.search(r"^[ck]", self.coordinate_type, re.I))

    def is_direct(self) -> bool:
        """Return True if DIRECT coordinate is set.

        Returns
        -------
        Boolean
            True if coordinate is direct (not cartesian)

        """
        return not self.is_cartesian()

    def __getitem__(self, item: int) -> NDArray[np.float64]:
        return self.positions[item]


class POSCAR(PosCarHead, PosCarPos):
    """Class for POSCAR (CONTCAR) format.

    This script does *NOT* support for constructing POSCAR
    from scratch. (Use ASE for this purpose.)

    It provides a way to slightly modify POSCAR or
    CONTCAR, which has already works well.

    Attributes
    ----------
    system_name, scaling_factor, cell_vecs

    """

    def __init__(self, arg: Sequence[str] | None = None) -> None:
        """Initialize.

        Parameters
        ----------
        arg: str
            POSCAR file name, or list of POSCAR text.

        """
        super().__init__()
        PosCarPos.__init__(self)
        if isinstance(arg, str | Path):
            poscar: tuple[str, ...] | list[str] = open_by_suffix(arg).readlines()
            self.load_array(poscar)
        if isinstance(arg, list | tuple):
            self.load_array(arg)

    def load_array(self, input_poscar: list[str] | tuple[str, ...]) -> None:
        """Parse POSCAR as list.

        Parameters
        ----------
        poscar: str, list, tuple
            POSCAR data

        """
        poscar = iter(map(str.rstrip, input_poscar))
        self.system_name = next(poscar)
        self.scaling_factor = float(next(poscar))
        self.cell_vecs[0] = [float(x) for x in next(poscar).split()]
        self.cell_vecs[1] = [float(x) for x in next(poscar).split()]
        self.cell_vecs[2] = [float(x) for x in next(poscar).split()]
        self.atom_types: list[str] = next(poscar).split()
        # parse POSCAR evenif the element names are not set.
        # At present, the String representation number
        #   are used for the  dummy name.
        if self.atom_types[0].isdigit():
            self.atomnums = [int(i) for i in self.atom_types]
        else:
            self.atomnums = [int(x) for x in next(poscar).split()]
        line7 = next(poscar)
        if re.search(r"^[\s]*Selective\b", line7, re.I):
            self.selective = True
            self.coordinate_type = next(poscar)
        else:
            self.selective = False
            self.coordinate_type = line7

        for line, _ in zip(poscar, self.site_label):
            tmp: list[str] = line.split()
            self.positions.append(np.asarray(np.array(tmp[:3]), dtype=np.float64))
            if self.selective:
                self.coordinate_changeflags.append(" ".join(tmp[3:]))

    def __iter__(self) -> Generator:
        yield from self.positions

    def sort(
        self,
        from_site: int = 0,
        to_site: int | None = None,
        axis: str | int = "z",
    ) -> None:
        """Sort positions attribute by coordinate.

        Parameters
        ----------
        from_site: int, default 0
            first index # for sort

        to_site: int, default None
            last index # for sort

        axis: str, default `z`
            Axis used for sort


        Notes
        -----
        The first site # is "0". It's the pythonic way.
        The element difference is **not** taken into account.

        """
        if to_site is None:
            to_site = sum(self.atomnums)
        if axis in ("x", "X", 0):
            axis = 0
        elif axis in ("y", "Y", 1):
            axis = 1
        elif axis in ("z", "Z", 2):
            axis = 2
        self.positions = (
            self.positions[:from_site]
            + sorted(
                self.positions[from_site:to_site],
                key=lambda sortaxis: sortaxis[axis],
            )
            + self.positions[to_site:]
        )

    def supercell(self, n_x: int, n_y: int, n_z: int) -> POSCAR:
        r"""Return the :math:`(n_x \\times n_y \\times n_z)` supercell.

        Parameters
        ----------
        n_x: int
            repeat number along x axis
        n_y: int
            repeat number along y axis
        n_z: int
            repeat number along z axis


        Returns
        -------
        POSCAR
            POSCAR object of the supercell

        """
        if (
            not isinstance(n_x, int)
            or not isinstance(n_y, int)
            or not isinstance(n_z, int)
        ):
            msg = "arguments must be positive integer"
            raise ValueError(msg)
        if n_x <= 0 or n_y <= 0 or n_z <= 0:
            msg = "arguments must be positive integer"
            raise ValueError(msg)
        #
        sposcar = copy.deepcopy(self)
        original_is_cartesian = sposcar.is_cartesian()
        if original_is_cartesian:
            sposcar.to_direct()
        sposcar.repack_in_cell()
        sposcar.cell_vecs[0] = sposcar.cell_vecs[0] * n_x
        sposcar.cell_vecs[1] = sposcar.cell_vecs[1] * n_y
        sposcar.cell_vecs[2] = sposcar.cell_vecs[2] * n_z
        sposcar.atomnums = [i * n_x * n_y * n_z for i in sposcar.atomnums]
        spositions: list[NDArray[np.float64]] = sposcar.positions
        sposcar.positions = []
        spositions = [
            np.array([x[0] / n_x, x[1] / n_y, x[2] / n_z]) for x in spositions
        ]
        for spos in spositions:
            for i_z in range(n_z):
                for i_y in range(n_y):
                    for i_x in range(n_x):
                        sposcar.positions.append(
                            np.array(
                                [
                                    spos[0] + i_x / n_x,
                                    spos[1] + i_y / n_y,
                                    spos[2] + i_z / n_z,
                                ],
                            ),
                        )
        sposcar.coordinate_changeflags = []
        for flags in self.coordinate_changeflags:
            for _ in range(n_x * n_y * n_z):
                sposcar.coordinate_changeflags.append(flags)
        sposcar.site_label = []
        if original_is_cartesian:
            sposcar.to_cartesian()
        return sposcar

    # class method? or independent function?
    def nearest(
        self,
        array: NDArray[np.float64],
        point: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Return the nearest position in the periodic space.

        Parameters
        ----------
        array: list of float  # << Check

        point: numpy.ndarray

        Returns
        -------
        numpy.ndarray

        """
        return min(array, key=lambda pos: np.linalg.norm(pos - point))

    # class method? or independent function?
    def make27candidate(self, position: ArrayLike) -> list[NDArray[np.float64]]:
        """Return 27 vectors set correspond the neiboring.

        Parameters
        ----------
        position: numpy.ndarray, list
            atom position defined in the coordinated by
            cell_vecs ( scaling facter is not accounted).

        Returns
        -------
        list

        """
        position = _vectorize(position)
        candidates27 = []
        if self.is_cartesian():
            for i, j, k in it.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1]):
                candidates27.append(
                    i * self.cell_vecs[0]
                    + j * self.cell_vecs[1]
                    + k * self.cell_vecs[2]
                    + position,
                )
        else:
            for i, j, k in it.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1]):
                candidates27.append(
                    i * np.array([1.0, 0.0, 0.0])
                    + j * np.array([0.0, 1.0, 0.0])
                    + k * np.array([0.0, 0.0, 1.0])
                    + position,
                )
        return candidates27

    def rotate_atom(
        self,
        site: int,
        axis_name: str,
        theta_deg: float,
        center: Sequence[float],
    ) -> None:
        """Rotate the atom.

        Parameters
        ----------
        site: int
            site # for rotation (The first atom is "0".).
        axis_name: str
            "X", "x", "Y", "y", "Z", or "z". Rotation axis.
        theta_deg: float
            Rotation angle (Degrees).
        center: numpy.ndarray, list, tuple
            center position for rotation.


        :Todo:  check the center in the Braves lattice.
                take into account the periodic boundary.

        """
        rotate_at: NDArray[np.float64] = _vectorize(center)
        if len(center) != 3:
            raise ValueError
        if not point_in_box(rotate_at / self.scaling_factor, self.cell_vecs):
            msg = "the center must be in the Braves lattice"
            raise ValueError(msg)
        if not isinstance(site, int):
            msg = "argument error in rotate_atom method"
            raise ValueError(msg)
        if not self.is_cartesian():
            self.to_cartesian()
        position: NDArray[np.float64] = self.positions[site]
        position -= rotate_at / self.scaling_factor
        rotate: dict[str, Callable[[float], NDArray[np.float64]]] = {
            "x": rotate_x,
            "y": rotate_y,
            "z": rotate_z,
        }
        position = rotate[axis_name.lower()](theta_deg).dot(position)
        position += rotate_at / self.scaling_factor
        self.positions[site] = position

    def rotate_atoms(
        self,
        site_list: Sequence[int],
        axis_name: str,
        theta_deg: float,
        center: Sequence[float],
    ) -> None:
        """Rotate atoms.

        Parameters
        ----------
        site_list:
            list array of site for rotation.
        axis_name:
            "X", "x", "Y", "y", "Z",or "z".  Rotation axis.
        theta_deg: float
            Rotation angle (Degrees).
        center: numpy.ndarray, list, tuple
            Position of rotation center

        """
        for site in site_list:
            self.rotate_atom(site, axis_name, theta_deg, center)

    def rotate_cell(self, theta_deg: float, axis_name: str = "Z") -> None:
        """Rotate unit-cell (rotation angle is set by degree).

        Parameters
        ----------
        theta_deg: float
            rotation angle
        axis_name: str
            axis name for rotation (x, y, or z)

        """
        original_is_cartesian = False
        if self.is_cartesian():
            original_is_cartesian = True
            self.to_direct()
        axis_name = axis_name.capitalize()
        rotate: dict[str, Callable[[float], NDArray[np.float64]]] = {
            "X": rotate_x,
            "Y": rotate_y,
            "Z": rotate_z,
        }
        self.cell_vecs = np.dot(rotate[axis_name](theta_deg), self.cell_vecs.T).T
        if original_is_cartesian:
            self.to_cartesian()

    def repack_in_cell(self) -> None:
        """Repack all atoms in the unit cell.

        No negative values in DIRECT coordinate.
        """
        original_is_cartesian = False
        if self.is_cartesian():
            original_is_cartesian = True
            self.to_direct()
        for pos in self.positions:
            for i in (0, 1, 2):
                while pos[i] < 0.0 or pos[i] > 1.0:
                    if pos[i] < 0.0:
                        pos[i] += 1.0
                    elif pos[i] > 1.0:
                        pos[i] -= 1.0
        if original_is_cartesian:
            self.to_cartesian()

    def __add__(self, other: POSCAR) -> POSCAR:
        """Add two poscar objects.

        Parameters
        ----------
        other: POSCAR

        Returns
        -------
        POSCAR

        Todo
        -----
        Check the lattice vectors, coordinate_type and so on.

        """
        if not isinstance(other, POSCAR):
            raise RuntimeError  # return NotImplemented
        dest_poscar = copy.deepcopy(self)
        if dest_poscar.scaling_factor != other.scaling_factor:
            msg = "scaling factor is different."
            raise ValueError(msg)
        if np.linalg.norm(dest_poscar.cell_vecs - other.cell_vecs) != 0:
            msg = "lattice vectors (cell matrix) are different."
            raise ValueError(msg)
        dest_poscar.atom_types.extend(other.atom_types)
        dest_poscar.atomnums.extend(other.atomnums)
        dest_poscar.positions.extend(other.positions)
        dest_poscar.coordinate_changeflags.extend(other.coordinate_changeflags)
        return dest_poscar

    def split(self, indexes: Sequence[int]) -> tuple[POSCAR, POSCAR]:
        """Split into two POSCAR object.

        Useful for differential charge distribution calculations.

        Parameters
        ----------
        indexes: tuple or list
            index array for a POSCAR file

        Returns
        -------
        one, other: POSCAR
            tuple of two POSCAR objects

        """
        one = copy.deepcopy(self)
        other = copy.deepcopy(self)
        one.positions = []
        other.positions = []
        one.coordinate_changeflags = []
        other.coordinate_changeflags = []
        logger.debug(f"one.positions: {one.positions}")
        atoms = tools.atom_types_atomnums_to_atoms(self.atom_types, self.atomnums)
        logger.debug(f"atoms: {atoms}")
        one_atoms: list[str] = []
        other_atoms: list[str] = []
        for i, (element, position, coordinate_flag) in enumerate(
            zip(atoms, self.positions, self.coordinate_changeflags),
        ):
            if i in indexes:
                one_atoms.append(element)
                one.positions.append(position)
                one.coordinate_changeflags.append(coordinate_flag)
            else:
                other_atoms.append(element)
                other.positions.append(position)
                other.coordinate_changeflags.append(coordinate_flag)
        logger.debug(f"one_atoms: {one_atoms}")
        one.atom_types, one.atomnums = tools.atoms_to_atom_types_atomnums(one_atoms)
        other.atom_types, other.atomnums = tools.atoms_to_atom_types_atomnums(
            other_atoms,
        )
        return one, other

    def merge(self, other: POSCAR) -> POSCAR:
        """Return POSCAR generated from two POSCARs.

        Even if the cell vectors and scaling factors are different,
        the 'merged' POSCAR is created.
        Use the cell vectors of the first POSCAR.

        Parameters
        ----------
        other: POSCAR
            POSCAR object

        Returns
        -------
        POSCAR
            added poscar object

        """
        if not isinstance(other, POSCAR):
            raise RuntimeError  ## return NotImplemented
        dest_poscar = copy.deepcopy(self)
        original_is_direct = False
        if dest_poscar.is_direct():
            original_is_direct = True
            dest_poscar.to_cartesian()
        other_poscar = copy.deepcopy(other)
        original_scaling_factor = dest_poscar.scaling_factor
        other_poscar.tune_scaling_factor(original_scaling_factor)
        other_poscar.to_cartesian()
        dest_poscar.atom_types.extend(other.atom_types)
        dest_poscar.atomnums.extend(other.atomnums)
        dest_poscar.positions.extend(other.positions)
        dest_poscar.coordinate_changeflags.extend(other.coordinate_changeflags)
        if original_is_direct:
            dest_poscar.to_direct()
        return dest_poscar

    def to_list(self) -> list[Any]:
        """Return POSCAR as list.

        Returns
        -------
        list
            a list representation of POSCAR.

        """
        out_list: list[Any] = []
        out_list.append(self.system_name)
        out_list.append(self.scaling_factor)
        out_list.append(self.cell_vecs[0])
        out_list.append(self.cell_vecs[1])
        out_list.append(self.cell_vecs[2])
        if not self.atom_types[0].isdigit():
            out_list.append(self.atom_types)
        out_list.append(self.atomnums)
        if self.selective:
            out_list.append("Selective Dynamics")
        out_list.append(self.coordinate_type)
        out_list.append(self.positions)
        out_list.append(self.coordinate_changeflags)
        out_list.append(self.site_label)
        return out_list

    def __str__(self) -> str:
        """Return as str.

        Returns
        -------
        str
            a string representation of POSCAR

        """
        tmp: list[str] = []
        tmp.append(self.system_name)
        tmp.append(str(self.scaling_factor))
        tmp.append("".join(f"   {i:20.17f}" for i in self.cell_vecs[0]))
        tmp.append("".join(f"   {i:20.17f}" for i in self.cell_vecs[1]))
        tmp.append("".join(f"   {i:20.17f}" for i in self.cell_vecs[2]))
        if not self.atom_types[0].isdigit():
            tmp.append(" " + " ".join(self.atom_types))
        tmp.append(" " + " ".join(str(i) for i in self.atomnums))
        if self.selective:
            tmp.append("Selective Dynamics")
        tmp.append(self.coordinate_type)
        for pos, t_or_f, atom in it.zip_longest(
            self.positions,
            self.coordinate_changeflags,
            self.site_label,
            fillvalue="",
        ):
            tmp.append(
                " ".join(f"  {i:20.17f}" for i in pos) + " " + t_or_f + " " + atom,
            )
        return "\n".join(tmp) + "\n"

    def str_short(self) -> str:
        """Return str object (short version).

        Returns
        -------
        str
            a string representation of POSCAR, with short (8) digit format.
            used in CHGCAR

        """
        tmp: list[str] = []
        tmp.append(self.system_name)
        tmp.append(f"  {self.scaling_factor:.14f}")
        for k in range(3):
            tmp.append("".join(f"{i:12.6f}" for i in self.cell_vecs[k]))
        if not self.atom_types[0].isdigit():
            tmp.append(" " + "".join([f"{i:>5}" for i in self.atom_types]))
        tmp.append("".join([f"{i:>6}" for i in self.atomnums]))
        tmp.append(self.coordinate_type)
        for pos in self.positions:
            tmp.append("".join(f"{i:10.6f}" for i in pos))
        return "\n".join(tmp) + "\n"

    def tune_scaling_factor(self, new_scaling_factor: float = 1.0) -> None:
        """Change scaling factor to new value.

        Parameters
        ----------
        new_scaling_factor: float

        Note
        -----
        **The Braves lattice are corrected (to be equal size)**

        Warning
        --------
        If you change the cell size, change scaling_factor attribute directly

        """
        self.cell_vecs *= self.scaling_factor / new_scaling_factor
        if self.is_cartesian():
            self.positions = [
                position * self.scaling_factor / new_scaling_factor
                for position in self.positions
            ]
        self.scaling_factor = new_scaling_factor

    def to_cartesian(self) -> None:
        """Change the coordinate to cartesian from direct."""
        if self.is_direct():
            self.coordinate_type = "Cartesian"
            mat = self.cell_vecs.transpose()
            self.positions = [mat.dot(v) for v in self.positions]

    def to_direct(self) -> None:
        """Change the coordinate to direct from cartesian."""
        if self.is_cartesian():
            self.coordinate_type = "Direct"
            mat = np.linalg.inv(np.transpose(self.cell_vecs))
            self.positions = [mat.dot(v) for v in self.positions]

    def guess_molecule(
        self,
        site_list: Sequence[int],
        center: Sequence[float] | None = None,
    ) -> None:
        """Arrange atom position to form a molecule.

        This method is effective to rotate a molecule.

        Parameters
        ----------
        site_list: list
            list of site number (the number begins with #0)
        center: list
            center position of "molecule" (Optional).

        Returns
        -------
        NDArray
            Array of Vector that represents "molecule".

        Note
        -----
        When the optional argument, center, is set, the atoms are
        re-arranged as to minimize the distance from this center.
        If not, atoms are re-arranged to minimize the total bonding
        length.  As the algorithm for mimizing the total length is
        not exhaustive, the resultant atom arrangement may different
        from what you expect, in spite of time-waste.  The center
        option is highly recommended to form a molecule.

        """
        molecule: list[NDArray[np.float64]] = [self.positions[j] for j in site_list]
        newposes = []
        for index, site in enumerate(site_list):
            target_atom = self.positions[site]
            atoms27 = self.make27candidate(target_atom)

            def func(pos: ArrayLike, center: Sequence[float] | None) -> float:
                molecule[index] = _vectorize(pos)
                if center is not None:  # bool([np.ndarray]) => Error
                    center_pos: NDArray[np.float64] = _vectorize(center)
                    return np.linalg.norm(_vectorize(pos) - center_pos)
                # fixme!! when the highest symmetry point
                # can be determined from the position list,
                # guess_molecule method does not require
                # the "center" option.
                # (molecule.
                #     product(molecule)).inject(0.0) do | s, vectors |
                s: float = 0.0
                for vectors in it.product(molecule, molecule):
                    s += np.linalg.norm(vectors[0] - vectors[1])
                return s

            newpos = min(atoms27, key=(lambda x: func(x, center)))
            newposes.append(newpos)
        for site, pos in zip(site_list, newposes):
            self.positions[site] = pos

    def translate(
        self,
        vector: Sequence[float],
        atomlist: Sequence[int],
    ) -> list[NDArray[np.float64]]:
        """Translate the selected atom(s) by vector.

        Parameters
        ----------
        vector: list, tuple, numpy.array
            translational vector (in Cartesian frame)
        atomlist: list
            list of the atom to be moved


        Note
        ------
        the first atom is "0", to follow the pythonic way.

        Returns
        -------
        NDArray
                position

        """
        if self.is_cartesian():
            vector_array: NDArray[np.float64] = _vectorize(vector)
            for i in atomlist:
                self.positions[i] = (
                    self.positions[i] + vector_array / self.scaling_factor
                )
        else:
            vector_array = _vectorize(vector)
            self.to_cartesian()
            for i in atomlist:
                self.positions[i] = (
                    self.positions[i] + vector_array / self.scaling_factor
                )
            self.to_direct()
        return self.positions

    @property
    def axes_lengths(self) -> tuple[float, float, float]:
        """Return cell axis lengths.

        Returns
        -------
        tuple
            cell axis length of x, y, and z

        """
        cell_x: float = np.linalg.norm(self.cell_vecs[0] * self.scaling_factor)
        cell_y: float = np.linalg.norm(self.cell_vecs[1] * self.scaling_factor)
        cell_z: float = np.linalg.norm(self.cell_vecs[2] * self.scaling_factor)
        return (cell_x, cell_y, cell_z)

    def translate_all(self, vector: Sequence[float]) -> None:
        """Translate **all** atoms by vector.

        Parameters
        ----------
        vector: list, numpy.array
            translational vector

        """
        atomrange = list(range(sum(self.atomnums)))
        self.translate(vector, atomrange)

    def save(self, filename: str) -> None:
        """Save POSCAR contents to the file named "filename".

        Parameters
        ----------
        filename: str
            File name for save

        """
        file: IO
        try:  # Version safety
            file = open(filename, mode="w", newline="\n")
        except TypeError:
            file = open(filename, mode="wb")
        with file:
            file.write(str(self))


def point_in_box(point: ArrayLike, cell_vecs: ArrayLike) -> bool:
    """Return True if point is located in the box.

    Parameters
    ----------
    point: numpy.ndarray, numpy.matrix, list, tuple
        vector representing the "point"
    cell_vecs: numpy.ndarray, numpy.matrix, list, tuple
        vectors defining the "box"

    Returns
    -------
    bool
        True if the point is interior in the "box" defined the cell_vecs.

    """
    if three_by_three(cell_vecs):
        thepoint: NDArray[np.float64] = np.array(point).flatten()
        cell_vectors: NDArray[np.float64] = np.array(cell_vecs)
        result: NDArray[np.float64] = np.dot(np.linalg.inv(cell_vectors.T), thepoint)
        return all((0 <= float(q) <= 1) for q in result)
    else:
        raise TypeError


def rotate_x(theta_deg: float) -> NDArray[np.float64]:
    """Rotation matrix around X-axis.

    Parameters
    ----------
    theta_deg: float
        Rotation angle (Degrees)

    Returns
    -------
    NDArray
        rotation matrix (Around X)

    Example
    ---------
    >>> rotate_x(60)
    array([[ 1.       ,  0.       ,  0.       ],
            [ 0.       ,  0.5      , -0.8660254],
            [ 0.       ,  0.8660254,  0.5      ]])

    """
    degree: float = np.pi / 180.0
    return np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, np.cos(theta_deg * degree), -np.sin(theta_deg * degree)],
            [0.0, np.sin(theta_deg * degree), np.cos(theta_deg * degree)],
        ],
    )


def rotate_y(theta_deg: float) -> NDArray[np.float64]:
    """Rotation matrix around Y-axis.

    Parameters
    ----------
    theta_deg: float
        Rotation angle

    Returns
    -------
    np.ndarray
        Rotation matrix (around Y)

    Example
    --------
    >>> rotate_y(60)
    array([[ 0.5      ,  0.       ,  0.8660254],
        [ 0.       ,  1.       ,  0.       ],
        [-0.8660254,  0.       ,  0.5      ]])

    """
    degree: float = np.pi / 180.0
    return np.array(
        [
            [np.cos(theta_deg * degree), 0.0, np.sin(theta_deg * degree)],
            [0.0, 1.0, 0.0],
            [-np.sin(theta_deg * degree), 0.0, np.cos(theta_deg * degree)],
        ],
    )


def rotate_z(theta_deg: float) -> NDArray[np.float64]:
    """Rotation matrix around Z-axis.

    Parameters
    ----------
    theta_deg: float
        Rotation angle

    Returns
    -------
    np.ndarray
        Rotation matrix (around Z)

    Example
    --------
    >>> rotate_z(60)
    array([[ 0.5      , -0.8660254,  0.       ],
        [ 0.8660254,  0.5      ,  0.       ],
        [ 0.       ,  0.       ,  1.       ]])

    """
    degree: float = np.pi / 180.0
    return np.array(
        [
            [np.cos(theta_deg * degree), -np.sin(theta_deg * degree), 0.0],
            [np.sin(theta_deg * degree), np.cos(theta_deg * degree), 0.0],
            [0.0, 0.0, 1.0],
        ],
    )


def three_by_three(vec: ArrayLike) -> bool:
    """Return True if vec can be converted into the 3x3 matrix.

    Parameters
    ----------
    vec: list-like
        Vector like object to check whether can be written by the 3x3 matrix

    Returns
    -------
    bool
        True if the vec can be written by 3x3 matrix

    """
    if not isinstance(vec, np.ndarray | np.matrix | list | tuple):
        return False
    if len(vec) != 3:
        return False
    return [3, 3, 3] == [len(_) for _ in vec]


def _vectorize(vector: ArrayLike) -> NDArray[np.float64]:
    """Return np.ndarray object.

    Parameters
    ----------
    vector: array-like
        array-like object to be wanted as the vector.

    Returns
    -------
    NDArray
        Numpy array (flatten)
    """
    if not isinstance(vector, np.ndarray | np.matrix | list | tuple):
        msg = "Cannot convert into vector."
        raise TypeError(msg)
    return np.array(vector).flatten()


# --------------------------

if __name__ == "__main__":
    import doctest

    doctest.testmod()
