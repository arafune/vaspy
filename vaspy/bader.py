# -*- coding: utf-8 -*-
"""Module for bader class.

Bader charge anaslysis performed by bader.
(http://theory.cm.utexas.edu/henkelman/code/bader/)

This program is not a part of vasp but it deeply connects with the vasp.
"""

from __future__ import division  # Version safety
from __future__ import print_function
from typing import IO, Any, List, Union
from pathlib import Path
from vaspy.tools import open_by_suffix


class BaderACF(object):
    r"""Class for storing bader charge analysis (ACF.dat).

    ACF.dat contains the coordinates of each atom, the charge
    associated with it according to Bader partitioning, percentage of
    the whole according to Bader partitioning and the minimum distance
    to the surface. This distance should be compared to maximum
    cut-off radius for the core region if pseudo potentials have been
    used.

    Attributes
    ----------
    natom: int
        Number of atoms
    positions: list
        Atom positions in :math:`\AA`
    charges: list
        charges
    vols: list
        volumes
    vaccharge: float
        vacuum charge
    vacvol: float
        vacuum volume
    nelectron: float
        number of electron

    """

    def __init__(self, filename: Union[str, Path] = "") -> None:
        """Initialize

        Parameters
        ----------
        filename : Union[str, Path], optional
            filename, by default ""
        """
        self.natom = 0
        self.positions: List[List[float]] = []
        self.charges: List[float] = []
        self.mindists: List[float] = []
        self.vols: List[float] = []
        self.vaccharge: float = 0.0
        self.vacvol: float = 0.0
        self.nelectron: float = 0.0
        if filename:
            self.parse(open_by_suffix(str(filename)))

    def parse(self, thefile: IO[Any]) -> None:
        """Parse AFC.dat.

        Parameters
        ----------
        thefile : IO[Any]
            'ACF.dat' file
        """
        # the first line is like:
        #     X     Y     Z     CHARGE      MIN DIST   ATOMIC VOL
        next(thefile)
        # the 2nd line is just "----------"
        separator = next(thefile)
        for line in thefile:
            if separator in line:
                break
            else:
                tmp = line.split()
                self.positions.append([float(pos) for pos in tmp[1:4]])
                self.charges.append(float(tmp[4]))
                self.mindists.append(float(tmp[5]))
                self.vols.append(float(tmp[6]))
        self.vaccharge = float(next(thefile).split()[-1])
        self.vacvol = float(next(thefile).split()[-1])
        self.nelectron = float(next(thefile).split()[-1])
        self.natom = len(self.positions)
        thefile.close()
