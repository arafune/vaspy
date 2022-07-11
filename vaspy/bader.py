# -*- coding: utf-8 -*-
"""Module for bader class.

Bader charge anaslysis performed by bader.
(http://theory.cm.utexas.edu/henkelman/code/bader/)

This program is not a part of vasp but it deeply connects with the vasp.
"""

from __future__ import annotations, division  # Version safety

from pathlib import Path
from typing import IO, Any

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
    n_atom: int
        Number of atoms
    positions: list
        Atom positions in :math:`\AA`
    charges: list
        charges
    vols: list
        volumes
    vac_charge: float
        vacuum charge
    vac_vol: float
        vacuum volume
    n_electron: float
        number of electron

    """

    def __init__(self, filename: str | Path = "") -> None:
        """Initialize

        Parameters
        ----------
        filename : str|Path
            filename, by default ""
        """
        self.n_atom = 0
        self.positions: list[list[float]] = []
        self.charges: list[float] = []
        self.mindists: list[float] = []
        self.vols: list[float] = []
        self.vac_charge: float = 0.0
        self.vac_vol: float = 0.0
        self.n_electron: float = 0.0
        if filename:
            self.parse(open_by_suffix(str(filename)))

    def parse(self, the_file: IO[Any]) -> None:
        """Parse AFC.dat.

        Parameters
        ----------
        the_file : IO[Any]
            'ACF.dat' file
        """
        # the first line is like:
        #     X     Y     Z     CHARGE      MIN DIST   ATOMIC VOL
        next(the_file)
        # the 2nd line is just "----------"
        separator = next(the_file)
        for line in the_file:
            if separator in line:
                break
            else:
                tmp = line.split()
                self.positions.append([float(pos) for pos in tmp[1:4]])
                self.charges.append(float(tmp[4]))
                self.mindists.append(float(tmp[5]))
                self.vols.append(float(tmp[6]))
        self.vac_charge = float(next(the_file).split()[-1])
        self.vac_vol = float(next(the_file).split()[-1])
        self.n_electron = float(next(the_file).split()[-1])
        self.n_atom = len(self.positions)
        the_file.close()
