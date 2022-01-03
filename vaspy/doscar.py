#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module provides DOSCAR and related classes.

From VASP webpage::

| The file DOSCAR contains the DOS and integrated DOS The units are
| "number of states/unit cell". For dynamic simulations and relaxations,
| an averaged DOS and an averaged integrated DOS is written to the file.
| For a description of how the averaging is done see 7.18, 7.32).
| The first few lines of the DOSCAR file are made up by a header,
| which is followed by NDOS lines holding three data
|
| energy     dos     integrated_dos
|
| For spin-polarized calculations each line holds five data
|
| energy     dos(up) dos(dwn)  integrated_dos(up) integrated_dos(dwn)

|  If RWIGS (Wigner Seitz radii, see section 7.29) is set in the INCAR file,
|  a l- and site-projected DOS is calculated and also written to the file
|  DOSCAR. One set of data is written for each ion,
|  each set of data holds NDOS lines with the following data
|
|   energy s-dos p-dos d-dos
|
|  and
|
|   energy s-dos(up) p-dos(up) d-dos(up) s-dos(dwn) p-dos(dwn) d-dos(dwn)
|
| for the non spin-polarized and spin polarized case respectively.
|  The units of the l- and site projected DOS are states/atom.
|  Please mind, that the site projected DOS is not evaluated
|  in the parallel version if NPAR is not equal to 1.

|  Mind: For relaxations the DOSCAR is usually useless.
|  If you want to get an accurate DOS for the final configuration
|  copy CONTCAR to POSCAR and make another static (ISTART=1; NSW=0)
|  calculation.
"""

from __future__ import annotations

import copy
import sys
from operator import add
import csv
from collections.abc import Iterable
from typing import Sequence, IO
from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.stderr.write("Install matplotlib, or you cannot use methods relating to draw\n")
from vaspy.tools import open_by_suffix


class DOSCAR(object):  # Version safety
    """Class for DOSCAR file.

    A container of DOS object

    Attributes
    ----------
    natom: int
        number of atoms

    tdos:

    pdoses:

    energies: tuple[float]
        Energy
    """

    def __init__(self, filename: str|Path = "") -> None:
        """Initialize

        Parameters
        ----------
        filename : str|Path, optional
            file name of "DOSCAR"
        """
        self.natom: int = 0
        self.tdos: TDOS|None = None
        self.pdoses: list[PDOS] = []
        self.energies: tuple[float, ...] = tuple()
        if filename:
            self.load_file(open_by_suffix(str(filename)))

    def load_file(self, thefile: IO[str]) -> None:
        """Parse DOSCAR file and store it in memory.

        Parameters
        ------------
        thefile: StringIO
            "DOSCAR" file

        """
        firstline = thefile.readline()
        self.natom = int(firstline[0:4])
        [thefile.readline() for _ in range(4)]
        header = thefile.readline()

        nedos: int = int(header[32:37])
        tmp: list[list[float]] = [
            [float(i) for i in next(thefile).split()] for _ in range(nedos)
        ]
        tmp_dos: list[tuple[float]] = [*zip(*tmp)]
        #
        self.energies = tmp_dos[0]
        if len(tmp_dos) == 3:
            self.tdos = TDOS((tmp_dos[1],))
        elif len(tmp_dos) == 5:
            self.tdos = TDOS((tmp_dos[1], tmp_dos[2]))
        else:
            raise RuntimeError("Runtime Error: Check DOS file")
        #

        try:
            line = next(thefile)
        except StopIteration:
            line = ""
        while line == header:
            tmp = [[float(i) for i in next(thefile).split()] for _ in range(nedos)]
            self.pdoses.append(PDOS([*zip(*tmp)][1:]))
            try:
                line = next(thefile)
            except StopIteration:
                line = ""

        thefile.close()

    def fermi_correction(self, fermi: float) -> None:
        """Correct energy by Fermi level.

        Parameters
        ----------
        fermi: float
            fermi level

        """
        self.energies = tuple([energy - fermi for energy in self.energies])


class DOS(Sequence):  # Version safety
    """Class for DOS.

    list object consisting two elements.
    The first element is the the energy.
    The latter element is list for the density.


    Attributes
    -----------
    dos: list
        the dos data.
        By default, the first column is the energy, the latter is the density.

    """

    def __init__(self, array: tuple[float]|None = None) -> None:
        """Initialize."""
        self.dos: list[tuple[float, ...]|list[float]] = []
        self.header: str|list[str] = ""
        if array:
            self.dos = [*zip(*array)]

    def __len__(self) -> int:
        """x.__len__() <=> len(x)."""
        return len(self.dos)

    def __getitem__(
        self, idx: int|slice
    ) -> tuple[float, ...]|list[tuple[float]]:
        return self.dos[idx]

    @property
    def T(self):
        return [*zip(*self.dos)]

    def _export_csv(
        self,
        filename: str,
        energy: list[float]|tuple[float, ...],
        header: list[str],
    ) -> None:
        """Export data to csv file"""
        assert len(energy) == len(self)
        assert len(header) == len(self[0]) + 1
        with open(filename, mode="w") as csvfile:
            writer = csv.writer(csvfile, delimiter="\t")
            writer.writerow(header)
            for e, d in zip(energy, self.dos):
                writer.writerow([e] + list(d))


class TDOS(DOS):
    """Class for total DOS.

    Attributes
    ----------
    header : list[str]

    dos:
    """

    def __init__(self, array: tuple[float]|None) -> None:
        """Initialize

        Parameters
        ----------
        array : tuple[float]|None
            DOS data
        """
        super().__init__(array)
        if len(self.dos[0]) == 1:
            self.header = ["Energy", "TDOS"]
        elif len(self.dos[0]) == 2:  # collinear spin
            self.header = ["Energy", "TDOS_up", "TDOS_down"]
            self.dos = [(d[0], -d[1]) for d in self.dos]

    def graphview(self) -> None:
        """Show graphview by matplotlib."""
        for density in self.dos[1:]:
            plt.plot(self.dos[0], density)
        plt.show()

    def export_csv(
        self,
        filename: str,
        energy: list[float]|tuple[float, ...],
    ) -> None:
        return super()._export_csv(filename, header=self.header, energy=energy)


class PDOS(DOS):
    """Class for partial DOS.

    Attributes
    ----------
    site: str
        name of the site.

    orbital_spin:
        name of the orbital with spin character.
        (If non-spin calculated results, :py:attr:`orbital_spin` is
        just orbital name)

    Parameters
    -----------
    array: list[tuple[float, ...]]
        DOS data
    site: str
        site name

    """

    orbitalnames = [
        "s",
        "py",
        "pz",
        "px",
        "dxy",
        "dyz",
        "dz2",
        "dxz",
        "dx2",
        "f-3",
        "f-2",
        "f-1",
        "f0",
        "f1",
        "f2",
        "f3",
    ]  # << see sphpro.F of vasp sourcev
    spins_soi = ("mT", "mX", "mY", "mZ")
    spins = ("up", "down")

    def __init__(
        self, array: tuple[float]|None = None, site: str = ""
    ) -> None:
        """Initialize."""
        super().__init__(array)
        self.site = site
        self.orbital_spin: list[str] = list()
        self.total: list[tuple[float, ...]] = []
        if array is not None:
            if len(self.dos[0]) == 9 or len(self.dos[0]) == 16:
                self.orbital_spin = self.orbitalnames
                for dos_at_energy in self.dos:
                    self.total.append((sum(dos_at_energy),))
            elif len(self.dos[0]) == 18 or len(self.dos[0]) == 32:  # Spin resolved
                self.orbital_spin = [
                    orb + "_" + spin for orb in self.orbitalnames for spin in self.spins
                ]
                # In collinear spin calculation, DOS of down-spin is
                # set by negative value.
                for i, dos_at_energy in enumerate(self.dos):
                    self.total.append(
                        (sum(dos_at_energy[0::2]), -sum(dos_at_energy[1::2]))
                    )
                    tmp = []
                    for j, energy in enumerate(dos_at_energy):

                        if j % 2 == 0:
                            tmp.append(energy)
                        else:
                            tmp.append(-energy)
                    self.dos[i] = tmp
            elif len(self.dos[0]) == 36:  # SOI
                self.orbital_spin = [
                    orb + "_" + spin
                    for orb in self.orbitalnames
                    for spin in self.spins_soi
                ]
                for dos_at_energy in self.dos:
                    self.total.append((sum(dos_at_energy),))
            else:
                print(len(self.dos[0]))
                raise RuntimeError("Check the DOSCAR file")

    def projected(self, orbital: str|int) -> tuple[float]|list[float]:
        if isinstance(orbital, int):
            idx = orbital
        else:
            idx = self.orbital_spin.index(orbital)
        return [d[idx] for d in self.dos]

    def graphview(self, *orbitalnames: str) -> None:
        """Show DOS graph by using matplotlib.  For 'just seeing' use."""
        try:
            alist: list[int] = [
                self.orbital_spin.index(orbname) for orbname in orbitalnames
            ]
        except ValueError:
            err = "Check argment of this function\n"
            err += "The following name(s) are accepted:\n"
            err += ", ".join(self.orbital_spin)
            raise ValueError(err)
        for orbital in alist:
            plt.plot(self.dos[0], self.dos[orbital + 1])
        plt.show()

    def export_csv(
        self, filename: str, energy: list[float]|tuple[float, ...]
    ) -> None:
        """Export data to file object (or file-like object) as csv format.

        Parameters
        ----------
        filename : str
            filename for output
        energy : list[float] |  tuple[float, ...]
            Energy data
        """
        header = ["Energy"]
        for i in self.orbital_spin:
            if self.site:
                header.append(self.site + "_" + i)
            else:
                header.append(i)
        super()._export_csv(filename, energy=energy, header=header)

    def plot_dos(
        self, orbitals: Sequence[str], fermi: float = 0.0
    ) -> None:  # Not implemented yet
        """Plot DOS spectra by matplotlib.pyplot.

        Parameters
        ----------
        orbitals: str
            orbital name

        fermi: float, optional (default is 0.0)

        Warning
        --------
        not implemented yet!!

        """
        pass

    def __add__(self, other: PDOS) -> PDOS:
        """Add two DOS objects.

        x.__add__(y) <-> x+y

        Parameters
        -----------
        other: PDOS
            len(other.energies) must be equal to len(self.energies).

        Returns
        -------
        PDOS

        """
        if not isinstance(other, PDOS):
            return NotImplemented
        if len(self.dos) == 0 and self.site == "":
            return copy.deepcopy(other)
        else:
            sum_pdos = PDOS()
            sum_pdos.site = self.site + other.site
            sum_pdos.orbital_spin = self.orbital_spin
            sum_pdos.dos = [list(map(add, x, y)) for x, y, in zip(self.dos, other.dos)]
            return sum_pdos
