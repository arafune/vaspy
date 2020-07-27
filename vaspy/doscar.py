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
import numpy as np
from typing import Iterator, Sequence, Union, Optional, IO, Tuple, List, Any


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

    nbands: int
        number of bands

        Notes
        -----
        VASP defaults is 301

    dos_container: list
        stores the DOS object. By default, the first item is the TDOS, and the
        latter items are PDOS.

    """

    def __init__(self, filename: Optional[str] = None) -> None:
        """Initialize."""
        self.natom: int = 0
        self.nbands: int = 0
        self.dos_container: List = list()

        if filename:
            self.load_file(open_by_suffix(filename))

    def load_file(self, thefile: Union[IO[str]]) -> None:
        """Parse DOSCAR file and store it in memory.

        Parameters
        ------------
        thefile: StringIO
            "DOSCAR" file

        """
        firstline = thefile.readline()
        self.natom = int(firstline[0:4])
        [thefile.readline() for i in range(4)]
        header = thefile.readline()
        self.nbands = int(header[32:37])
        tdos: np.ndarray = np.array(
            [next(thefile).split() for i in range(self.nbands)], dtype=np.float64
        )
        if tdos.shape[1] == 3:
            tdos = tdos[:, 0:2]
        elif tdos.shape[1] == 5:
            tdos = tdos[:, 0:3]
        else:
            raise RuntimeError
        self.dos_container = [tdos]
        try:
            line = next(thefile)
        except StopIteration:
            line = ""
        while line == header:
            self.dos_container.append(
                np.array(
                    [next(thefile).split() for i in range(self.nbands)],
                    dtype=np.float64,
                )
            )
            try:
                line = next(thefile)
            except StopIteration:
                line = ""
        thefile.close()


class DOS(object):  # Version safety
    """Class for DOS.

    List object consisting two elements.
    The first element is the the energy.
    The latter element is list for the density.


    Attributes
    -----------
    dos: list
        the dos data.
        By default, the first column is the energy, the latter is the density.

    """

    def __init__(self, array=None) -> None:
        """Initialize."""
        self.dos: np.ndarray = np.array([])
        if array is not None:
            self.dos = array.transpose()

    def __len__(self) -> int:
        """x.__len__() <=> len(x)."""
        return len(self.dos)

    def fermi_correction(self, fermi: float) -> None:
        """Correct energy by Fermi level.

        Parameters
        ----------
        fermi: float
            fermi level

        """
        self.dos[0] -= fermi

    def energies(self, i: Optional[str] = None):
        r"""Return the *i*-th energy of the object.

        Parameters
        ----------
        i : int, optional (default is all energy)
            index #

        Returns
        --------
        np.ndarray
            the energy value of the i-th point when i set. \
            If arg is null, return the all energies in DOS object.

        """
        if i is None:
            return self.dos[0]
        else:
            return self.dos[0][i]

    def export_csv(self, filename: str, header: Optional[str] = None) -> None:
        """Export data to file object (or file-like object) as csv format."""
        transposed_dos = self.dos.transpose()
        with open(filename, mode="wb") as fhandle:
            np.savetxt(
                fhandle, transposed_dos, header=header, delimiter="\t", newline="\n"
            )


class TDOS(DOS):
    """Class for total DOS.

    Parameters
    ----------
    array: np.array
        DOS data

    Attributes
    ----------
    header

    """

    def __init__(self, array) -> None:
        """Initialize."""
        super(TDOS, self).__init__(array)
        if len(self.dos) == 2:
            self.header: str = "Energy\tTDOS"
        else:
            self.header = "Energy\tTDOS_up\tTDOS_down"

    def export_csv(self, filename: str) -> None:
        """Export data to file object (or file-like object) as csv format."""
        header = self.header
        super(TDOS, self).export_csv(filename, header=header)

    def graphview(self) -> None:
        """Show graphview by matplotlib."""
        for density in self.dos[1:]:
            plt.plot(self.dos[0], density)
        plt.show()


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
    array: numpy.ndarray
        DOS data
    site: str
        site name

    """

    def __init__(self, array=None, site: Optional[str] = None) -> None:
        """Initialize."""
        super(PDOS, self).__init__(array)
        self.site = "" if site is None else site
        self.orbital_spin: List[str] = list()
        orbitalnames = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2"]
        # The above order is refered from sphpro.F of vasp source
        spins_soi = ("mT", "mX", "mY", "mZ")
        spins = ("up", "down")
        if array is not None:
            if len(self.dos) == 10:
                self.orbital_spin = orbitalnames
            elif len(self.dos) == 19:  # Spin resolved
                self.orbital_spin = [
                    orb + "_" + spin for orb in orbitalnames for spin in spins
                ]
                # In collinear spin calculation, DOS of down-spin is
                # set by negative value.
                for i in range(2, 19, 2):
                    self.dos[i] *= -1
            elif len(self.dos) == 37:  # SOI
                self.orbital_spin = [
                    orb + "_" + spin for orb in orbitalnames for spin in spins_soi
                ]
            else:
                raise ValueError("Check the DOS data")

    def graphview(self, *orbitalnames: str) -> None:
        """Show DOS graph by using matplotlib.  For 'just seeing' use."""
        try:
            alist: List[int] = [
                self.orbital_spin.index(orbname) for orbname in orbitalnames
            ]
        except ValueError:
            err = "Check argment of this function\n"
            err += "The following name(s) are accpted:\n"
            err += ", ".join(self.orbital_spin)
            raise ValueError(err)
        for orbital in alist:
            plt.plot(self.dos[0], self.dos[orbital + 1])
        plt.show()

    def export_csv(self, filename: str, site: Optional[str] = None) -> None:
        """Export data to file object (or file-like object) as csv format.

        Parameters
        ----------
        filename: str
            filename for output

        """
        tmp = ["Energy"]
        for i in self.orbital_spin:
            if site is not None:
                tmp.append(site + "_" + i)
            elif self.site == "":
                tmp.append(i)
            else:
                tmp.append(self.site + "_" + i)
        header = "\t".join(tmp)
        super(PDOS, self).export_csv(filename, header=header)

    def plot_dos(
        self, orbitals: Union[List[str], Tuple[str, ...]], fermi: float = 0.0
    ) -> None:  # Not implemented yet
        """Plot DOS spectra with matplotlib.pyplot.

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

    def __add__(self, other: "PDOS") -> "PDOS":
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
            energies = copy.deepcopy(self.dos[0])
            sum_pdos.dos = self.dos + other.dos
            sum_pdos.dos[0] = energies
            sum_pdos.site = self.site + other.site
            sum_pdos.orbital_spin = self.orbital_spin
            return sum_pdos
