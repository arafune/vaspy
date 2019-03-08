#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module provides DOSCAR and related classes.

From VASP webpage::

  The file DOSCAR contains the DOS and integrated DOS The units are
  "number of states/unit cell". For dynamic simulations and relaxations,
  an averaged DOS and an averaged integrated DOS is written to the file.
  For a description of how the averaging is done see 7.18, 7.32).
  The first few lines of the DOSCAR file are made up by a header,
  which is followed by NDOS lines holding three data

   energy     dos     integrated_dos

  For spin-polarized calculations each line holds five data

   energy     dos(up) dos(dwn)  integrated_dos(up) integrated_dos(dwn)

  If RWIGS (Wigner Seitz radii, see section 7.29) is set in the INCAR file,
  a l- and site-projected DOS is calculated and also written to the file
  DOSCAR. One set of data is written for each ion,
  each set of data holds NDOS lines with the following data

   energy s-dos p-dos d-dos

  and

   energy s-dos(up) p-dos(up) d-dos(up) s-dos(dwn) p-dos(dwn) d-dos(dwn)

  for the non spin-polarized and spin polarized case respectively.
  The units of the l- and site projected DOS are states/atom.
  Please mind, that the site projected DOS is not evaluated
  in the parallel version if NPAR is not equal to 1.

  Mind: For relaxations the DOSCAR is usually useless.
  If you want to get an accurate DOS for the final configuration
  copy CONTCAR to POSCAR and make another static (ISTART=1; NSW=0)
  calculation.
"""

from __future__ import division  # Version safety
from __future__ import print_function  # Version safety

import bz2
import copy
import os
import sys

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.stderr.write(
        'Install matplotlib, or you cannot use methods relating to draw\n')

# if _sys.version_info[0] >= 3:  # Version safety
#     from io import StringIO as _StringIO
# else:
#     from cStringIO import StringIO as _StringIO


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

    def __init__(self, filename=None):
        """Initialize."""
        self.natom = 0
        self.nbands = 0
        self.dos_container = list()

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
        tdos = np.array([next(thefile).split() for i in range(self.nbands)],
                        dtype=np.float64)
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
                np.array([next(thefile).split() for i in range(self.nbands)],
                         dtype=np.float64))
            try:
                line = next(thefile)
            except StopIteration:
                line = ""


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

    def __init__(self, array=None):
        """Initialize."""
        self.dos = np.array([])
        if array is not None:
            self.dos = array.transpose()

    def __len__(self):
        """x.__len__() <=> len(x)"""
        return len(self.dos)

    def fermi_correction(self, fermi):
        """Correct energy by Fermi level.

        Parameters
        ----------
        fermi: float
            fermi level

        """
        self.dos[0] -= fermi

    def energies(self, i=None):
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

    def export_csv(self, filename, header=None):
        """Export data to file object (or file-like object) as csv format."""
        transposed_dos = self.dos.transpose()
        if header is None:
            with open(filename, mode='wb') as fhandle:
                np.savetxt(
                    fhandle, transposed_dos, delimiter='\t', newline='\n')
        else:
            with open(filename, mode='wb') as fhandle:
                np.savetxt(
                    fhandle,
                    transposed_dos,
                    header=header,
                    delimiter='\t',
                    newline='\n')


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

    def __init__(self, array):
        """Initialize."""
        super(TDOS, self).__init__(array)
        if len(self.dos) == 2:
            self.header = "Energy\tTDOS"
        else:
            self.header = "Energy\tTDOS_up\tTDOS_down"

    def export_csv(self, filename):
        """Export data to file object (or file-like object) as csv format."""
        header = self.header
        super(TDOS, self).export_csv(filename, header=header)

    def graphview(self):
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

    def __init__(self, array=None, site=None):
        """Initialize."""
        super(PDOS, self).__init__(array)
        self.site = "" if site is None else site
        self.orbital_spin = list()
        orbitalnames = [
            "s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2"
        ]
        # The above order is refered from sphpro.F of vasp source
        spins_soi = ("mT", "mX", "mY", "mZ")
        spins = ("up", "down")
        if array is not None:
            flag = len(self.dos)
            if flag == 10:
                self.orbital_spin = orbitalnames
            elif flag == 19:  # Spin resolved
                self.orbital_spin = [
                    orb + "_" + spin for orb in orbitalnames for spin in spins
                ]
                # In collinear spin calculation, DOS of down-spin is
                # set by negative value.
                for i in range(2, 19, 2):
                    self.dos[i] *= -1
            elif flag == 37:  # SOI
                self.orbital_spin = [
                    orb + "_" + spin for orb in orbitalnames
                    for spin in spins_soi
                ]
            else:
                raise ValueError("Check the DOS data")

    def graphview(self, *orbitalnames):
        """Show DOS graph by using matplotlib.  For 'just seeing' use."""
        try:
            alist = [
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

    def export_csv(self, filename, site=None):
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

    def plot_dos(self, orbitals, fermi=0.0):  # Not implemented yet
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

    def __add__(self, other):
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
