#!/usr/bin/env python
# -*- coding: utf-8 -*-
''' ..py:module doscar

From VASP webpage:

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
in the parallel version if NPAR tex2html_wrap_inline5201 1.

Mind: For relaxations the DOSCAR is usually useless.
If you want to get an accurate DOS for the final configuration
copy CONTCAR to POSCAR and make another static (ISTART=1; NSW=0)
calculation.
'''

from __future__ import print_function  # Version safety
from __future__ import division  # Version safety
import copy
import numpy as np
import matplotlib.pyplot as plt
# if _sys.version_info[0] >= 3:  # Version safety
#     from io import StringIO as _StringIO
# else:
#     from cStringIO import StringIO as _StringIO


class DOSCAR(object):  # Version safety
    '''Class for DOSCAR file

    A container of DOS object

    :author: Ryuichi Arafune
    '''

    def __init__(self, arg=None):
        self.natom = 0
        self.nenergypnts = 0
        self.dos_container = list()

        if arg is not None:
            self.load_doscar_file(arg)

    def load_doscar_file(self, doscarfile):
        '''.. py:method load_doscar_file(doscarfile)

        Parse DOSCAR file and store it in memory

        :param doscar_file: file name of "DOSCAR"
        :type doscar_file: str
'''
        a_dos = list()
        thefile = open(doscarfile)
        line = thefile.readline()
        self.natom = int(line.split()[0])  # 1st line
        line = [thefile.readline() for i in range(4)]
        #for i in range(4):
        #    line = thefile.readline()      # 2-5 lines
        line = thefile.readline()          # 6 line
        try:
            self.nenergypnts = int(line.split()[2])
        except ValueError:
            self.nenergypnts = int(line[32:37])
        for dummy in range(self.nenergypnts):  # TDOS
            line = thefile.readline()
            data = line.split()
            if len(data) == 3:
                data = data[0:2]
            elif len(data) == 5:
                data = data[0:3]
            a_dos.append([np.float_(i) for i in data])
        self.dos_container.append(np.array(a_dos))
        #
        for dummy in range(self.natom):  # PDOS
            thefile.readline()
            a_dos = list()
            for dummy2 in range(self.nenergypnts):
                line = thefile.readline()
                data = line.split()
                a_dos.append([np.float_(i) for i in data])
            self.dos_container.append(np.array(a_dos))
        thefile.close()


class DOS(object):  # Version safety

    '''Class for DOS

    List object consisting two elements.
    The first element is the the energy.
    The latter element is list for the density.

    :attribute: dos
    '''

    def __init__(self, array=None):
        self.dos = np.array([])
        if array is not None:
            self.dos = array

    def __len__(self):
        """x.__len__() <=> len(x)"""
        return len(self.dos)

    def fermilevel_correction(self, fermi):
        '''..py:method:: fermilevel_correction(fermi)

        Fermi level Correction

        :param fermi: fermi level
        :type fermi: float
        '''
        tmp = self.dos.transpose()
        tmp[0] -= fermi
        self.dos = tmp.transpose()

    def energies(self, i=None):
        '''..py:method:: energies(i)

        Return the *i*-th energy of the object

        :param i: index #
        :type i: int
        :return: the energy value of the i-th point when i set.
                 If arg is null, return the all energies in DOS object.
        :rtype: np.ndarray
        '''
        tmp = self.dos.transpose()
        if i is None:
            return tmp[0]
        else:
            return self.dos[0][i]

    def export_csv(self, filename, header=None):
        '''..py:method:: export_csv(file[, header=header string])

        Export data to file object (or file-like object) as csv format.
        '''
        if header is None:
            with open(filename, mode='wb') as fhandle:
                np.savetxt(fhandle,
                           self.dos,
                           delimiter='\t', newline='\n')
        else:
            with open(filename, mode='wb') as fhandle:
                np.savetxt(fhandle,
                           self.dos, header=header,
                           delimiter='\t', newline='\n')

class TDOS(DOS):
    '''.. py:class:: TODS(array)

    Class for total DOS

    :param array: DOS data
    :type array: np.array

    .. py:attribute:: header
    '''


    def __init__(self, array):
        super(TDOS, self).__init__(array)
        if len(self.dos[0]) == 2:
            self.header = "Energy\tTDOS"
        else:
            self.header = "Energy\tTDOS_up\tTDOS_down"

    def export_csv(self, filename):
        '''..py:export_csv:: export_csv(filename)

        Export data to file object (or file-like object) as csv format.
        '''
        header = self.header
        super(TDOS, self).export_csv(filename, header=header)

    def graphview(self):
        ''' .. py:method:: graphview()

        Show graphview by using matplotlib'''
        data = self.dos.transpose()
        for density in data[1:]:
            plt.plot(data[0], density)
        plt.show()

class PDOS(DOS):
    '''.. py:class:: PDOS(array[, site])

    Class for partial DOS

    :param array: DOS data
    :type array: numpy.array
    :param site: site name
    :type site: str
    '''
    def __init__(self, array=None, site=None):
        super(PDOS, self).__init__(array)
        self.site = "" if site is None else site
        self.orbital_spin = list()
        orbitalnames = ["s", "py", "pz", "px", "dxy", "dyz", "dz2",
                        "dxz", "dx2"]
        # The above order is refered from sphpro.F of vasp source
        spins_soi = ("mT", "mX", "mY", "mZ")
        spins = ("up", "down")
        if array is not None:
            flag = len(self.dos[0])
            if flag == 10:
                self.orbital_spin = orbitalnames
            elif flag == 19:  # Spin resolved
                self.orbital_spin = [
                    orb + "_" + spin for orb in orbitalnames
                    for spin in spins]
                # In collinear spin calculation, DOS of down-spin is set by negative value.
                tmp = self.dos.transpose()
                for i in range(2, 19, 2):
                    tmp[i] *= -1
                self.dos = tmp.transpose()
            elif flag == 37:  # SOI
                self.orbital_spin = [
                    orb + "_" + spin for orb in orbitalnames
                    for spin in spins_soi]
            else:
                raise ValueError("Check the DOS data")

    def graphview(self, *orbitalnames):
        '''.. py:method:: graphview(**orbitalnames)

        Show DOS graph by using matplotlib.  For 'just seeing' use. '''
        try:
            alist = [self.orbital_spin.index(orbname) for orbname in orbitalnames]
        except ValueError:
            err = "Check argment of this function\n"
            err += "The following name(s) are accpted:\n"
            err += ", ".join(self.orbital_spin)
            raise ValueError(err)
        data = self.dos.transpose()
        for orbital in alist:
            plt.plot(data[0], data[orbital+1])
        plt.show()

    def export_csv(self, file, site=None):
        '''.. py:method:: export_csv(file[, site=site])

        Export data to file object (or file-like object) as csv format.
        '''
        tmp = ["Energy"]
        for i in self.orbital_spin:
            if site is not None:
                tmp.append(site + "_" + i)
            elif self.site == "":
                tmp.append(i)
            else:
                tmp.append(self.site + "_" + i)
        header = "\t".join(tmp)
        super(PDOS, self).export_csv(file, header=header)

    def plot_dos(self, orbitals, fermi=0.0):  # Not implemented yet
        '''..py:plot_dos(orbitals[, fermi=0.0])

        plot DOS spectra with matplotlib.pyplot

        :param orbitals: orbital name
        :type orbitals: str

        .. warning:: not implemented yet!!

        '''
        pass

    def __add__(self, other):
        '''x.__add__(y) <-> x+y

        :param addend: addend.energies.length must be equal to \
        self.energies.length.
        :type addend: PDOS
        :return: PDOS
        :rtype: PDOS
        '''
        if not isinstance(other, PDOS):
            return NotImplemented
        if self.dos == [] and self.site == "":
            return copy.deepcopy(other)
        else:
            sum_pdos = PDOS()
            dos = self.dos.transpose()
            otherdos = other.dos.transpose()
            energies = copy.deepcopy(dos[0])
            sum_pdos.dos = dos + otherdos
            sum_pdos.dos[0] = energies
            sum_pdos.dos = sum_pdos.dos.transpose()
            sum_pdos.site = self.site + other.site
            sum_pdos.orbital_spin = self.orbital_spin
            return sum_pdos
