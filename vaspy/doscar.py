#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function  # Version safety
from __future__ import division  # Version safety
import sys as _sys
import numpy as np
# import matplotlib.pyplot as plt
if _sys.version_info[0] >= 3:  # Version safety
    from io import StringIO as _StringIO
else:
    from cStringIO import StringIO as _StringIO


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
        for i in range(4):
            line = thefile.readline()      # 2-5 lines
        line = thefile.readline()          # 6 line
        try:
            self.nenergypnts = int(line.split()[2])
        except ValueError:
            self.nenergypnts = int(line[32:37])
        for energypnts in range(self.nenergypnts):  # TDOS
            line = thefile.readline()
            data = line.split()
            if len(data) == 3:
                data = data[0:2]
            elif len(data) == 5:
                data = data[0:3]
            a_dos.append([np.float_(i) for i in data])
        self.dos_container.append(np.array(a_dos))
        #
        for iatom in range(self.natom):  # PDOS
            thefile.readline()
            a_dos = list()
            for energypts in range(self.nenergypnts):
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
            return self.dos[i][0]

    def export_csv(self, filename, **kwargs):
        '''..py:method:: export_csv(file, **kwargs)

        Export data to file object (or file-like object) as csv format.
        kwargs are keyword options of csv.writer().

        See help(csv.writer) for detail.
        '''
        with open(filename, mode='wb') as fhandle:
            np.savetxt(fhandle,
                       self.dos, kwargs,
                       delimiter='\t', newline='\n')

    def __str__(self):
        """x.__str__() <=> str(x)

        :return: string representation of DOS object (tab deliminated).
        :rtype: str
        """
        with _StringIO() as stream:
            self.export_csv(stream, delimiter='\t', lineterminator='\n')
            return stream.getvalue()


class TDOS(DOS):
    '''Class for total DOS

    :author: Ryuichi Arafune
    '''

    def __init__(self, array):
        super(TDOS, self).__init__(array)
        if len(self.dos[0][1]) == 2:
            self.header = "#Energy\tTDOS"
        else:
            self.header = "#Energy\tTDOS_up\tTDOS_down"

    def export_csv(self, filename):
        '''..py:export_csv:: export_csv(file, kwargs)

        Export data to file object (or file-like object) as csv format.
        kwargs are keyword options of csv.writer().
        '''
        header = self.header.encode('utf-8')
        super(TDOS, self).export_csv(filename, header=header)


class PDOS(DOS):
    '''
    Class for partial DOS

    :author: Ryuichi Arafune
    '''
    def __init__(self, array=None, site=None):
        super(PDOS, self).__init__(array)
        self.site = "" if site is None else site
        self.orbital_spin = list()
        orbitalname = ["s", "py", "pz", "px", "dxy", "dyz", "dz2",
                       "dxz", "dx2"]
        # the above order is refered from sphpro.F of vasp source
        soi = ("mT", "mX", "mY", "mZ")
        spin = ("up", "down")
        spininfo = []
        if array is not None:
            flag = len(self.dos[0][1])
            if flag == 10:
                self.orbital_spin = orbitalname
            elif flag == 19:  # Spin resolved
                spininfo = ['up', 'down']
                self.orbital_spin = [
                    orb + "_" + spn for orb in orbitalname
                    for spn in spin]
            elif flag == 37:  # SOI
                spininfo = ['mT', 'mX', 'mY', 'mZ']
                self.orbital_spin = [
                    orb + "_" + spn for orb in orbitalname
                    for spn in soi]
            else:
                self.orbital_spin = []
        if len(spininfo) == 2:  # In collinear spin calculation.
                                # DOS of down-spin is set by
                                # negative value.
            tmp = self.dos.transpose()
            for i in range(2, 19, 2):
                tmp[i] *= -1
            self.dos = tmp.transpose()

    def export_csv(self, file, site=None):  # Not implemented yet
        """Export data to file object (or file-like object) as csv format.
        kwargs are keyword options of csv.writer().
        see help(csv.writer) for detail.
        """
#        csvwriter = _csv.writer(file, **kwargs)
        tmp = ["#energy"]
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

    def __add__(self, other):   # Not implemented yet
        """x.__add__(y) <-> x+y

        :param addend: addend.energies.length must be equal to
        self.energies.length.
        :type addend: PDOS
        :return: PDOS
        :rtype: PDOS
        """
        if not isinstance(other, PDOS):
            return NotImplemented
        if self.dos == [] and self.site == "":
            return other.deepcopy()
        else:
            sumPDOS = PDOS()
            for s, a in zip(self, other):
                # sumPDOS.push( [s[0], s[1].zip(a[1]).map{|i| i[0]+i[1]}])
                sumPDOS.append([s[0],
                                [sum(each) for each in zip(s[1], a[1])]])
            sumPDOS.site = self.site + other.site
            sumPDOS.orbital_spin = self.orbital_spin
            return sumPDOS

    def __str__(self, site=None):  # Not implemented yet
        """x.__str__() <-> str(x)

        Returns String representation of PDOS object.
        :param site: Site name to overwrite.
        :return: String representation of PDOS object.
        :rtype: str
        """
        with _StringIO() as stream:
            self.export_csv(stream, site=site,
                            delimiter='\t', lineterminator='\n')
            return stream.getvalue()


'''
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
