#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function  # Version safety
from __future__ import division  # Version safety
import copy as _copy
import csv as _csv
import sys as _sys
import numpy as np
# import matplotlib.pyplot as plt
from numbers import Number as _Number
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
        self.nAtom = 0
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
        aDOS = list()
        thefile = open(doscarfile)
        line = thefile.readline()
        self.nAtom = int(line.split()[0])  # 1st line
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
            aDOS.append([np.float_(i) for i in data])
        self.dos_container.append(np.array(aDOS))
        #
        for iatom in range(self.nAtom):  # PDOS
            thefile.readline()
            aDOS = list()
            for energypts in range(self.nenergypnts):
                line = thefile.readline()
                data = line.split()
                aDOS.append([np.float_(i) for i in data])
            self.dos_container.append(np.array(aDOS))
        thefile.close()


    def load_from_file(self, doscar_file):
        """Parse DOSCAR file to create dos_container

        :param doscar_file: file name of "DOSCAR"
        :type doscar_file: str
        :rtype: DOSCAR
        """
        with open(doscar_file) as the_doscar:
            separate_text = ""
            aDOS = list()
            for idx, line in enumerate(the_doscar, 1):
                if idx == 1:
                    self.nAtom = int(line.split()[0])
                elif 2 <= idx <= 5:
                    continue
                elif idx == 6:
                    try:
                        separate_text = line.rstrip('\n')
                        self.nenergypnts = int(separate_text.split()[2])
                    except ValueError:
                        self.nenergypnts = int(line[32:37])
                    continue
                else:  # idx >= 7
                    if separate_text in line:
                        self.dos_container.append(aDOS)
                        aDOS = list()
                    else:
                        aDOS.append([float(i) for i in line.split()])
            else:
                self.dos_container.append(aDOS)
                return self.dos_container.append(aDOS)


class DOS(object):  # Version safety

    '''Class for DOS

    List object consisting two elements.
    The first element is the the energy.
    The latter element is list for the density.

    :attribute: dos
    '''

    def __init__(self, array=None):
        self.dos = 0
        if array is not None:
            self.dos = array

    # def __deepcopy__(self, memo):
    #     "x.__deepcopy__() <-> copy.deepcopy(x)"
    #     dest = DOS()
    #     dest.dos = _copy.deepcopy(self.dos, memo)
    #     return dest

    # def deepcopy(self):
    #     '''call copy.deepcopy'''
    #     return _copy.deepcopy(self)

    # def __iter__(self):
    #     for ith_point in self.dos:
    #         yield ith_point

    # def append(self, dos_data):
    #     '''wrapper of list.append()

    #     :param dos_data: dos_data
    #     :type dos_data: np.array
    #     '''
    #     self.dos.append(_filter_dos_data(dos_data))

    # def pop(self, i=-1):
    #     '''Wrapper of list.pop()

    #     :return: return and remove the last element
    #     :rtype: np.array

    #     (the highest energy data) of the DOS object.
    #     '''
    #     return self.dos.pop(i)

    # self.unshift(X) => self.[0:0] = X
    # self.shift => self.pop(0)

    def __len__(self):
        """x.__len__() <=> len(x)"""
        return len(self.dos)

    # def __setitem__(self, idx, value):
    #     """x.__setitem__(i, y) <=> x[i]=y"""
    #     if isinstance(idx, int):
    #         value = _filter_dos_data(value)
    #     if isinstance(idx, slice):
    #         value = [_filter_dos_data(each) for each in value]
    #     self.dos.__setitem__(idx, value)

    # def __getitem__(self, idx):
    #     """x.__getitem__(i) <=> x[i]"""
    #     self.dos.__getitem__(idx)

    def fermilevel_correction(self, fermi):
        '''Fermi level Correction

        :param fermi: fermi level
        :type fermi: float
        '''
        tmp = self.dos.transpose()
        tmp[0] -= fermi
        self.dos = tmp.transpose()

    def energies(self, i=None):
        '''Return the *i*-th energy of the object

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

    def densities(self, i=None):
        '''Return the *i*-th density of the object

        :param i: index #
        :type i: int
        :return: the density of states for the i-th point when i set.
                 If arg is null, return the all energies in DOS object.
        :rtype: numpy.array
        '''
        if i is None:
            return [each[1] for each in self.dos]
        else:
            return self.dos[i][1]

    def export_csv(self, file, **kwargs):
        """Export data to file object (or file-like object) as csv format.
        kwargs are keyword options of csv.writer().

        See help(csv.writer) for detail.
        """
        csvwriter = _csv.writer(file, **kwargs)
        csvwriter.writerows([line[0]] + line[1] for line in self.dos)

    def __str__(self):
        """x.__str__() <=> str(x)

        :return: string representation of DOS object (tab deliminated).
        :rtype: str
        """
        with _StringIO() as stream:
            self.export_csv(stream, delimiter='\t', lineterminator='\n')
            return stream.getvalue()


def _filter_dos_data(data):
    listlikes = (list, tuple)
    if not isinstance(data, listlikes):
        raise RuntimeError("Invalid argument")
    if not isinstance(data[0], _Number):
        raise RuntimeError("non Numeric instance in header.")
    if not isinstance(data[1], listlikes):
        raise RuntimeError("non list/tuple instance in data.")
    return data


class TDOS(DOS):

    """Class for total DOS

    :author: Ryuichi Arafune
    """

    def __init__(self, array):
        super(TDOS, self).__init__(array)
        if len(self.dos[0][1]) == 2:
            self.header = ("TDOS", "intTDOS")
        else:
            self.header = ("TDOS_up", "TDOS_down",
                           "intTDOS_up", "intTDOS_down")

    def export_csv(self, file, **kwargs):
        """Export data to file object (or file-like object) as csv format.
        kwargs are keyword options of csv.writer().
        see help(csv.writer) for detail.
        """
        csvwriter = _csv.writer(file, **kwargs)
        csvwriter.writerow(["#energy"] + list(self.header))
        super(TDOS, self).export_csv(file, **kwargs)


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
        # this order above is refered from sphpro.F of vasp source
        soi = ("mT", "mX", "mY", "mZ")
        spin = ("up", "down")
        spininfo = []
        if array is not None:
            flag = len(self.dos[0][1])
            if flag == 9:
                self.orbital_spin = orbitalname
            elif flag == 18:  # Spin resolved
                spininfo = ['up', 'down']
                self.orbital_spin = [
                    orb + "_" + spn for orb in orbitalname
                    for spn in spin]
            elif flag == 36:  # SOI
                spininfo = ['mT', 'mX', 'mY', 'mZ']
                self.orbital_spin = [
                    orb + "_" + spn for orb in orbitalname
                    for spn in soi]
            else:
                self.orbital_spin = []

        if len(spininfo) == 2:  # In collinear spin calculation.
                                    # DOS of down-spin is set by
                                    # negative value.
            for line in self.dos:
                for index, density in enumerate(line[1]):
                    if index % 2 != 0:
                        line[1][index] = -density

    def __deepcopy__(self, memo):
        "x.__deepcopy__() <-> copy.deepcopy(x)"
        dest = PDOS()
        dest.dos = _copy.deepcopy(self.dos, memo)
        dest.site = _copy.deepcopy(self.site, memo)
        dest.orbital_spin = _copy.deepcopy(self.orbital_spin, memo)
        return dest

    def plot_dos(self, orbitals, fermi=0.0):
        """plot DOS spectra with matplotlib.pyplot

        :param orbitals: orbital name
        :type orbitals: str

        .. warning:: not implemented yet!!

        """
        pass

    def __add__(self, other):
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

    def export_csv(self, file, site=None, **kwargs):
        """Export data to file object (or file-like object) as csv format.
        kwargs are keyword options of csv.writer().
        see help(csv.writer) for detail.
        """
        csvwriter = _csv.writer(file, **kwargs)
        tmp = ["#energy"]
        for i in self.orbital_spin:
            if site is not None:
                tmp.append(site + "_" + i)
            elif self.site == "":
                tmp.append(i)
            else:
                tmp.append(self.site + "_" + i)
        csvwriter.writerow(tmp)

        super(PDOS, self).export_csv(file, **kwargs)

    def __str__(self, site=None):
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

# -----------------------------------------------------------------------

if __name__ == '__main__':
    import argparse
    import os as _os
    _path = _os.path
    mypath = _os.readlink(__file__) if _path.islink(__file__) else __file__
    _sys.path.append(_path.dirname(_path.abspath(mypath)))
    from outcar import OUTCAR as _OUTCAR
    import tools as _tools

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('doscar', metavar='DOSCAR_file')
    parser.add_argument('-o', '--outcar', metavar='OUTCAR_file',
                        nargs='?', default=NotImplemented,
                        help="""Use OUTCAR file.
If file path is not given,
try to open OUTCAR file
in the same directory to DOSCAR.""")
    # -o data/OUTCAR => args.outcar == "data/OUTCAR"
    # -o             => args.outcar == None
    # (not given)    => args.outcar == NotImplemented (default)
    parser.add_argument('-f', '--fermi', metavar='value', type=float,
                        default=0.0,
                        help="""Fermi level correction.
Energy shifts by this value.
if --outcar is set, this option is ignored.""")
    parser.add_argument('-s', '--site', metavar='atoms', action='append',
                        dest='atomset', type=_tools.parse_AtomselectionNum,
                        help="""atom # specified with range.
Use "-" or ","
(ex.) --site 1,2,7-9""")
    parser.add_argument('-a', '--as', metavar='name', action='append',
                        dest='atomsetname',
                        help="""the name of the range identified by --site.
(ex.) --as layer1
the name is used in the output filename.""")

    args = parser.parse_args()

    doscar = DOSCAR(args.doscar)
    atomlist = list()
    if args.outcar is not NotImplemented:
        if args.outcar is None:
            args.outcar = _path.join(_path.dirname(args.doscar), "OUTCAR")
        outcar = _OUTCAR(args.outcar)
        atomlist = outcar.atom_names
        args.fermi = outcar.fermi
    #
    if atomlist == []:
        atomlist.extend("atom" + str(i) for i in range(1, doscar.nAtom + 1))
    #
    # construct TDOS & PDOS objects
    #
    tmp = doscar.dos_container
    d = [TDOS(tmp.pop(0))]
    #
    d.extend(PDOS(*each) for each in zip(tmp, atomlist))  # tmp[1:] ?
    #
    if args.atomset is not None:
        if len(args.atomset) == len(args.atomsetname):
            sumPDOSs = list()
            for site, name in zip(args.atomset, args.atomsetname):
                each = PDOS()
                for atomNo in site:
                    # print(repr(each.dos), repr(each.site))
                    each += d[atomNo]
                each.site = name
                sumPDOSs.append(each)
            for summedPDOS in sumPDOSs:
                filename = summedPDOS.site + ".dat"
                try:  # Version safety
                    file = open(filename, mode='w', newline='')
                except TypeError:
                    file = open(filename, mode='wb')
                with file:
                    summedPDOS.export_csv(file, delimiter='\t')
                    # file.write(str(summedPDOS))
    #
    try:  # Version safety
        file = open("total.dat", mode='w', newline='')
    except TypeError:
        file = open("total.dat", mode='wb')
    with file:
        d[0].fermilevel_correction(args.fermi)
        d[0].export_csv(file, delimiter='\t')
        # file.write(str(d[0]))
    for i, n in zip(d[1:], atomlist):
        i.fermilevel_correction(args.fermi)
        if isinstance(i, PDOS) and i.site == "":
            i.site = n
        filename = n + "_dos.dat"
        try:  # Version safety
            file = open(filename, mode='w', newline='')
        except TypeError:
            file = open(filename, mode='wb')
        with file:
            i.export_csv(file, delimiter='\t')
            # file.write(str(i))


"""
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
"""
