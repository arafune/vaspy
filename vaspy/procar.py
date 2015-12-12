#!/usr/bin/env python
# -*- coding: utf-8 -*-
# translate from procar.rb of scRipt4VASP 2014/2/26 master branch

from __future__ import print_function  # Version safety
from __future__ import division        # Version safety
import re
import copy
import os
import sys
import csv
import functools as ft
if sys.version_info[0] >= 3:     # Version safety
    from io import StringIO
else:
    from cStringIO import StringIO
import numpy as np
try:
    from vaspy import tools
except ImportError:
    mypath = os.readlink(__file__) if os.path.islink(__file__) else __file__
    sys.path.append(os.path.dirname(os.path.abspath(mypath)))
    import tools


class PROCAR(object):  # Version safety

    '''Class for PROCAR

    .. py:class:: PROCAR(object)

    PROCAR consists of these lines.  Appear once per file.

    1. the first line

      :Example:   PROCAR lm decomposed + phase

    2. set number of k-points, bands and ions.
       (Appear once when spin-integrated, twice when spin-resolved.)

      :Example:

        # of k-points:   50         # of bands: 576         # of ions:  98

    3.  k-point character

      :Example:

        k-point    1 :    0.00000 0.00000 0.00000 weight = 0.02000000

      .. note:: that the first character is "blank".

    4. band character

      :Example:  band   1 # energy  -11.87868466 # occ.  2.00000000

    5. orbital contribution.

      :Example:
    1  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000

    :author: Ryuichi Arafune
    '''

    def __init__(self, arg=None, phase_read=False):
        self.__orbital = list()
        self.__phase = list()
        self.__spininfo = 0   # nospin: 1, spinresolved: 2, soi: 4
        self.__numk = 0
        self.__nBands = 0
        self.__nAtoms = 0
        self.__kvectors = list()
        self.__energies = list()

        if isinstance(arg, str):
            self.load_from_file(arg, phase_read)
        else:
            raise RuntimeError("the arg of PROCAR() should be the filename")

    @property
    def orbital(self):
        return self.__orbital

    @property
    def phase(self):
        return self.__phase

    @property
    def spininfo(self):
        '''spin-related information'''
        return self.__spininfo

    @property
    def numk(self):
        '''Number of k-points'''
        return self.__numk

    @property
    def nBands(self):
        '''Number of Bands'''
        return self.__nBands

    @property
    def nAtoms(self):
        '''Number of atoms'''
        return self.__nAtoms

    @property
    def kvectors(self):
        '''k-vectors'''
        return self.__kvectors

    @property
    def energies(self):
        '''Energies'''
        return self.__energies

    @property
    def orb_names(self):
        '''Name of orbitals'''
        return self.__orb_names

    def load_from_file(self, file, phase_read=False):
        '''.. py:method::  load_fromfile(file, [phase_read=False])

        A virtual parser of PROCAR

        :param str file: filename of *PROCAR* file
        :param Boolian phase_read: Switch for loading phase characters
        '''
        f = open(file)
        first_line = f.readline()
        if 'PROCAR lm decomposed + phase' not in first_line:
            close(f)
            raise RuntimeError("This PROCAR is not a proper format\n \
                                See INCAR in the calculation.\n")
        section = list()
        with f:
            for line in f:
                if line.isspace():
                    continue
                elif "k-points: " in line:
                    self.__numk, self.__nBands, self.__nAtoms = [
                        int(i) for i in re.split('\s|:',line) if i.isdigit()]
                elif "k-point " in line:
                    self.__kvectors.append(np.array(
                        [float(i) for i in line.split()[3:6]]))
                    section = []
                elif "band " in line:
                    self.__energies.append(float(line.split()[4]))
                    section = []
                elif "ion" in line:
                    if "tot" in line:
                        section = ['orbital']
                        self.__orb_names = tuple(line.split()[1:])
                    else:
                        section = ['phase']
                else:
                    if section == ['orbital']:
                        if "tot " in line[0:4]:
                            continue
                        tmp = [float(i) for i in line.split()[1:]]
                        self.__orbital.append(tmp)
                    elif section == ['phase']:
                        if not phase_read:
                            continue
                        tmp = [float(i) for i in line.split()[1:]]
                        self.__phase.append(tmp)

        self.__spininfo = (len(self.orbital) //
                           (self.numk * self.nBands * self.nAtoms))
        if len(self.orbital) % (self.numk * self.nBands * self.nAtoms) != 0:
            raise RuntimeError("PROCAR file may be broken")
        if self.spininfo == 1:
            self.__spininfo = ('',)    # standard
        elif self.spininfo == 2:
            self.__spininfo = ('_up', '_down')   # collinear
        elif self.spininfo == 4:
            self.__spininfo = ('_mT', '_mX', '_mY', '_mZ')  # non-collinear

    def __str__(self):
        '''x.__str__() <=> str(x)

        show the PROCAR character, not contents.
        '''
        template = '''The properties of this procar:
  # of k-points: {0.numk}
  # of bands: {0.nBands}
  # of ions: {0.nAtoms}
  # of kvectors: {1}
  # of energies: {2}
    ((# of k-points) * (# of bands) = {0.numk}*{0.nBands}={3})
  # of orbital component: {4}
    ((# of k-points) * (# of bands) * (# of ions) =
        {0.numk}*{0.nBands}*{0.nAtoms}={5})
  # of phase component: {6}
  Orbitals are: {0.orb_names}
  spininfo: {0.spininfo}
'''
        return template.format(self,
                               len(self.kvectors),
                               len(self.energies),
                               self.numk * self.nBands,
                               len(self.orbital),
                               self.numk * self.nBands * self.nAtoms,
                               len(self.phase))

    def __iter__(self):
        return iter(self.__orbital)

    def band(self):
        '''Return BandStructure object

        :rtype: BandStructure
        '''
        band = BandStructure()
        band.kvectors = self.kvectors[0:self.numk]
        band.nBands = self.nBands
        band.nAtoms = self.nAtoms
        band.spininfo = self.spininfo
        band.orb_names = list(self.orb_names)
        # BandStructure.isready() must be True
        band.orbitals = self.orbital
        band.phases = self.phase
        band.energies = self.energies
        return band


class BandStructure(object):

    ''' .. py:class:: BandStructure class

    The "Band structure" object deduced from PROCAR file.

    :class variables: kvectors, energies, states, spin, orb_names

    .. note:: Finally, Band, Orbital, State classes can be removed ?
    '''

    def __init__(self, arg=None):
        self.__nBands = 0
        self.__kvectors = list()
        self.__kdistance = list()
        self.__sitecomposed = []
        self.__orbitals = 0
        self.__phases = 0
        self.__energies = 0
        self.orb_names = ['s', 'py', 'pz', 'px',
                          'dxy', 'dyz', 'dz2', 'dxz', 'dx2',
                          'tot']

    def isready(self):
        '''Return True if numk, nBands, nAtoms, spininfo, and
        orb_names are set, otherwise raise ValueError.

        This method is used for check before when orbitals, phases, energies
        are set.'''
        if not hasattr(self, "numk"):
            raise ValueError("numk is not defined")
        if self.numk == 0:
            raise ValueError("numk is not defined")
        if self.nBands == 0:
            raise ValueError("nBands is not correctly set")
        if not hasattr(self, "nAtoms"):
            raise ValueError("nAtoms is not defined")
        if not hasattr(self, "spininfo"):
            raise ValueError("spininfo is not defined")
        elif not type(self.spininfo) == tuple:
            raise TypeError("spininfo type should be tuple")
        if not hasattr(self, "orb_names"):
            raise ValueError("orb_names is not defined")
        if 's' not in self.orb_names:
            raise ValueError("orb_names is not correctly set")
        return True

    @property
    def nBands(self):
        '''Number of bands'''
        return self.__nBands

    @nBands.setter
    def nBands(self, arg):
        self.__nBands = arg
        self.available_band = list(range(self.nBands))

    @property
    def orbitals(self):
        '''Number of bands'''
        return self.__orbitals

    @orbitals.setter
    def orbitals(self, arg):
        '''Setter for orbitals

        When standard (i.e. ISPIN = 0) or SOI, returns 4-rank tensor.
        When spin-resolved (i.e. ISPIN = 2 but w/o SOI),
        returns tuple that consists of 4-rank tensor (np.ndarray)
        arg is usually listlike object that can smoothly convert to ndarray.
        for testing, 4-rank tensor (ndarray) is accepted.
        (But this option is not so entirely tested.
        Do not use for real analysis)
        '''
        if type(arg) == np.ndarray and arg.ndim == 4:
            self.__orbitals = arg
        elif self.isready():
            if len(self.spininfo) == 1 or len(self.spininfo) == 4:
                self.__orbitals = \
                    np.array(arg).reshape(self.numk,
                                          self.nBands,
                                          self.nAtoms * len(self.spininfo),
                                          len(self.orb_names))
            elif len(self.spininfo) == 2:
                self.__orbitals = \
                    np.array(arg).reshape(2, self.numk,
                                          self.nBands,
                                          self.nAtoms,
                                          len(self.orb_names))
                self.__orbitals = [self.__orbitals[0], self.__orbitals[1]]

    @property
    def sitecomposed(self):
        return self.__sitecomposed

    @property
    def phases(self):
        '''phase data

        .. note:: At present I have no idea about this parameter.  How to use it?'''
        return self.__phases

    @phases.setter
    def phases(self, arg):
        '''phases

        Return is 4-rank tensor for the standard (i.e. ISPIN = 0) or SOI
        calculations.
        Return the tuple that consists of 4-rank tensor for the
        spin-resolved (i.e. ISPIN = 2 but w/o SOI),

        arg must be the list of the list.
        Two elements convert into the single complex ndarray.
        '''
        if len(arg) == 0:
            self.__phases = list()
            return self.__phases
        else:
            phase_re = np.array(arg[::2])
            phase_im = np.array(arg[1::2])
            phases = phase_re + phase_im * (0.0 + 1.0J)
        if self.isready():
            if len(self.spininfo) == 1 or len(self.spininfo) == 4:
                self.__phases = phases.reshape(self.numk, self.nBands,
                                               self.nAtoms,
                                               len(self.orb_names) - 1)
            elif len(self.spininfo) == 2:
                self.__phases = phases.reshape(2, self.numk, self.nBands,
                                               self.nAtoms,
                                               len(self.orb_names) - 1)
                self.__phases = (self.__phases[0], self.__phases[1])

    @property
    def kvectors(self):
        return self.__kvectors

    @kvectors.setter
    def kvectors(self, kvectors):
        if type(kvectors) != list:
            errmsg = 'kvectors must be an array of ndarray\n'
            raise TypeError(errmsg)
        self.__kvectors = kvectors
        self.numk = len(self.kvectors)
        self.__kdistance = list()
        for i, k in enumerate(self.kvectors):
            if i == 0:
                self.__kdistance.append(0.0)
            else:
                self.__kdistance.append(
                    self.kdistance[i - 1] +
                    np.linalg.norm(self.kvectors[i - 1] - k))

    @property
    def kdistance(self):
        '''Distance of k.  This is used for the horizontal axis in the band
structure drowing'''
        return self.__kdistance

    @property
    def energies(self):
        '''Energies'''
        return self.__energies

    @energies.setter
    def energies(self, arg):
        if self.isready():
            if len(self.spininfo) == 1 or len(self.spininfo) == 4:
                self.__energies = np.array(arg).reshape(self.numk, self.nBands)
            elif len(self.spininfo) == 2:
                self.__energies = (
                    np.array(
                        arg[:self.numk * len(self.available_band)]).reshape(
                            self.numk, self.nBands),
                    np.array(
                        arg[self.numk * len(self.available_band):]).reshape(
                            self.numk, self.nBands))
                self.__energies = np.array(self.__energies)
                
    def fermi_correction(self, fermi):
        self.__energies -= fermi
                        
    def compose_sites(self, arg):
        '''Make sitecomposed ndarray

        When sitecomposed ndarray has elements, the values remain.

        :param arg: site indexes to be composed. it contains unique numbers.
        :type arg: list, tuple, set
        '''
        # the element of site_number_list must be unique.
        site_numbers = tuple(set(arg))
        self.isready()  # if not ready, raise Error.
        if len(self.spininfo) == 1:
            cmporbs = np.array([[[np.sum(
                [y for x, y in enumerate(self.orbitals[i, j])
                 if x in site_numbers],
                axis=0)]
                for j in range(len(self.available_band))]
                for i in range(self.numk)])
            if self.__sitecomposed:
                self.__sitecomposed[0] = np.concatenate(
                    (self.__sitecomposed[0], cmporbs),
                    axis=2)
            else:
                self.__sitecomposed = [cmporbs]
        if len(self.spininfo) == 2:
            upspin_orbitals = self.orbitals[0]
            downspin_orbitals = self.orbitals[1]
            cmporbsUp = np.array([[[np.sum(
                [y for x, y in enumerate(upspin_orbitals[i, j])
                 if x in site_numbers],
                axis=0)]
                for j in range(len(self.available_band))]
                for i in range(self.numk)])
            cmporbsDown = np.array([[[np.sum(
                [y for x, y in enumerate(downspin_orbitals[i, j])
                 if x in site_numbers],
                axis=0)]
                for j in range(len(self.available_band))]
                for i in range(self.numk)])
            self.__orbitals[0] = np.concatenate((self.__orbitals[0],
                                                 cmporbsUp),
                                                axis=2)
            self.__orbitals[1] = np.concatenate((self.__orbitals[1],
                                                 cmporbsDown),
                                                axis=2)
            if self.__sitecomposed:
                self.__sitecomposed[0] = np.concatenate(
                    (self.__sitecomposed[0], cmporbsUp),
                    axis=2)
                self.__sitecomposed[1] = np.concatenate(
                    (self.__sitecomposed[1], cmporbsDown),
                    axis=2)
            else:
                self.__sitecomposed = [cmporbsUp, cmporbsDown]
        if len(self.spininfo) == 4:
            site_numbers_mT = tuple(x + self.nAtoms * 0 for x in site_numbers)
            site_numbers_mX = tuple(x + self.nAtoms * 1 for x in site_numbers)
            site_numbers_mY = tuple(x + self.nAtoms * 2 for x in site_numbers)
            site_numbers_mZ = tuple(x + self.nAtoms * 3 for x in site_numbers)
            #
            cmporbs_mT = np.array([[[np.sum(
                [y for x, y in enumerate(self.orbitals[i, j])
                 if x in site_numbers_mT],
                axis=0)]
                for j in range(len(self.available_band))]
                for i in range(self.numk)])
            cmporbs_mX = np.array([[[np.sum(
                [y for x, y in enumerate(self.orbitals[i, j])
                 if x in site_numbers_mX],
                axis=0)]
                for j in range(len(self.available_band))]
                for i in range(self.numk)])
            cmporbs_mY = np.array([[[np.sum(
                [y for x, y in enumerate(self.orbitals[i, j])
                 if x in site_numbers_mY],
                axis=0)]
                for j in range(len(self.available_band))]
                for i in range(self.numk)])
            cmporbs_mZ = np.array([[[np.sum(
                [y for x, y in enumerate(self.orbitals[i, j])
                 if x in site_numbers_mZ],
                axis=0)]
                for j in range(len(self.available_band))]
                for i in range(self.numk)])
            if self.__sitecomposed:
                self.__sitecomposed[0] = np.concatenate(
                    (self.__sitecomposed[0], cmporbs_mT),
                    axis=2)
                self.__sitecomposed[1] = np.concatenate(
                    (self.__sitecomposed[1], cmporbs_mX),
                    axis=2)
                self.__sitecomposed[2] = np.concatenate(
                    (self.__sitecomposed[2], cmporbs_mY),
                    axis=2)
                self.__sitecomposed[3] = np.concatenate(
                    (self.__sitecomposed[3], cmporbs_mZ),
                    axis=2)
            else:
                self.__sitecomposed = [cmporbs_mT, cmporbs_mX, cmporbs_mY,
                                       cmporbs_mZ]

    def check_orb_name(self, arg):
        '''Return arg without change if arg is a member of the 'orbital name'.
        i.e., if arg is an alias of the (more appropriate) orbital
        name, return it as is.  If arg is neither the appropriate
        orbital name nor the alias, raise ValueError.

        :param arg: the string to be checked as the orbital name
        :type arg: str
        :rtype: str

        '''
        translate_dict = {'pypx': 'pxpy', 'pzpx': 'pxpz', 'pzpy': 'pypz',
                          'pxpypz': 'p', 'pxpzpy': 'p', 'pypxpz': 'p',
                          'pypzpx': 'p', 'pzpxpy': 'p', 'pzpypx': 'p',
                          'spd': 'tot'}
        proper_orb_name_list = self.orb_names + [
            'sp', 'p', 'pxpy', 'pxpz', 'pypz', 'spd', 'd']
        if arg in translate_dict.keys():
            arg = translate_dict[arg]
        if arg in proper_orb_name_list:
            return arg
        else:
            errmsg = arg + ": (composed) orbital name was not defined."
            raise ValueError(errmsg)

    def get_orb_index(self, arg):
        '''Return the indexes corresponding orbitan name by tuple.

        This method returns the tuple of orbital number in self.orb_names.
        (i.e. self.orb_names.index(orbitalname).  If the orbital name has not
        been in self.orb_names (i.e. if the orbital name is not used in
        PROCAR file) but the orbital name is appropriate as the composed
        orbital ((ex.) sp, pxpy), returns the indexes of the orbitals to be
        composed as the tuple.

        :param arg: name of (composed) orbital
        :type arg: str
        :return: number corresponding to (composed) orbital name.
        :rtype: tuple
        '''
        orbname = self.check_orb_name(arg)
        if (orbname in self.orb_names and
                self.orb_names.index(orbname) <= self.orb_names.index('tot')):
            orb_indexes = (self.orb_names.index(orbname),)
        elif orbname == 'p':
            orb_indexes = (self.orb_names.index('px'),
                           self.orb_names.index('py'),
                           self.orb_names.index('pz'))
        elif orbname == 'd':
            orb_indexes = (self.orb_names.index('dxy'),
                           self.orb_names.index('dyz'),
                           self.orb_names.index('dx2'),
                           self.orb_names.index('dxz'),
                           self.orb_names.index('dz2'))
        elif orbname == 'sp':
            orb_indexes = (self.orb_names.index('s'),
                           self.orb_names.index('px'),
                           self.orb_names.index('py'),
                           self.orb_names.index('pz'))
        elif orbname == 'pxpy':
            orb_indexes = (self.orb_names.index('px'),
                           self.orb_names.index('py'))
        elif orbname == 'pypz':
            orb_indexes = (self.orb_names.index('py'),
                           self.orb_names.index('pz'))
        elif orbname == 'pxpz':
            orb_indexes = (self.orb_names.index('px'),
                           self.orb_names.index('pz'))
        else:
            err = str(orbname) + " is not a proper (composed) orbital name."
            raise RuntimeError(err)
        return orb_indexes

    def compose_orbital(self, arg):
        '''Add composed orbital contribution in each 'sites' stored in BandStructure.sitecomposed.

        :param arg: orbital names
        :type arg: str, list, tuple
        '''
        if not self.sitecomposed:
            err = "This method operates with on sitecomposed attribute,"
            err += " but it's null"
            raise RuntimeError(err)
        if type(arg) == str:
            if ':' in arg:
                arg = arg.split(':')
            else:
                arg = [arg]
        for orb in arg:
            orb = self.check_orb_name(orb)
            if orb in self.orb_names:
                continue
            else:
                orbindex = self.get_orb_index(orb)
                # calculate composed orbital...
                for l in range(len(self.sitecomposed)):
                    numk, numband, numsite, norbs = self.sitecomposed[l].shape
                    orbitalcomposed = np.array([[[[np.sum(
                        [y for x, y in enumerate(self.sitecomposed[l][i, j, k])
                         if x in orbindex])]
                        for k in range(numsite)]
                        for j in range(numband)]
                        for i in range(numk)])
                    self.sitecomposed[l] = np.concatenate(
                        (self.sitecomposed[l], orbitalcomposed),
                        axis=3)
                self.orb_names.append(orb)

    def del_band(self, band_indexes):
        '''not yet impremented'''
        if not self.sitecomposed:
            err = "This method operates with on sitecomposed attribute,"
            err += " but it's null"
            raise RuntimeError(err)

    def set_header(self, sitenames, orbnames):
        '''Return header of table

        :param sitenames: site names e.g., 'Ag', 'Graphene', '2ndLayer'...
        :type sitenames: tuple, list
        :param orbnames: orbital names e.g., (('s', 'pxpy',), ('pxpy', 'p'))
        :type orbnames: tuple, list
        '''
        if len(sitenames) != len(orbnames):
            raise ValueError("Length of sitenames and orbnames is not same.")
        numk, numband, numsite, norbs = self.sitecomposed[0].shape
        if numsite != len(sitenames):
            err = "Number of sitename is different with"
            err += "the size of sitecomposed."
            raise ValueError(err)
        if not self.sitecomposed:
            err = "This method operates with on sitecomposed attribute,"
            err += " but it's null"
            raise RuntimeError(err)
        if len(self.spininfo) == 1 or len(self.spininfo) == 4:
            header = ["#k", "energy"]
            for i, spin in enumerate(self.spininfo):
                for site, orbs in zip(sitenames, orbnames):
                    for orb in orbs:
                        header.append(site + "_" + orb + spin)
        elif len(self.spininfo) == 2:
            header = ["#k"]
            for i, spin in enumerate(self.spininfo):
                header.append("energy" + spin)
                for site, orbs in zip(sitenames, orbnames):
                    for orb in orbs:
                        header.append(site + "_" + orb + spin)
        return header

    def get_orbnums(self, orbnames):
        '''Return tuple whose size is same as that of arg, but
        the element is number determied from orb_names

        :param orbnames: orbital names e.g., (('s','pxpy','tot'),('s','p'))
        :type orbnames: list, tuple
        :rtype: tuple
        '''
        return tuple(tuple(self.orb_names.index(orb) for orb in orbs)
                     for orbs in orbnames)

    def list_sitecomposed_data(self, orbnames):
        '''Return list of sitecomposed attribute to 2D-list

        :param orbnames: orbital names  e.g., (('s','pxpy','p'),('s','pxpy','p'))
        :type orbnames: list, tuple
        '''
        orbnums = self.get_orbnums(orbnames)
        numk, numband, numsite, norbs = self.sitecomposed[0].shape
        table = list()
        if len(self.spininfo) == 1 or len(self.spininfo) == 4:
            for b in range(numband):
                for k in range(numk):
                    sitelist = list()
                    sitelist.append(self.kdistance[k])
                    sitelist.append(self.energies[k, b])
                    for sitecomposed in self.sitecomposed:
                        for site, norbs in zip(list(range(numsite)),
                                               orbnums):
                            for o in norbs:
                                sitelist.append(
                                    sitecomposed[k, b, site, o])

                    table.append(sitelist)
                table.append([])
        elif len(self.spininfo) == 2:
            for b in range(numband):
                for k in range(numk):
                    sitelist = list()
                    sitelist.append(self.kdistance[k])
                    sitelist.append(self.energies[0][k, b])
                    for site, norbs in zip(list(range(numsite)),
                                           orbnums):
                        for o in norbs:
                            sitelist.append(
                                self.sitecomposed[0][k, b, site, o])
                    sitelist.append(self.energies[1][k, b])
                    for site, norbs in zip(list(range(numsite)),
                                           orbnums):
                        for o in norbs:
                            sitelist.append(
                                self.sitecomposed[1][k, b, site, o])
                    table.append(sitelist)
                table.append([])
        else:
            raise RuntimeError("spininfo is incorrect")
        return table

    def get_sitecomposed_data(self, sitenames, orbnames):
        header = map(str, self.set_header(sitenames, orbnames))
        output = "\t".join(header) + "\n"
        lists = self.list_sitecomposed_data(orbnames)
        for l in lists:
            l = map(str, l)
            output += "\t".join(l) + "\n"
        return output
