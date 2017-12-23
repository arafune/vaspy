#!/usr/bin/env python
# -*- coding: utf-8 -*-
# translate from procar.rb of scRipt4VASP 2014/2/26 master branch
""".. py:module:: procar.py

This module provides PROCAR, BandWithProjection, EnergyBand,
Projection classes.

* PROCAR class is used to stored the PROCAR information in the memory.
* Band_with_Projection class is used in vaspy-procar.py script.
    *  This class is essentially used on ipython, but not so easy for use.
* EnergyBand class is used for drawing the energy band.
* Projection class is used for storing the orbital projection data
"""

from __future__ import print_function  # Version safety
from __future__ import division        # Version safety
import re
# import copy
import os
import bz2
# import csv
# import functools as ft
import numpy as np
from vaspy import eigenval


class PROCAR(eigenval.EIGENVAL):  # Version safety
    '''.. py:class:: PROCAR(PROCAR_file[, phase_read])

    Class for storing the data saved in PROCAR file.

    Parameters
    -----------

    PROCAR_file: str
        File name of "PROCAR".
    phase_read: boolean
        Set True is you read phase data.


    PROCAR consists of the following lines.  Appear once per file.

    1. The first line is used just as a comment.

      :Example:   PROCAR lm decomposed + phase

    2. Number of k-points, bands and ions.
       (Appear once when spin-integrated, twice when spin-resolved.)

      :Example:

        # of k-points:   50         # of bands: 576         # of ions:  98

    3.  k-point character

      :Example:

        k-point    1 :    0.00000 0.00000 0.00000 weight = 0.02000000

        .. Note::  That the first character must be "blank".

    4. Band character

      :Example:  band   1 # energy  -11.87868466 # occ.  2.00000000

    5. orbital contribution.

      :Example: 1 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000\
      0.000 0.000

    .. py:attribute:: orbital

      list object. Store orbital character

    .. py:attribute:: phase

      list (Not used at present)

    .. py:attribute:: spininfo

      tuple

        * nospin: ('',)
        * spinresolved: ('_up', '_down')
        * soi: ('_mT', '_mX', '_mY', '_mZ')

    '''

    def __init__(self, arg=None, phase_read=False):
        super(PROCAR, self).__init__(None)
        self.orbital = list()
        self.phase = list()
        self.orb_names = tuple()
        #
        if isinstance(arg, str):
            self.load_from_file(arg, phase_read)
        else:
            raise RuntimeError("the arg of PROCAR() should be the filename")

    def load_from_file(self, filename, phase_read=False):
        '''.. py:method::  load_from_file(filename, [phase_read=False])

        A virtual parser of PROCAR

        Parameters
        ----------

        filename: str
             Filename of *PROCAR* file
        phase_read: boolean
             Switch for loading phase characters
        '''
        if os.path.splitext(filename)[1] == ".bz2":
            try:
                procar_file = bz2.open(filename, mode='rt')
            except AttributeError:
                procar_file = bz2.BZ2File(filename, mode='r')
        else:
            procar_file = open(filename)
        first_line = next(procar_file)
        if 'PROCAR lm decomposed + phase' not in first_line:
            procar_file.close()
            raise RuntimeError("This PROCAR is not a proper format\n \
                                Check your INCAR the calculations.\n")
        section = list()
        with procar_file:
            for line in procar_file:
                if line.isspace():
                    continue
                elif "k-points: " in line:
                    self.numk, self.n_bands, self.n_atoms = [
                        int(i) for i in re.split(r'\s|:', line)
                        if i.isdigit()]
                elif "k-point " in line:
                    try:
                        self.kvectors.append(np.array(
                            [float(i) for i in line.split()[3:6]]))
                        section = []
                    except ValueError:
                        self.kvectors.append(np.array(
                            [np.float_(line[18:29]),
                             np.float_(line[29:40]),
                             np.float_(line[40:51])]))
                        section = []
                elif "band " in line:
                    self.energies.append(float(line.split()[4]))
                    section = []
                elif "ion" in line:
                    if "tot" in line:
                        section = ['orbital']
                        self.orb_names = tuple(line.split()[1:])
                    else:
                        section = ['phase']
                else:
                    if section == ['orbital']:
                        if "tot " in line[0:4]:
                            continue
                        self.orbital.append([float(i)
                                             for i in line.split()[1:]])
                    elif section == ['phase']:
                        if not phase_read:
                            continue
                        self.phase.append([float(i)
                                           for i in line.split()[1:]])
#
        self.spininfo = (len(self.orbital) //
                         (self.numk * self.n_bands * self.n_atoms))
        if len(self.orbital) % (self.numk * self.n_bands * self.n_atoms) != 0:
            raise RuntimeError("PROCAR file may be broken")
        if self.spininfo == 1:  # standard
            self.spininfo = ('',)
        elif self.spininfo == 2:   # collinear
            self.spininfo = ('_up', '_down')
            tmp = list(zip(*[iter(self.energies)]*self.numk*self.n_bands))
            self.energies = []
            for up, down in zip(tmp[0], tmp[1]):
                self.energies.append([up, down])
        elif self.spininfo == 4:  # non-collinear
            self.spininfo = ('_mT', '_mX', '_mY', '_mZ')

    def __str__(self):
        '''.. py:method:: __str__()

        __str__() <=> str(x)

        show the PROCAR character, not contents.
        '''
        template = '''The properties of this procar:
  # of k-points: {0.numk}
  # of bands: {0.n_bands}
  # of ions: {0.n_atoms}
  # of kvectors: {1}
  # of energies: {2}
    ((# of k-points) * (# of bands) = {0.numk}*{0.n_bands}={3})
  # of orbital component: {4}
    ((# of k-points) * (# of bands) * (# of ions) =
        {0.numk}*{0.n_bands}*{0.n_atoms}={5})
  # of phase component: {6}
  Orbitals are: {0.orb_names}
  spininfo: {0.spininfo}
'''
        return template.format(self,
                               len(self.kvectors),
                               len(self.energies),
                               self.numk * self.n_bands,
                               len(self.orbital),
                               self.numk * self.n_bands * self.n_atoms,
                               len(self.phase))

    def __iter__(self):
        return iter(self.orbital)

    def band(self, recvec=[[1.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0],
                           [0.0, 0.0, 1.0]]):
        '''.. py:method:: band()

        Return Band_with_projection object

        Parameters
        -----------

        recvec: array, numpy.ndarray
            reciprocal vector.

            .. Note:: Don't forget that the reciprocal vector
                      used in VASP need 2Pi to match
                      the conventional unit of the wavevector.

        Returns
        -------

        BandWithProjection
        '''
        recvecarray = np.array(recvec).T
        band = BandWithProjection()
        band.kvectors = [recvecarray.dot(kvector) for kvector in
                         self.kvectors[0:self.numk]]
        band.n_bands = self.n_bands
        band.n_atoms = self.n_atoms
        band.spininfo = self.spininfo
        band.orb_names = list(self.orb_names)
        # BandStructure.isready() must be True
        band.orbitals = self.orbital
        band.phases = self.phase
        band.energies = self.energies
        return band


class Projection(object):
    '''.. py:class:: Projection(projections, natom, numk, nbands[, soi])

    Orbital projection object for analyzing by using python.
'''

    def __init__(self, projection, natom=0, numk=0, nbands=0, soi=False):
        self.proj = np.array(projection)
        self.natom = natom
        self.numk = numk
        self.nbands = nbands
        self.soi = soi
        self.output_states = 0
        self.output_headers = []
        if soi:
            if 4 * natom * numk * nbands == len(projection):
                self.proj = self.proj.reshape(numk,
                                              nbands,
                                              natom,
                                              40).transpose()
            else:
                raise ValueError("Argments are mismatched.")
        else:
            if natom * numk * nbands == len(projection):
                self.proj = self.proj.reshape(numk,
                                              nbands,
                                              natom,
                                              10).transpose()
                # orbitals[orbindex][siteindex-1][bandindex][kindex]
            else:
                raise ValueError("Argments are mismatched.")

    def sum_states(self, states, axis=None):
        '''.. py:method:: sum_states(states[, axis=axis])

        Return summantion of states

        Parameters
        ----------

        states: tuple
             tuple of tuple of site index and orbital name
        axis: str
             quantization axis for SOI calculation. 'x', 'y' or 'z'.

        Note
        ------

        site index starts '1' not 0.

        Example
        --------

        ((1, "px"), (1, "py")) # to produce pxpy projection on the 1st element.

        ((1, "pz"), (43, "pz")) # to produce pz orbital projection \
        for 'surface' atom (Here, the 43 atoms is included in the unit cell \
        and 1st and 43the atoms are assumed to be identical.)
        '''
        result = 0
        orb_names = ['s', 'py', 'pz', 'px',
                     'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot']
        for a_state in states:
            if isinstance(a_state[1], int) and 0 <= a_state[1] <= 9:
                orbindex = a_state[1]
            elif a_state[1] in orb_names:
                orbindex = orb_names.index[a_state[1]]
            else:
                raise ValueError("Check your input for orbital name")
            if self.soi and (axis == 'x' or axis == 'X' or axis == 0):
                orbindex += 10
            elif self.soi and (axis == 'y' or axis == 'Y' or axis == 1):
                orbindex += 20
            elif self.soi and (axis == 'z' or axis == 'Z' or axis == 2):
                orbindex += 30
            result += self.proj[orbindex][a_state[0] - 1]
        return np.array([result])

    def add_output_states(self, name, state):
        '''.. py: method:: add_output_states(name, state)

        Construct states for output.

        Parameters
        -----------

        name: str
             name of the (summed) state.  Used for the header.
        state: numpy.ndarray
             orbital projection data.  From sum_states method.
        '''
        if name in self.output_headerlist:
            raise ValueError("Unique state nume is required")
        self.output_headers.append(name)
        try:
            self.output_states = np.concatenate((self.output_states, state))
        except (TypeError, ValueError):
            self.output_states = state
            self.output_headers = [name]


class BandWithProjection(object):
    '''.. py:class:: BandWithProjection()

    The "Band structure" object deduced from PROCAR file.

    This class provides the way Band data with orbital projection. And
    most of the functions are used in vaspy-procar.py script.  As this
    class is not so well sophisticated, the use is not easy.  Anyway,
    it works, however.

    :class variables: kvectors, energies, states, spin, orb_names

    Note
    -------

    Band, Orbital, State classes can be removed ?

    '''

    def __init__(self):
        self.__n_bands = 0
        self.__kvectors = list()
        self.kdistance = list()
        self.sitecomposed = []
        self.__orbitals = 0
        self.__phases = 0
        self.__energies = 0
        self.numk = 0
        self.n_atoms = 0
        self.spininfo = tuple()
        self.orb_names = ['s', 'py', 'pz', 'px',
                          'dxy', 'dyz', 'dz2', 'dxz', 'dx2',
                          'tot']

    def isready(self):
        '''.. py:method:: isready()

        Return True if numk, n_bands, n_atoms, spininfo, and
        orb_names are set, otherwise raise ValueError.

        Use for check before when orbitals, phases, energies
        are set.'''
        if not hasattr(self, "numk"):
            raise ValueError("numk is not defined")
        if self.numk == 0:
            raise ValueError("numk is not defined")
        if self.n_bands == 0:
            raise ValueError("n_bands is not correctly set")
        if not hasattr(self, "n_atoms"):
            raise ValueError("n_atoms is not defined")
        if not hasattr(self, "spininfo"):
            raise ValueError("spininfo is not defined")
        elif not isinstance(self.spininfo, tuple):
            raise TypeError("spininfo type should be tuple")
        if not hasattr(self, "orb_names"):
            raise ValueError("orb_names is not defined")
        if 's' not in self.orb_names:
            raise ValueError("orb_names is not correctly set")
        return True

    @property
    def n_bands(self):
        '''Number of bands'''
        return self.__n_bands

    @n_bands.setter
    def n_bands(self, arg):
        self.__n_bands = arg
        self.available_band = list(range(self.n_bands))

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
        if isinstance(arg, np.ndarray) and arg.ndim == 4:
            self.__orbitals = arg
        elif self.isready():
            if len(self.spininfo) == 1 or len(self.spininfo) == 4:
                self.__orbitals = \
                    np.array(arg).reshape(self.numk,
                                          self.n_bands,
                                          self.n_atoms * len(self.spininfo),
                                          len(self.orb_names))
            elif len(self.spininfo) == 2:
                self.__orbitals = \
                    np.array(arg).reshape(2, self.numk,
                                          self.n_bands,
                                          self.n_atoms,
                                          len(self.orb_names))
                self.__orbitals = [self.__orbitals[0], self.__orbitals[1]]

    @property
    def phases(self):
        '''Phase data

        Note
        ------

        At present I have no idea about this parameter. How to use it?
        '''
        return self.__phases

    @phases.setter
    def phases(self, arg):  # not checked
        '''Phases

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
                self.__phases = phases.reshape(self.numk, self.n_bands,
                                               self.n_atoms,
                                               len(self.orb_names) - 1)
            elif len(self.spininfo) == 2:
                self.__phases = phases.reshape(2, self.numk, self.n_bands,
                                               self.n_atoms,
                                               len(self.orb_names) - 1)
                self.__phases = (self.__phases[0], self.__phases[1])

    @property
    def kvectors(self):
        '''setter for kvector'''
        return self.__kvectors

    @kvectors.setter
    def kvectors(self, kvectors):
        if not isinstance(kvectors, list):
            errmsg = 'kvectors must be an array of ndarray\n'
            raise TypeError(errmsg)
        self.__kvectors = kvectors
        self.numk = len(self.kvectors)
        self.kdistance = list()
        for i, k in enumerate(self.kvectors):
            if i == 0:
                self.kdistance.append(0.0)
            else:
                self.kdistance.append(
                    self.kdistance[i - 1] +
                    np.linalg.norm(self.kvectors[i - 1] - k))

    @property    # checked
    def energies(self):
        '''Energies'''
        return self.__energies

    @energies.setter
    def energies(self, arg):
        '''Setter for energies property.
        '''
        if self.isready():
            if len(self.spininfo) == 1 or len(self.spininfo) == 4:
                self.__energies = np.array(arg).reshape(
                    self.numk, self.n_bands)
            elif len(self.spininfo) == 2:
                energies = np.array(arg).T.reshape(2, self.numk, self.n_bands)
                print(energies)
                self.__energies = energies

    def fermi_correction(self, fermi):
        '''.. py:method:: fermi_correction(fermi)
        Correct the Fermi level

        Parameters
        -----------

        fermi: float
             value of the Fermi level (from OUTCAR or vasprun.xml).
        '''
        self.__energies -= fermi

    def sum_site(self, arg):
        '''Make sitecomposed ndarray

        When sitecomposed ndarray has elements, the values remain.

        Parameters
        -----------
        arg: list, tuple, set
             site indexes to be aded. it contains unique numbers.
        '''
        # the element of site_number_list must be unique.
        site_numbers = tuple(set(arg))
        self.isready()  # if not ready, raise Error.
        if len(self.spininfo) == 1:
            cmporbs = np.array([[[np.sum(
                [y for x, y in enumerate(self.orbitals[i, j])
                 if x in site_numbers],
                axis=0)] for j in range(len(self.available_band))]
                                for i in range(self.numk)])
            if self.sitecomposed:
                self.sitecomposed[0] = np.concatenate(
                    (self.sitecomposed[0], cmporbs),
                    axis=2)
            else:
                self.sitecomposed = [cmporbs]
        if len(self.spininfo) == 2:
            upspin_orbitals = self.orbitals[0]
            downspin_orbitals = self.orbitals[1]
            cmporbs_up = np.array([[[np.sum(
                [y for x, y in enumerate(upspin_orbitals[i, j])
                 if x in site_numbers],
                axis=0)] for j in range(len(self.available_band))]
                                   for i in range(self.numk)])
            cmporbs_down = np.array([[[np.sum(
                [y for x, y in enumerate(downspin_orbitals[i, j])
                 if x in site_numbers],
                axis=0)] for j in range(len(self.available_band))]
                                     for i in range(self.numk)])
            self.__orbitals[0] = np.concatenate((self.__orbitals[0],
                                                 cmporbs_up),
                                                axis=2)
            self.__orbitals[1] = np.concatenate((self.__orbitals[1],
                                                 cmporbs_down),
                                                axis=2)
            if self.sitecomposed:
                self.sitecomposed[0] = np.concatenate(
                    (self.sitecomposed[0], cmporbs_up),
                    axis=2)
                self.sitecomposed[1] = np.concatenate(
                    (self.sitecomposed[1], cmporbs_down),
                    axis=2)
            else:
                self.sitecomposed = [cmporbs_up, cmporbs_down]
        if len(self.spininfo) == 4:
            site_numbers_mtotal = tuple(x + self.n_atoms *
                                        0 for x in site_numbers)
            site_numbers_mx = tuple(x + self.n_atoms *
                                    1 for x in site_numbers)
            site_numbers_my = tuple(x + self.n_atoms *
                                    2 for x in site_numbers)
            site_numbers_mz = tuple(x + self.n_atoms *
                                    3 for x in site_numbers)
            #
            cmporbs_mtotal = np.array([[[np.sum(
                [y for x, y in enumerate(self.orbitals[i, j])
                 if x in site_numbers_mtotal],
                axis=0)] for j in range(len(self.available_band))]
                                       for i in range(self.numk)])
            cmporbs_mx = np.array([[[np.sum(
                [y for x, y in enumerate(self.orbitals[i, j])
                 if x in site_numbers_mx],
                axis=0)] for j in range(len(self.available_band))]
                                   for i in range(self.numk)])
            cmporbs_my = np.array([[[np.sum(
                [y for x, y in enumerate(self.orbitals[i, j])
                 if x in site_numbers_my],
                axis=0)] for j in range(len(self.available_band))]
                                   for i in range(self.numk)])
            cmporbs_mz = np.array([[[np.sum(
                [y for x, y in enumerate(self.orbitals[i, j])
                 if x in site_numbers_mz],
                axis=0)] for j in range(len(self.available_band))]
                                   for i in range(self.numk)])
            if self.sitecomposed:
                self.sitecomposed[0] = np.concatenate(
                    (self.sitecomposed[0], cmporbs_mtotal),
                    axis=2)
                self.sitecomposed[1] = np.concatenate(
                    (self.sitecomposed[1], cmporbs_mx),
                    axis=2)
                self.sitecomposed[2] = np.concatenate(
                    (self.sitecomposed[2], cmporbs_my),
                    axis=2)
                self.sitecomposed[3] = np.concatenate(
                    (self.sitecomposed[3], cmporbs_mz),
                    axis=2)
            else:
                self.sitecomposed = [cmporbs_mtotal,
                                     cmporbs_mx, cmporbs_my, cmporbs_mz]

    def get_orb_index(self, arg):
        '''.. py:method::get_orb_index(arg)

        Return the indexes corresponding orbitan name by tuple.

        This method returns the tuple of orbital number in self.orb_names.
        (i.e. self.orb_names.index(orbitalname).  If the orbital name has not
        been in self.orb_names (i.e. if the orbital name is not used in
        PROCAR file) but the orbital name is appropriate as the composed
        orbital ((ex.) sp, pxpy), returns the indexes of the orbitals to be
        composed as the tuple.


        Parameters
        ----------

        arg: str
            name of (composed) orbital

        Returns
        --------

        tuple
            number corresponding to (composed) orbital name.
        '''
        orbname = check_orb_name(arg)
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

    def sum_orbital(self, arg):
        '''.. py:method:compose_orbital(arg)

        Add composed orbital contribution in
        each 'sites' stored in BandStructure.sitecomposed.

        Parameters
        -----------

        arg: str, list, tuple
             orbital names
        '''
        if not self.sitecomposed:
            err = "This method operates with on sitecomposed attribute,"
            err += " but it's null"
            raise RuntimeError(err)
        if isinstance(arg, str):
            if ':' in arg:
                arg = arg.split(':')
            else:
                arg = [arg]
        for orb in arg:
            orb = check_orb_name(orb)
            if orb in self.orb_names:
                continue
            else:
                orbindex = self.get_orb_index(orb)
                # calculate composed orbital...
                for l in range(len(self.sitecomposed)):
                    numk, numband, numsite, norbs = self.sitecomposed[l].shape
                    orbitalcomposed = np.array([[[[np.sum(
                        [y for x, y in enumerate(self.sitecomposed[l][i, j, k])
                         if x in orbindex])] for k in range(numsite)]
                                                 for j in range(numband)]
                                                for i in range(numk)])
                    self.sitecomposed[l] = np.concatenate(
                        (self.sitecomposed[l], orbitalcomposed),
                        axis=3)
                self.orb_names.append(orb)

    def del_band(self, band_indexes):  # not checked
        '''not yet impremented'''
        if not self.sitecomposed:
            err = "This method operates with on sitecomposed attribute,"
            err += " but it's null"
            raise RuntimeError(err)

    def set_header(self, sitenames, orbnames):
        '''Return header of table

        Parameters
        -----------

        sitenames: tuple, list
            site names e.g., 'Ag', 'Graphene', '2ndLayer'...
        orbnames: tuple, list
            orbital names e.g., (('s', 'pxpy',), ('pxpy', 'p'))

        Returns
        --------

        list
             header
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
        '''.. py:method::get_orbnums(orbnames)

        Return tuple whose size is same as that of arg, but
        the element is number determied from orb_names

        Parameters
        -----------

        orbnames: list, tuple
            orbital names e.g., (('s','pxpy','tot'),('s','p'))

        Returns
        ---------

        tuple
        '''
        return tuple(tuple(self.orb_names.index(orb) for orb in orbs)
                     for orbs in orbnames)

    def list_sitecomposed_data(self, orbnames):  # not checked
        '''.. py:method::list_sitecomposed_data(orbnames)

        Return list of sitecomposed attribute to 2D-list

        Parameters
        ----------

        orbnames: list, tuple
            orbital names  e.g., (('s','pxpy','p'),('s','pz','p'))
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
            for bandindex in range(numband):
                for kindex in range(numk):
                    sitelist = list()
                    sitelist.append(self.kdistance[kindex])
                    sitelist.append(self.energies[0][kindex, bandindex])
                    for site, norbs in zip(list(range(numsite)),
                                           orbnums):
                        for orb in norbs:
                            sitelist.append(
                                self.sitecomposed[0][kindex,
                                                     bandindex,
                                                     site, orb])
                    sitelist.append(self.energies[1][kindex, bandindex])
                    for site, norbs in zip(list(range(numsite)),
                                           orbnums):
                        for orb in norbs:
                            sitelist.append(
                                self.sitecomposed[1][kindex,
                                                     bandindex,
                                                     site, orb])
                    table.append(sitelist)
                table.append([])
        else:
            raise RuntimeError("spininfo is incorrect")
        return table

    def get_sitecomposed_data(self, sitenames, orbnames):  # not checked
        '''Return band structure with (gathered) orbital
        contributions as string'''
        header = map(str, self.set_header(sitenames, orbnames))
        output = "\t".join(header) + "\n"
        lists = self.list_sitecomposed_data(orbnames)
        for line in lists:
            line = map(str, line)
            output += "\t".join(line) + "\n"
        return output


def check_orb_name(arg):
    '''.. py:function::check_orb_name(arg)

    Return arg without change if arg is a member of the 'orbital name'.
    i.e., if arg is an alias of the (more appropriate) orbital
    name, return it as is.  If arg is neither the appropriate
    orbital name nor the alias, raise ValueError.

    Parameters
    ----------

    arg: str
        the string to be checked as the orbital name

    Returns
    --------

    str
    '''
    translate_dict = {'pypx': 'pxpy', 'pzpx': 'pxpz', 'pzpy': 'pypz',
                      'pxpypz': 'p', 'pxpzpy': 'p', 'pypxpz': 'p',
                      'pypzpx': 'p', 'pzpxpy': 'p', 'pzpypx': 'p',
                      'spd': 'tot'}
    proper_orb_name_list = ['s', 'py', 'pz', 'px',
                            'dxy', 'dyz', 'dz2', 'dxz', 'dx2',
                            'tot'] + ['sp', 'p',
                                      'pxpy', 'pxpz', 'pypz', 'spd', 'd']
    if arg in translate_dict.keys():
        arg = translate_dict[arg]
    if arg in proper_orb_name_list:
        return arg
    else:
        errmsg = arg + ": (composed) orbital name was not defined."
        raise ValueError(errmsg)
