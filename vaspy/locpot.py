#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
from __future__ import division
from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
try:
    from vaspy import poscar
    from vaspy import tools
except ImportError:
    mypath = os.readlink(__file__) if os.path.islink(__file__) else __file__
    sys.path.append(os.path.dirname(os.path.abspath(mypath)))
    import poscar
    import tools


class LOCPOT(poscar.POSCAR):

    '''Class for LOCPOT

    LOCPOT format is essentially same as that of CHGCAR

    :attribute: meshX, meshY, meshZ, potlist
'''

    def __init__(self, arg=None):
        super(LOCPOT, self).__init__(None)
        self.__meshX = 0
        self.__meshY = 0
        self.__meshZ = 0
        self.__potlist = []
        if arg:
            self.load_from_file(arg)

    def load_from_file(self, locpotfile):
        '''Parse LOCPOT file to construct LOCPOT object

        :param locpotfile: file name of 'LOCPOT'
        :type locpotfile: str
        :return: LOCPOT
        '''
        section = 'poscar'
        separator = None
        tmp = []
        with open(locpotfile) as f:
            for line in f:
                if section == 'poscar':
                    if line.isspace():
                        self.load_from_array(tmp)
                        section = 'define_separator'
                    else:
                        line = line.rstrip('\n')
                        tmp.append(line)
                elif section == 'define_separator':
                    separator = line if separator is None else separator
                    if self.meshX == self.meshY == self.meshZ == 0:
                        self.__meshX, self.__meshY, self.__meshZ = \
                            list(map(int, line.split()))
                        section = 'grid'
                elif section == 'grid':
                    line = line.rstrip('\n')
                    self.__potlist.extend(map(np.float64, line.split()))
            if len(self.__potlist) == self.meshX * self.meshY * self.meshZ:
                self.__potarray = np.array(self.__potlist).reshape(
                    self.meshZ, self.meshY, self.meshX)
            elif len(self.__potlist) == 2 * (self.meshX *
                                             self.meshY * self.meshZ) + 3 + sum(self.ionnums):  # LVHAR
                self.__potarray = np.array(
                    self.__potlist[:
                                   self.meshX *
                                   self.meshY *
                                   self.meshZ]).reshape(self.meshZ,
                                                        self.meshY,
                                                        self.meshX)
                self.__potarray2 = np.array(
                    self.__potlist[self.meshX *
                                   self.meshY *
                                   self.meshZ:
                                   self.meshX *
                                   self.meshY *
                                   self.meshZ + sum(self.ionnums)])
                self.__potarray3 = np.array(
                    self.__potlist[- self.meshX *
                                   self.meshY *
                                   self.meshZ:]).reshape(self.meshZ,
                                                         self.meshY,
                                                         self.meshX)

    @property
    def meshX(self):
        '''Number of mesh along the first axis of the cell'''
        return self.__meshX

    @property
    def meshY(self):
        '''Number of mesh along the second axis of the cell'''
        return self.__meshY

    @property
    def meshZ(self):
        '''Number of mesh along the second axis of the cell'''
        return self.__meshZ

    @property
    def potlist(self):
        '''Potential data (Array style)

        From Vasp Manual, the potential is written using the following commands in Fortran (as same as for charge density in CHGCAR)

        .. code-block :: fortran

           WRITE(IU, FORM) (((C(NX, NY, NZ), NX=1, NGXC), NY=1, NGYZ), NZ=1, NGZC)
'''
        return self.__potlist

    @property
    def potarray(self):
        '''Potential data (Array style)

        This data should be used until its meaning is clear.  I (RA)
        have not understand the meaning of the potarray 2 and 3...

        '''
        return self.__potarray

    @property  # <= Fixme!!!
    def potarray3(self):
        '''(maybe) Potential data (Array style)

        This data should NOT be used until its meaning is clear.  I
        (RA) have not understand the meaning of the potarray 2 and
        3...

        .. warning:: Do not use this attribute if you do not know what you treated

        '''
        return self.__potarray2

    @property  # <= Fixme!!!
    def potarray2(self):
        '''(maybe) Potential data (Array style)

        Usually this data should NOT be used.  I (RA) have not understand the meaning of the potarray 2 and 3...

        .. warning:: Do not use this attribute if you do not know what you treated
'''
        return self.__potarray3

    def get_mesh(self):
        '''Return mesh size

        :return: (meshX, meshY, meshZ)
        :rtype: tuple 
        '''
        return self.__meshX, self.__meshY, self.__meshZ

    def average_along_axis(self, axis_name, pottype='former'):
        '''Calculate average value of potential along 'axis'

        :param axis_name: 'X', 'Y', or 'Z'
        :type axis_name: str
        :return: averaged potential
        :rtype: 
        '''
        axis_name = axis_name.capitalize()
        if pottype == 'former':
            if axis_name == 'X':
                pot = np.average(np.average(
                    np.transpose(self.potarray, (2, 0, 1)), axis=2), axis=1)
            if axis_name == 'Y':
                pot = np.average(np.average(
                    np.transpose(self.potarray, (1, 0, 2)), axis=2), axis=1)
            if axis_name == 'Z':
                pot = np.average(np.average(self.potarray, axis=2), axis=1)
        elif pottype == 'latter':
            if axis_name == 'X':  # <= Fixme!!!
                pot = np.average(np.average(
                    np.transpose(self.potarray2, (2, 0, 1)), axis=2), axis=1)
            if axis_name == 'Y':  # <= Fixme!!!
                pot = np.average(np.average(
                    np.transpose(self.potarray2, (1, 0, 2)), axis=2), axis=1)
            if axis_name == 'Z':  # <= Fixme!!!
                pot = np.average(np.average(self.potarray2, axis=2), axis=1)
        else:
            raise RuntimeError("invalid pottype option")
        return pot

    def min_along_axis(self, axis_name, pottype='former'):
        '''Calculate minimum value of potential along 'axis'

        :param axis_name: 'X', 'Y', or 'Z'
        :type axis_name: str
        :return: min potential
        :rtype: 
        '''
        axis_name = axis_name.capitalize()
        if pottype == 'former':
            if axis_name == 'X':
                pot = np.min(np.min(
                    np.transpose(self.potarray, (2, 0, 1)), axis=2), axis=1)
            if axis_name == 'Y':
                pot = np.min(np.min(
                    np.transpose(self.potarray, (1, 0, 2)), axis=2), axis=1)
            if axis_name == 'Z':
                pot = np.min(np.min(self.potarray, axis=2), axis=1)
        elif pottype == 'latter':
            if axis_name == 'X':  # <= Fixme!!!
                pot = np.min(np.min(
                    np.transpose(self.potarray2, (2, 0, 1)), axis=2), axis=1)
            if axis_name == 'Y':  # <= Fixme!!!
                pot = np.min(np.min(
                    np.transpose(self.potarray2, (1, 0, 2)), axis=2), axis=1)
            if axis_name == 'Z':  # <= Fixme!!!
                pot = np.min(np.min(self.potarray2, axis=2), axis=1)
        else:
            raise RuntimeError("invalid pottype option")
        return pot

    def max_along_axis(self, axis_name, pottype='former'):
        '''Calculate maximum value of potential along 'axis'

        :param axis_name: 'X', 'Y', or 'Z'
        :type axis_name: str
        :return: maximum potential
        :rtype: np.ndarray
        '''
        axis_name = axis_name.capitalize()
        if pottype == 'former':
            if axis_name == 'X':
                pot = np.max(np.max(
                    np.transpose(self.potarray, (2, 0, 1)), axis=2), axis=1)
            if axis_name == 'Y':
                pot = np.max(np.max(
                    np.transpose(self.potarray, (1, 0, 2)), axis=2), axis=1)
            if axis_name == 'Z':
                pot = np.max(np.max(self.potarray, axis=2), axis=1)
        elif pottype == 'latter':
            if axis_name == 'X':  # <= Fixme!!!
                pot = np.max(np.max(
                    np.transpose(self.potarray2, (2, 0, 1)), axis=2), axis=1)
            if axis_name == 'Y':  # <= Fixme!!!
                pot = np.max(np.max(
                    np.transpose(self.potarray2, (1, 0, 2)), axis=2), axis=1)
            if axis_name == 'Z':  # <= Fixme!!!
                pot = np.max(np.max(self.potarray2, axis=2), axis=1)
        else:
            raise RuntimeError("invalid pottype option")
        return pot

    def median_along_axis(self, axis_name, pottype='former'):
        '''Calculate median value of potential along 'axis'

        :param axis_name: 'X', 'Y', or 'Z'
        :type axis_name: str
        :return: median potential
        :rtype: np.ndarray
        '''
        axis_name = axis_name.capitalize()
        if pottype == 'former':
            if axis_name == 'X':
                pot = np.median(np.median(
                    np.transpose(self.potarray, (2, 0, 1)), axis=2), axis=1)
            if axis_name == 'Y':
                pot = np.median(np.median(
                    np.transpose(self.potarray, (1, 0, 2)), axis=2), axis=1)
            if axis_name == 'Z':
                pot = np.median(np.median(self.potarray, axis=2), axis=1)
        elif pottype == 'latter':
            if axis_name == 'X':  # <= Fixme!!!
                pot = np.median(np.median(
                    np.transpose(self.potarray2, (2, 0, 1)), axis=2), axis=1)
            if axis_name == 'Y':  # <= Fixme!!!
                pot = np.median(np.median(
                    np.transpose(self.potarray2, (1, 0, 2)), axis=2), axis=1)
            if axis_name == 'Z':  # <= Fixme!!!
                pot = np.median(np.median(self.potarray2, axis=2), axis=1)
        else:
            raise RuntimeError("invalid pottype option")
        return pot

    def get_axes_lengthes(self):
        '''Return cell axis lengthes

        :return: cell axis length of x, y, and z
        :rtype: tuple
        '''
        x = self.latticeV1 * self.scaling_factor
        y = self.latticeV2 * self.scaling_factor
        z = self.latticeV3 * self.scaling_factor
        x = np.linalg.norm(x)
        y = np.linalg.norm(y)
        z = np.linalg.norm(z)
        return (x, y, z)

    def plot_potential_along_axis(self, axis_name, pottype='former'):
        '''Plot potential curve along the axis

        :param axis_name: the name of the axis (X, Y, or Z)
        :type axis_name: str
        '''
        axis_name = axis_name.capitalize()
        axesLength = self.get_axes_lengthes()
        if axis_name == 'X':
            X = np.linspace(0, axesLength[0], self.meshX)
            plt.clf()
            plt.xlim(xmax=axesLength[0])
        if axis_name == 'Y':
            X = np.linspace(0, axesLength[1], self.meshY)
            plt.clf()
            plt.xlim(xmax=axesLength[1])
        if axis_name == 'Z':
            X = np.linspace(0, axesLength[2], self.meshZ)
            plt.clf()
            plt.xlim(xmax=axesLength[2])
        Y1 = self.average_along_axis(axis_name, pottype)
        Y2 = self.max_along_axis(axis_name, pottype)
        Y3 = self.min_along_axis(axis_name, pottype)
        Y4 = self.median_along_axis(axis_name, pottype)
        plt.plot(X, Y1, label='average')
        plt.plot(X, Y2, label='max')
        plt.plot(X, Y3, label='min')
        plt.plot(X, Y4, label='median')
        xlabel = "Position along " + axis_name + "-axis (A)"
        title = "Local potential (" + axis_name + ")"
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel("Energy  ( eV )")
        plt.show()

# LVHAR-tag:  This tag is available in VASP.5.2.12 and newer version.
# It determines whether the total local potential (file LOCPOT ) contains
# the entire local potential (ionic plus Hartree plus exchange correlation)
# or the electrostatic contributions only (ionic plus Hartree). Note that
# in VASP.5.2.12, the default is to write the entire local potential,
# including the exchange correlation potential.
#
# Memo: I (RA) don't know why the potential data are  stored twice
# (exactly same data), when LVHAR = .TRUE. and LVTOT is not set.
