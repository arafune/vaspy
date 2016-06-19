#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
This module provides LOCPOT class
'''
#
from __future__ import division
from __future__ import print_function
import os
import sys
import bz2
import numpy as np
import matplotlib.pyplot as plt
try:
    from vaspy import poscar
    from vaspy import tools
except ImportError:
    MYPATH = os.readlink(__file__) if os.path.islink(__file__) else __file__
    sys.path.append(os.path.dirname(os.path.abspath(MYPATH)))
    import vaspy.poscar
    import vaspy.tools


class LOCPOT(poscar.POSCAR):
    '''.. py:class:: LOCPOT(locpotfile)

    Class for LOCPOT

    LOCPOT format is essentially same as that of CHGCAR

    Attributes
    -----------

    mesh_x, mesh_y, mesh_z, potlist
'''

    def __init__(self, arg=None):
        super(LOCPOT, self).__init__(None)
        self.mesh_x = 0
        self.mesh_y = 0
        self.mesh_z = 0
        self.potlist = []
        if arg:
            self.load_from_file(arg)

    def load_from_file(self, locpotfile):
        '''.. py:method:: load_from_file(locpotfile)

        Parse LOCPOT file to construct LOCPOT object

        Parameters
        ----------
        locpotfile: str
             file name of 'LOCPOT'

        Returns
        --------

        LOCPOT
        '''
        section = 'poscar'
        separator = None
        tmp = []
        if os.path.splitext(locpotfile)[1] == '.bz2':
            locpot_file = bz2.BZ2File(locpotfile, mode='r')
        else:
            locpot_file = open(locpotfile)
        with locpot_file:
            for line in locpot_file:
                if section == 'poscar':
                    if line.isspace():
                        self.load_from_array(tmp)
                        section = 'define_separator'
                    else:
                        line = line.rstrip('\n')
                        tmp.append(line)
                elif section == 'define_separator':
                    separator = line if separator is None else separator
                    if self.mesh_x == self.mesh_y == self.mesh_z == 0:
                        self.mesh_x, self.mesh_y, self.mesh_z = \
                                [int(i) for i in line.split()]
#                            list(map(int, line.split()))
                        section = 'grid'
                elif section == 'grid':
                    line = line.rstrip('\n')
                    self.potlist.extend([np.float64(i) for i in
                                         line.split()])
#                        map(np.float64, line.split()))
            if len(self.potlist) == self.mesh_x * self.mesh_y * self.mesh_z:
                self.potarray = np.array(self.potlist).reshape(
                    self.mesh_z, self.mesh_y, self.mesh_x)
            elif len(self.potlist) == 2 * (self.mesh_x *
                                           self.mesh_y *
                                           self.mesh_z) + 3 + sum(self.ionnums):
                self.potarray = np.array(
                    self.potlist[:
                                 self.mesh_x *
                                 self.mesh_y *
                                 self.mesh_z]).reshape(self.mesh_z,
                                                       self.mesh_y,
                                                       self.mesh_x)
                self.potarray2 = np.array(
                    self.potlist[self.mesh_x *
                                 self.mesh_y *
                                 self.mesh_z:
                                 self.mesh_x *
                                 self.mesh_y *
                                 self.mesh_z + sum(self.ionnums)])
                self.potarray3 = np.array(
                    self.potlist[- self.mesh_x *
                                 self.mesh_y *
                                 self.mesh_z:]).reshape(self.mesh_z,
                                                        self.mesh_y,
                                                        self.mesh_x)

    def get_mesh(self):
        '''.. py:method:: get_mesh()

        Return mesh size

        Returns
        -------
        tuple

            (mesh_x, mesh_y, mesh_z)
        '''
        return self.mesh_x, self.mesh_y, self.mesh_z

    def average_along_axis(self, axis_name, pottype='former'):
        '''.. py:method:: average_along_axis(axis_name, pottype)

        Calculate average value of potential along 'axis'

        Parameters
        ----------

        axis_name: str
             'X', 'Y', or 'Z'
        pottype: str
             'pottype' (very optional)

        Returns
        -------
        numpy.ndarray
            averaged potential
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
        '''.. py:method:: min_along_axis(axis_name, pottype)

        Calculate minimum value of potential along 'axis'

        Parameters
        -----------

        axis_name: str
            'X', 'Y', or 'Z'
        pottype: str
            'pottype' (very optional)

        Returns
        -------
        numpy.ndarray
            min potential
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
        '''.. py:method:: man_along_axis(axis_name, pottype)

        Calculate maximum value of potential along 'axis'

        Parameters
        -----------

        axis_name: str
            'X', 'Y', or 'Z'
        pottype: str
            'pottype' (very optional)

        Returns
        -------
        numpy.ndarray
            max potential
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
        '''.. py:method:: median_along_axis(axis_name, pottype)

        Calculate median value of potential along 'axis'

        Parameters
        -----------

        axis_name: str
            'X', 'Y', or 'Z'
        pottype: str
            'pottype' (very optional)

        Returns
        -------
        numpy.ndarray
            median potential
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
        '''.. py:method:: get_axes_lengthes()

        Return cell axis lengthes

        Returns
        -------

        tuple
            cell axis length of x, y, and z
        '''
        cell_x = np.linalg.norm(self.cell_vecs[0] * self.scaling_factor)
        cell_y = np.linalg.norm(self.cell_vecs[1] * self.scaling_factor)
        cell_z = np.linalg.norm(self.cell_vecs[2] * self.scaling_factor)
        return (cell_x, cell_y, cell_z)

    def plot_potential_along_axis(self, axis_name, pottype='former'):
        '''.. py:method:: plot_potential_along_axis(axis_name, pottype)

        Plot potential curve along the axis

        Parameters
        ----------

        axis_name: str
            the name of the axis (X, Y, or Z)
        potential: str
            'pottype' (very optional)
        '''
        axis_name = axis_name.capitalize()
        axes_length = self.get_axes_lengthes()
        if axis_name == 'X':
            horizontal_axis = np.linspace(0, axes_length[0], self.mesh_x)
            plt.clf()
            plt.xlim(xmax=axes_length[0])
        if axis_name == 'Y':
            horizontal_axis = np.linspace(0, axes_length[1], self.mesh_y)
            plt.clf()
            plt.xlim(xmax=axes_length[1])
        if axis_name == 'Z':
            horizontal_axis = np.linspace(0, axes_length[2], self.mesh_z)
            plt.clf()
            plt.xlim(xmax=axes_length[2])
        y_average = self.average_along_axis(axis_name, pottype)
        y_max = self.max_along_axis(axis_name, pottype)
        y_min = self.min_along_axis(axis_name, pottype)
        y_median = self.median_along_axis(axis_name, pottype)
        plt.plot(horizontal_axis, y_average, label='average')
        plt.plot(horizontal_axis, y_max, label='max')
        plt.plot(horizontal_axis, y_min, label='min')
        plt.plot(horizontal_axis, y_median, label='median')
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
