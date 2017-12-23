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
import numpy as np
import matplotlib.pyplot as plt
try:
    from vaspy import mesh3d
except ImportError:
    MYPATH = os.readlink(__file__) if os.path.islink(__file__) else __file__
    sys.path.append(os.path.dirname(os.path.abspath(MYPATH)))
    import mesh3d


class LOCPOT(mesh3d.VASPGrid):
    '''.. py:class:: LOCPOT(locpotfile)

    Class for LOCPOT

    LOCPOT format is essentially same as that of CHGCAR

'''

    def __init__(self, arg=None):
        super(LOCPOT, self).__init__(None)
        if arg:
            self.load_from_file(arg)

    def plot_potential_along_axis(self, axis_name, type=0):
        '''.. py:method:: plot_potential_along_axis(axis_name, pottype)

        Plot potential curve along the axis

        Parameters
        ----------

        axis_name: str
            the name of the axis (X, Y, or Z)
        type: int
            'select potential type' (very optional)
        '''
        axis_name = axis_name.capitalize()
        axes_length = self.get_axes_lengthes()
        if axis_name == 'X':
            horizontal_axis = np.linspace(0, axes_length[0], self.meshsize[0])
            plt.clf()
            plt.xlim(xmax=axes_length[0])
        if axis_name == 'Y':
            horizontal_axis = np.linspace(0, axes_length[1], self.meshsize[1])
            plt.clf()
            plt.xlim(xmax=axes_length[1])
        if axis_name == 'Z':
            horizontal_axis = np.linspace(0, axes_length[2], self.meshsize[2])
            plt.clf()
            plt.xlim(xmax=axes_length[2])
        y_average = self.average_along_axis(axis_name, type)
        y_max = self.max_along_axis(axis_name, type)
        y_min = self.min_along_axis(axis_name, type)
        y_median = self.median_along_axis(axis_name, type)
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
