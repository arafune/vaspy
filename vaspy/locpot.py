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

    '''class for LOCPOT format.

LOCPOT format is essentially same as CHGCAR but simpler.
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
        '''
        :param locpotfile: LOCPOT file name
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
                    self.__potlist.extend(map(float, line.split()))
            self.__potarray = np.array(self.__potlist).reshape(
                self.meshZ, self.meshY, self.meshX)
            
    @property
    def meshX(self):
        return self.__meshX

    @property
    def meshY(self):
        return self.__meshY

    @property
    def meshZ(self):
        return self.__meshZ

    @property
    def potlist(self):
        return self.__potlist

    @property
    def potarray(self):
        return self.__potarray

    def get_mesh(self):
        return self.__meshX, self.__meshY, self.__meshZ 

    def average_along_axis(self, axis_name):
        axis_name = axis_name.capitalize()
        if axis_name == 'X':
            pot = np.sum(np.sum(
                np.transpose(self.potarray, (2,0,1)), axis=2), axis=1)
            pot /= (self.meshY * self.meshZ)
        if axis_name == 'Y':
            pot= np.sum(np.sum(
                np.transpose(self.potarray, (1,0,2)), axis=2), axis=1)
            pot /= (self.meshX * self.meshZ)
        if axis_name == 'Z':
            pot = np.sum(np.sum(self.potarray, axis=2), axis=1)
            pot /= (self.meshX * self.meshY)
        return pot

    def max_along_axis(self, axis_name):
        axis_name = axis_name.capitalize()
        if axis_name == 'X':
            pot = np.min(np.min(
                np.transpose(self.potarray, (2,0,1)), axis=2), axis=1)
        if axis_name == 'Y':
            pot= np.min(np.min(
                np.transpose(self.potarray, (1,0,2)), axis=2), axis=1)
        if axis_name == 'Z':
            pot = np.min(np.min(self.potarray, axis=2), axis=1)
        return pot

    def min_along_axis(self, axis_name):
        axis_name = axis_name.capitalize()
        if axis_name == 'X':
            pot = np.max(np.max(
                np.transpose(self.potarray, (2,0,1)), axis=2), axis=1)
        if axis_name == 'Y':
            pot= np.max(np.max(
                np.transpose(self.potarray, (1,0,2)), axis=2), axis=1)
        if axis_name == 'Z':
            pot = np.max(np.max(self.potarray, axis=2), axis=1)
        return pot

    def median_along_axis(self, axis_name):
        axis_name = axis_name.capitalize()
        if axis_name == 'X':
            pot = np.median(np.median(
                np.transpose(self.potarray, (2,0,1)), axis=2), axis=1)
        if axis_name == 'Y':
            pot= np.median(np.median(
                np.transpose(self.potarray, (1,0,2)), axis=2), axis=1)
        if axis_name == 'Z':
            pot = np.median(np.median(self.potarray, axis=2), axis=1)
        return pot


    def get_axes_lengthes(self):
        x = self.latticeV1 * self.scaling_factor
        y = self.latticeV2 * self.scaling_factor
        z = self.latticeV3 * self.scaling_factor
        x = np.linalg.norm(x)
        y = np.linalg.norm(y)
        z = np.linalg.norm(z)
        return (x, y, z)

    def plot_potential_along_axis(self, axis_name):
        axis_name = axis_name.capitalize()
        axesLength = self.get_axes_lengthes()
        if axis_name == 'X':
            X = np.linspace(0, axesLength[0], self.meshX)
        if axis_name == 'Y':
            X = np.linspace(0, axesLength[1], self.meshY)
        if axis_name == 'Z':
            X = np.linspace(0, axesLength[2], self.meshZ)
        Y1 = self.average_along_axis(axis_name)
        Y2 = self.max_along_axis(axis_name)
        Y3 = self.min_along_axis(axis_name)
        Y4 = self.median_along_axis(axis_name)
        plt.plot(X, Y1, Y2, Y3, Y4)

    
