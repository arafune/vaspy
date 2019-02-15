# -*- coding: utf-8 -*-
'''This module provide XDATCAR class
'''

import bz2
import os

import numpy as np

import vaspy.poscar


class XDATCAR(vaspy.poscar.POSCAR_HEAD):
    '''class for XDATCAR format

    Attributes
    ----------

    configurations: list
'''

    def __init__(self, filename=None):
        '''
        Parameters
        ----------

        arg: str
            XDATCAR file name
'''
        super(XDATCAR, self).__init__()
        self.configurations = []
        if filename:
            if os.path.splitext(filename)[1] == ".bz2":
                try:
                    thefile = bz2.open(filename, mode='rt')
                except AttributeError:
                    thefile = bz2.BZ2File(filename, mode='r')
            else:
                thefile = open(filename)
            self.load_file(thefile)

    def load_file(self, thefile):
        '''A virtual parser of PROCAR

        Parameters
        ----------

        thefile: StringIO
            'XDATCAR' file
'''
        self.system_name = next(thefile).strip()
        self.scaling_factor = float(next(thefile).strip())
        self.cell_vecs[0] = np.array([float(x) for x in next(thefile).split()])
        self.cell_vecs[1] = np.array([float(x) for x in next(thefile).split()])
        self.cell_vecs[2] = np.array([float(x) for x in next(thefile).split()])
        self.iontypes = next(thefile).split()
        self.ionnums = [int(x) for x in next(thefile).split()]
        positions = []
        for line in thefile:
            if 'Direct configuration=' in line:
                if positions:
                    self.configurations.append(positions)
                    positions = []
            else:
                position = np.array([float(x) for x in line.strip().split()])
                positions.append(position)
        self.configurations.append(positions)

    def __str__(self):
        '''
        Returns
        -------

        str
            a string representation of XDATCAR
'''
        tmp = self.system_name + '\n'
        tmp += '        {}\n'.format(self.scaling_factor)
        for i in range(3):
            tmp += '      {:#.6f}   {:#.6f}    {:6f}\n'.format(
                self.cell_vecs[i][0], self.cell_vecs[i][1],
                self.cell_vecs[i][2])
        for element in self.iontypes:
            tmp += '    {}'.format(element)
        tmp += '\n'
        for ionnum in self.ionnums:
            tmp += '    {}'.format(ionnum)
        tmp += '\n'
        for frame_index, positions in enumerate(self.configurations):
            tmp += 'Direct configuration=    {}\n'.format(frame_index + 1)
            for position in positions:
                tmp += '    {:#.6f}    {:#.6f}    {:6f}\n'.format(
                    position[0], position[1], position[2])
        return tmp
