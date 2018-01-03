# -*- coding: utf-8 -*-
'''
Module for bader class

Bader charge anaslysis performed by bader
http://theory.cm.utexas.edu/henkelman/code/bader/

This program is not a part of vasp but it deeply connects with the vasp.
'''

from __future__ import division  # Version safety
from __future__ import print_function
import os.path
import bz2
import numpy as np


class BaderACF(object):
    '''class for storing bader charge analysis (ACF.dat)

    ACF.dat contains the coordinates of each atom, the charge
    associated with it according to Bader partitioning, percentage of
    the whole according to Bader partitioning and the minimum distance
    to the surface. This distance should be compared to maximum
    cut-off radius for the core region if pseudo potentials have been
    used.


    Properties
    ----------
    natom: int
       Number of atoms
    positions: list
       Atom positions in AA
    charges: list
       charges
    vols: list
       volumes
    vaccharge: float
       vacuum charge
    vacvol: float
       vacuum volume
    nelectron: float
       number of electron
    '''
    def __init__(self, filename=None):
        self.natom = 0
        self.positions = []
        self.charges = []
        self.mindists = []
        self.vols = []
        self.vaccharge = 0
        self.vacvol = 0
        self.nelectron = 0
        if filename:
            if os.path.splitext(filename)[1] == '.bz2':
                try:
                    self.baderfile = bz2.open(filename, mode='rt')
                except AttributeError:
                    self.baderfile = bz2.BZ2File(filename, mode='r')
            else:
                self.baderfile = open(filename)
            self.parse()

    def parse(self):
        '''Parse AFC.dat'''
        # the first line is like:
        #     X     Y     Z     CHARGE      MIN DIST   ATOMIC VOL
        _ = next(self.baderfile)
        # the 2nd line is just "----------"
        _ = next(self.baderfile)
        for line in self.baderfile:
            if '------------------------------------' in line:
                break
            else:
                tmp = line.split()
                self.positions.append([float(pos) for pos in tmp[1:4]])
                self.charges.append(float(tmp[4]))
                self.mindists.append(float(tmp[5]))
                self.vols.append(float(tmp[6]))
        self.vaccharge = float(next(self.baderfile).split()[-1])
        self.vacvol = float(next(self.baderfile).split()[-1])
        self.nelectron = float(next(self.baderfile).split()[-1])
        self.natom = len(self.positions)
        self.baderfile.close()
