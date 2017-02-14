'''.. py:module:: eigenmval.py

This module provides EIGENVAL.
'''


from __future__ import print_function
from __future__ import division
import os
import bz2
import numpy as np
import matplotlib.pyplot as plt


class EIGENVAL(object):
    '''.. py:class:: EIGENVAL(EIGENVAL_file)

    Class for storing the data of EIGENVAL file.

    Parameters
    -----------

    EIGENVAL_file: str
        File name of 'EIGENVAL'

    Attributes
    ----------
    n_atom: int
        Number of atoms
    numk: int
        Number of k vectors
    spininfo: int
        No_spin or non-collinear:1 collinear spin: 2
    kvectors: list
        List of kvalues, the length of kvectors must be equal to numk.
    energies: list
        Energy values
    '''

    def __init__(self, arg=None):
        self.n_atoms = 0
        self.numk = 0
        self.kvectors = list()
        self.n_bands = 0
        self.energies = list()
        self.spininfo = 0
        #
        if isinstance(arg, str):
            self.load_from_file(arg)

    def load_from_file(self, filename):
        '''..py:method:: load_from_file(filename)

        A virtual parser of EIGENVAL

        parameters
        -----------

        filename: str
           Filename of EIGENVAL file. bziped file is also readable.
        '''
        if os.path.splitext(filename)[1] == '.bz2':
            try:
                thefile = bz2.open(filename, mode='rt')
            except AttributeError:
                thefile = bz2.BZ2File(filename, mode='r')
        else:
            thefile = open(filename)
        with thefile:
            self.n_atoms, _, _, self.spininfo = [int(i) for i in
                                                 next(thefile).strip().split()]
            next(thefile)
            next(thefile)
            next(thefile)
            next(thefile)
            _, self.numk, self.n_bands = [int(i) for i in
                                          next(thefile).strip().split()]
            for ki in range(self.numk):
                # the first line in the sigleset begin with the blank line
                next(thefile)
                self.kvectors.append(np.array(
                    [float(i) for i in next(thefile).strip().split()[0:3]]))
                for bi in range(self.n_bands):
                    if self.spininfo == 1:
                        self.energies.append(float(
                            next(thefile).strip().split()[1]))
                    else:
                        self.energies.append(
                            np.array([float(i) for i in
                                      next(thefile).strip().split()[1:3]]))
        self.energies = np.array(self.energies)

    def onlyband(self, recvec=[[1.0, 0.0, 0.0],
                               [0.0, 1.0, 0.0],
                               [0.0, 0.0, 1.0]]):
        '''.. py:method:: onlyband(recvec)

        Return Band_with_projection object

        Parameters
        -----------

        recvec: array, numpy.ndarray
            reciprocal vector

            .. Note:: Don't forget that the reciprocal vector used
                      in VASP needs 2PI to match  the conventional
                      unit of the wavevector.

        Returns
        ---------

        EnergyBand
'''
        recvecarray = np.array(recvec).T
        physical_kvector = [recvecarray.dot(kvector) for kvector in
                            self.kvectors[0:self.numk]]
        return EnergyBand(physical_kvector, self.energies, self.spininfo)

    def __str__(self):
        '''..py:method:: __str__()

        __str__() <=> str(x)

        Show the EIGENVAL character, not contents itself.
        '''
        template = '''The parameter of EIGENVALUE
# of atoms: {0.n_atoms}
# of kpoints: {0.numk}
# of bands: {0.n_bands}
'''
        return template.format(self)


class EnergyBand(object):
    '''.. py:class:: EnergyBand(kvectors, energies, spininfo)

    Simple band structure object for analyzing by using ipython.

    Parameters
    -----------

    kvectors: numpy.ndarray
         1D array data of k-vectors.
    energies: numpy.ndarray
         1D array data of energies
    spininfo: int, tuple
         Spin type.  1 or ("",) means No-spin.  2 or ('_up', '_down')
         means collinear spin, 4 or ('_mT', '_mX', '_mY', '_mZ') means
         collinear spin. This class does not distinguish  non-collinear spin
         and No-spin.
'''

    def __init__(self, kvectors, energies, spininfo=1):
        self.kvectors = np.array(kvectors)
        self.kdistances = np.cumsum(
            np.linalg.norm(
                np.concatenate(
                    (np.array([[0, 0, 0]]),
                     np.diff(kvectors, axis=0))), axis=1))
        self.numk = len(self.kvectors)
        self.nbands = len(energies) // len(kvectors)
        self.spininfo = spininfo
        if self.spininfo == 1:  # standard
            self.spininfo = ('',)
        elif self.spininfo == 2 or len(self.spininfo) == 2:   # collinear
            self.spininfo = ('_up', '_down')
            self.nbands = self.nbands // 2
        elif self.spininfo == 4:  # non-collinear
            self.spininfo = ('_mT', '_mX', '_mY', '_mZ')
        if spininfo == 2 or spininfo == ('_up', '_down'):
            self.energies = np.array(energies).reshape(
                (2, self.numk, self.nbands))
        else:
            self.energies = np.array(energies).reshape(
                (self.numk, self.nbands))

    def fermi_correction(self, fermi):
        '''.. py:method:: fermi_correction(fermi)

        Correct the Fermi level

        Parameters
        ----------

        fermi: float
             value of the Fermi level.
'''
        self.energies -= fermi

    def __str__(self):
        '''.. py:method:: __str__()

        Returns
        --------

        str
            a string represntation of EnergyBand.  Useful for gnuplot and Igor.
        '''
        if self.spininfo == 2 or len(self.spininfo) == 2:
            output = '#k\tEnergy_up\tEnergy_down\n'
            for k_i in range(self.numk):
                for k, up, down in zip(self.kdistances,
                                       self.energies[0][k_i],
                                       self.energies[1][k_i]):
                    output += '{0:.9e}\t{1:.9e}\t{2:.9e}\n'.format(k, up, down)
                output += '\n'
        else:
            output = '#k\tEnergy\n'
            for k_i in range(self.numk):
                for k, en in zip(self.kdistances, self.energies[k_i]):
                    output += '{0:.9e}\t{1:.9e}\n'.format(k, en)
                output += '\n'
        return output

    def showband(self, yrange=None, spin=None):  # How to set default value?
        '''.. py:method:: showband(yrange)

        Draw band structure by using maptlotlib.
        For 'just seeing' use.

        Parameters
        ----------

        yrange: tuple
             Minimum and maximum value of the y-axis.
             If not specified, use the matplotlib default value.
'''
        if self.spininfo == 2 and spin == 'up':
            energies = np.swapaxes(self.energies[0], 1, 0)
        elif self.spininfo == 2 and spin == 'down':
            energies = np.swapaxes(self.energies[1], 1, 0)
        else:
            energies = np.swapaxes(self.energies, 1, 0)
        for i in range(0, energies.shape[0]):
            plt.plot(self.kdistances,
                     energies[i],
                     color='blue')
        if yrange is not None:
            plt.ylim([yrange[0], yrange[1]])
        plt.xlim([self.kdistances[0],
                  self.kdistances[-1]])
        plt.ylabel("Energy (eV)")
        plt.show()
