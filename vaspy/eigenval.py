'''
This module provides EIGENVAL.
'''


from __future__ import print_function
from __future__ import division
import os
import bz2
import numpy as np
try:
    import matplotlib.pyplot as plt
except ImportError:
    print('Install matplotlib, or you cannot use methods relating to draw\n')


class EnergyBand(object):
    '''
    Simple band structure object for analyzing by using ipython.

    Attributes
    ----------

    kvecs: numpy.array
        kvectors
    kdistances: numpy.array
        kdisance
    numk: int
        number of kpoints
    nbands: int
        number of bands
    energies: numpy.array
        energies[spin_i, k_i, band_i], where spin_i, k_i, and band_i are spin-,
        k- and band-index, respectively.
    label: list of str
        used as a label (data 'title' such as '#k', 'Energy') in str format


    Parameters
    ----------

    kvecs: numpy.ndarray
         1D array data of k-vectors.
    energies: numpy.ndarray
         1D array data of energies
    spininfo: int, tuple, optional
         Spin type.  1 or ("",) mean No-spin.  2 or ('_up', '_down')
         mean collinear spin, 4 or ('_mT', '_mX', '_mY', '_mZ') mean
         collinear spin. This class does not distinguish non-collinear spin
         and No-spin.  (default is 1)
    '''

    def __init__(self, kvecs=(), energies=(), spininfo=1):
        self.kvecs = np.array(kvecs)
        self.numk = len(self.kvecs)
        try:
            self.nbands = len(energies) // len(kvecs)
        except ZeroDivisionError:
            self.nbands = 0
        self.energies = energies
        self.spininfo = spininfo
        if self.spininfo == 1:  # standard
            self.spininfo = ('',)
        elif self.spininfo == 2 or len(self.spininfo) == 2:   # collinear
            self.spininfo = ('_up', '_down')
        elif self.spininfo == 4:  # non-collinear
            self.spininfo = ('_mT', '_mX', '_mY', '_mZ')
        self.label = ['#k']

    @property
    def kdistances(self):
        '''Return kdistances'''
        return np.cumsum(np.linalg.norm(
            np.concatenate((np.array([[0, 0, 0]]),
                            np.diff(self.kvecs, axis=0))), axis=1))

    def __str__(self):  # << FIXME
        '''
        Returns
        --------

        str
            a string represntation of EnergyBand.
            Useful for gnuplot and Igor.
        '''
        output = self.label[0]
        for label in self.label[1]:
            output += '\t'+label
        output += '\n'
        for band_i in range(self.energies.shape[2]):
            for k, energies in zip(self.kdistances, self.energies[:, :, band_i].T):
                output += str(k)
                for energy in energies:
                    output += "\t{0:.8e}".format(energy)
                output += "\n"
            output += "\n"
        return output

    def figure(self, color='blue', spin_i=0):
        '''
        Return Axes object of the energy band.

        Parameters
        ----------

        color: str, optional (default is 'blue')
            color of the band line

        spin_i: spin_index
            default is 0

        Returns
        --------
        matplotlib.pyplot.Axes

        Example
        -------

        Here is a typical code::

           fig = plt.figure()
           ax = band.figure(color='blue')
           ax.set_ylabel('Energy  ( eV )')
           ax.set_ylim(-5, 5)
           ax.set_xlim(0, 4)
           plt.show()
        '''
        _ = [plt.plot(self.kdistances, self.energies[spin_i, :, band_i],
                      color=color)
             for band_i in range(self.energies.shape[2])]
        return plt.gca()

    def show(self, yrange=None, spin_i=0):  # How to set default value?
        '''
        Draw band structure by using maptlotlib.
        For 'just seeing' use.

        Parameters
        ----------

        yrange: tuple, optional  (default: all range)
             Minimum and maximum value of the y-axis.
             If not specified, use the matplotlib default value.

        spin_i: int  (default is 0 for no spin or 'up' spin)
             Spin index. For spin-polarized collinear band
        '''
        for band_i in range(self.energies.shape[2]):
            plt.plot(self.kdistances, self.energies[spin_i, :, band_i],
                     color='blue')
        if yrange is not None:
            plt.ylim([yrange[0], yrange[1]])
        plt.xlim([self.kdistances[0],
                  self.kdistances[-1]])
        plt.ylabel(self.label[spin_i+1] + ' (eV)')
        plt.show()


    def to_physical_kvector(self, recvec=((1.0, 0.0, 0.0),
                                          (0.0, 1.0, 0.0),
                                          (0.0, 0.0, 1.0))):
        '''Change kvec unit to inverse AA

        Parameters
        -----------
        recvec: array, numpy.ndarray, optional (default is the unit vector)
            reciprocal vector

            .. Note:: Don't forget that the reciprocal vector used
                      in VASP needs 2PI to match  the conventional
                      unit of the wavevector.
        '''
        recvec = np.array(recvec)
        self.kvecs = np.array(
            [recvec.dot(kvecs) for kvecs in self.kvecs])


class EIGENVAL(EnergyBand):
    '''
    Class for storing the data of EIGENVAL file.

    Parameters
    -----------

    filename: str
        File name of 'EIGENVAL'

    Attributes
    ----------
    natom: int
        Number of atoms
    '''

    def __init__(self, filename=None):
        super(EIGENVAL, self).__init__()
        self.natom = 0
        #
        if filename:
            if os.path.splitext(filename)[1] == '.bz2':
                try:
                    self.thefile = bz2.open(filename, mode='rt')
                except AttributeError:
                    self.thefile = bz2.BZ2File(filename, mode='r')
            else:
                self.thefile = open(filename)
            self.load_file()

    def load_file(self):
        '''
        A virtual parser of EIGENVAL
        '''
        self.natom, _, _, self.spininfo = [int(i) for i in
                                           next(self.thefile).split()]
        if self.spininfo == 2:
            self.label.extend(['Energy_up', 'Energy_down'])
        else:
            self.label.append('Energy')
        next(self.thefile)
        next(self.thefile)
        next(self.thefile)
        next(self.thefile)
        _, self.numk, self.nbands = [int(i) for i in
                                     next(self.thefile).split()]
        self.kvecs = []
        self.energies = []
        for _ in range(self.numk):
            # the first line in the sigleset begins with the blank
            next(self.thefile)
            self.kvecs.append(
                [float(i) for i in next(self.thefile).split()[0:3]])
            for _ in range(self.nbands):
                self.energies.append(
                    [float(i) for i in
                     next(self.thefile).split()[1:self.spininfo+1]])
        self.kvecs = np.array(self.kvecs)
        self.energies = np.array(
            self.energies).T.reshape(self.spininfo, self.numk, self.nbands)
        self.thefile.close()
