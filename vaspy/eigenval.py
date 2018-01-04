'''
This module provides EIGENVAL.
'''


from __future__ import print_function
from __future__ import division
import os
import sys
import bz2
import numpy as np
try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.stderr.write(
        'Install matplotlib, or you cannot use methods relating to draw')


class EIGENVAL(object):
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
    numk: int
        Number of k vectors
    nbands: int
        Number of bands
    spininfo: int
        No_spin or non-collinear:1 collinear spin: 2
    kvecs[ki]: list or numpy array
        List of kvalues, the length of kvectors must be equal to numk.
    energies[bi+ki*nbands]: list or numpy.array
        Energy values (two-value array for spin-polarized eigenvalu)
    '''

    def __init__(self, filename=None):
        self.natom = 0
        self.numk = 0
        self.kvecs = list()
        self.nbands = 0
        self.energies = list()
        self.spininfo = 0
        #
        if filename:
            if os.path.splitext(filename)[1] == '.bz2':
                try:
                    thefile = bz2.open(filename, mode='rt')
                except AttributeError:
                    thefile = bz2.BZ2File(filename, mode='r')
            else:
                thefile = open(filename)
            self.load_file(thefile)

    def load_file(self, thefile):
        '''
        A virtual parser of EIGENVAL

        thefile: StringIO
            'EIGENVAL' file
        '''
        self.natom, _, _, self.spininfo = [int(i) for i in
                                           next(thefile).split()]
        next(thefile)
        next(thefile)
        next(thefile)
        next(thefile)
        _, self.numk, self.nbands = [int(i) for i in
                                     next(thefile).split()]
        for _ in range(self.numk):
            # the first line in the sigleset begins with the blank
            next(thefile)
            self.kvecs.append(np.array(
                [float(i) for i in next(thefile).split()[0:3]]))
            for _ in range(self.nbands):
                if self.spininfo == 1:
                    self.energies.append(float(
                        next(thefile).split()[1]))
                else:
                    self.energies.append(
                        np.array([float(i) for i in
                                  next(thefile).split()[1:3]]))
        self.energies = np.array(self.energies)
        thefile.close()

    def to_band(self, recvec=((1.0, 0.0, 0.0),
                              (0.0, 1.0, 0.0),
                              (0.0, 0.0, 1.0))):
        '''
        Return EnergyBand object

        Parameters
        -----------

        recvec: array, numpy.ndarray, optional (default is the unit vector)
            reciprocal vector

            .. Note:: Don't forget that the reciprocal vector used
                      in VASP needs 2PI to match  the conventional
                      unit of the wavevector.

        Returns
        ---------

        vaspy.eigenval.EnergyBand
        '''
        recvecarray = np.array(recvec).T
        kvector_physical = [recvecarray.dot(kvector) for kvector in
                            self.kvecs[0:self.numk]]
        return EnergyBand(kvector_physical, self.energies, self.spininfo)

    def __str__(self):
        '''..py:method:: __str__()

        __str__() <=> str(x)

        Show the EIGENVAL character, not contents itself.
        '''
        template = '''The parameter of EIGENVALUE
# of atoms: {0.natom}
# of kpoints: {0.numk}
# of bands: {0.nbands}
'''
        return template.format(self)


class EnergyBand(object):
    '''
    Simple band structure object for analyzing by using ipython.

    Attributes
    ----------

    kvecs: numpy array
        kvectors
    kdistances: nparray
        kdisance
    numk: int   # should be __numk ?
        number of kpoints
    nbands: int  # 
        number of bands

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

    def __init__(self, kvecs, energies, spininfo=1):
        self.kvecs = np.array(kvecs)
        self.kdistances = np.cumsum(
            np.linalg.norm(
                np.concatenate(
                    (np.array([[0, 0, 0]]),
                     np.diff(kvecs, axis=0))), axis=1))
        self.numk = len(self.kvecs)
        self.nbands = len(energies) // len(kvecs)
        self.spininfo = spininfo
        if self.spininfo == 1:  # standard
            self.spininfo = ('',)
        elif self.spininfo == 2 or len(self.spininfo) == 2:   # collinear
            self.spininfo = ('_up', '_down')
        elif self.spininfo == 4:  # non-collinear
            self.spininfo = ('_mT', '_mX', '_mY', '_mZ')
        if spininfo == 2 or spininfo == ('_up', '_down'):
            self.energies = np.array(energies).reshape(
                (self.numk, self.nbands, 2))
        else:
            self.energies = np.array(energies).reshape(
                (self.numk, self.nbands))

    def fermi_correction(self, fermi):
        '''
        Correct the Fermi level

        Parameters
        ----------

        fermi: float
             value of the Fermi level.
        '''
        self.energies -= fermi

    def __str__(self):
        '''
        Returns
        --------

        str
            a string represntation of EnergyBand.
            Useful for gnuplot and Igor.
        '''
        energies = np.swapaxes(self.energies, 1, 0)
        if self.spininfo == 2 or len(self.spininfo) == 2:
            output = '#k\tEnergy_up\tEnergy_down\n'
            for band_i in range(self.nbands):
                for k_i, energy in zip(self.kdistances, energies[band_i]):
                    output += '{0:.9e}\t{1:.9e}\t{2:.9e}\n'.format(k_i,
                                                                   energy[0],
                                                                   energy[1])
                output += '\n'
        else:
            output = '#k\tEnergy\n'
            for band_i in range(self.nbands):
                for k_i, enenergy in zip(self.kdistances, energies[band_i]):
                    output += '{0:.9e}\t{1:.9e}\n'.format(k_i, enenergy)
                output += '\n'
        return output

    def figure(self, color='blue', spin=None):
        '''
        Return Axes object of the energy band.

        Parameters
        ----------

        color: str, optional (default is 'blue')
            color of the band line

        spin: str
            up or down

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
        energies = np.swapaxes(self.energies, 1, 0)
        draws = []
        if self.spininfo == 2 or len(self.spininfo) == 2:
            if spin == 'down':
                for band_i in range(0, self.nbands):
                    draws.append(
                        plt.plot(self.kdistances, energies[band_i].T[1],
                                 color=color))
            else:
                for band_i in range(0, self.nbands):
                    draws.append(
                        plt.plot(self.kdistances, energies[band_i].T[0],
                                 color=color))
        else:
            for band_i in range(0, self.nbands):
                draws.append(
                    plt.plot(self.kdistances, energies[band_i],
                             color=color))
        return plt.gca()

    def show(self, yrange=None, spin=None):  # How to set default value?
        '''
        Draw band structure by using maptlotlib.
        For 'just seeing' use.

        Parameters
        ----------

        yrange: tuple, optional  (default: all range)
             Minimum and maximum value of the y-axis.
             If not specified, use the matplotlib default value.

        spin: str  (default is no spin or 'up' spin)
             Spin direction for spin-polarized collinear band
        '''
        energies = np.swapaxes(self.energies, 1, 0)
        if self.spininfo == 2 or len(self.spininfo) == 2:
            if spin == 'down':
                for band_i in range(0, self.nbands):
                    plt.plot(self.kdistances, energies[band_i].T[1],
                             color='blue')
            else:
                for band_i in range(0, self.nbands):
                    plt.plot(self.kdistances, energies[band_i].T[0],
                             color='blue')
        else:
            for band_i in range(0, self.nbands):
                plt.plot(self.kdistances, energies[band_i],
                         color='blue')
        if yrange is not None:
            plt.ylim([yrange[0], yrange[1]])
        plt.xlim([self.kdistances[0],
                  self.kdistances[-1]])
        plt.ylabel("Energy (eV)")
        plt.show()
