'''.. py:module:: eigenmval.py

This module provides EIGENVAL.
'''


from __future__ import print_function
from __future__ import division
import os
import bz2
import numpy as np
import vaspy.procar as procar # for EnergyBand class


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
        self.spininfo = 0  #          
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
            self.n_atoms, _, _, self.spininfo = [int(i) for i in next(thefile).strip().split()]
            next(thefile)
            next(thefile)
            next(thefile)
            next(thefile)
            _, self.numk , self.n_bands = [int(i) for i in next(thefile).strip().split()]
            for ki in range(self.numk):
                next(thefile)  # the first line in the sigleset begin with the blank line
                self.kvectors.append(np.array([float(i) for
                                              i in next(thefile).strip().split()[0:3]]))
                for bi in range(self.n_bands):
                    if self.spininfo == 1:
                        self.energies.append(float(next(thefile).strip().split()[1]))
                    else:
                        self.energies.append(np.array(
                            [float(i) for i in next(thefile).strip().split()[1:3]]))
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
        return procar.EnergyBand(physical_kvector, self.energies, self.spininfo)    

    def __str__(self):
        '''..py:method:: __str__()

        __str__() <=> str(x)

        Show the EIGENVAL character, not contents itself.
        '''
        template='''The parameter of EIGENVALUE
# of atoms: {0.n_atoms}
# of kpoints: {0.numk}
# of bands: {0.n_bands}
'''
        return template.format(self)
