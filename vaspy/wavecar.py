# -*- coding: utf-8 -*-
'''.. py:module:: wavecar
Module for WAVECAR class
'''

from __future__ import division, print_function  # Version safety
import os
import re
import numpy as np


class WAVECAR(object):
    '''.. py:class:: WAVECAR(WAVECAR_file)

    class for storing the data of WAVECAR file.

    Parameters
    ----------

    WAEVCAR_file: str
        File name of 'WAVECAR'

    Attributes
    ------------
    

    recl: numpy.int
        Record length
    nspin: numpy.int
        Number of spins
    rtag: numpy.int
        Tag for precsion in WAVECAR
    nkpts: numpy.int
        Number of k-points
    nbands: numpy.int
        Number of bands
    encut: numpy.float
        Cut off energy in eV unit.
    realcell: numpy.array
        Vectors for the unit cell in real space
    rcpcell: numpy.array
        Vectors for the unit cell in reciprocal space
        
    '''

    def __init__(self, filename='WAVECAR'):
        '''initialization of WAVECAR class'''
        self.filename =filename
        self.wfc = open(self.filename, 'rb')
        #
        # read the basic information
        self.readheader()
        # read the band information
        self.readband()
        self.wfc.close()

    def readheader(self):
        '''.. py:method:: readheader()

        Read the information of the system stored in the first
        two record in WAVECAR file

        rec1: recl, nspin, rtag
        rec2: nkpts, nbands ,encut ((cell(i, j) i=1, 3), j=1, 3)
        '''
        self.wfc.seek(0)
        self.recl, self.nspin, self.rtag = np.array(
            np.fromfile(self.wfc, dtype=np.float, count=3),
            dtype=int)
        if self.rtag == 45200:
            self.wfprec = np.complex64
        elif self.rtag == 45210:
            self.wfprec = np.complex128
        else:
            raise ValueError('Invalid TAG value: {}'.format(self.rtag)) 
        self.wfc.seek(self.recl)
        dump = np.fromfile(self.wfc, dtype=np.float, count=12)

        self.nkpts = int(dump[0])
        self.nbands = int(dump[1])
        self.encut = dump[2]
        self.realcell = dump[3:].reshape((3,3))
        self.volume = np.linalg.det(self.realcell)
        self.rcpcell = np.linalg.inv(self.realcell).T
                           

    def readband(self):
        '''.. py:method:: readband()
'''
        pass


    def __str__(self):
        ''' .. py:method:: __str__()
        Print out the system parameters
        '''
        string = "record length  =       {0}  spins =           {1}  prec flag        {2}".format(self.recl, self.nspin, self.rtag)
        string += "\nno. k points =          {0}".format(self.nkpts)
        string += "\nno. bands =          {0}".format(self.nbands)
        string +="\nreal space lattice vectors:"
        for i in range(3):
            string += "\na"+str(i+1)
            string += " = {0}    {1}    {2}".format(self.realcell[i][0],
                                                    self.realcell[i][1],
                                                    self.realcell[i][2])
        string += "\n\n"
        string += "\nvolume unit cell =   {0}".format(self.volume)
        string += "\nReciprocal lattice vectors:"
        for i in range(3):
            string += "\nb" + str(i+1)
            string += " = {0}    {1}    {2}".format(self.rcpcell[i][0],
                                                    self.rcpcell[i][1],
                                                    self.rcpcell[i][2])
        # string +="\nreciprocal lattice vector magnitudes:"
        return string
