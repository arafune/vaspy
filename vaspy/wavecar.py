# -*- coding: utf-8 -*-
'''.. py:module:: wavecar
Module for WAVECAR class
'''

from __future__ import division, print_function  # Version safety
import numpy as np
import vaspy.const as const
from scipy.fftpack import ifftn


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
    nplws: numpy.int
        Number of plane waves
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
    kvecs: numpy.array
        kvector
    bands: numpy.array
        Energy
    occs: numpy.array
        Occupation
    '''

    def __init__(self, filename='WAVECAR'):
        '''initialization of WAVECAR class'''
        self.filename = filename
        self.wfc = open(self.filename, 'rb')
        #
        # read the basic information
        self.readheader()
        # read the band information
        self.readband()

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
        self.realcell = dump[3:].reshape((3, 3))
        self.volume = np.linalg.det(self.realcell)
        self.rcpcell = np.linalg.inv(self.realcell).T
        unit_cell_vector_magnitude = np.linalg.norm(self.realcell, axis=1)
        cutof = np.ceil(
            np.sqrt(self.encut / const.RytoeV) / (2*np.pi / (
                unit_cell_vector_magnitude / const.au_to_A)))
        self.ngrid = np.array(2 * cutof + 1, dtype=int)
        #         self.ngrid = np.array(cutof, dtype=int) でよい？

    def readband(self):
        '''.. py:method:: readband()

        Read the information about the band from WAVECAR file

        The infomation obtained by this method is as follows:

        * Number of plane waves (nplws)
        * A integer set for k-vectors
        * energy of the band (as a function of spin-, k-, and band index)
        * occupation  (as a function of spin-, k-, and band index)

        '''
        self.nplws = np.zeros(self.nkpts, dtype=int)
        self.kvecs = np.zeros((self.nkpts, 3), dtype=float)
        self.bands = np.zeros((self.nspin, self.nkpts, self.nbands),
                              dtype=float)
        self.occs = np.zeros((self.nspin, self.nkpts, self.nbands),
                             dtype=float)
        for spin_i in range(self.nspin):
            for k_i in range(self.nkpts):
                pos = 2 + spin_i * self.nkpts * (self.nbands + 1)
                pos += k_i * (self.nbands + 1) + 1 - 1
                # the last '1' corresponds iband=1
                self.wfc.seek(pos * self.recl)
                dump = np.fromfile(self.wfc, dtype=np.float,
                                   count=4+3*self.nbands)
                if spin_i == 0:
                    self.nplws[k_i] = int(dump[0])
                    self.kvecs[k_i] = dump[1:4]
                dump = dump[4:].reshape((-1, 3))
                self.bands[spin_i, k_i, :] = dump[:, 0]
                self.occs[spin_i, k_i, :] = dump[:, 2]
        self.kpath = np.concatenate(([0, ],
                                     np.cumsum(
                                         np.linalg.norm(
                                             np.dot(
                                                 np.diff(self.kvecs, axis=0),
                                                 self.rcpcell),
                                             axis=1))))
        if self.nkpts == 1:
            self.kpath = None
        return self.kpath, self.bands

    def gvectors(self, k_i=0):
        '''
        G-vectors :math:`G` is determined by the following condition:
            :math:`(G+k)^2 / 2 < E_{cut}`

        Parameters
        ------------
        k_i: int
           k index

        Returns
        ---------
        numpy.ndarray
           G vectors
        '''

        kvec = self.kvecs[k_i]
        fx = [ii if ii < self.ngrid[0] / 2 + 1
              else ii - self.ngrid[0]
              for ii in range(self.ngrid[0])]
        fy = [ii if ii < self.ngrid[1] / 2 + 1
              else ii - self.ngrid[1]
              for ii in range(self.ngrid[1])]
        fz = [ii if ii < self.ngrid[2] / 2 + 1
              else ii - self.ngrid[2]
              for ii in range(self.ngrid[2])]
        kgrid = np.array([(fx[ix], fy[iy], fz[iz])
                          for iz in range(self.ngrid[2])
                          for iy in range(self.ngrid[1])
                          for ix in range(self.ngrid[0])], dtype=float)
        energy_k = const.HSQDTM * np.linalg.norm(
            np.dot(kgrid + kvec[np.newaxis, :], 2*np.pi*self.rcpcell),
            axis=1)**2
        g_vec = kgrid[np.where(energy_k < self.encut)[0]]
        return np.asarray(g_vec, dtype=int)

    def bandcoeff(self, spin_i=0, k_i=0, band_i=0, norm=False):
        '''Read the coefficient of the planewave of the KS
        states specified by the `spin_i`, `k_i` and `band_i`

        Parameters
        ----------
        spin_i: int
           spin index (0 or 1)
        k_i: int
           k index. Starts with 0
        band_i: int
            band index. starts with 0
        norm: bool
            If true the Band coeffients are normliazed
        '''
        pos = 2 + spin_i * self.nkpts * (self.nbands + 1)
        pos += k_i * (self.nbands + 1) + band_i + 1
        self.wfc.seek(pos * self.recl)
        nplw = self.nplws[k_i]
        dump = np.fromfile(self.wfc, dtype=self.wfprec, count=nplw)
        cg = np.array(dump, dtype=np.complex128)
        if norm:
            cg /= np.linalg.norm(cg)
        return cg

    def realspace_wfc(self, spin_i=0, k_i=0, band_i=0,
                      gvec=None, ngrid=None, norm=False):
        '''.. py:method:: realspace_wfc(spin_i, k_i, band_i, gvec, ngrid, norm)

        Calculate the pseudo-wavefunction of the KS states in
        the real space by using FFT transformation of the reciprocal
        space planewave coefficients.

        The 3D FE grid size is detemined by ngrid, which defaults
        to self.ngrid if it is not provided.  GVectors of the KS
        states is used to put 1D plane wave coefficient back to 3D
        grid.

        Parameters
        ----------
        spin_i: int
           spin index (0 or 1)
        k_i: int
           k index. Starts with 0. default is 0
        band_i: int
            band index. starts with 0. default is 0.
        norm: bool
            If true the Band coeffients are normliazed
        gvec: numpy.array
            G-vector for calculation. If not set, use gvectors(k_i)
        ngrid: numpy.array
            Ngrid for calculation. If not set, use self.ngrid.
        '''
        if ngrid is None:
            ngrid = self.ngrid.copy()
        else:
            ngrid = np.array(ngrid, dtype=int)
        if gvec is None:
            gvec = self.gvectors(k_i)
        phi_k = np.zeros(ngrid, dtype=np.complex128)
        gvec %= ngrid[np.newaxis, :]
        phi_k[gvec[:, 0],
              gvec[:, 1],
              gvec[:, 2]] = self.bandcoeff(spin_i,
                                           k_i,
                                           band_i,
                                           norm)
        return ifftn(phi_k)

    def __str__(self):
        ''' .. py:method:: __str__()
        Print out the system parameters
        '''
        the1stline = "record length  =       {0}  "
        the1stline += "spins =           {1}  "
        the1stline += "prec flag        {2}"
        string = the1stline.format(self.recl, self.nspin, self.rtag)
        string += "\nno. k points =          {0}".format(self.nkpts)
        string += "\nno. bands =          {0}".format(self.nbands)
        string += "\nreal space lattice vectors:"
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
