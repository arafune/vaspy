# -*- coding: utf-8 -*-
'''
Module for WAVECAR class
'''

from __future__ import division, print_function  # Version safety
import numpy as np
from scipy.fftpack import ifftn
import vaspy.mesh3d as mesh3d
import vaspy.poscar as poscar

Ry_in_eV = 13.605826
au_in_AA = 0.529177249


class WAVECAR(object):
    ''' 
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
    nplwvs: numpy.int
        Number of plane waves.
    numk: numpy.int
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
        self.wfc = open(filename, 'rb')
        #
        self.header()
        self.band()

    def header(self):
        ''' 
        Read the information of the system stored in the first
        two record in WAVECAR file

        rec1: recl, nspin, rtag
        rec2: numk, nbands ,encut, ((cell(i, j) i=1, 3), j=1, 3), efermi
        '''
        self.wfc.seek(0)
        self.recl, self.nspin, self.rtag = np.array(
            np.fromfile(self.wfc, dtype=np.float, count=3),
            dtype=int)
        self.wfc.seek(self.recl)
        #        print(self.wfc.tell())
        #
        dump = np.fromfile(self.wfc, dtype=np.float, count=13)
        #
        self.numk = int(dump[0])
        self.nbands = int(dump[1])
        self.encut = dump[2]
        self.realcell = dump[3:12].reshape((3, 3))
        self.efermi = dump[12]
        #        print(self.wfc.tell())
        self.volume = np.linalg.det(self.realcell)
        self.rcpcell = np.linalg.inv(self.realcell).T
        unit_cell_vector_magnitude = np.linalg.norm(self.realcell, axis=1)
        cutoff = np.ceil(
            np.sqrt(self.encut / Ry_in_eV) / (2*np.pi / (
                unit_cell_vector_magnitude / au_in_AA)))
        # FFT Minimum grid size
        self.ngrid = np.array(2 * cutoff + 1, dtype=int)
        #         self.ngrid = np.array(cutof, dtype=int) でよい？

    @property
    def prec(self):
        '''Return precision determined from self.rtag'''
        if self.rtag == 45200:
            return np.complex64
        elif self.rtag == 45210:
            return np.complex128
        else:
            raise ValueError('Invalid TAG value: {}'.format(self.rtag))

    def band(self):
        ''' 
        Read the information about the band from WAVECAR file

        The infomation obtained by this method is as follows:

        * Number of plane waves (nplwv)
        * A integer set for k-vectors
        * energy of the band (as a function of spin-, k-, and band index)
        * occupation  (as a function of spin-, k-, and band index)

        '''
        self.kvecs = np.zeros((self.numk, 3), dtype=float)
        self.bands = np.zeros((self.nspin, self.numk, self.nbands),
                              dtype=float)
        self.nplwvs = np.zeros(self.numk, dtype=int)
        self.occs = np.zeros((self.nspin, self.numk, self.nbands),
                             dtype=float)
        for spin_i in range(self.nspin):
            for k_i in range(self.numk):
                pos = 2 + spin_i * self.numk * (self.nbands + 1)
                pos += k_i * (self.nbands + 1)
                self.wfc.seek(pos * self.recl)
#                print(self.wfc.tell())
                dump = np.fromfile(self.wfc, dtype=np.float,
                                   count=4+3*self.nbands)
                if spin_i == 0:
                    self.nplwvs[k_i] = int(dump[0])
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
        if self.numk == 1:
            self.kpath = None
        return self.kpath, self.bands

    def gvectors(self, k_i=0):
        r'''
        G-vectors :math:`G` is determined by the following condition:
        :math:`\frac{(G+k)^2}{ 2} < E_{cut}`

        note: hbar2over2m is :math:`\frac{\hbar^2}{2m}`

        * 0.529177249 is au unit in AA
        * 13.605826 is Ry unit in eV

        Parameters
        ------------
        k_i: int, optional
           k index :math:`k_i` (the default value is 0).

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
        hbar2over2m = 13.605826 * 0.529177249 * 0.529177249
        energy_k = hbar2over2m * np.linalg.norm(
            np.dot(kgrid + kvec[np.newaxis, :], 2*np.pi*self.rcpcell),
            axis=1)**2
        g_vec = kgrid[np.where(energy_k < self.encut)[0]]
        return np.asarray(g_vec, dtype=int)

    def bandcoeff(self, spin_i=0, k_i=0, band_i=0, norm=False):
        ''' 
        Read the coefficient of the planewave of the KS
        states specified by the `spin_i`, `k_i` and `band_i`

        Parameters
        ----------
        spin_i: int, optional
           spin index :math:`s_i` (0 or 1) (default value is 0)
        k_i: int, optioanl
           k index :math:`k_i`. Starts with 0 (default value is 0)
        band_i: int, optioanl
            band index :math:`b_i`. starts with 0 (default value is 0)
        norm: bool, optioanl
            If true the Band coeffients are normliazed (default is false)
        '''
        irec = 3 + spin_i * self.numk * (self.nbands + 1)
        irec += k_i * (self.nbands + 1) + band_i
        self.wfc.seek(irec * self.recl)
#        print(self.wfc.tell())
        nplw = self.nplwvs[k_i]
        dump = np.fromfile(self.wfc, dtype=self.prec, count=nplw)
        print(self.wfc.tell())        
        cg = np.array(dump, dtype=np.complex128)
        if norm:
            cg /= np.linalg.norm(cg)
        return cg

    def realspace_wfc(self, spin_i=0, k_i=0, band_i=0,
                      gvec=None, ngrid=None, norm=False,
                      poscar=poscar.POSCAR()):
        r''' 
        Calculate the pseudo-wavefunction of the KS states in
        the real space by using FFT transformation of the reciprocal
        space planewave coefficients. 

        The 3D FE grid size is detemined by ngrid, which defaults
        to self.ngrid if it is not provided.  GVectors of the KS
        states is used to put 1D plane wave coefficient back to 3D
        grid.

        Parameters
        -----------
        spin_i: int
           spin index (0 or 1)
        k_i: int
           k index :math:`k_i`. Starts with 0. default is 0
        band_i: int
            band index :math:`b_i`. starts with 0. default is 0.
        norm: bool
            If true the Band coeffients are normliazed
        gvec: numpy.array
            G-vector for calculation. If not set, use gvectors(k_i)
        ngrid: numpy.array
            Ngrid for calculation. If not set, use self.ngrid.
        poscar: vaspy.poscar, optional
            POSCAR object (defalut is blank POSCAR object)

        Returns
        -----------

        phi_r.T: numpy.array
            If poscar is not specified, for Collinear-wavecar file.
            phi_r is mesh data for the wavefunction in the real space.
            .T is due to the original data is made by fortran program
            (i.e. 'VASP', of course.

        phi_r_up.T, phi_r_down.T: tuple of numpy.array
            If poscar is not specified, for SOI-wavecar file.
            phi_r_up corresponds to the 'up' spin spinor wavefunction.
            phi_r_down corresponds to the 'down' spin spinor wavefunction.

        vaspgrid: VASPGrid
            Returns VASPGrid object, if poscar is specified.  The former frame
            represents the real part of the wavefunction at :math:`k_i` and
            :math:`b_i` in the real space, the latter frame the imaginary
            part. On the other hand, the SOI-wavecar has 4 frames
            The first and second are for the "up" spin, and the third
            and fourth are "down" spin. (Judging SOI by
            gvectors(k_i).shape[0] :math:`\neq` bandcoeff(k_i).size)
        '''
        if ngrid is None:
            ngrid = self.ngrid.copy()
        else:
            ngrid = np.array(ngrid, dtype=int)
        if gvec is None:
            gvec = self.gvectors(k_i)
        gvec %= ngrid[np.newaxis, :]
        try:  # Collininear
            phi_k = np.zeros(ngrid, dtype=np.complex128)
            phi_k[gvec[:, 0], gvec[:, 1], gvec[:, 2]] = self.bandcoeff(spin_i,
                                                                       k_i,
                                                                       band_i,
                                                                       norm)
            phi_r = ifftn(phi_k)
            if poscar.scaling_factor == 0.:
                return phi_r.T
            else:
                vaspgrid = mesh3d.VASPGrid()
                vaspgrid.poscar = poscar
                vaspgrid.grid.shape = ngrid
                # checking consistency between POSCAR and WAVECAR
                np.testing.assert_array_almost_equal(
                    poscar.scaling_factor * poscar.cell_vecs,
                    self.realcell)
                re = np.real(phi_r.flatten('F'))
                im = np.imag(phi_r.flatten('F'))
                vaspgrid.grid.data = np.concatenate((re, im))
                return vaspgrid
        except ValueError:   # SOI
            phi_k_up = np.zeros(ngrid, dtype=np.complex128)
            phi_k_down = np.zeros(ngrid, dtype=np.complex128)
            bandcoeff = self.bandcoeff(spin_i, k_i, band_i, norm)
            bandcoeff_up = bandcoeff.reshape(2, bandcoeff.size//2)[0]
            bandcoeff_down = bandcoeff.reshape(2, bandcoeff.size//2)[1]
            phi_k_up[gvec[:, 0], gvec[:, 1], gvec[:, 2]] = bandcoeff_up
            phi_k_down[gvec[:, 0], gvec[:, 1], gvec[:, 2]] = bandcoeff_down
            phi_r_up = ifftn(phi_k_up)
            phi_r_down = ifftn(phi_k_down)
            if poscar.scaling_factor == 0.:
                return phi_r_up.T, phi_r_down.T
            else:
                vaspgrid = mesh3d.VASPGrid()
                vaspgrid.poscar = poscar
                vaspgrid.grid.shape = phi_k_up.shape
                # checking consistency between POSCAR and WAVECAR
                np.testing.assert_array_almost_equal(
                    poscar.scaling_factor * poscar.cell_vecs,
                    self.realcell)
                up_re = np.real(phi_r.up.flatten('F'))
                up_im = np.imag(phi_r_up.flatten('F'))
                down_re = np.real(phi_r_down.flatten('F'))
                down_im = np.imag(phi_r_down.flatten('F'))
                vaspgrid.grid.data = np.concatenate((up_re, up_im,
                                                     down_re, down_im))
                return vaspgrid

    def __str__(self):
        ''' 
        Print out the system parameters
        '''
        the1stline = "record length  =       {0}  "
        the1stline += "spins =           {1}  "
        the1stline += "prec flag        {2}"
        string = the1stline.format(self.recl, self.nspin, self.rtag)
        string += "\nno. k points =          {0}".format(self.numk)
        string += "\nno. bands =          {0}".format(self.nbands)
        string += "\nreal space lattice vectors:"
        for i in range(3):
            string += "\na"+str(i+1)
            string += " = {0}    {1}    {2}".format(self.realcell[i][0],
                                                    self.realcell[i][1],
                                                    self.realcell[i][2])
        string += "\n"
        string += "\nvolume unit cell =   {0}".format(self.volume)
        string += "\nReciprocal lattice vectors:"
        for i in range(3):
            string += "\nb" + str(i+1)
            string += " = {0}    {1}    {2}".format(self.rcpcell[i][0],
                                                    self.rcpcell[i][1],
                                                    self.rcpcell[i][2])
        # string +="\nreciprocal lattice vector magnitudes:"
        return string
