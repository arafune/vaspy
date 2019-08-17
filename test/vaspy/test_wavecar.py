#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Test for WAVECAR class'''
import os

import numpy as np
# import tempfile
from nose.tools import assert_almost_equal, eq_, ok_, with_setup

import vaspy.poscar as poscar
import vaspy.wavecar as wavecar


class TestHatomWavecar(object):
    '''Class for Test WAVECAR module by using H_gamma.wavecar'''

    def setup(self):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + 'H_gamma.wavecar'
        self.h = wavecar.WAVECAR(data_file)

    @with_setup(setup=setup)
    def test_wavecar_header_and_band(self):
        '''Test for H atom WAVECAR header'''
        ok_(self.h.gamma)
        eq_(1, self.h.nspin)
        eq_(np.complex64, self.h.prec)
        eq_(1, self.h.numk)
        eq_(27, self.h.nbands)
        eq_(600, self.h.encut)
        np.testing.assert_array_equal(
            [[25., 0., 0.], [0., 25., 0.], [0., 0., 25.]], self.h.realcell)
        assert_almost_equal(25**3, self.h.volume)
        np.testing.assert_array_almost_equal([101, 101, 101], self.h.ngrid)
        #
        eq_(260834, self.h.nplwvs[0])
        np.testing.assert_array_equal((260834, 3), self.h.gvectors().shape)


class TestCOWavecar(object):
    '''Class for Test WAVECAR module by using WAVECAR.CO.wavecar'''

    def setup(self):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + 'WAVECAR.CO.wavecar'
        self.co = wavecar.WAVECAR(data_file)

    @with_setup(setup=setup)
    def test_wavecar_header(self):
        '''Test for CO molecule WAVECAR property

        rtag, wfprec
        '''
        eq_(711680, self.co.recl)  # record length
        eq_(2, self.co.nspin)  # spin
        eq_(45200, self.co.rtag)  # precision flag
        eq_(np.complex64, self.co.prec)

        eq_(1, self.co.numk)
        eq_(54, self.co.nbands)
        eq_(400, self.co.encut)
        assert_almost_equal(-8.05748789, self.co.efermi)
        np.testing.assert_array_equal(
            [[17., 0., 0.], [0., 17., 0.], [0., 0., 17.]], self.co.realcell)
        assert_almost_equal(4913, self.co.volume)
        np.testing.assert_array_almost_equal([57, 57, 57], self.co.ngrid)

    @with_setup(setup=setup)
    def test_wavecar_band(self):
        '''test for CO WAVECAR band'''
        kpath, kbands = self.co.band()
        eq_(None, kpath)
        eq_((self.co.nspin, self.co.numk, self.co.nbands), kbands.shape)
        np.testing.assert_array_almost_equal([
            -29.49120151, -14.03737717, -11.88106463, -11.88106223, -9.00976389
        ], kbands[0, 0, 0:5])  # The value can be taken from EIGENVAL

    @with_setup(setup=setup)
    def test_wavecar_str(self):
        """Test for __str__ special method."""
        output = self.co.__str__()
        teststr = "record length  =       711680  spins =           2  "
        teststr += "prec flag        45200\nno. k points =          1"
        teststr += "\nno. bands =          54\nreal space lattice vectors:\n"
        teststr += "a1 = 17.0    0.0    0.0\na2 = 0.0    17.0    0.0\n"
        teststr += "a3 = 0.0    0.0    17.0\n\nvolume unit cell =   "
        teststr += "4913.0000000000055\nReciprocal lattice vectors:\n"
        teststr += "b1 = 0.058823529411764705    0.0    0.0"
        teststr += "\nb2 = 0.0    0.058823529411764705"
        teststr += "    0.0\nb3 = 0.0    0.0    0.058823529411764705"
        eq_(teststr, output)


class TestGrapheneWavecar(object):
    '''Class for Test WAVECAR module by using Graphene.wavecar'''

    def setup(self):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + 'Graphene.wavecar'
        self.gr = wavecar.WAVECAR(data_file)
        self.gr_poscar = poscar.POSCAR(datadir + 'POSCAR.Graphene')

    def test_wavecar_header(self):
        '''Test for Graphene WAVECAR property'''
        eq_(30864, self.gr.recl)  # record length
        eq_(1, self.gr.nspin)  # spin
        eq_(45200, self.gr.rtag)  # precision flag
        eq_(np.complex64, self.gr.prec)
        #
        eq_(240, self.gr.numk)
        eq_(27, self.gr.nbands)
        eq_(550, self.gr.encut)
        assert_almost_equal(-2.9500542715, self.gr.efermi)
        np.testing.assert_array_almost_equal(
            [[2.139081750, -1.235000000, 0.], [2.139081750, 1.235000000, 0.],
             [0., 0., 24.7]], self.gr.realcell)
        # volume of cell in OUTCAR is not enough
        assert_almost_equal(130.50323848575005, self.gr.volume)
        # FIXME!!: Where these value come from ?
        np.testing.assert_array_almost_equal([11, 11, 97], self.gr.ngrid)

    @with_setup(setup=setup)
    def test_wavecar_band(self):
        '''Test for Graphene wavecar band'''
        kpath, kbands = self.gr.band()
        # from OUTCAR
        # maximum and minimum number of plane-waves per node :  3857  3785
        eq_(3857, self.gr.nplwvs.max())
        eq_(3785, self.gr.nplwvs.min())
        eq_((self.gr.numk, ), kpath.shape)  # gr.numk = 240
        eq_((self.gr.nspin, self.gr.numk, self.gr.nbands), kbands.shape)
        np.testing.assert_array_almost_equal(
            [-22.516876, -10.623282, -6.106901, -6.094072, 0.245639, 1.006991],
            kbands[0, 0, 0:6])  # The value can be taken from EIGENVAL

    @with_setup(setup=setup)
    def test_realsapece_wfc(self):
        '''Test for generation real space wfc (Graphene)'''
        np.testing.assert_array_almost_equal([
            0.00013770 + 0.0001743j, 0.00014605 + 0.00018611j, 0.00017262 +
            0.00022051j, 0.00021561 + 0.00027499j, 0.00026360 + 0.00033486j
        ],
                                             self.gr.realspace_wfc()[0][0][:5])
        vaspgrid = self.gr.realspace_wfc(poscar=self.gr_poscar)
        np.testing.assert_array_almost_equal(
            [0.00013770, 0.00014605, 0.00017262, 0.00021561, 0.00026360],
            vaspgrid.grid.data[:5])
        eq_(2, vaspgrid.grid.nframe)


class test_RestoreGammaGrid(object):
    #
    def test_check_symmetry(self):
        '''Test for check_symmetry'''
        ok_(not wavecar.check_symmetry(
            np.arange(3 * 3 * 3).reshape((3, 3, 3))))
        okgrid = np.array([[[0, 1, 1], [3, 4, 7], [3, 7, 4]],
                           [[9, 10, 11], [12, 13, 14], [15, 16, 17]],
                           [[9, 11, 10], [15, 17, 16], [12, 14, 13]]])
        ok_(wavecar.check_symmetry(okgrid))

    def test_RestorGammaGrid(self):
        '''Test for RestoreGammaGrid function'''
        grid332 = np.arange(3 * 3 * 3).reshape((3, 3, 3))[:, :, :2]
        grid332[2][0][0] = np.conjugate(grid332[1][0][0])
        grid332[2][2][0] = np.conjugate(grid332[1][1][0])
        grid332[1][2][0] = np.conjugate(grid332[2][1][0])
        grid332[0][2][0] = np.conjugate(grid332[0][1][0])
        result = wavecar.restore_gamma_grid(grid332, para=True)
        ok_(wavecar.check_symmetry(result))
        #
        grid553 = np.arange(5 * 5 * 5).reshape((5, 5, 5))[:, :, :3]
        #
        grid553[3, 0, 0] = np.conjugate(grid553[2, 0, 0])
        grid553[4, 0, 0] = np.conjugate(grid553[1, 0, 0])
        #
        grid553[0, 3, 0] = np.conjugate(grid553[0, 2, 0])
        grid553[1, 3, 0] = np.conjugate(grid553[4, 2, 0])
        grid553[2, 3, 0] = np.conjugate(grid553[3, 2, 0])
        grid553[3, 3, 0] = np.conjugate(grid553[2, 2, 0])
        grid553[4, 3, 0] = np.conjugate(grid553[1, 2, 0])
        #
        grid553[0, 4, 0] = np.conjugate(grid553[0, 1, 0])
        grid553[1, 4, 0] = np.conjugate(grid553[4, 1, 0])
        grid553[2, 4, 0] = np.conjugate(grid553[3, 1, 0])
        grid553[3, 4, 0] = np.conjugate(grid553[2, 1, 0])
        grid553[4, 4, 0] = np.conjugate(grid553[1, 1, 0])
        result = wavecar.restore_gamma_grid(grid553, para=True)
        ok_(wavecar.check_symmetry(result))
        #
        grid355 = np.arange(5 * 5 * 5).reshape((5, 5, 5))[:3, :, :]
        #
        grid355[0, 0, 3] = np.conjugate(grid355[0, 0, 2])
        grid355[0, 0, 4] = np.conjugate(grid355[0, 0, 1])
        #
        grid355[0, 3, 0] = np.conjugate(grid355[0, 2, 0])
        grid355[0, 3, 1] = np.conjugate(grid355[0, 2, 4])
        grid355[0, 3, 2] = np.conjugate(grid355[0, 2, 3])
        grid355[0, 3, 3] = np.conjugate(grid355[0, 2, 2])
        grid355[0, 3, 4] = np.conjugate(grid355[0, 2, 1])
        #
        grid355[0, 4, 0] = np.conjugate(grid355[0, 1, 0])
        grid355[0, 4, 1] = np.conjugate(grid355[0, 1, 4])
        grid355[0, 4, 2] = np.conjugate(grid355[0, 1, 3])
        grid355[0, 4, 3] = np.conjugate(grid355[0, 1, 2])
        grid355[0, 4, 4] = np.conjugate(grid355[0, 1, 1])
        result = wavecar.restore_gamma_grid(grid355, para=False)
        ok_(wavecar.check_symmetry(result))


class TestCobaltWavecar(object):
    '''Class for Test WAVECAR module by using Co.wavecar (SOI)'''

    def setup(self):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + 'Co.wavecar'
        self.co = wavecar.WAVECAR(data_file)
        self.co_poscar = poscar.POSCAR(datadir + 'Co.POSCAR')

    @with_setup(setup=setup)
    def test_wavecar_header(self):
        '''Test for Cobalt property'''
        eq_(7968, self.co.recl)  # record length
        eq_(1, self.co.nspin)  # spin
        eq_(45200, self.co.rtag)  # precision flag
        eq_(np.complex64, self.co.prec)
        #
        eq_(9, self.co.numk)
        eq_(54, self.co.nbands)
        eq_(400, self.co.encut)
        assert_almost_equal(-0.84118407515959326, self.co.efermi)
        np.testing.assert_array_almost_equal(
            np.array([[1.0, 0.0, 0.], [0.5, 0.8660254, 0.], [0., 0., 2.0]]) *
            2.501, self.co.realcell)
        # volume of cell in OUTCAR is not enough
        assert_almost_equal(27.095782694613046, self.co.volume)
        # FIXME!!: Where these value come from ?
        #  Maximum number of reciprocal cells 2x 2x 4 (in OUTCAR)
        np.testing.assert_array_almost_equal([11, 11, 19], self.co.ngrid)

    @with_setup(setup=setup)
    def test_wavecar_band(self):
        '''Test for Co wavecar band'''
        kpath, kbands = self.co.band()
        # from OUTCAR
        eq_(996, self.co.nplwvs.max())
        eq_(958, self.co.nplwvs.min())
        eq_((self.co.numk, ), kpath.shape)  # co.numk = 240
        eq_((self.co.nspin, self.co.numk, self.co.nbands), kbands.shape)
        np.testing.assert_array_almost_equal(
            [-6.532492, -5.858599, -4.037263, -3.418892, -3.265558, -2.642231],
            kbands[0, 0, 0:6])  # The value can be taken from EIGENVAL

    @with_setup(setup=setup)
    def test_realsapece_wfc(self):
        '''Test for generation real space wfc (Cobalt)'''
        np.testing.assert_array_almost_equal(
            [
                -7.84915157e-05 - 5.61362047e-05j, -7.88077954e-05 -
                5.63624713e-05j, -7.92001487e-05 - 5.66431905e-05j,
                -7.92675956e-05 - 5.66915450e-05j,
                -7.90423140e-05 - 5.65305135e-05j
            ],
            self.co.realspace_wfc()[0][0][0][:5])
        vaspgrid = self.co.realspace_wfc(poscar=self.co_poscar)
        np.testing.assert_array_almost_equal([
            -7.84915157e-05, -7.88077954e-05, -7.92001487e-05, -7.92675956e-05,
            -7.90423140e-05
        ], vaspgrid.grid.data[:5])
        eq_(4, vaspgrid.grid.nframe)
