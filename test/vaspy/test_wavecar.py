#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Test for WAVECAR class'''
# import unittest
import os
# import tempfile
from nose.tools import eq_, ok_
from nose.tools import with_setup, raises, assert_almost_equal
import numpy as np
import vaspy.wavecar as wavecar
import vaspy.poscar as poscar

class TestCOWavecar(object):
    '''Class for Test WAVECAR module by using WAVECAR.CO.wavecar'''
    def setup(self):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + 'WAVECAR.CO.wavecar'
        self.co = wavecar.WAVECAR(data_file)
        
    @with_setup(setup=setup)
    def test_wavecar_readheader(self):
        '''Test for CO molecule WAVECAR property

        rtag, wfprec
        '''
        eq_(711680, self.co.recl)  # record length
        eq_(2, self.co.nspin)      # spin
        eq_(45200, self.co.rtag)   # precision flag
        eq_(np.complex64, self.co.prec)
        
        eq_(1, self.co.numk)
        eq_(54, self.co.nbands)
        eq_(400, self.co.encut)
        np.testing.assert_array_equal(
            [[17.,   0.,   0.],
             [  0.,  17.,   0.],
             [  0.,   0.,  17.]], self.co.realcell)
        assert_almost_equal(4913, self.co.volume)
        np.testing.assert_array_almost_equal([57, 57, 57], self.co.ngrid)
        
    @with_setup(setup=setup)
    def test_wavecar_band(self):
        '''test for CO WAVECAR band'''
        kpath, kbands = self.co.band()
        eq_(None, kpath)
        eq_((self.co.nspin, self.co.numk, self.co.nbands), kbands.shape)
        np.testing.assert_array_almost_equal(
            [-29.49120151, -14.03737717, -11.88106463, -11.88106223,  -9.00976389],
            kbands[0, 0, 0:5])  # The value can be taken from EIGENVAL

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
        eq_(1, self.gr.nspin)      # spin
        eq_(45200, self.gr.rtag)   # precision flag
        eq_(np.complex64, self.gr.prec)
        
        eq_(240, self.gr.numk)
        eq_(27, self.gr.nbands)
        eq_(550, self.gr.encut)
        np.testing.assert_array_almost_equal(
            [[2.139081750,   -1.235000000,   0.],
             [2.139081750,   1.235000000,   0.],
             [  0.,   0.,  24.7]], self.gr.realcell)
        # volume of cell in OUTCAR is not enough
        assert_almost_equal(130.50323848575005, self.gr.volume)
        # FIXME!!: Where these value come from ? 
        np.testing.assert_array_almost_equal([11, 11, 97], self.gr.ngrid) 
        
    @with_setup(setup=setup)
    def test_wavecar_band(self):
        '''Test for Graphene wavecar band'''
        kpath, kbands = self.gr.band()
        # from OUTCAR
        # maximum and minimum number of plane-waves per node :      3857     3785
        eq_(3857, self.gr.nplwvs.max())
        eq_(3785, self.gr.nplwvs.min())
        eq_((self.gr.numk,), kpath.shape)  # gr.numk = 240
        eq_((self.gr.nspin, self.gr.numk, self.gr.nbands), kbands.shape)
        np.testing.assert_array_almost_equal(
            [-22.516876, -10.623282, -6.106901, -6.094072, 0.245639, 1.006991],
            kbands[0, 0, 0:6])  # The value can be taken from EIGENVAL

    @with_setup(setup=setup)
    def test_realsapece_wfc(self):
        '''Test for generation real space wfc (Graphene)'''
        np.testing.assert_array_almost_equal(
            [0.00013770+0.0001743j, 0.00014605+0.00018611j,
             0.00017262+0.00022051j, 0.00021561+0.00027499j,
             0.00026360+0.00033486j], 
         self.gr.realspace_wfc()[0][0][:5])
        vaspgrid = self.gr.realspace_wfc(poscar=self.gr_poscar)
        np.testing.assert_array_almost_equal(
            [0.00013770, 0.00014605, 0.00017262, 0.00021561,
             0.00026360], vaspgrid.grid.data[:5])
        eq_(2, vaspgrid.grid.num_frame)
