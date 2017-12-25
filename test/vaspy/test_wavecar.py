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

class TestCOWavecar(object):
    '''Class for Test WAVECAR module by using WAVECAR.CO.wavecar'''
    def setup(self):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + 'WAVECAR.CO.wavecar'
        self.co = wavecar.WAVECAR(data_file)
        
    @with_setup(setup=setup)
    def test_wavecar_readheader(self):
        '''Test for WAVECAR property 

        rtag, wfprec
        '''
        eq_(711680, self.co.recl)  # record length
        eq_(2, self.co.nspin)      # spin
        eq_(45200, self.co.rtag)   # precision flag
        eq_(np.complex64, self.co.wfprec)
        #
        eq_(1, self.co.nkpts)
        eq_(54, self.co.nbands)
        eq_(400, self.co.encut)
        np.testing.assert_array_equal(
            [[17.,   0.,   0.],
             [  0.,  17.,   0.],
             [  0.,   0.,  17.]], self.co.realcell)
        assert_almost_equal(4913, self.co.volume)
        np.testing.assert_array_almost_equal([57, 57, 57], self.co.ngrid)
        
    @with_setup(setup=setup)
    def test_wavecar_readband(self):
        kpath, bands = co.readband()
        
