# -*- coding: utf-8 -+-
'''Unit test for OUTCAR class'''
import os

import numpy as np
from nose.tools import eq_, ok_, with_setup

import vaspy.outcar


class TestOUTCAR(object):
    '''class for test of vaspy.outcar module
    '''

    def setup(self):
        '''Reading OUTCAR file for test'''
        datadir = os.path.abspath(os.path.dirname(__file__)) + '/data/'
        datafile = datadir + 'OUTCAR'
        self.outcar = vaspy.outcar.OUTCAR(datafile)

    @with_setup(setup=setup)
    def test_read_basic_properties(self):
        '''Test for OUTCAR reading basic properties'''
        eq_(-0.7681, self.outcar.fermi)
        eq_(54, self.outcar.nbands)
        #
        # reciprocal vectors
        np.testing.assert_almost_equal(
            np.array([0.389375981, -0.224806327, 0.000000000]),
            self.outcar.recvec[0])
        np.testing.assert_almost_equal(
            np.array([0.000000000, 0.449612655, 0.000000000]),
            self.outcar.recvec[1])
        np.testing.assert_almost_equal(
            np.array([0.000000000, 0.000000000, 0.023484312]),
            self.outcar.recvec[2])
        # kvecs and weights
        eq_(33, len(self.outcar.kvecs))
        eq_(33, len(self.outcar.weights))
        np.testing.assert_array_almost_equal(
            [[0.000000, 0.000000, 0.000000], [0.058824, 0.000000, 0.000000],
             [0.117647, 0.000000, 0.000000], [0.176471, -0.000000, 0.000000],
             [0.235294, 0.000000, 0.000000], [0.294118, 0.000000, 0.000000],
             [0.352941, -0.000000, 0.000000]], self.outcar.kvecs[:7])
        np.testing.assert_array_almost_equal([
            1.000000, 6.000000, 6.000000, 6.000000, 6.000000, 6.000000,
            6.000000
        ], self.outcar.weights[:7])
        # total charges
        np.testing.assert_array_almost_equal(
            [[0.454, 0.282, 9.198], [0.446, 0.366, 9.192], [
                0.447, 0.372, 9.185
            ], [0.451, 0.377, 9.184], [0.451, 0.381, 9.185],
             [0.450, 0.376, 9.192], [0.456, 0.288, 9.200]],
            self.outcar.total_charges[-1])
        eq_(11, len(self.outcar.total_charges))
        eq_(7, len(self.outcar.total_charges[0]))
        eq_(3, len(self.outcar.total_charges[0][0]))
