#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
from nose.tools import eq_
from nose.tools import ok_
from nose.tools import with_setup
from nose.tools import assert_equal
from nose.tools import raises
import vaspy.eigenval as eigenval
import numpy as np

class TestEIGENVAL(object):
    '''Class for EIGENVAL class test'''

    def setup(self):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        self.eigenval_spin = eigenval.EIGENVAL(datadir + 'EIGENVAL.spin')
        self.eigenval_soi = eigenval.EIGENVAL(datadir + 'EIGENVAL.soi')
        
    @with_setup(setup=setup)
    def test_check_basic_parameters(self):
        '''Checking the basic parameters stored'''
        eq_(344, self.eigenval_spin.n_atoms)
        eq_(1890, self.eigenval_spin.n_bands)
        eq_(1, self.eigenval_spin.numk)
        eq_((1890, 2), self.eigenval_spin.energies.shape)
        np.testing.assert_array_equal([0., 0., 0.],
                                      self.eigenval_spin.kvectors[0])
        np.testing.assert_array_equal([-22.893072, -22.890912],
                                      self.eigenval_spin.energies[0])
        np.testing.assert_array_equal([  5.806516,   5.813531],
                                      self.eigenval_spin.energies[-1])
        #
        eq_(1, self.eigenval_soi.n_atoms)
        eq_(54, self.eigenval_soi.n_bands)
        eq_(625, self.eigenval_soi.numk)
        np.testing.assert_array_equal(-9.548122,
                                      self.eigenval_soi.energies[0])
        np.testing.assert_array_almost_equal(26.193038,
                                      self.eigenval_soi.energies[-1])

