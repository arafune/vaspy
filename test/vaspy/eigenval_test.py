# -*- coding: utf-8 -*-
'''Test for EIGENVAL class'''
import os

import numpy as np
from nose.tools import assert_equal, eq_, ok_, raises, with_setup

import vaspy.eigenval as eigenval


class TestEIGENVAL(object):
    '''Class for EIGENVAL class test'''

    def setup(self):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        self.eigenval_spin = eigenval.EIGENVAL(datadir + 'EIGENVAL.spin')
        self.eigenval_soi = eigenval.EIGENVAL(datadir + 'EIGENVAL.soi')

    @with_setup(setup=setup)
    def test_check_basic_parameters(self):
        '''Checking the basic parameters stored'''
        eq_(344, self.eigenval_spin.natom)
        eq_(2, self.eigenval_spin.nspin)
        eq_(1890, self.eigenval_spin.nbands)
        eq_(1, self.eigenval_spin.numk)
        eq_((2, 1, 1890), self.eigenval_spin.energies.shape)
        np.testing.assert_array_equal([0., 0., 0.],
                                      self.eigenval_spin.kvecs[0])
        np.testing.assert_array_almost_equal([[
            -22.893072, -22.891405, -22.841659, -22.832997, -22.761726,
            -22.737184
        ]], self.eigenval_spin.energies[0, :, :6])
        np.testing.assert_array_equal([5.806516],
                                      self.eigenval_spin.energies[0, :, -1])
        np.testing.assert_array_equal([5.813531],
                                      self.eigenval_spin.energies[1, :, -1])
        #
        eq_(1, self.eigenval_soi.natom)
        eq_(1, self.eigenval_soi.nspin)
        eq_(54, self.eigenval_soi.nbands)
        eq_(625, self.eigenval_soi.numk)
        np.testing.assert_array_equal(
            [-9.548122, -9.491253, -9.322813, -9.050432, -8.691472, -8.285128],
            self.eigenval_soi.energies[0, :, 0][0:6])
        np.testing.assert_array_almost_equal(
            26.193038, self.eigenval_soi.energies[0, :, -1][-1])

    @with_setup(setup=setup)
    def test_make_label(self):
        '''test for make_label function'''
        labels_soi = self.eigenval_soi.make_label('k', 'energy')
        eq_('#k', labels_soi[0])
        eq_('Energy', labels_soi[1])
        labels_spin = self.eigenval_spin.make_label('k', 'energy')
        eq_('#k', labels_spin[0])
        eq_('Energy_up', labels_spin[1])
        eq_('Energy_down', labels_spin[2])
