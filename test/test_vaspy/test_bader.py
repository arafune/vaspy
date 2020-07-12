# -*- coding: utf-8 -+-
'''Unit test for OUTCAR class'''
import os

import numpy as np
from nose.tools import eq_, ok_, with_setup

import vaspy.bader


class TestBader(object):
    """Class for test of vaspy.outcar module."""

    def setup(self):
        '''Reading ACF.chg.dat for test'''
        datadir = os.path.abspath(os.path.dirname(__file__)) + '/data/'
        datafile = datadir + 'ACF.chg.dat'
        self.acf = vaspy.bader.BaderACF(datafile)

    @with_setup(setup=setup)
    def test_read_basic_properties(self):
        """Test for reading basic properties."""
        """
        VACUUM CHARGE:               1.2829
        VACUUM VOLUME:             157.7708
        NUMBER OF ELECTRONS:      2828.0000
        """
        eq_(1.2829, self.acf.vaccharge)
        eq_(157.7708, self.acf.vacvol)
        eq_(2828.000, self.acf.nelectron)
        eq_(344, self.acf.natom)
        np.testing.assert_array_almost_equal([6.745000, -37.491966, 11.674620],
                                             self.acf.positions[0])
        eq_(3.526327, self.acf.charges[0])
        eq_(0.423581, self.acf.mindists[0])
        eq_(17.540723, self.acf.vols[0])
