#! /usr/bin/env python
# -*- coding: utf-8 -+-

import os
from nose.tools import with_setup
from nose.tools import eq_
from nose.tools import ok_
import numpy as np
import vaspy.outcar

class TestOUTCAR(object):
    '''class for test of vaspy.outcar module
    '''
    def setup(self):
        datadir = os.path.abspath(os.path.dirname(__file__)) + '/data/'
        datafile = datadir + 'OUTCAR'
        self.outcar = vaspy.outcar.OUTCAR(datafile)


    @with_setup(setup=setup)
    def test_fermi(self):
        eq_(-0.7681, self.outcar.fermi)

    def test_recvec(self):
        np.testing.assert_almost_equal(
            np.array([0.389375981, -0.224806327, 0.000000000]),
            self.outcar.recvec[0])
        np.testing.assert_almost_equal(
            np.array([0.000000000, 0.449612655, 0.000000000]),
            self.outcar.recvec[1])
        np.testing.assert_almost_equal(
            np.array([0.000000000, 0.000000000, 0.023484312]),
            self.outcar.recvec[2])
