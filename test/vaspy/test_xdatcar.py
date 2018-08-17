#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
# import tempfile
from nose.tools import eq_, ok_
from nose.tools import with_setup, assert_equal, raises
from nose.tools import assert_false
import numpy as np
import vaspy.xdatcar as xdatcar


class TestXDATCAR(object):
    '''Class for Test of Vsim_asc

'''
    def setup(self):
        '''XDATCAR'''
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        datafile = datadir + "XDATCAR.0.bz2"
        self.xdatcar_test = xdatcar.XDATCAR(datafile)

    @with_setup(setup=setup)
    def test_(self):
        eq_('MoS2', self.xdatcar_test.system_name)
        eq_(1.0, self.xdatcar_test.scaling_factor)
        np.testing.assert_allclose([3.184000, 0.000000, 0.000000],
                                   self.xdatcar_test.cell_vecs[0])
        np.testing.assert_allclose([1.592000, 2.757425, 0.000000],
                                   self.xdatcar_test.cell_vecs[1])
        np.testing.assert_allclose([0, 0, 38.0],
                                   self.xdatcar_test.cell_vecs[2])
        eq_(['Mo', 'S'], self.xdatcar_test.iontypes)
        eq_([1, 2], self.xdatcar_test.ionnums)
        eq_(5, len(self.xdatcar_test.configurations))
        #
        np.testing.assert_allclose([[0.33333333,  0.33333333,  0.24999815],
                                    [0.66666667,  0.66666667,  0.29113380],
                                    [0.66666667,  0.66666667,  0.20886348]],
                                    self.xdatcar_test.configurations[0])
        np.testing.assert_allclose([[   0.33333333,  0.33333333,  0.24999815],
                                    [   0.66666667,  0.66666667,  0.29113854],
                                    [0.66666667,  0.66666667,  0.20885874]],
                                    self.xdatcar_test.configurations[4])
        
        # #  mode is '3', second atom
        # np.testing.assert_allclose([0.000000 + -0.188937j,
        #                             -0.188937 + 0.000000j,
        #                             -0.000000 + 0.000000],
        #                            self.xdatcar_test.d_vectors[3][1])
        # np.testing.assert_allclose([-0.112619 - 0.065020j,
        #                             0.065021 - 0.112619j,
        #                             0.0 - 0.0j],
        #                            self.xdatcar_testy.d_vectors[5][1])



   
