#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
# import tempfile
from nose.tools import eq_, ok_
from nose.tools import with_setup, assert_equal, raises
from nose.tools import assert_false
import numpy as np
import vaspy.vsim_asc as vsim_asc


class TestVsimASCII(object):
    '''Class for Test of Vsim_asc

'''
    def setup(self):
        '''VSIM'''
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        datafile = datadir + "monolayer_hBN_phonon.ascii"
        self.hBN = vsim_asc.VSIM_ASC(datafile)

    @with_setup(setup=setup)
    def test_(self):
        eq_('Phonopy generated file for v_sim 3.6', self.hBN.system_name)
        np.testing.assert_allclose([2.484999131, 0, 0],
                                   self.hBN.lattice_vector[0])
        np.testing.assert_allclose([1.242498697, 2.152072878, 0],
                                   self.hBN.lattice_vector[1])
        np.testing.assert_allclose([0, 0, 24.85], 
                                   self.hBN.lattice_vector[2])
        eq_(['B', 'N'], self.hBN.ions)
        #
        #
        np.testing.assert_allclose([9.354567, 17.815346, 26.162733,
                                   32.356041, 35.832499, 39.191721],
                                   self.hBN.freqs)
        #  mode is '3', second atom
        np.testing.assert_allclose([0.000000 + -0.188937j,
                                    -0.188937 + 0.000000j,
                                    -0.000000 + 0.000000],
                                   self.hBN.d_vectors[3][1])
        np.testing.assert_allclose([-0.112619 -0.065020j,
                                    0.065021 - 0.112619j,
                                    0.0 -0.0j],
                                   self.hBN.d_vectors[5][1])

