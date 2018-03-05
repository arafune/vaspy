#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
from nose.tools import eq_
from nose.tools import ok_
from nose.tools import with_setup
from nose.tools import assert_equal
from nose.tools import raises
import vaspy.chgcar
import numpy as np


class TestCHGCAR(object):
    '''Class for test of CHGCAR module

    Use CHGCAR for test data'''
    def setup(self):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file_spin = datadir + 'CHGCAR_spin'
        data_file_spin_2 = datadir + 'CHGCAR_spin_2'
        self.chgcar_spin = vaspy.chgcar.CHGCAR(data_file_spin)
        self.chgcar_spin_2 = vaspy.chgcar.CHGCAR(data_file_spin_2)

    @with_setup(setup=setup)
    def test_CHGsum(self):
        '''test for CHGCAR __add__ method'''
        add = self.chgcar_spin + self.chgcar_spin_2
        np.testing.assert_array_almost_equal(
            [9.33086564, 7.45015796, 3.32200981, -0.02888063, -0.83756195],
            add.grid.data[0:5])
        np.testing.assert_array_equal([1, 1], add.poscar.ionnums)
        np.testing.assert_array_equal(['Co', 'Co'], add.poscar.iontypes)

    @with_setup(setup=setup)
    def test_CHGmerge(self):
        '''test for CHGCAR merge method'''
        merge = self.chgcar_spin.merge(self.chgcar_spin_2)
        np.testing.assert_array_almost_equal(
            [9.33086564, 7.45015796, 3.32200981, -0.02888063, -0.83756195],
            merge.grid.data[0:5])
        np.testing.assert_array_equal([1], merge.poscar.ionnums)
        np.testing.assert_array_equal(['Co'], merge.poscar.iontypes)

    @with_setup(setup=setup)
    def testCHGdiff(self):
        '''test for CHGCAR __sub__ method'''
        sub = self.chgcar_spin - self.chgcar_spin_2
        np.testing.assert_array_almost_equal(
            [-4.41988373, -3.52902219, -1.57358359, 0.0136803, 0.39673987],
            sub.grid.data[0:5])
        np.testing.assert_array_equal([1], sub.poscar.ionnums)
        np.testing.assert_array_equal(['Co'], sub.poscar.iontypes)

    @with_setup(setup=setup)
    def test_mag_spin(self):
        '''test for CHGCAR magnetization
        '''
        magCHG = self.chgcar_spin.magnetization()
        eq_(4*4*6, magCHG.grid.data.size)
        np.testing.assert_array_almost_equal(
            [2.57589571650E+00, 2.09247290400E+00, 1.03267062280E+00,
             1.74930752210E-01, -2.99056784380E-02],
            magCHG.grid.data[0:5])
        eq_(["up-down"], magCHG.spin)

    @with_setup(setup=setup)
    def test_mag_majority(self):
        '''CHGCAR majority spin'''
        majorspin = self.chgcar_spin.majorityspin()
        eq_(4*4*6, majorspin.grid.data.size)
        np.testing.assert_array_almost_equal(
            [2.51569334, 2.02652039, 0.95344186,
             0.08366529, -0.12515836],
            majorspin.grid.data[0:5])

    @with_setup(setup=setup)
    def test_mag_minority(self):
        '''CHGCAR minority spin'''
        minorspin = self.chgcar_spin.minorityspin()
        eq_(4*4*6, minorspin.grid.data.size)
        np.testing.assert_array_almost_equal(
            [-0.06020238, -0.06595251, -0.07922876,
             -0.09126546, -0.09525268],
            minorspin.grid.data[0:5])

    @with_setup
    def spindirec_CHGCAR_SOI(self):
        pass
