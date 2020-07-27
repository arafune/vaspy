#!/usr/bin/python
# -*- coding: utf-8 -*-
import os

import numpy as np

import vaspy
import vaspy.chgcar

datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"

class TestCHGCAR(object):
    """Class for test of CHGCAR module.

    Use CHGCAR for test data."""

    def setup_method(self):
        data_file_spin = datadir + "CHGCAR_spin"
        data_file_spin_2 = datadir + "CHGCAR_spin_2"
        self.chgcar_spin = vaspy.load(data_file_spin)
        self.chgcar_spin_2 = vaspy.chgcar.CHGCAR(data_file_spin_2)

    def test_CHGsum(self):
        """Test for CHGCAR __add__ method."""
        add = self.chgcar_spin + self.chgcar_spin_2
        np.testing.assert_array_almost_equal(
            [9.33086564, 7.45015796, 3.32200981, -0.02888063, -0.83756195],
            add.grid.data[0:5],
        )
        np.testing.assert_array_equal([1, 1], add.poscar.atomnums)
        np.testing.assert_array_equal(["Co", "Co"], add.poscar.atomtypes)

    def test_CHGmerge(self):
        """Test for CHGCAR merge method."""
        merge = self.chgcar_spin.merge(self.chgcar_spin_2)
        np.testing.assert_array_almost_equal(
            [9.33086564, 7.45015796, 3.32200981, -0.02888063, -0.83756195],
            merge.grid.data[0:5],
        )
        np.testing.assert_array_equal([1], merge.poscar.atomnums)
        np.testing.assert_array_equal(["Co"], merge.poscar.atomtypes)

    def testCHGdiff(self):
        """Test for CHGCAR __sub__ method."""
        sub = self.chgcar_spin - self.chgcar_spin_2
        np.testing.assert_array_almost_equal(
            [-4.41988373, -3.52902219, -1.57358359, 0.0136803, 0.39673987],
            sub.grid.data[0:5],
        )
        np.testing.assert_array_equal([1], sub.poscar.atomnums)
        np.testing.assert_array_equal(["Co"], sub.poscar.atomtypes)

    def test_mag_spin(self):
        """Test for CHGCAR magnetization."""
        magCHG = self.chgcar_spin.magnetization()
        assert 4 * 4 * 6 == magCHG.grid.data.size
        np.testing.assert_array_almost_equal(
            [
                2.57589571650e00,
                2.09247290400e00,
                1.03267062280e00,
                1.74930752210e-01,
                -2.99056784380e-02,
            ],
            magCHG.grid.data[0:5],
        )
        assert ["up-down"] == magCHG.spin

    def test_mag_majority(self):
        """Test for CHGCAR majority spin."""
        majorspin = self.chgcar_spin.majorityspin()
        assert 4 * 4 * 6 == majorspin.grid.data.size
        np.testing.assert_array_almost_equal(
            [2.51569334, 2.02652039, 0.95344186, 0.08366529, -0.12515836],
            majorspin.grid.data[0:5],
        )

    def test_mag_minority(self):
        """Test for CHGCAR minority spin."""
        minorspin = self.chgcar_spin.minorityspin()
        assert 4 * 4 * 6 == minorspin.grid.data.size
        np.testing.assert_array_almost_equal(
            [-0.06020238, -0.06595251, -0.07922876, -0.09126546, -0.09525268],
            minorspin.grid.data[0:5],
        )

    def spindirec_CHGCAR_SOI(self):
        pass
