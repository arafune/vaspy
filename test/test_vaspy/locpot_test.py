#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os

import numpy as np
from nose.tools import assert_equal, assert_false, eq_, ok_, raises
from nose.tools import assert_raises

from vaspy.locpot import LOCPOT


class TestLOCPOT(object):
    def setup_method(self, method):
        """Setup: load data."""
        currentpath = os.path.abspath(os.path.dirname(__file__))
        fullpath_1 = currentpath + "/data/LOCPOT.dummy"
        fullpath_2 = currentpath + "/data/LOCPOT.dummy2"
        self.testlocpot1 = LOCPOT(fullpath_1)
        self.testlocpot2 = LOCPOT(fullpath_2)

    def test_locpot_poscar_part(self):
        """Test whether LOCPOT correctly treats POSCAR."""
        assert_equal("hBN-Cu", self.testlocpot1.poscar.system_name)
        len_x, len_y, len_z = self.testlocpot2.poscar.axes_lengthes
        assert_equal(
            (len_x, len_y, len_z),
            (6.7395854049999997, 6.7395858129882784, 36.271284721999997),
        )

    def test__str__(self):
        """Test for __str__() special method."""
        output = self.testlocpot2.__str__()
        teststr = "hBN-Cu\n  0.99700000000000\n    6.759865    0.000000"
        ok_(teststr in output)
        eq_(3363, len(output))

    def test___str__grid(self):
        """Test for __str__() special method in grid object."""
        output = self.testlocpot2.grid.__str__()
        teststr = "\n  3  4  6\n   0.00000000000E+00"
        ok_(teststr in output)
        eq_(1466, len(output))

    def test_grid_slice_z(self):
        """Test for grid slice (z)"""
        z0 = self.testlocpot2.grid.slice(0, "z")
        z1 = self.testlocpot2.grid.slice(1, "z")
        testarray2D = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]])
        np.testing.assert_equal(testarray2D, z0)
        np.testing.assert_equal(testarray2D + 12, z1)

    def test_grid_slice_y(self):
        """Test for grid slice (y)"""
        y0 = self.testlocpot2.grid.slice(0, "y")
        y1 = self.testlocpot2.grid.slice(1, "y")
        testarray2D = np.array(
            [
                [0, 1, 2],
                [12, 13, 14],
                [24, 25, 26],
                [36, 37, 38],
                [48, 49, 50],
                [60, 61, 62],
            ]
        )
        np.testing.assert_equal(testarray2D, y0)
        np.testing.assert_equal(testarray2D + 3, y1)

    def test_grid_slice_x(self):
        """Test for grid slice (x)"""
        x0 = self.testlocpot2.grid.slice(0, "x")
        x1 = self.testlocpot2.grid.slice(1, "x")
        testarray2D = np.array(
            [
                [0, 3, 6, 9],
                [12, 15, 18, 21],
                [24, 27, 30, 33],
                [36, 39, 42, 45],
                [48, 51, 54, 57],
                [60, 63, 66, 69],
            ]
        )
        np.testing.assert_equal(testarray2D, x0)
        np.testing.assert_equal(testarray2D + 1, x1)

    def test_grid_integrate_z(self):  # <<< FIXME
        """Test for grid integrate (z)"""
        z_0 = self.testlocpot2.grid.integrate("z", 0, 1)
        z_all = self.testlocpot2.grid.integrate("z")
        np.testing.assert_equal(self.testlocpot2.grid.slice(0, "z"), z_0)
        np.testing.assert_equal(
            self.testlocpot2.grid.slice(0, "z")
            + self.testlocpot2.grid.slice(1, "z")
            + self.testlocpot2.grid.slice(2, "z")
            + self.testlocpot2.grid.slice(3, "z")
            + self.testlocpot2.grid.slice(4, "z")
            + self.testlocpot2.grid.slice(5, "z"),
            z_all,
        )

    #        np.testing.assert_equal(testarray2D + 12, z1)

    def test_grid_integrate_y(self):
        """Test for grid integrate (y)"""
        y_0 = self.testlocpot2.grid.integrate("y", 0, 1)
        y_all = self.testlocpot2.grid.integrate("y")
        np.testing.assert_equal(self.testlocpot2.grid.slice(0, "y"), y_0)
        np.testing.assert_equal(
            self.testlocpot2.grid.slice(0, "y")
            + self.testlocpot2.grid.slice(1, "y")
            + self.testlocpot2.grid.slice(2, "y")
            + self.testlocpot2.grid.slice(3, "y"),
            y_all,
        )

    def test_grid_integrate_x(self):  # <<< FIXME
        """Test for grid integrate (x)"""
        x_0 = self.testlocpot2.grid.integrate("x", 0, 1)
        x_all = self.testlocpot2.grid.integrate("x")
        np.testing.assert_equal(self.testlocpot2.grid.slice(0, "x"), x_0)
        np.testing.assert_equal(
            self.testlocpot2.grid.slice(0, "x")
            + self.testlocpot2.grid.slice(1, "x")
            + self.testlocpot2.grid.slice(2, "x"),
            x_all,
        )

    def test_average_along(self):
        """Test LOCPOT averaging."""
        # along Z-axis
        np.testing.assert_equal(
            np.array([5.5, 17.5, 29.5, 41.5, 53.5, 65.5]),
            self.testlocpot2.grid.average_along_axis("z"),
        )
        # along X-axis
        np.testing.assert_almost_equal(
            np.array([34.5, 35.5, 36.5]), self.testlocpot2.grid.average_along_axis("x")
        )
        # along Y-axis
        np.testing.assert_equal(
            np.array([31.0, 34.0, 37.0, 40.0]),
            self.testlocpot2.grid.average_along_axis("Y"),
        )

    def test_min_along(self):
        """Test for LOCPOT min."""
        np.testing.assert_equal(
            np.array([0.0, 1.0, 2.0]), self.testlocpot2.grid.min_along_axis("x")
        )
        np.testing.assert_equal(
            np.array([0.0, 3.0, 6.0, 9.0]), self.testlocpot2.grid.min_along_axis("y")
        )
        np.testing.assert_equal(
            np.array([0.0, 12.0, 24.0, 36.0, 48.0, 60.0]),
            self.testlocpot2.grid.min_along_axis("z"),
        )

    def test_max_along(self):
        """Test for LOCPOT max."""
        np.testing.assert_equal(
            np.array([69.0, 70.0, 71.0]), self.testlocpot2.grid.max_along_axis("x")
        )
        np.testing.assert_equal(
            np.array([62.0, 65.0, 68.0, 71.0]),
            self.testlocpot2.grid.max_along_axis("y"),
        )
        np.testing.assert_equal(
            np.array([11.0, 23.0, 35.0, 47.0, 59.0, 71.0]),
            self.testlocpot2.grid.max_along_axis("z"),
        )

    def test_median_along(self):
        """Test for LOCPOT median."""
        np.testing.assert_equal(
            np.array([34.5, 35.5, 36.5]), self.testlocpot2.grid.median_along_axis("x")
        )
        np.testing.assert_equal(
            np.array([31.0, 34.0, 37.0, 40.0]),
            self.testlocpot2.grid.median_along_axis("y"),
        )
        np.testing.assert_equal(
            np.array([5.5, 17.5, 29.5, 41.5, 53.5, 65.5]),
            self.testlocpot2.grid.median_along_axis("z"),
        )

    @raises(TypeError, AssertionError, ValueError)
    def test_exception_in_locpot(self):
        """Test for Exception."""
        self.testlocpot2.grid.median_along_axis("k")
        self.testlocpot2.grid.max_along_axis("k")
        self.testlocpot2.grid.min_along_axis("k")
        self.testlocpot2.grid.average_along_axis("k")
