"""Unit test for OUTCAR class."""
import os

import numpy as np
import pytest

import vaspy


class TestOUTCAR:
    """class for test of vaspy.outcar module."""

    @pytest.fixture()
    def outcar0(self):
        """Reading OUTCAR file for test."""
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        datafile = datadir + "OUTCAR"
        return vaspy.load(datafile)

    def test_read_basic_properties(self, outcar0: vaspy.outcar.OUTCAR):
        """Test for OUTCAR reading basic properties."""
        assert outcar0.fermi == -0.7681
        assert outcar0.n_bands == 54
        #
        # reciprocal vectors
        np.testing.assert_almost_equal(
            np.array([0.389375981, -0.224806327, 0.000000000]), outcar0.recvec[0],
        )
        np.testing.assert_almost_equal(
            np.array([0.000000000, 0.449612655, 0.000000000]), outcar0.recvec[1],
        )
        np.testing.assert_almost_equal(
            np.array([0.000000000, 0.000000000, 0.023484312]), outcar0.recvec[2],
        )
        # k_vectors and weights
        assert len(outcar0.k_vectors) == 33
        assert len(outcar0.weights) == 33
        np.testing.assert_array_almost_equal(
            [
                [0.000000, 0.000000, 0.000000],
                [0.058824, 0.000000, 0.000000],
                [0.117647, 0.000000, 0.000000],
                [0.176471, -0.000000, 0.000000],
                [0.235294, 0.000000, 0.000000],
                [0.294118, 0.000000, 0.000000],
                [0.352941, -0.000000, 0.000000],
            ],
            outcar0.k_vectors[:7],
        )
        np.testing.assert_array_almost_equal(
            [1.000000, 6.000000, 6.000000, 6.000000, 6.000000, 6.000000, 6.000000],
            outcar0.weights[:7],
        )
        # total charges
        np.testing.assert_array_almost_equal(
            [
                [0.454, 0.282, 9.198],
                [0.446, 0.366, 9.192],
                [0.447, 0.372, 9.185],
                [0.451, 0.377, 9.184],
                [0.451, 0.381, 9.185],
                [0.450, 0.376, 9.192],
                [0.456, 0.288, 9.200],
            ],
            outcar0.total_charges[-1],
        )
        assert len(outcar0.total_charges) == 11
        assert len(outcar0.total_charges[0]) == 7
        assert len(outcar0.total_charges[0][0]) == 3
        np.testing.assert_almost_equal(-7.39238998, outcar0.totens[-1])
