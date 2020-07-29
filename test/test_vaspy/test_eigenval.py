# -*- coding: utf-8 -*-
"""Test for EIGENVAL class."""
import os
from typing import Callable

import pytest
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
import vaspy
import vaspy.eigenval as eigenval


datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"


@pytest.fixture
def cobalt() -> eigenval.EIGENVAL:
    cobalt = vaspy.load(datadir + "EIGENVAL.Co-collinear", mode="EIGENVAL")
    return cobalt


class TestEIGENVAL(object):
    """Class for EIGENVAL class test."""

    def setup_method(self) -> None:
        self.eigenval_spin = eigenval.EIGENVAL(datadir + "EIGENVAL.spin")
        self.eigenval_soi = eigenval.EIGENVAL(datadir + "EIGENVAL.soi")

    def test_check_basic_parameters(self) -> None:
        """Check the basic parameters stored."""
        assert 344 == self.eigenval_spin.natom
        assert 2 == self.eigenval_spin.nspin
        assert 1890 == self.eigenval_spin.nbands
        assert 1 == self.eigenval_spin.numk
        assert (2, 1, 1890) == self.eigenval_spin.energies.shape
        assert_array_equal([0.0, 0.0, 0.0], self.eigenval_spin.kvecs[0])
        assert_array_almost_equal(
            [[-22.893072, -22.891405, -22.841659, -22.832997, -22.761726, -22.737184]],
            self.eigenval_spin.energies[0, :, :6],
        )
        assert_array_equal([5.806516], self.eigenval_spin.energies[0, :, -1])
        assert_array_equal([5.813531], self.eigenval_spin.energies[1, :, -1])
        #
        assert 1 == self.eigenval_soi.natom
        assert 1 == self.eigenval_soi.nspin
        assert 54 == self.eigenval_soi.nbands
        assert 625 == self.eigenval_soi.numk
        assert_array_equal(
            [-9.548122, -9.491253, -9.322813, -9.050432, -8.691472, -8.285128],
            self.eigenval_soi.energies[0, :, 0][0:6],
        )
        assert_array_almost_equal(26.193038, self.eigenval_soi.energies[0, :, -1][-1])

    def test_to_physical_kvector(self) -> None:
        """Test for to_physical_kvector."""
        assert_array_almost_equal([[0.0, 0.0, 0.0]], self.eigenval_spin.kvecs)
        assert_array_almost_equal(
            [[0.0, 0.0, 0.0], [4.0e-02, 0, 0.0], [8.0e-02, 0, 0]],
            self.eigenval_soi.kvecs[0:3],
        )
        self.eigenval_soi.to_physical_kvector(
            recvec=np.array(((2, 0, 0), (0, 2, 0), (0, 0, 2)))
        )
        assert_array_almost_equal(
            [[0.0, 0.0, 0.0], [8.0e-02, 0, 0.0], [16.0e-02, 0, 0]],
            self.eigenval_soi.kvecs[0:3],
        )

    def test_make_label(self) -> None:
        """Test for make_label function."""
        labels_soi = self.eigenval_soi.make_label("k", "energy")
        assert "#k" == labels_soi[0]
        assert "Energy" == labels_soi[1]
        labels_spin = self.eigenval_spin.make_label("k", "energy")
        assert "#k" == labels_spin[0]
        assert "Energy_up" == labels_spin[1]
        assert "Energy_down" == labels_spin[2]

    def test__getitem_spi(self) -> None:
        assert_array_equal(self.eigenval_spin[0][0], [0.0, 0.0, 0.0])
        assert_array_equal(self.eigenval_soi[0][1][0], [-9.548122])
        assert_array_equal(
            self.eigenval_soi[0][1][:5],
            [[-9.548122], [-9.31537], [-6.490189], [-6.186601], [-6.034935],],
        )
        assert_array_equal(
            [x[0] for x in self.eigenval_soi[:3]],
            [
                [0.000000e00, 0.000000e00, 0.000000e00],
                [4.0e-2, 2.161984e-18, 0.0],
                [8.000000e-02, 4.323969e-18, 0.000000e00],
            ],
        )
        assert_array_equal(
            self.eigenval_spin[0][1][:5],
            [
                [-22.893072, -22.890912],
                [-22.891405, -22.889365],
                [-22.841659, -22.843597],
                [-22.832997, -22.830941],
                [-22.761726, -22.760568],
            ],
        )

    def test_cobalt(self, cobalt: eigenval.EIGENVAL):
        """Test EIGENVAL by using EIGENVAL.Co-collinear

        Collinear magnetism

        """

        assert_array_equal(cobalt[0][0], [0, 0, 0])
        assert_array_equal(cobalt[0][1][0], [-10.128206, -9.857677])
        assert_array_equal(cobalt[1][1][0], [-10.071635, -9.799951])
