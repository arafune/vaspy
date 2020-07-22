# -*- coding: utf-8 -*-
"""Test for EIGENVAL class."""
import os

import numpy as np
from numpy.core.numeric import array_equal
from numpy.testing import assert_array_equal, assert_array_almost_equal
import vaspy.eigenval as eigenval


class TestEIGENVAL(object):
    """Class for EIGENVAL class test."""

    def setup_method(self) -> None:
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        self.eigenval_spin = eigenval.EIGENVAL(datadir + "EIGENVAL.spin")
        self.eigenval_soi = eigenval.EIGENVAL(datadir + "EIGENVAL.soi")
        self.cobalt_col = eigenval.EIGENVAL(datadir + "EIGENVAL.Co-collinear")

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
        np.testing.assert_array_almost_equal(
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
            self.eigenval_soi[0][1],
            [
                -9.548122,
                -9.31537,
                -6.490189,
                -6.186601,
                -6.034935,
                -5.296314,
                -5.206645,
                -4.810655,
                -4.443255,
                -4.306386,
                -3.301869,
                -3.295532,
                -3.105542,
                -3.024424,
                0.677385,
                0.788049,
                1.180932,
                1.401581,
                1.422082,
                1.604341,
                2.671816,
                2.855335,
                3.087653,
                3.232587,
                5.089931,
                5.271577,
                5.587036,
                5.709463,
                8.358908,
                8.522839,
                8.866504,
                9.009401,
                12.386523,
                12.575418,
                12.932918,
                13.054403,
                17.193601,
                17.376781,
                17.718552,
                17.840641,
                20.020032,
                20.379719,
                22.706773,
                22.879311,
                23.217494,
                23.332574,
                23.864843,
                24.073247,
                24.475428,
                24.673471,
                25.073306,
                25.122272,
                25.260491,
                25.311589,
            ],
        )
        assert_array_equal(
            self.eigenval_soi[:3][0],
            [
                [0.000000e00, 0.000000e00, 0.000000e00],
                [4.0e-2, 2.161984e-18, 0.0],
                [8.000000e-02, 4.323969e-18, 0.000000e00],
            ],
        )
        assert_array_equal(self.cobalt_col[0][0], [0, 0, 0])
        assert_array_equal(self.cobalt_col[0][1][0], [-10.128206, -9.857677])
        assert_array_equal(self.cobalt_col[1][1][0], [-10.071635, -9.799951])
