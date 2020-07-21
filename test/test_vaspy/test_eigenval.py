# -*- coding: utf-8 -*-
"""Test for EIGENVAL class."""
import os

import numpy as np

import vaspy.eigenval as eigenval


class TestEIGENVAL(object):
    """Class for EIGENVAL class test."""

    def setup_method(self) -> None:
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        self.eigenval_spin = eigenval.EIGENVAL(datadir + "EIGENVAL.spin")
        self.eigenval_soi = eigenval.EIGENVAL(datadir + "EIGENVAL.soi")

    def test_check_basic_parameters(self) -> None:
        """Check the basic parameters stored."""
        assert 344 == self.eigenval_spin.natom
        assert 2 == self.eigenval_spin.nspin
        assert 1890 == self.eigenval_spin.nbands
        assert 1 == self.eigenval_spin.numk
        assert (2, 1, 1890) == self.eigenval_spin.energies.shape
        np.testing.assert_array_equal([0.0, 0.0, 0.0], self.eigenval_spin.kvecs[0])
        np.testing.assert_array_almost_equal(
            [[-22.893072, -22.891405, -22.841659, -22.832997, -22.761726, -22.737184]],
            self.eigenval_spin.energies[0, :, :6],
        )
        np.testing.assert_array_equal([5.806516], self.eigenval_spin.energies[0, :, -1])
        np.testing.assert_array_equal([5.813531], self.eigenval_spin.energies[1, :, -1])
        #
        assert 1 == self.eigenval_soi.natom
        assert 1 == self.eigenval_soi.nspin
        assert 54 == self.eigenval_soi.nbands
        assert 625 == self.eigenval_soi.numk
        np.testing.assert_array_equal(
            [-9.548122, -9.491253, -9.322813, -9.050432, -8.691472, -8.285128],
            self.eigenval_soi.energies[0, :, 0][0:6],
        )
        np.testing.assert_array_almost_equal(
            26.193038, self.eigenval_soi.energies[0, :, -1][-1]
        )

    def test_to_physical_kvector(self) -> None:
        """Test for to_physical_kvector."""
        np.testing.assert_array_almost_equal(
            [[0.0, 0.0, 0.0]], self.eigenval_spin.kvecs
        )
        np.testing.assert_array_almost_equal(
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
