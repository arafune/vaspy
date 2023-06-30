#! /usr/bin/env python
import os

import numpy as np

from vaspy import xdatcar


class TestXDATCAR:
    """Class for Test of Vsim_asc."""

    def setup_method(self, method):
        """XDATCAR."""
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        datafile = datadir + "XDATCAR.0.bz2"
        self.xdatcar_test = xdatcar.XDATCAR(datafile)

    def test_(self):
        assert self.xdatcar_test.system_name == "MoS2"
        assert self.xdatcar_test.scaling_factor == 1.0
        np.testing.assert_allclose(
            [3.184000, 0.000000, 0.000000], self.xdatcar_test.cell_vecs[0],
        )
        np.testing.assert_allclose(
            [1.592000, 2.757425, 0.000000], self.xdatcar_test.cell_vecs[1],
        )
        np.testing.assert_allclose([0, 0, 38.0], self.xdatcar_test.cell_vecs[2])
        assert ["Mo", "S"] == self.xdatcar_test.atom_types
        assert [1, 2] == self.xdatcar_test.atomnums
        assert len(self.xdatcar_test.configurations) == 5
        #
        np.testing.assert_allclose(
            [
                [0.33333333, 0.33333333, 0.24999815],
                [0.66666667, 0.66666667, 0.29113380],
                [0.66666667, 0.66666667, 0.20886348],
            ],
            self.xdatcar_test.configurations[0],
        )
        np.testing.assert_allclose(
            [
                [0.33333333, 0.33333333, 0.24999815],
                [0.66666667, 0.66666667, 0.29113854],
                [0.66666667, 0.66666667, 0.20885874],
            ],
            self.xdatcar_test.configurations[4],
        )

        # #  mode is '3', second atom
        # np.testing.assert_allclose([0.000000 + -0.188937j,
        #                             -0.188937 + 0.000000j,
        #                             -0.000000 + 0.000000],
        #                            self.xdatcar_test.d_vectors[3][1])
        # np.testing.assert_allclose([-0.112619 - 0.065020j,
        #                             0.065021 - 0.112619j,
        #                             0.0 - 0.0j],
        #                            self.xdatcar_testy.d_vectors[5][1])
