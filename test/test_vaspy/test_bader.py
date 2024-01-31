"""Unit test for OUTCAR class."""

import os

import numpy as np

import vaspy


class TestBader:
    """Class for test of vaspy.outcar module."""

    def setup_method(self):
        """Reading ACF.chg.dat for test."""
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        datafile = datadir + "ACF.chg.dat"
        self.acf = vaspy.load(datafile)

    def test_read_basic_properties(self):
        """Test for reading basic properties."""
        """
        VACUUM CHARGE:               1.2829
        VACUUM VOLUME:             157.7708
        NUMBER OF ELECTRONS:      2828.0000
        """
        assert self.acf.vac_charge == 1.2829
        assert self.acf.vac_vol == 157.7708
        assert self.acf.n_electron == 2828.000
        assert 344, self.acf.n_atom
        np.testing.assert_array_almost_equal(
            [6.745000, -37.491966, 11.674620],
            self.acf.positions[0],
        )
        assert self.acf.charges[0] == 3.526327
        assert self.acf.mindists[0] == 0.423581
        assert self.acf.vols[0] == 17.540723
