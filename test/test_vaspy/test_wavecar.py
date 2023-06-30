#!/usr/bin/env python
"""Test for WAVECAR class."""
import os

import numpy as np
from numpy.testing import assert_almost_equal

import vaspy
from vaspy import poscar, wavecar


class TestHatomWavecar:
    """Class for Test WAVECAR module by using H_gamma.wavecar."""

    def setup_method(self, method):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + "H_gamma.wavecar"
        self.h = vaspy.load(data_file)

    def test_wavecar_header_and_band(self):
        """Test for H atom WAVECAR header."""
        self.h.gamma
        self.h.check_DwNGZHalf()
        assert self.h.n_spin == 1
        assert np.complex64 == self.h.prec
        assert self.h.num_k == 1
        assert self.h.n_bands == 27
        assert self.h.encut == 600
        np.testing.assert_array_equal(
            [[25.0, 0.0, 0.0], [0.0, 25.0, 0.0], [0.0, 0.0, 25.0]], self.h.realcell,
        )
        assert_almost_equal(25**3, self.h.volume)
        np.testing.assert_array_almost_equal([101, 101, 101], self.h.ngrid)
        #
        self.h.nplwvs[0] == 260834
        np.testing.assert_array_equal((260834, 3), self.h.gvectors().shape)
        np.testing.assert_array_almost_equal(
            [0, 0, 0, 0, 0], self.h.realspace_wfc()[0][0][:5],
        )


class TestCOWavecar:
    """Class for Test WAVECAR module by using WAVECAR.CO.wavecar."""

    def setup_method(self):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + "WAVECAR.CO.wavecar"
        self.co = wavecar.WAVECAR(data_file)

    def test_wavecar_header(self):
        """Test for CO molecule WAVECAR property.

        rtag, wfprec
        """
        assert self.co.recl == 711680  # record length
        assert self.co.n_spin == 2  # spin
        assert self.co.rtag == 45200  # precision flag
        assert np.complex64 == self.co.prec

        assert self.co.num_k == 1
        assert self.co.n_bands == 54
        assert self.co.encut == 400
        assert_almost_equal(-8.05748789, self.co.efermi)
        np.testing.assert_array_equal(
            [[17.0, 0.0, 0.0], [0.0, 17.0, 0.0], [0.0, 0.0, 17.0]], self.co.realcell,
        )
        assert_almost_equal(4913, self.co.volume)
        np.testing.assert_array_almost_equal([57, 57, 57], self.co.ngrid)

    def test_wavecar_band(self):
        """Test for CO WAVECAR band."""
        self.co.band()
        assert (self.co.n_spin, self.co.num_k, self.co.n_bands) == self.co.bands.shape
        np.testing.assert_array_almost_equal(
            [-29.49120151, -14.03737717, -11.88106463, -11.88106223, -9.00976389],
            self.co.bands[0, 0, 0:5],
        )  # The value can be taken from EIGENVAL

    def test_wavecar_str(self):
        """Test for __str__ special method."""
        output = self.co.__str__()
        teststr = "record length  =       711680  spins =           2  "
        teststr += "prec flag        45200\nno. k points =          1"
        teststr += "\nno. bands =          54\nreal space lattice vectors:\n"
        teststr += "a1 = 17.0    0.0    0.0\na2 = 0.0    17.0    0.0\n"
        teststr += "a3 = 0.0    0.0    17.0\n\nvolume unit cell =   "
        teststr += "4913.0000000000055\nReciprocal lattice vectors:\n"
        teststr += "b1 = 0.058823529411764705    0.0    0.0"
        teststr += "\nb2 = 0.0    0.058823529411764705"
        teststr += "    0.0\nb3 = 0.0    0.0    0.058823529411764705"
        assert teststr == output


class TestGrapheneWavecar:
    """Class for Test WAVECAR module by using Graphene.wavecar."""

    def setup_method(self):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + "Graphene.wavecar"
        self.gr = wavecar.WAVECAR(data_file)
        self.gr_poscar = poscar.POSCAR(datadir + "POSCAR.Graphene")

    def test_wavecar_header(self):
        """Test for Graphene WAVECAR property."""
        assert self.gr.recl == 30864  # record length
        assert self.gr.n_spin == 1  # spin
        assert self.gr.rtag == 45200  # precision flag
        assert np.complex64 == self.gr.prec
        assert self.gr.num_k == 240
        assert self.gr.n_bands == 27
        assert self.gr.encut == 550
        assert_almost_equal(-2.9500542715, self.gr.efermi)
        np.testing.assert_array_almost_equal(
            [
                [2.139081750, -1.235000000, 0.0],
                [2.139081750, 1.235000000, 0.0],
                [0.0, 0.0, 24.7],
            ],
            self.gr.realcell,
        )
        # volume of cell in OUTCAR is not enough
        assert_almost_equal(130.50323848575005, self.gr.volume)
        # FIXME!!: Where these value come from ?
        np.testing.assert_array_almost_equal([11, 11, 97], self.gr.ngrid)

    def test_wavecar_band(self):
        """Test for Graphene wavecar band."""
        self.gr.band()
        # from OUTCAR
        # maximum and minimum number of plane-waves per node :  3857  3785
        assert self.gr.nplwvs.max() == 3857
        assert self.gr.nplwvs.min() == 3785
        assert (self.gr.num_k,) == self.gr.kpath.shape  # gr.num_k = 240
        assert (self.gr.n_spin, self.gr.num_k, self.gr.n_bands) == self.gr.bands.shape
        np.testing.assert_array_almost_equal(
            [-22.516876, -10.623282, -6.106901, -6.094072, 0.245639, 1.006991],
            self.gr.bands[0, 0, 0:6],
        )  # The value can be taken from EIGENVAL

    def test_realsapece_wfc(self):
        """Test for generation real space wfc (Graphene)."""
        np.testing.assert_array_almost_equal(
            [
                0.00013770 + 0.0001743j,
                0.00014605 + 0.00018611j,
                0.00017262 + 0.00022051j,
                0.00021561 + 0.00027499j,
                0.00026360 + 0.00033486j,
            ],
            self.gr.realspace_wfc()[0][0][:5],
        )
        vaspgrid = self.gr.realspace_wfc(poscar=self.gr_poscar)
        np.testing.assert_array_almost_equal(
            [0.00013770, 0.00014605, 0.00017262, 0.00021561, 0.00026360],
            vaspgrid.grid.data[:5],
        )
        assert vaspgrid.grid.n_frame == 2


class test_RestoreGammaGrid:
    #
    def test_check_symmetry(self):
        """Test for check_symmetry."""
        assert not wavecar.check_symmetry(np.arange(3 * 3 * 3).reshape((3, 3, 3)))
        k_grid = np.array(
            [
                [[0, 1, 1], [3, 4, 7], [3, 7, 4]],
                [[9, 10, 11], [12, 13, 14], [15, 16, 17]],
                [[9, 11, 10], [15, 17, 16], [12, 14, 13]],
            ],
        )
        assert wavecar.check_symmetry(k_grid)

    def test_RestorGammaGrid(self):
        """Test for RestoreGammaGrid function."""
        grid332 = np.arange(3 * 3 * 3).reshape((3, 3, 3))[:, :, :2]
        grid332[2][0][0] = np.conjugate(grid332[1][0][0])
        grid332[2][2][0] = np.conjugate(grid332[1][1][0])
        grid332[1][2][0] = np.conjugate(grid332[2][1][0])
        grid332[0][2][0] = np.conjugate(grid332[0][1][0])
        result = wavecar.restore_gamma_grid(grid332, para=True)
        assert wavecar.check_symmetry(result)
        #
        grid553 = np.arange(5 * 5 * 5).reshape((5, 5, 5))[:, :, :3]
        #
        grid553[3, 0, 0] = np.conjugate(grid553[2, 0, 0])
        grid553[4, 0, 0] = np.conjugate(grid553[1, 0, 0])
        #
        grid553[0, 3, 0] = np.conjugate(grid553[0, 2, 0])
        grid553[1, 3, 0] = np.conjugate(grid553[4, 2, 0])
        grid553[2, 3, 0] = np.conjugate(grid553[3, 2, 0])
        grid553[3, 3, 0] = np.conjugate(grid553[2, 2, 0])
        grid553[4, 3, 0] = np.conjugate(grid553[1, 2, 0])
        #
        grid553[0, 4, 0] = np.conjugate(grid553[0, 1, 0])
        grid553[1, 4, 0] = np.conjugate(grid553[4, 1, 0])
        grid553[2, 4, 0] = np.conjugate(grid553[3, 1, 0])
        grid553[3, 4, 0] = np.conjugate(grid553[2, 1, 0])
        grid553[4, 4, 0] = np.conjugate(grid553[1, 1, 0])
        result = wavecar.restore_gamma_grid(grid553, para=True)
        assert wavecar.check_symmetry(result)
        #
        grid355 = np.arange(5 * 5 * 5).reshape((5, 5, 5))[:3, :, :]
        #
        grid355[0, 0, 3] = np.conjugate(grid355[0, 0, 2])
        grid355[0, 0, 4] = np.conjugate(grid355[0, 0, 1])
        #
        grid355[0, 3, 0] = np.conjugate(grid355[0, 2, 0])
        grid355[0, 3, 1] = np.conjugate(grid355[0, 2, 4])
        grid355[0, 3, 2] = np.conjugate(grid355[0, 2, 3])
        grid355[0, 3, 3] = np.conjugate(grid355[0, 2, 2])
        grid355[0, 3, 4] = np.conjugate(grid355[0, 2, 1])
        #
        grid355[0, 4, 0] = np.conjugate(grid355[0, 1, 0])
        grid355[0, 4, 1] = np.conjugate(grid355[0, 1, 4])
        grid355[0, 4, 2] = np.conjugate(grid355[0, 1, 3])
        grid355[0, 4, 3] = np.conjugate(grid355[0, 1, 2])
        grid355[0, 4, 4] = np.conjugate(grid355[0, 1, 1])
        result = wavecar.restore_gamma_grid(grid355, para=False)
        assert wavecar.check_symmetry(result)


class TestCobaltWavecar:
    """Class for Test WAVECAR module by using Co.wavecar (SOI)."""

    def setup_method(self, method):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + "Co.wavecar"
        self.co = wavecar.WAVECAR(data_file)
        self.co_poscar = poscar.POSCAR(datadir + "Co.POSCAR")

    def test_wavecar_header(self):
        """Test for Cobalt property."""
        assert self.co.recl == 7968  # record length
        assert self.co.n_spin == 1  # spin
        assert self.co.rtag == 45200  # precision flag
        assert np.complex64 == self.co.prec
        assert self.co.num_k == 9
        assert self.co.n_bands == 54
        assert self.co.encut == 400
        assert_almost_equal(-0.84118407515959326, self.co.efermi)
        np.testing.assert_array_almost_equal(
            np.array([[1.0, 0.0, 0.0], [0.5, 0.8660254, 0.0], [0.0, 0.0, 2.0]]) * 2.501,
            self.co.realcell,
        )
        # volume of cell in OUTCAR is not enough
        assert_almost_equal(27.095782694613046, self.co.volume)
        # FIXME!!: Where these value come from ?
        #  Maximum number of reciprocal cells 2x 2x 4 (in OUTCAR)
        np.testing.assert_array_almost_equal([11, 11, 19], self.co.ngrid)

    def test_wavecar_band(self):
        """Test for Co wavecar band."""
        self.co.band()
        # from OUTCAR
        assert self.co.nplwvs.max() == 996
        assert self.co.nplwvs.min() == 958
        assert (self.co.num_k,) == self.co.kpath.shape  # co.num_k = 240
        assert (self.co.n_spin, self.co.num_k, self.co.n_bands) == self.co.bands.shape
        np.testing.assert_array_almost_equal(
            [-6.532492, -5.858599, -4.037263, -3.418892, -3.265558, -2.642231],
            self.co.bands[0, 0, 0:6],
        )  # The value can be taken from EIGENVAL

    def test_realsapece_wfc(self):
        """Test for generation real space wfc (Cobalt)."""
        np.testing.assert_array_almost_equal(
            [
                -7.84915157e-05 - 5.61362047e-05j,
                -7.88077954e-05 - 5.63624713e-05j,
                -7.92001487e-05 - 5.66431905e-05j,
                -7.92675956e-05 - 5.66915450e-05j,
                -7.90423140e-05 - 5.65305135e-05j,
            ],
            self.co.realspace_wfc()[0][0][0][:5],
        )
        vaspgrid = self.co.realspace_wfc(poscar=self.co_poscar)
        np.testing.assert_array_almost_equal(
            [
                -7.84915157e-05,
                -7.88077954e-05,
                -7.92001487e-05,
                -7.92675956e-05,
                -7.90423140e-05,
            ],
            vaspgrid.grid.data[:5],
        )
        assert vaspgrid.grid.n_frame == 4
