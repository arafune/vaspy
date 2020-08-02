#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test for PROCAR class"""
import os

import numpy as np

# import tempfile


import vaspy.procar as procar


class TestSinglePROCAR(object):
    """Class for Test of PROCAR module.

    Use PROCAR_single for test data (dummy PROCAR)

    """

    def setup_method(self, method):
        """PROCAR object by using file of "PROCAR_single" in data directory."""
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + "PROCAR_single"
        self.singleprocar = procar.PROCAR(data_file)

    def test_sigleprocar_firstcheck(self):
        """Load test for PROCAR_single."""
        assert [""] == self.singleprocar.label["spin"]
        assert 3 == self.singleprocar.natom
        assert 1 == self.singleprocar.nbands
        assert 1 == self.singleprocar.numk
        np.testing.assert_array_equal([[[-15.0]]], self.singleprocar.energies)
        np.testing.assert_array_equal(
            ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2", "tot"],
            self.singleprocar.label["orbital"],
        )
        np.testing.assert_array_equal([[0.0, 0.0, 0.0]], self.singleprocar.kvecs)

    def test_projection1(self):
        """Test Projection class."""
        # testarray is taken from PROCAR_single
        assert 5 == self.singleprocar.proj.ndim
        np.testing.assert_array_equal((1, 1, 1, 3, 10), self.singleprocar.proj.shape)
        np.testing.assert_array_equal(
            [
                0.0000,
                0.0001,
                0.0002,
                0.0003,
                0.0004,
                0.0005,
                0.0006,
                0.0007,
                0.0008,
                0.0036,
            ],
            self.singleprocar.proj[0, 0, 0, 0],
        )
        np.testing.assert_array_equal(
            [
                0.0010,
                0.0011,
                0.0012,
                0.0013,
                0.0014,
                0.0015,
                0.0016,
                0.0017,
                0.0018,
                0.0126,
            ],
            self.singleprocar.proj[0, 0, 0, 1],
        )
        np.testing.assert_array_equal(
            [
                0.0020,
                0.0021,
                0.0022,
                0.0023,
                0.0024,
                0.0025,
                0.0026,
                0.0027,
                0.0028,
                0.0216,
            ],
            self.singleprocar.proj[0, 0, 0, 2],
        )

    def test_singleprocar_fermi_correction(self):
        """Test for fermi_correction."""
        self.singleprocar.fermi_correction(1.0)
        np.testing.assert_array_equal(self.singleprocar.energies, [[[-16.0]]])

    def test_singleprocar_sum_site(self):
        """Test for append_sumsite  (single)."""
        self.singleprocar.append_sumsite((0, 2), "zero-two")
        assert self.singleprocar.label["site"][-1] == "zero-two"
        np.testing.assert_array_almost_equal(
            self.singleprocar.proj[0, 0, 0, 3],
            [
                0.0020,
                0.0022,
                0.0024,
                0.0026,
                0.0028,
                0.0030,
                0.0032,
                0.0034,
                0.0036,
                0.0252,
            ],
        )

    def test_singleprocar_sum_orbital(self):
        """Test for append_sumorbital  (single)."""
        self.singleprocar.append_sumsite((0, 2), "zero-two")
        self.singleprocar.append_sumorbital((1, 3), "pxpy")
        self.singleprocar.label["orbital"][-1] == "pxpy"
        np.testing.assert_array_almost_equal(
            self.singleprocar.proj[0, 0, 0, 3],
            [
                0.0020,
                0.0022,
                0.0024,
                0.0026,
                0.0028,
                0.0030,
                0.0032,
                0.0034,
                0.0036,
                0.0252,
                0.0048,
            ],
        )

    def test_make_label(self):
        """Test for make label from PROCAR_single w/o append_sum* ."""
        label = self.singleprocar.make_label((0, 2), ((0, 3, 1), (4, 7)))
        assert label == ["#k", "Energy", "0_s", "0_px", "0_py", "2_dxy", "2_dxz"]

    def test_singleprocar_sum_site_orbital(self):
        """Testing: sum by site and then sum by orbital (PROCAR_single).

        Also test for orbital_index

        """
        self.singleprocar.append_sumsite((0, 2), site_name="zero_two")
        p = self.singleprocar.orbital_index("p")
        self.singleprocar.append_sumorbital(p, "p")
        pxpy = self.singleprocar.orbital_index("pxpy")
        self.singleprocar.append_sumorbital(pxpy, "pxpy")
        d = self.singleprocar.orbital_index("d")
        self.singleprocar.append_sumorbital(d, "d")
        np.testing.assert_allclose(
            self.singleprocar.proj[0, 0, 0, 3],
            [
                0.0020,
                0.0022,
                0.0024,
                0.0026,
                0.0028,
                0.0030,
                0.0032,
                0.0034,
                0.0036,
                0.0252,
                0.0072,
                0.0048,
                0.0160,
            ],
        )
        assert self.singleprocar.label["orbital"] == [
            "s",
            "py",
            "pz",
            "px",
            "dxy",
            "dyz",
            "dz2",
            "dxz",
            "dx2",
            "tot",
            "p",
            "pxpy",
            "d",
        ]

        assert self.singleprocar.label["site"][-1] == "zero_two"
        assert 4 == len(self.singleprocar.label["site"])

    def test_singleprocar_setheader(self):
        """Test for Band_with_projection.set_header."""
        self.singleprocar.append_sumsite((0, 2), "zero_two")
        for orbital in ("p", "pxpy", "d"):
            self.singleprocar.append_sumorbital(
                self.singleprocar.orbital_index(orbital), orbital
            )
        assert self.singleprocar.make_label((3,), ((10, 11, 12),)) == [
            "#k",
            "Energy",
            "zero_two_p",
            "zero_two_pxpy",
            "zero_two_d",
        ]

        # same as above

        assert self.singleprocar.make_label(
            (3,), ((self.singleprocar.orbital_index(o) for o in ("p", "pxpy", "d")),),
        ) == ["#k", "Energy", "zero_two_p", "zero_two_pxpy", "zero_two_d"]

    def test_text_sheet(self):
        """Test for text output used for Igor or gnuplot."""
        assert (
            self.singleprocar.text_sheet()
            == """#k	Energy
 0.00000000e+00	-1.50000000e+01

"""
        )

        self.singleprocar.append_sumsite((0, 2), "zero_two")
        for orbital in ("p", "pxpy", "d"):
            self.singleprocar.append_sumorbital(
                self.singleprocar.orbital_index(orbital), orbital
            )
        assert self.singleprocar.make_label((3,), ((10, 11, 12),)) == [
            "#k",
            "Energy",
            "zero_two_p",
            "zero_two_pxpy",
            "zero_two_d",
        ]
        assert (
            self.singleprocar.text_sheet((3,), ((10, 11, 12),))
            == """#k	Energy	zero_two_p	zero_two_pxpy	zero_two_d
 0.00000000e+00	-1.50000000e+01	7.20000000e-03	4.80000000e-03	1.60000000e-02

"""
        )


# # ------------------------------


class TestSpinPolarizedPROCAR(object):
    """Class for Test of PROCAR module.

    Use PROCAR_spin_dummy for test data (dummy PROCAR)"""

    def setup_method(self, method):
        """PROCAR object by using file of "PROCAR_single" in data directory."""
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + "PROCAR_spin_dummy"
        self.spinprocar = procar.PROCAR(data_file)

    def test_spinprocar_firstcheck(self):
        """Load test for PROCAR_spin_dummy."""
        assert ["_up", "_down"] == self.spinprocar.label["spin"]
        assert 3 == self.spinprocar.natom
        assert 4 == self.spinprocar.nbands
        assert 3 == self.spinprocar.numk
        np.testing.assert_array_almost_equal(
            np.array(
                [
                    [
                        [-10.0, -5.0, 0.0, 5.0],
                        [-7.0, -4.0, -1.0, 4.0],
                        [-6.0, -1.0, -3.0, 0.0],
                    ],
                    [
                        [-10.5, -5.5, -0.5, -5.5],
                        [-7.5, -4.5, -1.5, -4.5],
                        [-6.5, -1.5, -3.5, -0.5],
                    ],
                ]
            ),
            self.spinprocar.energies,
        )
        np.testing.assert_array_equal(
            ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2", "tot"],
            self.spinprocar.label["orbital"],
        )
        np.testing.assert_array_equal(
            [[0.0, 0.0, 0.0], [0.25, 0.25, 0.00], [0.50, 0.50, 0.00]],
            self.spinprocar.kvecs,
        )
        np.testing.assert_array_equal(720, self.spinprocar.proj.size)
        np.testing.assert_array_equal(
            [
                0.0000,
                0.0001,
                0.0002,
                0.0003,
                0.0004,
                0.0005,
                0.0006,
                0.0007,
                0.0008,
                0.0036,
            ],
            self.spinprocar.proj[0, 0, 0, 0],
        )
        np.testing.assert_array_equal(
            [
                0.0010,
                0.0011,
                0.0012,
                0.0013,
                0.0014,
                0.0015,
                0.0016,
                0.0017,
                0.0018,
                0.0126,
            ],
            self.spinprocar.proj[0, 0, 0, 1],
        )
        np.testing.assert_array_equal(
            [
                0.0020,
                0.0021,
                0.0022,
                0.0023,
                0.0024,
                0.0025,
                0.0026,
                0.0027,
                0.0028,
                0.0216,
            ],
            self.spinprocar.proj[0, 0, 0, 2],
        )

    def test_make_label(self):
        """test for make label from PROCAR_spin w/o append_sum*."""
        label = self.spinprocar.make_label((0, 2), ((0, 3, 1), (4, 7)))
        assert label == [
            "#k",
            "Energy_up",
            "Energy_down",
            "0_up_s",
            "0_down_s",
            "0_up_px",
            "0_down_px",
            "0_up_py",
            "0_down_py",
            "2_up_dxy",
            "2_down_dxy",
            "2_up_dxz",
            "2_down_dxz",
        ]

    def test_spinprocar_band(self):
        """Band_with_projection object test generated by PROCAR_spin_dummy."""
        np.testing.assert_array_almost_equal(
            self.spinprocar.kdistances, [0.0, 0.353553, 0.707107]
        )
        """test for Band_with_projection.energies setter (SPIN)."""
        assert self.spinprocar.energies.shape == (2, 3, 4)
        np.testing.assert_array_equal(
            self.spinprocar.energies,
            np.array(
                [
                    [
                        [-10.0, -5.0, 0.0, 5.0],
                        [-7.0, -4.0, -1.0, 4.0],
                        [-6.0, -1.0, -3.0, 0.0],
                    ],
                    [
                        [-10.5, -5.5, -0.5, -5.5],
                        [-7.5, -4.5, -1.5, -4.5],
                        [-6.5, -1.5, -3.5, -0.5],
                    ],
                ]
            ),
        )

    def test_spinprocar_fermi_correction(self):
        """test for fermi_correction (SPIN).
        """
        self.spinprocar.fermi_correction(1.0)
        np.testing.assert_array_equal(
            self.spinprocar.energies,
            np.array(
                [
                    [
                        [-11.0, -6.0, -1.0, 4.0],
                        [-8.0, -5.0, -2.0, 3.0],
                        [-7.0, -2.0, -4.0, -1.0],
                    ],
                    [
                        [-11.5, -6.5, -1.5, -6.5],
                        [-8.5, -5.5, -2.5, -5.5],
                        [-7.5, -2.5, -4.5, -1.5],
                    ],
                ]
            ),
        )

    def test_spinprocar_band_orbitalread(self):
        """test for Band_with_projection.orbitals setter (SPIN)."""
        assert self.spinprocar.proj[0].shape == (3, 4, 3, 10)
        assert self.spinprocar.proj[1].shape == (3, 4, 3, 10)
        # for up spin, ik = 1, ib = 0, iatom = 2
        #            ->  (In PROCAR, k# =2, band# = 1, atom# = 3)
        np.testing.assert_array_almost_equal(
            self.spinprocar.proj[0][1][0][2],
            [
                0.0080,
                0.0081,
                0.0082,
                0.0083,
                0.0084,
                0.0085,
                0.0086,
                0.0087,
                0.0088,
                0.0756,
            ],
        )

    def test_spinprocar_band_sum_orbital1(self):
        """test for sum_orbital (SPIN).

        raise RuntimeError when no item in sitecomposed."""
        self.spinprocar.append_sumsite((0, 2), "zero_two")
        for orbital in ("p", "pxpy", "d"):
            self.spinprocar.append_sumorbital(
                self.spinprocar.orbital_index(orbital), orbital
            )
        np.testing.assert_allclose(
            self.spinprocar.proj[0][0][0][3],
            [
                0.0020,
                0.0022,
                0.0024,
                0.0026,
                0.0028,
                0.0030,
                0.0032,
                0.0034,
                0.0036,
                0.0252,
                0.0072,
                0.0048,
                0.0160,
            ],
        )
        np.testing.assert_allclose(
            self.spinprocar.proj[1][0][0][3],
            [
                2.0020,
                2.0022,
                2.0024,
                2.0026,
                2.0028,
                2.0030,
                2.0032,
                2.0034,
                2.0036,
                18.0252,
                6.0072,
                4.0048,
                10.0160,
            ],
        )
        assert self.spinprocar.label["orbital"] == [
            "s",
            "py",
            "pz",
            "px",
            "dxy",
            "dyz",
            "dz2",
            "dxz",
            "dx2",
            "tot",
            "p",
            "pxpy",
            "d",
        ]

    def test_spinprocar_make_label(self):
        """test for Band_with_projection.set_header  (SPIN).
        """
        self.spinprocar.append_sumsite((0, 2), "test")
        for orbital in ("p", "pxpy", "d"):
            self.spinprocar.append_sumorbital(
                self.spinprocar.orbital_index(orbital), orbital
            )
        assert self.spinprocar.make_label((3,), (((10,), (11,), (12,)),)) == [
            "#k",
            "Energy_up",
            "Energy_down",
            "test_up_p",
            "test_down_p",
            "test_up_pxpy",
            "test_down_pxpy",
            "test_up_d",
            "test_down_d",
        ]

    def test_text_sheet(self):
        """test for simple band data output (Spin polarized)."""
        assert (
            self.spinprocar.text_sheet()
            == """#k	Energy_up	Energy_down
 0.00000000e+00	-1.00000000e+01	-1.05000000e+01
 3.53553391e-01	-7.00000000e+00	-7.50000000e+00
 7.07106781e-01	-6.00000000e+00	-6.50000000e+00

 0.00000000e+00	-5.00000000e+00	-5.50000000e+00
 3.53553391e-01	-4.00000000e+00	-4.50000000e+00
 7.07106781e-01	-1.00000000e+00	-1.50000000e+00

 0.00000000e+00	0.00000000e+00	-5.00000000e-01
 3.53553391e-01	-1.00000000e+00	-1.50000000e+00
 7.07106781e-01	-3.00000000e+00	-3.50000000e+00

 0.00000000e+00	5.00000000e+00	-5.50000000e+00
 3.53553391e-01	4.00000000e+00	-4.50000000e+00
 7.07106781e-01	0.00000000e+00	-5.00000000e-01

"""
        )


# # -------------------------


class TestSOIPROCAR(object):
    """Class for Test of PROCAR module.

    Use PROCAR_SOI_dummy for test data (dummy PROCAR)"""

    def setup(self):
        """PROCAR object by using file of "PROCAR_single" in data directory"""
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + "PROCAR_SOI_dummy"
        self.soiprocar = procar.PROCAR(data_file)

    def test_soiprocar_firstcheck(self):
        """Load test for PROCAR_SOI_dummy."""
        assert ["_mT", "_mX", "_mY", "_mZ"] == self.soiprocar.label["spin"]
        assert 3 == self.soiprocar.natom
        assert 2 == self.soiprocar.nbands
        assert 3 == self.soiprocar.numk
        np.testing.assert_array_equal(
            np.array([[[-10.0, -5.0], [-7.0, -4.0], [-6.0, -1.0]]]),
            self.soiprocar.energies,
        )
        np.testing.assert_array_equal(
            ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2", "tot"],
            self.soiprocar.label["orbital"],
        )
        np.testing.assert_array_equal(
            [[0.0, 0.0, 0.0], [0.25, 0.25, 0.00], [0.50, 0.50, 0.00]],
            self.soiprocar.kvecs,
        )
        # 72 = natom * nbands * numk * 4
        np.testing.assert_array_equal(720, self.soiprocar.proj.size)
        np.testing.assert_array_equal(
            [
                0.0000,
                0.0001,
                0.0002,
                0.0003,
                0.0004,
                0.0005,
                0.0006,
                0.0007,
                0.0008,
                0.0036,
            ],
            self.soiprocar.proj[0, 0, 0, 0],
        )
        np.testing.assert_array_equal(
            [
                0.0010,
                0.0011,
                0.0012,
                0.0013,
                0.0014,
                0.0015,
                0.0016,
                0.0017,
                0.0018,
                0.0126,
            ],
            self.soiprocar.proj[0, 0, 0, 1],
        )
        np.testing.assert_array_equal(
            [
                0.0020,
                0.0021,
                0.0022,
                0.0023,
                0.0024,
                0.0025,
                0.0026,
                0.0027,
                0.0028,
                0.0216,
            ],
            self.soiprocar.proj[0, 0, 0, 2],
        )

    def test_make_label(self):
        """test for make label from PROCAR_soi w/o append_sum*
        """
        label = self.soiprocar.make_label((0, 2), ((0, 3, 1), (4, 7)))
        assert label == [
            "#k",
            "Energy",
            "0_mT_s",
            "0_mX_s",
            "0_mY_s",
            "0_mZ_s",
            "0_mT_px",
            "0_mX_px",
            "0_mY_px",
            "0_mZ_px",
            "0_mT_py",
            "0_mX_py",
            "0_mY_py",
            "0_mZ_py",
            "2_mT_dxy",
            "2_mX_dxy",
            "2_mY_dxy",
            "2_mZ_dxy",
            "2_mT_dxz",
            "2_mX_dxz",
            "2_mY_dxz",
            "2_mZ_dxz",
        ]

    #            ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot'],

    def test_soiprocar_band(self):
        """Band_with_projection object test generated from PROCAR_soi_dummy"""
        np.testing.assert_array_almost_equal(
            self.soiprocar.kdistances, [0.0, 0.353553, 0.707107]
        )
        """test for Band_with_projection.energies setter (SOI)"""
        assert self.soiprocar.energies.shape == (
            1,
            self.soiprocar.numk,
            self.soiprocar.nbands,
        )
        np.testing.assert_array_equal(
            self.soiprocar.energies, [[[-10, -5], [-7, -4], [-6, -1]]]
        )
        """test for Band_with_projection.fermi_correction
        """
        self.soiprocar.fermi_correction(1.0)
        np.testing.assert_array_equal(
            self.soiprocar.energies, [[[-11, -6], [-8, -5], [-7, -2]]]
        )

    def test_soiprocar_band_orbitalread(self):
        """test for Band_with_projection.orbitals setter (SOI)"""
        assert self.soiprocar.proj.shape == (
            4,
            self.soiprocar.numk,
            self.soiprocar.nbands,
            self.soiprocar.natom,
            10,
        )

        # for ik = 0, ib = 1, atom=2, spin=mY,
        #                      (k#=1,  band#=2, atom#=3, spin=mY)
        np.testing.assert_array_equal(
            self.soiprocar.proj[2][0][1][2],
            #                   ^ This two means "mY"
            [
                2.0050,
                2.0051,
                2.0052,
                2.0053,
                2.0054,
                2.0055,
                2.0056,
                2.0057,
                2.0058,
                2.0486,
            ],
        )
        # for ik = 1, ib = 0, atom=1, spin=mZ,
        #                      (k#=2,  band#=1, atom#=2, spin=mZ)
        np.testing.assert_array_equal(
            self.soiprocar.proj[3][1][0][1],
            #                   ^ This two means "mY"
            [
                3.0070,
                3.0071,
                3.0072,
                3.0073,
                3.0074,
                3.0075,
                3.0076,
                3.0077,
                3.0078,
                3.0666,
            ],
        )

    def test_soiprocar_band_sum_site(self):
        """test for Band_with_projection.sum_site (SOI)"""
        self.soiprocar.append_sumsite((0, 2), "test")
        assert self.soiprocar.label["site"].index("test") == 3
        np.testing.assert_allclose(
            self.soiprocar.proj[0][0][0][self.soiprocar.label["site"].index("test")],
            [
                0.0020,
                0.0022,
                0.0024,
                0.0026,
                0.0028,
                0.0030,
                0.0032,
                0.0034,
                0.0036,
                0.0252,
            ],
        )

    def test_soiprocar_band_sum_orbital1(self):
        """test for Band_with_projection.sum_orbital (SOI)

        raise RuntimeError when no item in sitecomposed
        """
        self.soiprocar.append_sumsite((0, 2), "test")
        for orbital in ("p", "pxpy", "d"):
            self.soiprocar.append_sumorbital(
                self.soiprocar.orbital_index(orbital), orbital
            )
        np.testing.assert_allclose(
            self.soiprocar.proj[0][0][0][3],  # mT
            [
                0.0020,
                0.0022,
                0.0024,
                0.0026,
                0.0028,
                0.0030,
                0.0032,
                0.0034,
                0.0036,
                0.0252,
                0.0072,
                0.0048,
                0.0160,
            ],
        )
        np.testing.assert_allclose(
            self.soiprocar.proj[2][0][0][3],  # mY
            [
                4.002,
                4.0022,
                4.0024,
                4.0026,
                4.0028,
                4.0030,
                4.0032,
                4.0034,
                4.0036,
                4.0252,
                12.0072,
                8.0048,
                20.016,
            ],
        )
        assert self.soiprocar.label["orbital"] == [
            "s",
            "py",
            "pz",
            "px",
            "dxy",
            "dyz",
            "dz2",
            "dxz",
            "dx2",
            "tot",
            "p",
            "pxpy",
            "d",
        ]

    def test_soiprocar_make_label(self):
        """test for make_label  (SOI)
        """
        self.soiprocar.append_sumsite((0, 2), "test")
        for orbital in ("p", "pxpy", "d"):
            self.soiprocar.append_sumorbital(
                self.soiprocar.orbital_index(orbital), orbital
            )
        for orbital in ("p", "pxpy", "d"):
            self.soiprocar.append_sumorbital(
                self.soiprocar.label["orbital"].index(orbital), orbital
            )
        assert self.soiprocar.make_label((3,), ((10, 11, 12),)) == [
            "#k",
            "Energy",
            "test_mT_p",
            "test_mX_p",
            "test_mY_p",
            "test_mZ_p",
            "test_mT_pxpy",
            "test_mX_pxpy",
            "test_mY_pxpy",
            "test_mZ_pxpy",
            "test_mT_d",
            "test_mX_d",
            "test_mY_d",
            "test_mZ_d",
        ]

    def test_text_sheet(self):
        """test for simple band data output (SOI)"""
        assert (
            self.soiprocar.text_sheet()
            == """#k	Energy
 0.00000000e+00	-1.00000000e+01
 3.53553391e-01	-7.00000000e+00
 7.07106781e-01	-6.00000000e+00

 0.00000000e+00	-5.00000000e+00
 3.53553391e-01	-4.00000000e+00
 7.07106781e-01	-1.00000000e+00

"""
        )


class test_functions_in_procarpy(object):
    def test_tiny_check_function(self):
        """Test for shortfunction"""
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_single = open(datadir + "PROCAR_single")
        data_spin = open(datadir + "PROCAR_spin_dummy")

        result_single = procar.tiny_check(data_single)
        assert 1 == result_single[0]  # numk
        assert 1 == result_single[1]  # nbands
        assert 3 == result_single[2]  # natom
        assert [
            "s",
            "py",
            "pz",
            "px",
            "dxy",
            "dyz",
            "dz2",
            "dxz",
            "dx2",
            "tot",
        ] == result_single[3]
        assert result_single[4]  # collinear
        result_spin = procar.tiny_check(data_spin)
        assert 3 == result_spin[0]  # numk
        assert 4 == result_spin[1]  # nbands
        assert 3 == result_spin[2]  # natom
        assert [
            "s",
            "py",
            "pz",
            "px",
            "dxy",
            "dyz",
            "dz2",
            "dxz",
            "dx2",
            "tot",
        ] == result_spin[3]

        assert result_spin[4]  # collinear
        data_soi = open(datadir + "PROCAR_SOI_dummy")
        result_soi = procar.tiny_check(data_soi)
        assert 3 == result_soi[0]  # numk
        assert 2 == result_soi[1]  # nbands
        assert 3 == result_soi[2]  # natom
        assert [
            "s",
            "py",
            "pz",
            "px",
            "dxy",
            "dyz",
            "dz2",
            "dxz",
            "dx2",
            "tot",
        ] == result_soi[3]

        assert not (result_soi[4])  # collinear

