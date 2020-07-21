#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function  # Version safety

import os
from pathlib import Path
import tempfile

import numpy as np
from numpy.testing import assert_array_equal

import vaspy.poscar


class TestPOSCAR(object):
    def setup_method(self) -> None:
        global test_poscar_string
        filePOSCAR = tempfile.mkstemp()
        f = open(filePOSCAR[1], "w")
        f.write(test_poscar_string)
        f.close()
        self.testposcar = vaspy.poscar.POSCAR(filePOSCAR[1])
        self.blancposcar = vaspy.poscar.POSCAR()
        os.remove(filePOSCAR[1])
        datadir = Path(__file__).parent / "data"
        self.fepc = vaspy.poscar.POSCAR(str(datadir / "POSCAR.FePc"))

    def test_fundamentals(self):
        """Test for POSCAR class: fundamental data read"""
        # a=os.getcwd()  # return the directory where nose execute.
        assert "NiC4S4" == self.testposcar.system_name
        np.testing.assert_allclose(
            np.array([0.866025404, -0.5, 0.0]), self.testposcar.cell_vecs[0]
        )
        np.testing.assert_allclose(
            np.array([0.866025404, 0.5, 0.0]), self.testposcar.cell_vecs[1]
        )
        self.testposcar.cell_vecs[2] = (1, 0, 0)
        np.testing.assert_allclose(np.array([1, 0, 0]), self.testposcar.cell_vecs[2])
        assert [
            "#0:Ni1",
            "#1:Ni2",
            "#2:Ni3",
            "#3:C1",
            "#4:C2",
            "#5:C3",
            "#6:C4",
            "#7:C5",
            "#8:C6",
            "#9:C7",
            "#10:C8",
            "#11:C9",
            "#12:C10",
            "#13:C11",
            "#14:C12",
            "#15:S1",
            "#16:S2",
            "#17:S3",
            "#18:S4",
            "#19:S5",
            "#20:S6",
            "#21:S7",
            "#22:S8",
            "#23:S9",
            "#24:S10",
            "#25:S11",
            "#26:S12",
        ] == self.testposcar.site_label

    def test_point_in_box(self):
        assert not (
            vaspy.poscar.point_in_box((0.5, 0.5, 0.2), self.testposcar.cell_vecs)
        )
        assert vaspy.poscar.point_in_box((0.5, 0.1, 0.2), self.testposcar.cell_vecs)
        assert vaspy.poscar.point_in_box(
            (0.5, 0.5, 0.2), ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        )

    def test_is_cartesian(self):
        assert not (self.testposcar.is_cartesian())
        self.testposcar.to_cartesian()
        assert self.testposcar.is_cartesian()

    def test_is_direct(self):
        assert self.testposcar.is_direct()
        self.testposcar.to_cartesian()
        assert not self.testposcar.is_direct()

    def test_is_selective(self):
        assert self.testposcar.selective
        assert not self.blancposcar.selective

    def test_pos(self):
        np.testing.assert_array_equal(
            np.array([0.0, 0.5, 0.5]), self.testposcar.positions[2]
        )
        np.testing.assert_allclose(
            np.array([0.23764, 0.429027113, 0.5]),
            self.testposcar.positions[3],
            rtol=1e-04,
        )

    def test_tune_scaling_factor(self):
        self.testposcar.tune_scaling_factor(1.0)
        np.testing.assert_allclose(
            np.array([12.66995166052, -7.315, 0.0]),
            self.testposcar.cell_vecs[0],
            rtol=1e-07,
        )

    def test_tune_scaling_factor_withCartesian(self):
        self.testposcar.to_cartesian()
        self.testposcar.tune_scaling_factor(1.0)
        np.testing.assert_allclose(
            np.array([12.66995166052, -7.315, 0.0]),
            self.testposcar.cell_vecs[0],
            rtol=1e-07,
        )
        np.testing.assert_allclose(
            np.array([12.66995166052, 0.0, 7.5000000001850005]),
            self.testposcar.positions[0],
            rtol=1e-06,
        )
        np.testing.assert_allclose(
            np.array([8.4466344319, -1.4000000014, 7.500000000185]),
            self.testposcar.positions[6],
            rtol=1e-06,
        )

    def test_to_cartesian(self):
        self.testposcar.to_cartesian()
        np.testing.assert_allclose(
            np.array([0.494477, -0.047847, 0.512645]),
            self.testposcar.positions[5],
            rtol=1e-05,
        )
        self.testposcar.to_direct()

    def test_to_direct(self):
        tmp = self.testposcar.positions[5]
        self.testposcar.to_direct()
        np.testing.assert_array_equal(tmp, self.testposcar.positions[5])

    def test_to_list(self):
        tmp = self.testposcar.to_list()
        assert "NiC4S4" == tmp[0]
        assert 14.63 == tmp[1]
        np.testing.assert_allclose([0.8660254, -0.5, 0], tmp[2])
        np.testing.assert_allclose([0.8660254, 0.5, 0], tmp[3])
        np.testing.assert_allclose([0, 0, 1.0252904990], tmp[4])
        assert ["Ni", "C", "S"] == tmp[5]
        assert [3, 12, 12] == tmp[6]
        assert "Selective Dynamics" == tmp[7]
        assert "Direct" == tmp[8]
        np.testing.assert_allclose([0.5, 0.5, 0.5], tmp[9][0])
        np.testing.assert_allclose([0.5, 0.0, 0.5], tmp[9][1])
        np.testing.assert_allclose([0.0, 0.5, 0.5], tmp[9][2])
        np.testing.assert_allclose([0.237639553, 0.429027113, 0.5], tmp[9][3])
        # ...
        assert "T T T" == tmp[10][0]
        assert "T T T" == tmp[10][1]
        assert "T F T" == tmp[10][2]
        #
        assert "#0:Ni1" == tmp[11][0]
        assert "#1:Ni2" == tmp[11][1]
        assert "#2:Ni3" == tmp[11][2]
        assert "#3:C1" == tmp[11][3]

    def test_to_str(self):
        global tmpstr_original
        global tmpstr_after_rotate
        assert tmpstr_original == self.testposcar.__str__()
        self.testposcar.rotate_atom(0, "z", 90, (0, 0, 0))
        self.testposcar.to_direct()

    def test_poscar_supercell1(self):
        """Tests for poscar supercell method."""
        supercell = self.testposcar.supercell(3, 2, 1)
        np.testing.assert_allclose(
            np.array([2.59807621, -1.5, 0.0]), supercell.cell_vecs[0]
        )
        np.testing.assert_allclose(
            np.array([1.73205081, 1.0, 0.0]), supercell.cell_vecs[1]
        )
        np.testing.assert_allclose(
            np.array([0.0, 0.0, 1.02529049]), supercell.cell_vecs[2]
        )
        assert "NiC4S4" == supercell.system_name
        assert ["Ni", "C", "S"] == supercell.atomtypes
        assert [18, 72, 72] == supercell.atomnums
        supercell = self.testposcar.supercell(1, 1, 1)
        np.testing.assert_allclose(self.testposcar.positions[0], supercell.positions[0])
        supercell = self.testposcar.supercell(3, 2, 1)
        np.testing.assert_allclose(
            supercell.positions[0],
            np.array(
                [
                    self.testposcar.positions[0][0] / 3,
                    self.testposcar.positions[0][1] / 2,
                    self.testposcar.positions[0][2] / 1,
                ]
            ),
        )
        np.testing.assert_allclose(
            supercell.positions[1],
            np.array(
                [
                    self.testposcar.positions[0][0] / 3 + 1 * (1 / 3),
                    self.testposcar.positions[0][1] / 2,
                    self.testposcar.positions[0][2] / 1,
                ]
            ),
        )
        supercell = self.testposcar.supercell(3, 2, 1)
        assert 6 * len(self.testposcar.positions) == len(supercell.positions)

    def test_poscar_split(self):
        """Test for POSCAR.split."""
        one, other = self.testposcar.split([2, 3, 4, 5, 6])
        assert 5 == len(one.positions)
        assert ["Ni", "C"] == one.atomtypes
        assert [1, 4] == one.atomnums
        assert one.selective
        assert other.selective

    def test_nearest(self):
        pass

    def test_make27candidate(self):
        pass

    def test_guess_molecule(self):
        pass

    def test_atom_rotate(self):
        pass

    def test_plus(self):
        pass

    def test__getitem(self) -> None:
        print(self.fepc)
        assert_array_equal(
            self.fepc[0], [0.1362036989956494, 0.0552692292217493, 0.0000000000000000]
        )
        assert self.fepc[0][1] == 0.0552692292217493
        assert_array_equal(
            self.fepc[0:3],
            [
                [0.1362036989956494, 0.0552692292217493, 0.0000000000000000],
                [0.2050800706767476, 0.0348757867600226, 0.0000000000000000],
                [0.2642020159712474, 0.0709072792541197, 0.0000000000000000],
            ],
        )


test_poscar_string = """NiC4S4
14.63
0.866025404	-0.5	0
0.866025404	0.5	0
0	0	1.025290499
Ni	C	S
3	12	12
Selective Dynamics
Direct
0.5	0.5	0.5     T  T  T
0.5	0       0.5	T  T  T
0	0.5	0.5     T  F  T
0.237639553	0.429027113	0.5  T  T  T
0.237639553	0.333333333	0.5  T  T  T
0.333333333	0.237639553	0.5  T  T  T
0.429027113	0.237639553	0.5  T  T  T
0.429027113	0.333333333	0.5  T  T  T
0.333333333	0.429027113	0.5  T  T  T
-0.237639553	-0.429027113	0.5  T  T  T
-0.237639553	-0.333333333	0.5  T  T  T
-0.333333333	-0.237639553	0.5  T  T  T
-0.429027113	-0.237639553	0.5  T  T  T
-0.429027113	-0.333333333	0.5  T  T  T
-0.333333333	-0.429027113	0.5  T  T  T
0.11323764	0.553429027	0.5  T  T  T
0.11323764	0.333333333	0.5  T  T  T
0.333333333	0.11323764	0.5  T  T  T
0.553429027	0.11323764	0.5  T  T  T
0.553429027	0.333333333	0.5  T  T  T
0.333333333	0.553429027	0.5  T  T  T
-0.11323764	-0.553429027	0.5  T  T  T
-0.11323764	-0.333333333	0.5  T  T  T
-0.333333333	-0.11323764	0.5  T  T  T
-0.553429027	-0.11323764	0.5  T  T  T
-0.553429027	-0.333333333	0.5  T  T  T
-0.333333333	-0.553429027	0.5  T  T  T
"""

tmpstr_original = """NiC4S4
14.63
    0.86602540400000005   -0.50000000000000000    0.00000000000000000
    0.86602540400000005    0.50000000000000000    0.00000000000000000
    0.00000000000000000    0.00000000000000000    1.02529049900000002
 Ni C S
 3 12 12
Selective Dynamics
Direct
   0.50000000000000000    0.50000000000000000    0.50000000000000000 T T T #0:Ni1
   0.50000000000000000    0.00000000000000000    0.50000000000000000 T T T #1:Ni2
   0.00000000000000000    0.50000000000000000    0.50000000000000000 T F T #2:Ni3
   0.23763955300000000    0.42902711300000002    0.50000000000000000 T T T #3:C1
   0.23763955300000000    0.33333333300000001    0.50000000000000000 T T T #4:C2
   0.33333333300000001    0.23763955300000000    0.50000000000000000 T T T #5:C3
   0.42902711300000002    0.23763955300000000    0.50000000000000000 T T T #6:C4
   0.42902711300000002    0.33333333300000001    0.50000000000000000 T T T #7:C5
   0.33333333300000001    0.42902711300000002    0.50000000000000000 T T T #8:C6
  -0.23763955300000000   -0.42902711300000002    0.50000000000000000 T T T #9:C7
  -0.23763955300000000   -0.33333333300000001    0.50000000000000000 T T T #10:C8
  -0.33333333300000001   -0.23763955300000000    0.50000000000000000 T T T #11:C9
  -0.42902711300000002   -0.23763955300000000    0.50000000000000000 T T T #12:C10
  -0.42902711300000002   -0.33333333300000001    0.50000000000000000 T T T #13:C11
  -0.33333333300000001   -0.42902711300000002    0.50000000000000000 T T T #14:C12
   0.11323764000000000    0.55342902699999996    0.50000000000000000 T T T #15:S1
   0.11323764000000000    0.33333333300000001    0.50000000000000000 T T T #16:S2
   0.33333333300000001    0.11323764000000000    0.50000000000000000 T T T #17:S3
   0.55342902699999996    0.11323764000000000    0.50000000000000000 T T T #18:S4
   0.55342902699999996    0.33333333300000001    0.50000000000000000 T T T #19:S5
   0.33333333300000001    0.55342902699999996    0.50000000000000000 T T T #20:S6
  -0.11323764000000000   -0.55342902699999996    0.50000000000000000 T T T #21:S7
  -0.11323764000000000   -0.33333333300000001    0.50000000000000000 T T T #22:S8
  -0.33333333300000001   -0.11323764000000000    0.50000000000000000 T T T #23:S9
  -0.55342902699999996   -0.11323764000000000    0.50000000000000000 T T T #24:S10
  -0.55342902699999996   -0.33333333300000001    0.50000000000000000 T T T #25:S11
  -0.33333333300000001   -0.55342902699999996    0.50000000000000000 T T T #26:S12
"""

tmpstr_after_rotate = """NiC4S4
14.63
    0.86602540400000005   -0.50000000000000000    0.00000000000000000
    0.86602540400000005    0.50000000000000000    0.00000000000000000
    0.00000000000000000    0.00000000000000000    1.02529049900000002
 Ni C S
 3 12 12
Selective Dynamics
Direct
  -0.86602540400000005    0.86602540400000005    0.50000000000000000 T T T #0:Ni1
   0.50000000000000000    0.00000000000000000    0.50000000000000000 T T T #1:Ni2
   0.00000000000000000    0.50000000000000000    0.50000000000000000 T F T #2:Ni3
   0.23763955299999995    0.42902711299999996    0.50000000000000000 T T T #3:C1
   0.23763955299999998    0.33333333300000001    0.50000000000000000 T T T #4:C2
   0.33333333300000001    0.23763955299999998    0.50000000000000000 T T T #5:C3
   0.42902711300000002    0.23763955300000000    0.50000000000000000 T T T #6:C4
   0.42902711299999996    0.33333333300000001    0.50000000000000000 T T T #7:C5
   0.33333333300000001    0.42902711299999996    0.50000000000000000 T T T #8:C6
  -0.23763955299999995   -0.42902711299999996    0.50000000000000000 T T T #9:C7
  -0.23763955299999998   -0.33333333300000001    0.50000000000000000 T T T #10:C8
  -0.33333333300000001   -0.23763955299999998    0.50000000000000000 T T T #11:C9
  -0.42902711300000002   -0.23763955300000000    0.50000000000000000 T T T #12:C10
  -0.42902711299999996   -0.33333333300000001    0.50000000000000000 T T T #13:C11
  -0.33333333300000001   -0.42902711299999996    0.50000000000000000 T T T #14:C12
   0.11323763999999997    0.55342902699999996    0.50000000000000000 T T T #15:S1
   0.11323763999999999    0.33333333300000001    0.50000000000000000 T T T #16:S2
   0.33333333300000001    0.11323763999999999    0.50000000000000000 T T T #17:S3
   0.55342902699999996    0.11323764000000003    0.50000000000000000 T T T #18:S4
   0.55342902699999996    0.33333333299999995    0.50000000000000000 T T T #19:S5
   0.33333333299999995    0.55342902699999996    0.50000000000000000 T T T #20:S6
  -0.11323763999999997   -0.55342902699999996    0.50000000000000000 T T T #21:S7
  -0.11323763999999999   -0.33333333300000001    0.50000000000000000 T T T #22:S8
  -0.33333333300000001   -0.11323763999999999    0.50000000000000000 T T T #23:S9
  -0.55342902699999996   -0.11323764000000003    0.50000000000000000 T T T #24:S10
  -0.55342902699999996   -0.33333333299999995    0.50000000000000000 T T T #25:S11
  -0.33333333299999995   -0.55342902699999996    0.50000000000000000 T T T #26:S12
"""
