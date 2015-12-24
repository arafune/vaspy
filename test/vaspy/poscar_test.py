#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import os
import tempfile
import numpy as np
from vaspy.poscar import POSCAR


class TestPOSCAR(unittest.TestCase):

    def setUp(self):
        global test_poscar_string
        filePOSCAR = tempfile.mkstemp()
        f = open(filePOSCAR[1], 'w')
        f.write(test_poscar_string)
        f.close()
        self.testposcar = POSCAR(filePOSCAR[1])
        self.blancposcar = POSCAR()
        os.remove(filePOSCAR[1])

    def test_system_name(self):
        # a=os.getcwd()  # return the directory where nose execute.
        self.assertEqual('NiC4S4', self.testposcar.system_name)

    def test_lattice_vec1(self):
        np.testing.assert_array_equal(np.array([0.866025404, -0.5, 0.]),
                                      self.testposcar.lattice_vec1)

    def test_lattice_vec2(self):
        np.testing.assert_array_equal(np.array([0.866025404, 0.5, 0.]),
                                      self.testposcar.lattice_vec2)

    def test_lattice_vec3(self):
        np.testing.assert_array_equal(np.array([0.0, 0.0, 1.025290499]),
                                      self.testposcar.lattice_vec3)

    def test_atom_identifer(self):
        self.assertEqual(['#1:Ni1', '#2:Ni2', '#3:Ni3',
                          '#4:C1', '#5:C2', '#6:C3', '#7:C4',
                          '#8:C5', '#9:C6', '#10:C7', '#11:C8',
                          '#12:C9', '#13:C10', '#14:C11', '#15:C12',
                          '#16:S1', '#17:S2', '#18:S3', '#19:S4',
                          '#20:S5', '#21:S6', '#22:S7', '#23:S8',
                          '#24:S9', '#25:S10', '#26:S11', '#27:S12'],
                         self.testposcar.atom_identifer)

    def test_is_cartesian(self):
        self.assertFalse(self.testposcar.is_cartesian)
        self.testposcar.to_cartesian()
        self.assertTrue(self.testposcar.is_cartesian)

    def test_is_direct(self):
        self.assertTrue(self.testposcar.is_direct)

    def test_is_selective(self):
        self.assertTrue(self.testposcar.is_selective)
        self.assertFalse(self.blancposcar.is_selective)

    def test_pos(self):
        np.testing.assert_array_equal(np.array([0., 0.5, 0.5]),
                                      self.testposcar.pos(3))
        np.testing.assert_allclose(np.array([0.23764, 0.429027113, 0.5]),
                                   self.testposcar.pos(4),
                                   rtol=1e-04)

    def test_pos_replace(self):
        np.testing.assert_array_equal(
            np.array([0.333333333, 0.237639553, 0.5]),
            self.testposcar.pos(6))

        self.assertRaises(RuntimeError,
                          self.testposcar.pos_replace,
                          6,
                          [0, 0, 0])
        self.testposcar.to_cartesian()
        self.testposcar.pos_replace(6, [0, 0, 0])
        np.testing.assert_array_equal(np.array([0, 0, 0]),
                                      self.testposcar.pos(6))
        self.testposcar.pos_replace(7, [1, 3, 4])
        np.testing.assert_array_equal(np.array([1, 3, 4]),
                                      self.testposcar.pos(7))

    def test_pos_raise_value_error(self):
        self.assertRaises(ValueError, self.testposcar.pos,-1)
        self.assertRaises(ValueError, self.testposcar.pos,0)
        
    def test_tune_scaling_factor(self):
        self.testposcar.tune_scaling_factor(1.0)
        np.testing.assert_allclose(
            np.array([12.66995166052, -7.315, 0.0]),
            self.testposcar.lattice_vec1,
            rtol=1e-07)

    def test_tune_scaling_factor_withCartesian(self):
        self.testposcar.to_cartesian()
        self.testposcar.tune_scaling_factor(1.0)
        np.testing.assert_allclose(
            np.array([12.66995166052, -7.315, 0.0]),
            self.testposcar.lattice_vec1,
            rtol=1e-07)
        np.testing.assert_allclose(
            np.array([12.66995166052, 0.0, 7.5000000001850005]),
            self.testposcar.pos(1),
            rtol=1e-06)
        np.testing.assert_allclose(
            np.array([8.4466344319, -1.4000000014, 7.500000000185]),
            self.testposcar.pos(7),
            rtol=1e-06)

    def test_to_cartesian(self):
        tmp = self.testposcar.pos(6)
        self.testposcar.to_cartesian()
        np.testing.assert_allclose(
            np.array([0.494477, -0.047847, 0.512645]),
            self.testposcar.pos(6),
            rtol=1e-05)
        self.testposcar.to_direct()
        np.testing.assert_allclose(
            tmp,
            self.testposcar.pos(6),
            rtol=1e-05)

    def test_to_direct(self):
        tmp = self.testposcar.pos(6)
        self.testposcar.to_direct()
        np.testing.assert_array_equal(tmp, self.testposcar.pos(6))

    def test_to_list(self):
        tmp = self.testposcar.to_list()
        self.assertEqual("NiC4S4", tmp[0])
        self.assertEqual(14.63, tmp[1])
        np.testing.assert_allclose([0.8660254, -0.5, 0], tmp[2])
        np.testing.assert_allclose([0.8660254, 0.5, 0], tmp[3])
        np.testing.assert_allclose([0, 0, 1.0252904990], tmp[4])
        self.assertEqual(['Ni', 'C', 'S'], tmp[5])
        self.assertEqual([3, 12, 12], tmp[6])
        self.assertEqual('Selective Dynamics', tmp[7])
        self.assertEqual('Direct', tmp[8])
        np.testing.assert_allclose([0.5, 0.5, 0.5], tmp[9][0])
        np.testing.assert_allclose([0.5, 0., 0.5], tmp[9][1])
        np.testing.assert_allclose([0., 0.5, 0.5], tmp[9][2])
        np.testing.assert_allclose([0.237639553, 0.429027113, 0.5],
                                   tmp[9][3])
        # ...
        self.assertEqual('T T T', tmp[10][0])
        self.assertEqual('T T T', tmp[10][1])
        self.assertEqual('T F T', tmp[10][2])
        #
        self.assertEqual('#1:Ni1', tmp[11][0])
        self.assertEqual('#2:Ni2', tmp[11][1])
        self.assertEqual('#3:Ni3', tmp[11][2])
        self.assertEqual('#4:C1', tmp[11][3])

    def test_to_str(self):
        global tmpstr_original
        global tmpstr_after_rotate
        self.assertEqual(tmpstr_original, self.testposcar.__str__())
        self.testposcar.atom_rotate(1, "z", 90, (0, 0, 0))
        self.testposcar.to_direct()
        self.assertEqual(tmpstr_after_rotate,
                         self.testposcar.__str__())


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
   0.50000000000000000    0.50000000000000000    0.50000000000000000 T T T #1:Ni1
   0.50000000000000000    0.00000000000000000    0.50000000000000000 T T T #2:Ni2
   0.00000000000000000    0.50000000000000000    0.50000000000000000 T F T #3:Ni3
   0.23763955300000000    0.42902711300000002    0.50000000000000000 T T T #4:C1
   0.23763955300000000    0.33333333300000001    0.50000000000000000 T T T #5:C2
   0.33333333300000001    0.23763955300000000    0.50000000000000000 T T T #6:C3
   0.42902711300000002    0.23763955300000000    0.50000000000000000 T T T #7:C4
   0.42902711300000002    0.33333333300000001    0.50000000000000000 T T T #8:C5
   0.33333333300000001    0.42902711300000002    0.50000000000000000 T T T #9:C6
  -0.23763955300000000   -0.42902711300000002    0.50000000000000000 T T T #10:C7
  -0.23763955300000000   -0.33333333300000001    0.50000000000000000 T T T #11:C8
  -0.33333333300000001   -0.23763955300000000    0.50000000000000000 T T T #12:C9
  -0.42902711300000002   -0.23763955300000000    0.50000000000000000 T T T #13:C10
  -0.42902711300000002   -0.33333333300000001    0.50000000000000000 T T T #14:C11
  -0.33333333300000001   -0.42902711300000002    0.50000000000000000 T T T #15:C12
   0.11323764000000000    0.55342902699999996    0.50000000000000000 T T T #16:S1
   0.11323764000000000    0.33333333300000001    0.50000000000000000 T T T #17:S2
   0.33333333300000001    0.11323764000000000    0.50000000000000000 T T T #18:S3
   0.55342902699999996    0.11323764000000000    0.50000000000000000 T T T #19:S4
   0.55342902699999996    0.33333333300000001    0.50000000000000000 T T T #20:S5
   0.33333333300000001    0.55342902699999996    0.50000000000000000 T T T #21:S6
  -0.11323764000000000   -0.55342902699999996    0.50000000000000000 T T T #22:S7
  -0.11323764000000000   -0.33333333300000001    0.50000000000000000 T T T #23:S8
  -0.33333333300000001   -0.11323764000000000    0.50000000000000000 T T T #24:S9
  -0.55342902699999996   -0.11323764000000000    0.50000000000000000 T T T #25:S10
  -0.55342902699999996   -0.33333333300000001    0.50000000000000000 T T T #26:S11
  -0.33333333300000001   -0.55342902699999996    0.50000000000000000 T T T #27:S12
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
  -0.86602540400000005    0.86602540400000005    0.50000000000000000 T T T #1:Ni1
   0.50000000000000000    0.00000000000000000    0.50000000000000000 T T T #2:Ni2
   0.00000000000000000    0.50000000000000000    0.50000000000000000 T F T #3:Ni3
   0.23763955299999995    0.42902711299999996    0.50000000000000000 T T T #4:C1
   0.23763955299999998    0.33333333300000001    0.50000000000000000 T T T #5:C2
   0.33333333300000001    0.23763955299999998    0.50000000000000000 T T T #6:C3
   0.42902711300000002    0.23763955300000000    0.50000000000000000 T T T #7:C4
   0.42902711299999996    0.33333333300000001    0.50000000000000000 T T T #8:C5
   0.33333333300000001    0.42902711299999996    0.50000000000000000 T T T #9:C6
  -0.23763955299999995   -0.42902711299999996    0.50000000000000000 T T T #10:C7
  -0.23763955299999998   -0.33333333300000001    0.50000000000000000 T T T #11:C8
  -0.33333333300000001   -0.23763955299999998    0.50000000000000000 T T T #12:C9
  -0.42902711300000002   -0.23763955300000000    0.50000000000000000 T T T #13:C10
  -0.42902711299999996   -0.33333333300000001    0.50000000000000000 T T T #14:C11
  -0.33333333300000001   -0.42902711299999996    0.50000000000000000 T T T #15:C12
   0.11323763999999997    0.55342902699999996    0.50000000000000000 T T T #16:S1
   0.11323763999999999    0.33333333300000001    0.50000000000000000 T T T #17:S2
   0.33333333300000001    0.11323763999999999    0.50000000000000000 T T T #18:S3
   0.55342902699999996    0.11323764000000003    0.50000000000000000 T T T #19:S4
   0.55342902699999996    0.33333333299999995    0.50000000000000000 T T T #20:S5
   0.33333333299999995    0.55342902699999996    0.50000000000000000 T T T #21:S6
  -0.11323763999999997   -0.55342902699999996    0.50000000000000000 T T T #22:S7
  -0.11323763999999999   -0.33333333300000001    0.50000000000000000 T T T #23:S8
  -0.33333333300000001   -0.11323763999999999    0.50000000000000000 T T T #24:S9
  -0.55342902699999996   -0.11323764000000003    0.50000000000000000 T T T #25:S10
  -0.55342902699999996   -0.33333333299999995    0.50000000000000000 T T T #26:S11
  -0.33333333299999995   -0.55342902699999996    0.50000000000000000 T T T #27:S12
"""

if __name__ == '__main__':
    unittest.main()
