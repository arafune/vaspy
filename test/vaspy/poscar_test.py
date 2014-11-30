#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import os
import  numpy as np
from vaspy.poscar import POSCAR

class TestPOSCAR(unittest.TestCase):
    def setUp(self):
#       self.testposcar = POSCAR()
        test_poscar_string="""NiC4S4
14.63		
0.866025404	-0.5	0
0.866025404	0.5	0
0	0	1.025290499
Ni	C	S
3	12	12
Selective Dynamics
Direct
0.5	0.5	0.5
0.5	0       0.5	
0	0.5	0.5
0.237639553	0.429027113	0.5
0.237639553	0.333333333	0.5
0.333333333	0.237639553	0.5
0.429027113	0.237639553	0.5
0.429027113	0.333333333	0.5
0.333333333	0.429027113	0.5
-0.237639553	-0.429027113	0.5
-0.237639553	-0.333333333	0.5
-0.333333333	-0.237639553	0.5
-0.429027113	-0.237639553	0.5
-0.429027113	-0.333333333	0.5
-0.333333333	-0.429027113	0.5
0.11323764	0.553429027	0.5
0.11323764	0.333333333	0.5
0.333333333	0.11323764	0.5
0.553429027	0.11323764	0.5
0.553429027	0.333333333	0.5
0.333333333	0.553429027	0.5
-0.11323764	-0.553429027	0.5
-0.11323764	-0.333333333	0.5
-0.333333333	-0.11323764	0.5
-0.553429027	-0.11323764	0.5
-0.553429027	-0.333333333	0.5
-0.333333333	-0.553429027	0.5
"""
        f = open('testPOSCAR', 'w') 
        f.write(test_poscar_string) 
        f.close() 
        self.testposcar = POSCAR("./testPOSCAR")
##        
##  
##      
    def system_name_test(self):
        # a=os.getcwd()  # return the directory where nose execute.
        self.assertEqual('NiC4S4',self.testposcar.system_name)
    def latticeV1_test(self):
        np.testing.assert_array_equal(np.array([[0.866025404, -0.5, 0.]]),
                                      self.testposcar.latticeV1)
    def latticeV2_test(self):
        np.testing.assert_array_equal(np.array([[0.866025404, 0.5, 0.]]),
                                      self.testposcar.latticeV2)
    def latticeV3_test(self):
        np.testing.assert_array_equal(np.array([[0.0, 0.0, 1.025290499]]),
                                      self.testposcar.latticeV3)

    def atom_identifer_test(self):
        self.assertEqual(['#1:Ni1', '#2:Ni2', '#3:Ni3',
                          '#4:C1', '#5:C2', '#6:C3', '#7:C4',
                          '#8:C5', '#9:C6', '#10:C7', '#11:C8',
                          '#12:C9', '#13:C10', '#14:C11', '#15:C12',
                          '#16:S1', '#17:S2', '#18:S3', '#19:S4',
                          '#20:S5', '#21:S6', '#22:S7', '#23:S8',
                          '#24:S9', '#25:S10', '#26:S11', '#27:S12'],
                         self.testposcar.atom_identifer)

    def is_cartesian_test(self):
        self.assertFalse(self.testposcar.is_cartesian)

    def is_direct_test(self):
        self.assertTrue(self.testposcar.is_direct)

    def is_selective_test(self):
        self.assertTrue(self.testposcar.is_selective)

    def pos_test(self):
        pass



if __name__ == '__main__':
    unittest.main()
    
