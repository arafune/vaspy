#!/usr/env/bin python
# -*- coding: utf-8 -*-

import unittest
import  numpy as np
from vaspy.poscar import POSCAR

class TestPOSCAR(unittest.TestCase):
    def setUp(self):
        self.p = POSCAR()
        self.poscar = POSCAR()
        self.p._POSCAR__latticeV1 = np.array([[0.0, 0.0, 0.1]])

##        
##  for object "p"
##      
    def latticeV1_test(self):
        np.testing.assert_array_equal(np.array([[0., 0., 0.1]]),
                                      self.p.latticeV1)
    def latticeV2_test(self):
        np.testing.assert_array_equal(np.array([[0., 0., 0.]]),
                                      self.p.latticeV2)
    def latticeV3_test(self):
        np.testing.assert_array_equal(np.array([[0., 0., 0.]]),
                                      self.p.latticeV3)
    

if __name__ == '__main__':
    unittest.main()
    
