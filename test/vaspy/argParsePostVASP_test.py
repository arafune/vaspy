#!/usr/env/bin python
# -*- coding: utf-8 -*-

import unittest
from vaspy.argParsePostVASP import APPV

class TestAPPV(unittest.TestCase):
    def setUp(self):
        self.v=APPV()
    def parse_Atomselection_test(self):
        self.assertEqual(["1", "10", "11", "12", "13", "14", "15", "2", "3", "4", "5", "8", "9"], self.v.parse_Atomselection("1-5,8,8,9-15,10"))
    def parse_AtomselectionNum_test(self):
        self.assertEqual([1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15 ], self.v.parse_AtomselectionNum("1-5,8,8,9-15,10"))



#parse_AtomselectionNum_test
if __name__ == '__main__':
    unittest.main()
    
