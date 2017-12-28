#! /usr/bin/env python
# -*- coding: utf-8 -*-
import unittest
import os
from vaspy.locpot import LOCPOT
import numpy as np


class TestLOCPOT(unittest.TestCase):

    def setUp(self):
        currentpath = (os.path.abspath(os.path.dirname(__file__)))
        fullpath_1 = currentpath + "/data/LOCPOT.dummy"
        fullpath_2 = currentpath + "/data/LOCPOT.dummy2"
        self.testlocpot1 = LOCPOT(fullpath_1)
        self.testlocpot2 = LOCPOT(fullpath_2)

    def test_locpot_poscar_part(self):
        '''test whether LOCPOT correctly treats POSCAR'''
        self.assertEqual('hBN-Cu', self.testlocpot1.poscar.system_name)
        x, y, z = self.testlocpot2.poscar.axes_lengthes
        self.assertEqual((x, y, z),
                         (6.7395854049999997,
                          6.7395858129882784,
                          36.271284721999997))

    def test_average_along(self):
        '''test LOCPOT averaging'''
        # along Z-axis
        np.testing.assert_equal(np.array([5.5, 17.5, 29.5, 41.5, 53.5, 65.5]),
                                self.testlocpot2.grid.average_along_axis('z'))
        # along X-axis
        np.testing.assert_almost_equal(
            np.array([34.5, 35.5, 36.5]),
            self.testlocpot2.grid.average_along_axis('x'))
        # along Y-axis
        np.testing.assert_equal(np.array([31., 34., 37., 40.]),
                                self.testlocpot2.grid.average_along_axis('Y'))

if __name__ == '__main__':
    unittest.main()
