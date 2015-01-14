#! /usr/bin/env python
# -*- coding: utf-8 -*-
import unittest
import os
from vaspy.locpot import LOCPOT
import numpy as np


class TestLOCPOT(unittest.TestCase):

    def setUp(self):
        currentpath = (os.path.abspath(os.path.dirname(__file__)))
        fullpath_1 = currentpath + "/LOCPOT.dummy"
        fullpath_2 = currentpath + "/LOCPOT.dummy2"
        self.testlocpot1 = LOCPOT(fullpath_1)
        self.testlocpot2 = LOCPOT(fullpath_2)

    def test_system_name(self):
        self.assertEqual('hBN-Cu', self.testlocpot1.system_name)

    def test_average_along_axisZ(self):
        np.testing.assert_equal(np.array([5.5, 17.5, 29.5, 41.5, 53.5, 65.5]),
                                self.testlocpot2.average_along_axis('z'))

    def test_average_along_axisX(self):
        np.testing.assert_almost_equal(
            np.array([34.5,  35.5,  36.5]),
            self.testlocpot2.average_along_axis('x'))

    def test_average_along_axisY(self):
        np.testing.assert_equal(np.array([31.,  34.,  37.,  40.]),
                                self.testlocpot2.average_along_axis('Y'))

    def test_axis_length(self):
        x, y, z = self.testlocpot2.get_axes_lengthes()
        self.assertEqual((x, y, z),
                         (6.7395854049999997,
                          6.7395858129882784,
                          36.271284721999997))

if __name__ == '__main__':
    unittest.main()
