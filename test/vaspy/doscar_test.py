#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import os
import numpy as np
from vaspy.doscar import DOSCAR, DOS, TDOS, PDOS


class TestDOSCAR(unittest.TestCase):

    def setUp(self):
        currentpath = (os.path.abspath(os.path.dirname(__file__)))
        testDOS = currentpath + "/DOSCAR_dummy"
        self.testDOSCAR = DOSCAR(testDOS)
    pass


class TestDOS(unittest.TestCase):
    pass


class TestTDOS(unittest.TestCase):
    pass


class TestPDOS(unittest.TestCase):
    pass


if __name__ == '__main__':
    unittest.main()
