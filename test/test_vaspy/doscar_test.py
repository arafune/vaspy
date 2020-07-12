#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import numpy as np
from nose.tools import assert_equal, assert_false, eq_, ok_, raises

from vaspy.doscar import DOS, DOSCAR, PDOS, TDOS


class TestDOSCAR(object):
    def setUp(self):
        currentpath = os.path.abspath(os.path.dirname(__file__))
        testDOS = currentpath + "/DOSCAR_dummy"
        self.testDOSCAR = DOSCAR(testDOS)


class TestDOS(object):
    pass


class TestTDOS(object):
    pass


class TestPDOS(object):
    pass
