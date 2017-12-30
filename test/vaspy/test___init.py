#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This module provides test routines for __init__.py for vaspy module
'''
import os
from nose.tools import ok_, eq_
from nose.tools import with_setup

import vaspy

class TestInitModule(object):
    '''Class for Test of vaspy.module

    The functions should be tested are in __init__.py
    '''
    def setup(self):
        self.datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"

    def test_load_in_init(self):
        '''Test for __init__.py'''
        poscar = self.datadir + 'POSCAR_dummy'
        procar_single = self.datadir + 'PROCAR_single'
        procar_spin = self.datadir + 'PROCAR_spin_dummy'
        procar_soi = self.datadir + 'PROCAR_soi_dummy'
        doscar = self.datadir + 'DOSCAR_dummy'
        ok_(isinstance(vaspy.load(poscar), vaspy.poscar.POSCAR))
        ok_(isinstance(vaspy.load(procar_single), vaspy.procar.PROCAR))
        ok_(isinstance(vaspy.load(procar_spin), vaspy.procar.PROCAR))
        ok_(isinstance(vaspy.load(procar_soi), vaspy.procar.PROCAR))
        ok_(isinstance(vaspy.load(doscar), vaspy.doscar.DOSCAR))
