#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from nose.tools import eq_, ok_, with_setup

import vaspy.poscar


class TestPOSCAR_Head(object):
    def setup(self):
        self.poscar_head = vaspy.poscar.POSCAR_HEAD()
        self.poscar_head.atomtypes = ['Ag', 'Si']
        self.poscar_head.atomnums = [3, 5]
        self.poscar_head.system_name = 'testPOSCAR'

    @with_setup(setup=setup)
    def test_poscar_head(self):
        eq_('testPOSCAR', self.poscar_head.system_name)
        eq_('#0:Ag1', self.poscar_head.site_label[0])
