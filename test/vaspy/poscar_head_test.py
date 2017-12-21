#!/usr/bin/env python
# -*- coding: utf-8 -*-

from nose.tools import with_setup
from nose.tools import eq_
from nose.tools import ok_
import numpy as np
import vaspy.poscar


class TestPOSCAR_Head(object):
    def setup(self):
        self.poscar_head=vaspy.poscar.POSCAR_HEAD()
        self.poscar_head.iontypes=['Ag', 'Si']
        self.poscar_head.ionnums=[3, 5]
        self.poscar_head.system_name='testPOSCAR'

    @with_setup(setup=setup)
    def test_poscar_head(self):
        eq_('testPOSCAR', self.poscar_head.system_name)
        eq_('#0:Ag1', self.poscar_head.atom_identifer[0])
