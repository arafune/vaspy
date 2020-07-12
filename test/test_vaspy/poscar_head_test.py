#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

# from nose.tools import with_setup

from vaspy import poscar


class Test_POSCAR_Head(object):
    def setup_method(self, method):
        self.poscar_head = poscar.POSCAR_HEAD()
        self.poscar_head.atomtypes = ["Ag", "Si"]
        self.poscar_head.atomnums = [3, 5]
        self.poscar_head.system_name = "testPOSCAR"

    # @with_setup(setup=setup)
    def test_poscar_head(self):
        assert "testPOSCAR" == self.poscar_head.system_name
        assert "#0:Ag1" == self.poscar_head.site_label[0]
