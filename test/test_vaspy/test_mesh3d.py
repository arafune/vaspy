#!/usr/bin/python
# -*- coding: utf-8 -*-
import os


import vaspy.mesh3d as mesh3d


class TestMesh3d(object):
    """Class for EIGENVAL class test"""

    def setup_method(self):
        os.path.abspath(os.path.dirname(__file__)) + "/data/"

    #        self.eigenval_soi = eigenval.EIGENVAL(datadir + 'EIGENVAL.soi')

    def test_check_basic_parameters(self):
        """Checking the basic parameters stored"""
        meshdata = list(range(3 * 4 * 5))
        self.simple_grid = mesh3d.Grid3D((3, 4, 5), meshdata)
