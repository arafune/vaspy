#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
from nose.tools import eq_
from nose.tools import ok_
from nose.tools import with_setup
from nose.tools import assert_equal
from nose.tools import raises
import vaspy.mesh3d as mesh3d
import numpy as np

class TestMesh3d(object):
    '''Class for EIGENVAL class test'''

    def setup(self):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
#        self.eigenval_soi = eigenval.EIGENVAL(datadir + 'EIGENVAL.soi')
        
    @with_setup(setup=setup)
    def test_check_basic_parameters(self):
        '''Checking the basic parameters stored'''
        meshdata = list (range(3*4*5))
        self.simple_grid = mesh3d.Grid3D((3,4,5), meshdata)
        
