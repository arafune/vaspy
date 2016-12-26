#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
from nose.tools import eq_
from nose.tools import ok_
from nose.tools import with_setup
from nose.tools import assert_equal
from nose.tools import raises
import vaspy.chgcar
import numpy as np


class TestCHGCAR(object):
    '''Class for test of CHGCAR module

    Use CHGCAR for test data'''
    def setup(self):
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file1 = datadir + 'CHGCAR_1'
        data_file2 = datadir + 'CHGCAR_2'
        data_file_spin = datadir + 'CHGCAR_spin'
        data_file_soi =  datadir + 'CHGCAR_soi'
        #self.chgcar1 = vaspy.chgcar.CHGCAR(data_file1)
        #self.chgcar2 = vaspy.chgcar.CHGCAR(data_file2)
        #self.chgcar_spin = vasp.chgcar.CHGCAR(data_file_spin)
        #self.chgcar_soi = vasp.chgcar.CHGCAR(data_file_soi)

    @with_setup(setup=setup)
    def testCHGsum(self):
        pass

    @with_setup(setup=setup)
    def testCHGdiff(self):
        pass

    @with_setup(setup=setup)
    def magCHG(self):
        pass

    @with_setup
    def spindirec_CHGCAR_SOI(self):
        pass
