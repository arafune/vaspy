#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
from __future__ import division
from __future__ import print_function
import os
import sys
import matplotlib.pyplot as plt
try:
    from vaspy import poscar
    from vaspy import tools
except ImportError:
    mypath = os.readlink(__file__) if os.path.islink(__file__) else __file__
    sys.path.append(os.path.dirname(os.path.abspath(mypath)))
    import poscar
    import tools


class LOCPOOT(poscar.POSCAR):

    '''class for LOCPOT format.

LOCPOT format is essentially same as CHGCAR but simpler.
'''

    def __init__(self, arg=None):
        pass

    def load_from_file(self, locpotfile):
        pass

    def average_along_axis(self, axis_name):
        axis_name = axis_name.capitalize()
        pass

    def max_along_axis(self, axis_name):
        axis_name = axis_name.capitalize()
        pass

    def min_along_axis(self, axis_name):
        axis_name = axis_name.capitalize()
        pass

    def mid_along_axis(self, axis_name):
        axis_name = axis_name.capitalize()
        pass

    def plot_potential_along_axis(self, axis_name):
        axis_name = axis_name.capitalize()
        pass
