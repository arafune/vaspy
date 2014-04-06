# -*- coding: utf-8 -*-

from distutils.core import setup
import sys

reqpkgs = ['numpy']
if sys.hexversion <= 0x20600f0:
    reqpkgs.append('argparse')

setup(name='vaspy',
      version='0.1.0a',
      py_modules=[
          'vaspy.poscar',
          'vaspy.chgcar',
          'vaspy.doscar',
          'vaspy.outcar',
          'vaspy.procar',
          'vaspy.tools',
          ],
      requires=reqpkgs
      )
