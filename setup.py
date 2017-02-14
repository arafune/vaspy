# -*- coding: utf-8 -*-
'''Setup script for vaspy
'''
from distutils.core import setup
import sys

reqpkgs = ['numpy', 'matplotlib']
if sys.hexversion < 0x20700f0 or 0x30000f0 <= sys.hexversion < 0x30200f0:
    reqpkgs.append('argparse')

setup(name='VASPy',
      version='0.3.1',
      author='Mao Kanno',
      author_email='mkanno.t.k@13.alumni.u-tokyo.ac.jp',
      maintainer='Ryuichi Arafune',
      maintainer_email='ARAFUNE.Ryuichi@nims.go.jp',
      description='VASP postprocessing scripts',
      py_modules=[
          'vaspy.poscar',
          'vaspy.chgcar',
          'vaspy.doscar',
          'vaspy.outcar',
          'vaspy.procar',
          'vaspy.locpot',
          'vaspy.tools',
          'vaspy.mesh3d',
      ],
      scripts=['scripts/vaspy-chgcar.py',
               'scripts/vaspy-doscar.py',
               'scripts/vaspy-outcar.py',
               'scripts/vaspy-poscar.py',
               'scripts/vaspy-procar.py'],
      requires=reqpkgs,
      data_files=[('etc', ['rmvaspy.py'])])
