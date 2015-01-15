# -*- coding: utf-8 -*-

from distutils.core import setup
import sys

reqpkgs = ['numpy']
if sys.hexversion <= 0x20600f0 or 0x30000f0 <= sys.hexversion < 0x30200f0:
    reqpkgs.append('argparse')

setup(name='VASPy',
      version='1.1.0b',
      author='Mao Kanno',
      author_email='0643839020@mail.ecc.u-tokyo.ac.jp',
      description='VASP postprocessing scripts',
      py_modules=[
          'vaspy.poscar',
          'vaspy.chgcar',
          'vaspy.doscar',
          'vaspy.outcar',
          'vaspy.procar',
          'vaspy.locpot',
          'vaspy.tools',
      ],
      scripts=['scripts/vaspy-chgcar.py',
               'scripts/vaspy-doscar.py',
               'scripts/vaspy-outcar.py',
               'scripts/vaspy-poscar.py',
               'scripts/vaspy-procar.py',
               ],
      requires=reqpkgs
      )
