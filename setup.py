# -*- coding: utf-8 -*-

from distutils.core import setup
import sys

reqpkgs = ['numpy']
if sys.hexversion <= 0x20600f0 or 0x30000f0 <= sys.hexversion < 0x30200f0:
    reqpkgs.append('argparse')

setup(name='vaspy',
      version='1.0.1',
      author='Mao Kanno',
      author_email='0643839020@mail.ecc.u-tokyo.ac.jp',
      description='vasp postprocessing scripts',
      py_modules=[
          'vaspy.poscar',
          'vaspy.chgcar',
          'vaspy.doscar',
          'vaspy.outcar',
          'vaspy.procar',
          'vaspy.tools',
          ],
      scripts=['scripts/chgcar.py',
               'scripts/doscar.py',
               'scripts/outcar.py',
               'scripts/poscar.py',
               'scripts/procar.py',
               ],
      requires=reqpkgs
      )
