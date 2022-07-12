# -*- coding: utf-8 -*-
"""Setup script for vaspy
"""
import sys
from distutils.core import setup

from vaspy import __version__

reqpkgs = ["numpy", "matplotlib"]
if sys.hexversion < 0x20700F0 or 0x30000F0 <= sys.hexversion < 0x30200F0:
    reqpkgs.append("argparse")

setup(
    name="VASPy",
    version=__version__,
    author="Mao Kanno",
    author_email="mkanno.t.k@13.alumni.u-tokyo.ac.jp",
    maintainer="Ryuichi Arafune",
    maintainer_email="ARAFUNE.Ryuichi@nims.go.jp",
    description="VASP postprocessing scripts",
    py_modules=[
        "vaspy.bader",
        "vaspy.chgcar",
        "vaspy.const",
        "vaspy.doscar",
        "vaspy.eigenval",
        "vaspy.incar",
        "vaspy.locpot",
        "vaspy.mesh3d",
        "vaspy.outcar",
        "vaspy.poscar",
        "vaspy.procar",
        "vaspy.tools",
        "vaspy.utility",
        "vaspy.vsim_asc",
        "vaspy.wavecar",
        "vaspy.xdatcar",
    ],
    scripts=[
        "scripts/vaspy-chgcar.py",
        "scripts/vaspy-doscar.py",
        "scripts/vaspy-incar.py",
        "scripts/vaspy-outcar.py",
        "scripts/vaspy-poscar.py",
        "scripts/vaspy-procar.py",
        "scripts/vaspy-vsim.py",
        "scripts/vaspy-wavecar.py",
        "scripts/vaspy-xdatcar.py",
    ],
    requires=reqpkgs,
    data_files=[("etc", ["rmvaspy.py"])],
)
