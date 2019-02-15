# -*- coding: utf-8 -*-
"""
VASPy
======

modules for VASP pre/post-process
"""
import logging
import os.path
import re

from vaspy import *

logging.disable(logging.DEBUG)
__all__ = [
    'chgcar', 'doscar', 'outcar', 'poscar', 'procar', 'locpot', 'eigenval',
    'wavecar', 'mesh3d', 'const', 'vsim_asc', 'tools'
]


def load(filename, mode=None, additional=None):
    '''
    load file.  Guess the file type by the filename.
    Use `mode` option to specify the file type.

    Parameters
    ----------

    filename: str
        filename

    mode: str, optional
       optional argument mode

    Notes
    -----
        * 'mode' is the optional argument, because this function judges
    the file mode mainly from its name
        * mode is poscar, outcar, chgcar, procar (case insensitive)
    '''
    from . import poscar, outcar, chgcar, doscar, locpot, procar
    from . import eigenval, wavecar

    filenamebase = os.path.basename(filename).lower()
    if isinstance(mode, str):
        mode = mode.lower()
    if re.search(r'poscar|contcar', filenamebase) or mode == 'poscar':
        return poscar.POSCAR(filename)
    elif re.search(r'outcar', filenamebase) or mode == 'outcar':
        return outcar.OUTCAR(filename)
    elif re.search(r'chgcar|parchg', filenamebase) or mode == 'chgcar':
        return chgcar.CHGCAR(filename, pickleddata=additional)
    elif re.search(r'procar', filenamebase) or mode == 'procar':
        return procar.PROCAR(filename)
    elif re.search(r'locpot', filenamebase) or mode == 'locpot':
        return locpot.LOCPOT(filename, pickleddata=additional)
    elif re.search(r'doscar', filenamebase) or mode == 'doscar':
        return doscar.DOSCAR(filename)
    elif re.search(r'vasp', filenamebase):
        return poscar.POSCAR(filename)
    elif re.search(r'eigenval', filenamebase) or mode == 'eigenval':
        return eigenval.EIGENVAL(filename)
    elif re.search(r'wavecar', filenamebase) or mode == 'wavecar':
        return wavecar.WAVECAR(filename)
    else:
        raise RuntimeError("The loding mode cannot be identified!  Set 'mode'")
