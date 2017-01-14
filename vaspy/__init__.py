# -*- coding: utf-8 -*-
"""
   VASPy
   ======

   modules for VASP pre/post-process
"""
import os.path
import re
__all__ = ['chgcar',
           'doscar',
           'outcar',
           'poscar',
           'procar',
           'locpot',
           'tools',
           ]

from vaspy import *


def load(filename, mode=None):
    ''' .. py:function:: load(filename, mode)

    load file

    Parameters
    ----------

    filename: str
        filename

    mode: str
       optional argument mode

    Notes
    -----
        * 'mode' is the optional argument, because this function judges
    the file mode mainly from its name
        * mode is poscar, outcar, chgcar, procar (case insensitive)
    '''
    from . import poscar, outcar, chgcar, doscar, locpot, procar
    filenamebase = os.path.basename(filename).lower()
    if isinstance(mode, str):
        mode = mode.lower()
    if re.search(r'poscar|contcar', filenamebase) or mode == 'poscar':
        output = poscar.POSCAR(filename)
    elif re.search(r'outcar', filenamebase) or mode == 'outcar':
        output = outcar.OUTCAR(filename)
    elif re.search(r'chgcar|parchg', filenamebase) or mode == 'chgcar':
        output = chgcar.CHGCAR(filename)
    elif re.search(r'procar', filenamebase) or mode == 'procar':
        output = procar.PROCAR(filename)
    elif re.search(r'locpot', filenamebase) or mode == 'locpot':
        output = locpot.LOCPOT(filename)
    elif re.search(r'doscar', filenamebase) or mode == 'doscar':
        output = doscar.DOSCAR(filename)
    elif re.search(r'vasp', filenamebase):
        output = poscar.POSCAR(filename)
    else:
        raise RuntimeError("The loding mode cannot be identified!  Set 'mode'")
    return output
