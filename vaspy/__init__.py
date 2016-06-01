# -*- coding: utf-8 -*-
"""
   VASPy
   ======

    modules made for VASP pre\postprocess
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

def load(filename, type=None):
    ''' .. py:function:: load(filename, type)

    load file

    :param str filename: filename

    :param str type: optional argument type

        * This is the optional argument, because this function judges
    the file type mainly from its name
        * type is poscar, outcar, chgcar, procar (case insensitive)
    '''
    from . import poscar, outcar, chgcar, doscar, locpot, procar
    filenamebase = os.path.basename(filename).lower()
    if isinstance(type, str):
        type = type.lower()
    if re.search(r'poscar|vasp', filenamebase) or type == 'poscar':
        output = poscar.POSCAR(filename)
    elif re.search(r'outcar', filenamebase) or type == 'outcar':
        output = outcar.OUTCAR(filename)
    elif re.search(r'chgcar', filenamebase) or type == 'chgcar':
        output = chgcar.CHGCAR(filename)
    elif re.search(r'procar', filenamebase) or type == 'procar':
        output = procar.PROCAR(filename)
    elif re.search(r'locpot', filenamebase) or type == 'locpot':
        output = locpot.LOCPOT(filename)
    elif re.search(r'doscar', filenamebase) or type == 'doscar':
        output = doscar.DOSCAR(filename)
    else:
        raise RuntimeError("The file type cannot be identified!  Set 'type'")        
    return output
        
        
