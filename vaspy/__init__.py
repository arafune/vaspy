# -*- coding: utf-8 -*-
"""VASPy.

modules for VASP pre/post-process
"""
import os.path
import re

from vaspy import (
    bader,
    chgcar,
    doscar,
    outcar,
    poscar,
    procar,
    locpot,
    eigenval,
    wavecar,
    mesh3d,
    const,
    vsim_asc,
    tools,
    incar,
)
from typing import Optional, Union, Any

__all__: list[str] = [
    "bader",
    "chgcar",
    "doscar",
    "outcar",
    "poscar",
    "procar",
    "locpot",
    "eigenval",
    "wavecar",
    "mesh3d",
    "const",
    "vsim_asc",
    "tools",
    "load",
    "utility",
    "incar",
]

__version__: str = "0.6.0"


def load(
    filename: str, mode: Optional[str] = None, additional: Union[bool, str, None] = None
) -> Any:
    """Load files.

    Guess the file type by the filename.
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
            the file mode mainly from its name.

        * mode is poscar, outcar, chgcar, procar, locpot, doscar,
            eigenval, wavecar (case insensitive).

    """

    filenamebase = os.path.basename(filename).lower()
    if isinstance(mode, str):
        mode = mode.lower()
    if re.search(r"poscar|contcar", filenamebase) or mode == "poscar":
        return poscar.POSCAR(filename)
    elif re.search(r"incar", filenamebase) or mode == "incar":
        return incar.Incar(filename)
    elif re.search(r"outcar", filenamebase) or mode == "outcar":
        return outcar.OUTCAR(filename)
    elif re.search(r"chgcar|parchg", filenamebase) or mode == "chgcar":
        return chgcar.CHGCAR(filename, pickles=additional)
    elif re.search(r"procar", filenamebase) or mode == "procar":
        return procar.PROCAR(filename)
    elif re.search(r"locpot", filenamebase) or mode == "locpot":
        return locpot.LOCPOT(filename, pickles=additional)
    elif re.search(r"doscar", filenamebase) or mode == "doscar":
        return doscar.DOSCAR(filename)
    elif re.search(r"vasp", filenamebase):
        return poscar.POSCAR(filename)
    elif re.search(r"eigenval", filenamebase) or mode == "eigenval":
        return eigenval.EIGENVAL(filename)
    elif re.search(r"wavecar", filenamebase) or mode == "wavecar":
        return wavecar.WAVECAR(filename)
    elif re.search(r"acf", filenamebase) or mode == "bader":
        return bader.BaderACF(filename)
    else:
        raise RuntimeError("The loding mode cannot be identified!  Set 'mode'")
