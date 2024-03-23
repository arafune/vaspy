"""VASPy.

modules for VASP pre/post-process
"""

from __future__ import annotations

import os.path
from pathlib import Path
import re
from typing import TypeAlias

from vaspy import (
    bader,
    chgcar,
    const,
    doscar,
    eigenval,
    incar,
    locpot,
    mesh3d,
    outcar,
    poscar,
    procar,
    tools,
    vsim_asc,
    wavecar,
)

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

VASP: TypeAlias = (
    poscar.POSCAR
    | incar.Incar
    | outcar.OUTCAR
    | chgcar.CHGCAR
    | procar.PROCAR
    | locpot.LOCPOT
    | doscar.DOSCAR
    | eigenval.EIGENVAL
    | bader.BaderACF
    | wavecar.WAVECAR
)


def load(filename: str, mode: str = "", additional: str = "") -> VASP:
    """Load files.

    Guess the file type by the filename.
    Use `mode` option to specify the file type.

    Parameters
    ----------
    filename: str
        filename
    mode: str, optional
        optional argument mode
    additional: str, optional
        additional arg to load the VASP file.


    Notes
    -----
        * 'mode' is the optional argument, because this function judges
            the file mode mainly from its name.

        * mode is poscar, outcar, chgcar, procar, locpot, doscar,
            eigenval, wavecar (case insensitive).

    """
    filenamebase = Path(filename).name.lower()
    mode = mode.lower()
    if re.search(r"poscar|contcar", filenamebase) or mode == "poscar":
        return poscar.POSCAR(filename)
    if re.search(r"incar", filenamebase) or mode == "incar":
        return incar.Incar(filename)
    if re.search(r"outcar", filenamebase) or mode == "outcar":
        return outcar.OUTCAR(filename)
    if re.search(r"chgcar|parchg", filenamebase) or mode == "chgcar":
        return chgcar.CHGCAR(filename, pickles=additional)
    if re.search(r"procar", filenamebase) or mode == "procar":
        return procar.PROCAR(filename)
    if re.search(r"locpot", filenamebase) or mode == "locpot":
        return locpot.LOCPOT(filename, pickles=additional)
    if re.search(r"doscar", filenamebase) or mode == "doscar":
        return doscar.DOSCAR(filename)
    if re.search(r"vasp", filenamebase):
        return poscar.POSCAR(filename)
    if re.search(r"eigenval", filenamebase) or mode == "eigenval":
        return eigenval.EIGENVAL(filename)
    if re.search(r"wavecar", filenamebase) or mode == "wavecar":
        return wavecar.WAVECAR(filename)
    if re.search(r"acf", filenamebase) or mode == "bader":
        return bader.BaderACF(filename)
    msg = "The loading mode cannot be identified!  Set 'mode'"
    raise RuntimeError(msg)
