# -*- coding: utf-8 -*-
"""Module for CHGCAR class.

translate from chgcar.rb in scRipt4VASP, 2014/2/26 master branch
"""
from __future__ import annotations

import copy
from typing import Optional, IO, Union, List
from pathlib import Path
from vaspy.mesh3d import VASPGrid
from vaspy.tools import open_by_suffix


class CHGCAR(VASPGrid):
    """Class for CHGCAR.

    An example of the first few lines of CHGCAR. ::

        hBN-Cu                            # 1st line poscar.POSCAR[0]
        1.00000000000000                  # 2nd line poscar.POSCAR[1]
        6.762964    0.000000    0.000000  # 3rd line poscar.POSCAR[2]
        3.381482    5.856898    0.000000  # 4th line poscar.POSCAR[3]
        0.000000    0.000000   29.004836  # 5th line poscar.POSCAR[4]
        B    Cu   N    Si                 # 6th line poscar.POSCAR[5]
        7    21     7     6               # 7th line poscar.POSCAR[6]
        Direct                            # 8th line poscar.POSCAR[7]
        0.047680  0.261795  0.361962      # 9th line poscar.POSCAR[8]
        ...
                                            # the single blanck line
        240   240   288                     # number of gridmesh
        0.0000 0.0005 0.0002 0.0020 0.0001  # five columns in each line
        0.0030 0.0025 0.0001 0.0023 0.0003  #  ...
        ...                                 #  ...


    Attributes
    ----------
    spin: int or list
        Represents spin character

    Notes
    -----
    Current verstion ignores "augmentation occupacies".

    """

    def __init__(
        self, filename: Union[str, Path] = "", pickles: Optional[str] = None
    ) -> None:
        """Initialize."""
        super(CHGCAR, self).__init__(None)
        self.spin: List[str] = [""]
        if filename:
            self.load_file(open_by_suffix(str(filename)), pickles)

    def load_file(self, thefile: IO[str], pickles: Optional[str] = None) -> None:
        """Parse CHGCAR file to construct CHGCAR object.

        Parameters
        ----------
        thefile: IO
            CHGCAR file

        """
        super(CHGCAR, self).load_file(thefile, pickles)
        if self.grid.nframe == 1:
            self.spin = [""]
        elif self.grid.nframe == 2:
            self.spin = ["up+down", "up-down"]
        elif self.grid.nframe == 4:
            self.spin = ["mT", "mX", "mY", "mZ"]
        else:
            raise RuntimeError("CHGCAR is correct?")
        thefile.close()

    def magnetization(self, direction: Optional[str] = None) -> CHGCAR:
        """Return CHGCAR for magnetization.

        For collinear spin-polarized calculations
        (``ISPIN=2`` but ``LSORBIT=.FALSE.``),
        two sets of data are found in CHGCAR file. The first set
        is the total charge density (spin-up plus spin-down),
        the second one the magnetization density (spin-up minus spin-down).

        For non-collinear spin-polarized calculations
        (``ISPIN=2`` and ``LSORBIT=.TRUE.``),
        CHGCAR file stores the total charge density and the
        magnetisation density in the x, y and z direction in this order.

        For collinear spinpolarized calculation the argument does
        not make a sense.  For non-collinear CHGCAR, direction
        should be one of 'x', 'y', 'z' and 't'

        Parameters
        ----------
        direction: str, optional
            specify x, y, z or t in noncollinear calculation.
            't' means the total.

        Returns
        -------
        CHGCAR
            CHGCAR of the spin-distribution

        """
        if len(self.spin) == 1:
            raise RuntimeError("This CHGCAR is not spinresolved version")
        dest = copy.deepcopy(self)
        if len(self.spin) == 2:
            dest.grid = dest.grid.frame(1)
            dest.spin = ["up-down"]
        elif len(self.spin) == 4:
            if direction is None or direction == "t":
                dest.grid = dest.grid.frame(0)
                dest.spin = ["mT"]
            if direction == "x":
                dest.grid = dest.grid.frame(1)
                dest.spin = ["mX"]
            elif direction == "y":
                dest.grid = dest.grid.frame(2)
                dest.spin = ["mY"]
            elif direction == "z":
                dest.grid = dest.grid.frame(3)
                dest.spin = ["mZ"]
        dest.grid.data = dest.grid.data.flatten()
        return dest

    def majorityspin(self) -> CHGCAR:
        """Return CHGCAR for majority spin.

        This method is for CHGCAR given by ``ISPIN=2`` but not-SOI
        calculations.

        Returns
        -------
        vaspy.chgcar.CHGCAR
            CHGCAR for the majority spin charge

        """
        assert len(self.spin) == 2, "This CHGCAR is not spinresolved version"
        dest = copy.deepcopy(self)
        tmp = dest.grid.data.reshape(2, self.grid.size)
        dest.grid.data = (tmp[0] + tmp[1]) / 2
        dest.spin = ["up"]
        return dest

    def minorityspin(self) -> CHGCAR:
        """Return CHGCAR for minority spin.

        This method is for CHGCAR given by ``ISPIN=2`` but not-SOI
        calculations.

        Returns
        ---------
        vaspy.chgcar.CHGCAR
            CHGCAR for the minority spin charge

        """
        assert len(self.spin) == 2, "This CHGCAR is not spinresolved version"
        dest = copy.deepcopy(self)
        tmp = dest.grid.data.reshape(2, self.grid.size)
        dest.grid.data = (tmp[0] - tmp[1]) / 2
        dest.spin = ["down"]
        return dest
