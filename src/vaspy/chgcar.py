"""Module for CHGCAR class.

translate from chgcar.rb in scRipt4VASP, 2014/2/26 master branch
"""

from __future__ import annotations

import copy
from typing import IO, TYPE_CHECKING

from vaspy.mesh3d import VASPGrid
from vaspy.tools import open_by_suffix

if TYPE_CHECKING:
    from pathlib import Path


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
                                            # the single blanc line
        240   240   288                     # number of grid mesh
        0.0000 0.0005 0.0002 0.0020 0.0001  # five columns in each line
        0.0030 0.0025 0.0001 0.0023 0.0003  #  ...
        ...                                 #  ...


    Attributes
    ----------
    spin: int or list
        Represents spin character

    Notes
    -----
    Current version ignores "augmentation occupancies".

    """

    def __init__(self, filename: str | Path = "", pickles: str = "") -> None:
        """Initialize."""
        super().__init__()
        self.spin: list[str] = [""]
        if filename:
            self.load_file(open_by_suffix(str(filename)), pickles)

    def load_file(self, the_file: IO[str], pickles: str = "") -> None:
        """Parse CHGCAR file to construct CHGCAR object.

        Parameters
        ----------
        the_file: IO
            CHGCAR file

        pickles: str
            filename of griddata stored by np.save or np.savez
        """
        super().load_file(the_file, pickles)
        if self.grid.n_frame == 1:
            self.spin = [""]
        elif self.grid.n_frame == 2:  # noqa: PLR2004
            self.spin = ["up+down", "up-down"]
        elif self.grid.n_frame == 4:  # noqa: PLR2004
            self.spin = ["mT", "mX", "mY", "mZ"]
        else:
            msg = "CHGCAR is correct?"
            raise RuntimeError(msg)
        the_file.close()

    def magnetization(self, direction: str = "") -> CHGCAR:
        """Return CHGCAR for magnetization.

        For collinear spin-polarized calculations
        (``ISPIN=2`` but ``LSORBIT=.FALSE.``),
        two sets of data are found in CHGCAR file. The first set
        is the total charge density (spin-up plus spin-down),
        the second one the magnetization density (spin-up minus spin-down).

        For non-collinear spin-polarized calculations
        (``ISPIN=2`` and ``LSORBIT=.TRUE.``),
        CHGCAR file stores the total charge density and the
        magnetization density in the x, y and z direction in this order.

        For collinear spin polarized calculation the argument does
        not make a sense.  For non-collinear CHGCAR, direction
        should be one of 'x', 'y', 'z' and 't'

        Parameters
        ----------
        direction: str, optional
            specify x, y, z or t in non-collinear calculation.
            't' means the total.

        Returns
        -------
        CHGCAR
            CHGCAR of the spin-distribution

        """
        if len(self.spin) == 1:
            msg = "This CHGCAR is not spin resolved version"
            raise RuntimeError(msg)
        dest = copy.deepcopy(self)
        if len(self.spin) == 2:  # noqa: PLR2004
            dest.grid = dest.grid.frame(1)
            dest.spin = ["up-down"]
        elif len(self.spin) == 4:  # noqa: PLR2004
            if direction == "t":
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

    def majority_spin(self) -> CHGCAR:
        """Return CHGCAR for majority spin.

        This method is for CHGCAR given by ``ISPIN=2`` but not-SOI
        calculations.

        Returns
        -------
        vaspy.chgcar.CHGCAR
            CHGCAR for the majority spin charge

        """
        assert (
            len(self.spin) == 2  # noqa: PLR2004
        ), "This CHGCAR is not spin resolved version"
        dest = copy.deepcopy(self)
        tmp = dest.grid.data.reshape(2, self.grid.size)
        dest.grid.data = (tmp[0] + tmp[1]) / 2
        dest.spin = ["up"]
        return dest

    def minority_spin(self) -> CHGCAR:
        """Return CHGCAR for minority spin.

        This method is for CHGCAR given by ``ISPIN=2`` but not-SOI
        calculations.

        Returns
        -------
        vaspy.chgcar.CHGCAR
            CHGCAR for the minority spin charge

        """
        assert (
            len(self.spin) == 2  # noqa: PLR2004
        ), "This CHGCAR is not spin resolved version"
        dest = copy.deepcopy(self)
        tmp = dest.grid.data.reshape(2, self.grid.size)
        dest.grid.data = (tmp[0] - tmp[1]) / 2
        dest.spin = ["down"]
        return dest
