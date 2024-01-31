"""XDATCAR class."""

from __future__ import annotations

from typing import IO, TYPE_CHECKING

import numpy as np

from vaspy.poscar import PosCarHead
from vaspy.tools import open_by_suffix

if TYPE_CHECKING:
    from pathlib import Path


class XDATCAR(PosCarHead):
    """Class for XDATCAR.

    Attributes
    ----------
    configurations: list

    """

    def __init__(self, filename: str | Path = "") -> None:
        """Initialize.

        Parameters
        ----------
        filename: str | Path
            XDATCAR file name

        """
        super().__init__()
        self.configurations = []
        if filename:
            self.load_file(open_by_suffix(str(filename)))

    def load_file(self, the_file: IO[str]) -> None:
        """Parse PROCAR.

        Parameters
        ----------
        the_file: StringIO
            'XDATCAR' file

        """
        self.system_name = next(the_file).strip()
        self.scaling_factor = float(next(the_file).strip())
        self.cell_vecs[0] = np.array([float(x) for x in next(the_file).split()])
        self.cell_vecs[1] = np.array([float(x) for x in next(the_file).split()])
        self.cell_vecs[2] = np.array([float(x) for x in next(the_file).split()])
        self.atom_types = next(the_file).split()
        self.atomnums = [int(x) for x in next(the_file).split()]
        positions = []
        for line in the_file:
            if "Direct configuration=" in line:
                if positions:
                    self.configurations.append(positions)
                    positions = []
            else:
                position = np.array([float(x) for x in line.strip().split()])
                positions.append(position)
        self.configurations.append(positions)
        the_file.close()

    def __str__(self) -> str:
        """Return as str.

        Returns
        -------
        str
            a string representation of XDATCAR

        """
        tmp = self.system_name + "\n"
        tmp += f"        {self.scaling_factor}\n"
        for i in range(3):
            tmp += "      {:#.6f}   {:#.6f}    {:6f}\n".format(
                self.cell_vecs[i][0],
                self.cell_vecs[i][1],
                self.cell_vecs[i][2],
            )
        for element in self.atom_types:
            tmp += f"    {element}"
        tmp += "\n"
        for atomnum in self.atomnums:
            tmp += f"    {atomnum}"
        tmp += "\n"
        for frame_index, positions in enumerate(self.configurations):
            tmp += f"Direct configuration=    {frame_index + 1}\n"
            for position in positions:
                tmp += "    {:#.6f}    {:#.6f}    {:6f}\n".format(
                    position[0],
                    position[1],
                    position[2],
                )
        return tmp
