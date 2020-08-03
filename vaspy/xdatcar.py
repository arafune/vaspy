# -*- coding: utf-8 -*-
"""This module provide XDATCAR class."""

import numpy as np

import vaspy.poscar
from vaspy.tools import open_by_suffix
from pathlib import Path
from typing import IO, Optional, Union


class XDATCAR(vaspy.poscar.POSCAR_HEAD):
    """Class for XDATCAR.

    Attributes
    ----------
    configurations: list

    """

    def __init__(self, filename: Union[str, Path, None] = None) -> None:
        """Initialize.

        Parameters
        ----------
        arg: str
            XDATCAR file name

        """
        super(XDATCAR, self).__init__()
        self.configurations = []
        if filename:
            self.load_file(open_by_suffix(str(filename)))

    def load_file(self, thefile: IO[str]) -> None:
        """Parse PROCAR.

        Parameters
        ----------
        thefile: StringIO
            'XDATCAR' file

        """
        self.system_name = next(thefile).strip()
        self.scaling_factor = float(next(thefile).strip())
        self.cell_vecs[0] = np.array([float(x) for x in next(thefile).split()])
        self.cell_vecs[1] = np.array([float(x) for x in next(thefile).split()])
        self.cell_vecs[2] = np.array([float(x) for x in next(thefile).split()])
        self.atomtypes = next(thefile).split()
        self.atomnums = [int(x) for x in next(thefile).split()]
        positions = []
        for line in thefile:
            if "Direct configuration=" in line:
                if positions:
                    self.configurations.append(positions)
                    positions = []
            else:
                position = np.array([float(x) for x in line.strip().split()])
                positions.append(position)
        self.configurations.append(positions)
        thefile.close()

    def __str__(self) -> str:
        """Return as str.

        Returns
        -------
        str
            a string representation of XDATCAR

        """
        tmp = self.system_name + "\n"
        tmp += "        {}\n".format(self.scaling_factor)
        for i in range(3):
            tmp += "      {:#.6f}   {:#.6f}    {:6f}\n".format(
                self.cell_vecs[i][0], self.cell_vecs[i][1], self.cell_vecs[i][2]
            )
        for element in self.atomtypes:
            tmp += "    {}".format(element)
        tmp += "\n"
        for atomnum in self.atomnums:
            tmp += "    {}".format(atomnum)
        tmp += "\n"
        for frame_index, positions in enumerate(self.configurations):
            tmp += "Direct configuration=    {}\n".format(frame_index + 1)
            for position in positions:
                tmp += "    {:#.6f}    {:#.6f}    {:6f}\n".format(
                    position[0], position[1], position[2]
                )
        return tmp
