#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module provides LOCPOT class."""
#
import sys
from typing import Optional

import numpy as np

from vaspy import mesh3d

try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.stderr.write("Install matplotlib, or you cannot use methods relating to draw\n")


class LOCPOT(mesh3d.VASPGrid):
    """Class for LOCPOT.

    LOCPOT format is essentially same as that of CHGCAR

    """

    def __init__(
        self, filename: Optional[str] = None, pickles: Optional[str] = None
    ) -> None:
        """Initialize."""
        super(LOCPOT, self).__init__(filename, pickles)

    def plot_potential_along_axis(
        self, axis_name: str, frame: int = 0, save: Optional[str] = None
    ) -> None:  # FIXME!!
        """Plot potential curve along the axis.

        Parameters
        ----------
        axis_name: str
            the name of the axis (X, Y, or Z)
        frame: int, optional  (default is 0)
            'select frame potential' (very optional)

        """
        axis_name = axis_name.capitalize()
        axes_length = self.poscar.axes_lengthes
        if axis_name == "X":
            horizontal_axis = np.linspace(0, axes_length[0], self.grid.shape[0])
            plt.clf()
            plt.xlim(xmax=axes_length[0])
        elif axis_name == "Y":
            horizontal_axis = np.linspace(0, axes_length[1], self.grid.shape[1])
            plt.clf()
            plt.xlim(xmax=axes_length[1])
        elif axis_name == "Z":
            horizontal_axis = np.linspace(0, axes_length[2], self.grid.shape[2])
            plt.clf()
            plt.xlim(xmax=axes_length[2])
        else:
            raise ValueError("Wrong axis name")
        y_average = self.grid.average_along_axis(axis_name, frame)
        y_max = self.grid.max_along_axis(axis_name, frame)
        y_min = self.grid.min_along_axis(axis_name, frame)
        y_median = self.grid.median_along_axis(axis_name, frame)
        plt.plot(horizontal_axis, y_average, label="average")
        plt.plot(horizontal_axis, y_max, label="max")
        plt.plot(horizontal_axis, y_min, label="min")
        plt.plot(horizontal_axis, y_median, label="median")
        xlabel = "Position along " + axis_name + "-axis (A)"
        title = "Local potential (" + axis_name + ")"
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel("Energy  ( eV )")
        if save:
            plt.savefig(save, dpi=600, format="png")
        plt.show()


# LVHAR-tag:  This tag is available in VASP.5.2.12 and newer version.
# It determines whether the total local potential (file LOCPOT ) contains
# the entire local potential (ionic plus Hartree plus exchange correlation)
# or the electrostatic contributions only (ionic plus Hartree). Note that
# in VASP.5.2.12, the default is to write the entire local potential,
# including the exchange correlation potential.
#
# Memo: I (RA) don't know why the potential data are stored twice
# (exactly same data), when LVHAR = .TRUE. and LVTOT is not set.
