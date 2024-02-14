#!/usr/bin/env python
"""LOCPOT class."""
#
from typing import TYPE_CHECKING, Literal

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from vaspy import mesh3d

if TYPE_CHECKING:
    from numpy.typing import NDArray


class LOCPOT(mesh3d.VASPGrid):
    """Class for LOCPOT.

    LOCPOT format is essentially same as that of CHGCAR

    """

    def __init__(self, filename: str = "", pickles: str = "") -> None:
        """Initialize."""
        super().__init__(filename, pickles)

    def plot_potential_along_axis(
        self,
        axis_name: Literal["X", "Y", "Z", "x", "y", "z"],
        ax: Axes | None = None,
        frame: int = 0,
        save: str = "",
        **kwargs,
    ) -> Axes:
        """Plot potential curve along the axis.

        Parameters
        ----------
        axis_name: str
            the name of the axis (X, Y, or Z)
        frame: int, optional  (default is 0)
            'select frame potential' (very optional)
        ax: Axes
            Matplotlib axes object.
        save: str | Path
            string or Path object to save the plot image.
        kwargs: Incomplete
            pass to plt.subplots

        """
        axes_length = self.poscar.axes_lengths

        if ax is None:
            _, ax = plt.subplots(figsize=kwargs.pop("figsize", (11, 5)))
        assert isinstance(ax, Axes)
        if axis_name in ("X", "x"):
            horizontal_axis: NDArray[np.float64] = np.linspace(
                0,
                axes_length[0],
                self.grid.shape[0],
            )
            ax.set_xlim(xmax=axes_length[0])
        elif axis_name in ("Y", "y"):
            horizontal_axis = np.linspace(0, axes_length[1], self.grid.shape[1])
            ax.set_xlim(xmax=axes_length[1])
        elif axis_name in ("Z", "z"):
            horizontal_axis = np.linspace(0, axes_length[2], self.grid.shape[2])
            ax.set_xlim(xmax=axes_length[2])
        else:
            msg = "Wrong axis name"
            raise ValueError(msg)
        y_average = self.grid.average_along_axis(axis_name, frame)
        y_max = self.grid.max_along_axis(axis_name, frame)
        y_min = self.grid.min_along_axis(axis_name, frame)
        y_median = self.grid.median_along_axis(axis_name, frame)
        ax.plot(horizontal_axis, y_average, label="average")
        ax.plot(horizontal_axis, y_max, label="max")
        ax.plot(horizontal_axis, y_min, label="min")
        ax.plot(horizontal_axis, y_median, label="median")
        x_label = "Position along " + axis_name + "-axis (A)"
        title = "Local potential (" + axis_name + ")"
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
        ax.set_title(title)
        ax.set_xlabel(x_label)
        ax.set_ylabel("Energy  ( eV )")
        if save:
            plt.savefig(save, dpi=600, format="png")
        return ax


# LVHAR-tag:  This tag is available in VASP.5.2.12 and newer version.
# It determines whether the total local potential (file LOCPOT ) contains
# the entire local potential (ionic plus Hartree plus exchange correlation)
# or the electrostatic contributions only (ionic plus Hartree). Note that
# in VASP.5.2.12, the default is to write the entire local potential,
# including the exchange correlation potential.
#
# Memo: I (RA) don't know why the potential data are stored twice
# (exactly same data), when LVHAR = .TRUE. and LVTOT is not set.
