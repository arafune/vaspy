"""mesh3D Module to provide class VASPGRID FFT-grid NG(X,Y,Z)F.

That is this is class VASPGRID is the parent class of CHGCAR,
LOCPOT, and ELFCAR.  The ELFCAR class has not yet implemented yet, though.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import IO, TYPE_CHECKING, Literal

import numpy as np

from vaspy import poscar, tools
from vaspy.tools import open_by_suffix

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import NDArray


class VASPGrid:
    """Class for VaspGrid used in CHGCAR, LOCPOT, ELFCAR, etc.

    General format of the file uses VaspGrid format::

        hBN-Cu                            # 1st line poscar.POSCAR[0]
        1.00000000000000                  # 2nd line poscar.POSCAR[1]
        6.762964    0.000000    0.000000  # 3rd line poscar.POSCAR[2]
        3.381482    5.856898    0.000000  # 4th line poscar.POSCAR[3]
        0.000000    0.000000   29.004836  # 5th line poscar.POSCAR[4]
        B    Cu   N    Si                 # 6th line poscar.POSCAR[5]
        7    21     7     6               # 7th line poscar.POSCAR[6]
        Direct                            # 8th line poscar.POSCAR[7]
        0.047680  0.261795  0.361962      # 9th line poscar.POSCAR[8]
        ....
                                            # the single blanc line
        240   240   288                     # number of grid mesh
        0.0000 0.0005 0.0002 0.0020 0.0001  # five columns in each line
        0.0030 0.0025 0.0001 0.0023 0.0003  #  ...
        ...                                 #  ...

    In CHGCAR file, information about "augmentation occupancies" is
    filled after ech charge, magnetization values. Currently,
    "augmentation occupancies" is ignored.

    In LOCPOT file, additional information (probably) about the ion
    itself by usually the value of unity is filled.  i.e. 0.1000E+1
    appears after each potential value by the times of
    self.ionnums. Currently this information is ignored.

    Attributes
    ----------
    poscar, grid, additional

    """

    def __init__(self, filename: str | Path = "", pickles: str = "") -> None:
        """Initialize."""
        self.poscar = poscar.POSCAR()
        self.grid = Grid3D()
        self.additional: list[str] = []
        if filename:
            self.load_file(open_by_suffix(str(filename)), pickles)

    def load_file(self, the_file: IO[str], pickles: str = "") -> None:
        """Construct the object from the file.

        Parameters
        ----------
        the_file: StringIO
            file

        pickles: str
            filename of griddata stored by np.save or np.savez

        """
        separator: str = ""
        tmp: list[str] = []
        griddata = ""
        # read POSCAR part
        line: str = the_file.readline()
        while not line.isspace():
            tmp.append(line.strip("\n"))
            line = the_file.readline()
        self.poscar.load_array(tmp)
        # read grid size and use it as separator
        separator = the_file.readline()
        self.grid.shape = tuple([int(string) for string in separator.split()])
        # Volumetric data
        if pickles:
            if Path(pickles).suffix == ".npy":
                self.grid.data = np.load(pickles)
            elif Path(pickles).suffix == ".npz":
                self.grid.data = np.load(pickles)["arr_0"]
            else:
                msg = "Check the volmetric data type"
                raise TypeError(msg)
            the_file.close()
            return
        griddata += next(the_file).replace("***********", "Nan")
        if self.grid.size % len(griddata.split()) == 0:
            lines_for_mesh = self.grid.size // len(griddata.split())
        else:
            lines_for_mesh = self.grid.size // len(griddata.split()) + 1
        for _ in range(lines_for_mesh - 1):  # read the first frame
            griddata += next(the_file).replace("***********", "Nan")
        section = "grid"
        for line in the_file:
            if section == "aug":
                if separator in line:
                    for _ in range(lines_for_mesh):
                        griddata += next(the_file).replace("***********", "Nan")
                    section = "grid"
                elif "augmentation occupancies " in line:
                    pass  # Used for CHGCAR, not LOCPOT. not implemented
                else:
                    pass  # Used for CHGCAR, not LOCPOT. not implemented
            elif section == "grid":
                if "augmentation occupancies " in line:
                    section = "aug"
                elif separator in line:
                    for _ in range(lines_for_mesh):
                        griddata += next(the_file).replace("***********", "Nan")
                else:
                    # for unused data stored in LOCPOT
                    self.additional.extend(line.split())
        self.grid.data = np.fromstring(griddata, dtype=float, sep=" ")
        the_file.close()

    def __str__(self) -> str:
        """Return th string representation.

        Returns
        -------
        str
            a string representation of VASPGrid object

        """
        poscarstr = self.poscar.str_short()
        meshstr = self.grid.__str__()
        return poscarstr + meshstr + "\n"

    def save(self, filename: str) -> None:
        """Save object as the same file-style.

        Parameters
        ----------
        filename: str
            file name

        """
        the_file: IO[str]
        with Path(filename).open(mode="w", newline="\n") as the_file:
            the_file.write(str(self))

    def frame(self, frame_i: int) -> VASPGrid:
        """Return VASPGrid object for only frame_i th frame.

        Parameters
        ----------
        frame_i: int
            frame index

        """
        output_vaspgrid = VASPGrid()
        output_vaspgrid.poscar = self.poscar
        output_vaspgrid.grid = self.grid.frame(frame_i)
        return output_vaspgrid

    def merge(self, other: VASPGrid) -> VASPGrid:
        """Add density data.

        Add two VASPGrid object, but POSCAR part remains same as original.
        Use this method to calculate Bader analysis, for example.

        Parameters
        ----------
        other: VASPGrid
            Addition VaspGrid object

        Returns
        -------
        VASPGrid
            Rusultant by summing two grid values

        """
        add_grid = copy.deepcopy(self)
        try:
            add_grid.grid.data = self.grid.data + other.grid.data
        except ValueError as v_err:
            msg = "The mesh shapes are different each other"
            raise RuntimeError(msg) from v_err
        return add_grid

    def __add__(self, other: VASPGrid) -> VASPGrid:
        """Add both density and atom position.

        x.__add__(y) <=> x + y

        Parameters
        ----------
        other: VASPGrid
            Addition VaspGrid object

        Returns
        -------
        Grid3D
            Rusultant by summing two grid values and poscar is also added.

        """
        add_grid = copy.deepcopy(self)
        add_grid.poscar = self.poscar + other.poscar
        try:
            add_grid.grid.data = self.grid.data + other.grid.data
        except ValueError as v_err:
            msg = "The mesh shapes are different each other"
            raise RuntimeError(msg) from v_err
        return add_grid

    def __sub__(self, other: VASPGrid) -> VASPGrid:
        """Subtract the density.

        x.__sub__(y) <=> x - y

        Parameters
        ----------
        other: VASPGrid
            difference VASPGrid object

        Returns
        -------
        Grid3D
            Resultant by difference between two objects.

        Note:
        ----
        The resultant grid data is the difference between two objects,
        of course. On the other hand, the atom position information
        unchange by this method.  Use the 'minuend' object.  The atom
        information in subrtrahend object is totally ignored.

        """
        diff_grid = copy.deepcopy(self)
        try:
            diff_grid.grid.data = self.grid.data - other.grid.data
        except ValueError as v_err:
            msg = "The mesh shapes are different each other"
            raise RuntimeError(msg) from v_err
        return diff_grid


class Grid3D:
    """Class for NG(X,Y,Z)F in VASP.

    This class is used chg_array in CHGCAR, Potential in LOCPOT,
    electron localization function (ELF) in ELFCAR

    Attributes
    ----------
    size: tuple
        number of mesh in the single frame
    n_frame: int
        number of frames
    shape: tuple
        shape[0], shape[1], shape[2]
    data: NDArray
        1D-list or 1D-numpy array.
        The length of grid is shape[0] * shape[1] * shape[2]

    """

    def __init__(
        self,
        shape: tuple[int, ...] | NDArray[np.int_] = (0, 0, 0),
        data: Sequence[float] | None = None,
    ) -> None:
        """Initialize."""
        self.shape = shape
        if data is None:
            self.data: NDArray[np.float64] = np.array([])
        else:
            self.data = np.array(data)

    @property
    def size(self) -> int:
        """Return the number of meshes in the frame."""
        return self.shape[0] * self.shape[1] * self.shape[2]

    @property
    def n_frame(self) -> int:
        """Return the number of grid frames."""
        return divmod(self.data.size, self.size)[0]

    def frame(self, frame_i: int) -> Grid3D:
        """Return the i-th frame.

        Parameters
        ----------
        frame_i:int
            frame index

        """
        assert frame_i < self.n_frame
        dest = copy.deepcopy(self)
        dest.data = self.data.reshape(self.n_frame, self.size)[frame_i]
        return dest

    def slice(
        self,
        position: int,
        axis: Literal["x", "y", "z"] = "z",
        frame_i: int = 0,
    ) -> NDArray[np.float64]:
        """Parameters
        ----------
        axis: str
            'x', 'y', or 'z'.  Case insensitive.
        position: int
            position for slice
        frame_i: int
            frame index (0-3)

        Return:
        ------
        NDArray
            2D numpy array that sliced from 3D mesh data.

        """
        griddata = self.data[frame_i * self.size : (frame_i + 1) * self.size]
        if axis in ("x", "X"):
            return griddata.reshape(self.shape[2], self.shape[1], self.shape[0])[
                :,
                :,
                position,
            ]
        if axis in ("y", "Y"):
            return griddata.reshape(self.shape[2], self.shape[1], self.shape[0])[
                :,
                position,
                :,
            ]
        if axis in ("z", "Z"):
            return griddata.reshape(self.shape[2], self.shape[1], self.shape[0])[
                position,
                :,
                :,
            ]
        msg = 'axis must be "x", "y" or "z".'
        raise RuntimeError(msg)

    def integrate(
        self,
        axis: Literal["X", "Y", "Z", "x", "y", "z"],
        from_coor: int | None = None,
        to_coor: int | None = None,
        frame_i: int = 0,
    ) -> NDArray[np.float64]:
        """Return 2D data integrated occupacy along the 'axis'.

        Integration range can be specified by from_coor and to_coor.
        If not specified, from 0 to end.

        Parameters
        ----------
        axis: str
            'x', 'y', or 'z'.  Case insensitive
        from_coor: int
            'from' value of range of interval integration
        to_coor: int
            'to' value of range interval integration
        frame_i: int
            frame index

        Return
        ------
        NDArray
            2D numpy array that integrated from 3D mesh data

        """
        griddata = self.data[frame_i * self.size : (frame_i + 1) * self.size]
        if axis in ("x", "X"):
            return np.sum(
                griddata.reshape(self.shape[2], self.shape[1], self.shape[0])[
                    :,
                    :,
                    from_coor:to_coor,
                ],
                axis=2,
            )
        if axis in ("y", "Y"):
            return np.sum(
                griddata.reshape(self.shape[2], self.shape[1], self.shape[0])[
                    :,
                    from_coor:to_coor,
                    :,
                ],
                axis=1,
            )
        if axis in ("z", "Z"):
            return np.sum(
                griddata.reshape(self.shape[2], self.shape[1], self.shape[0])[
                    from_coor:to_coor,
                    :,
                    :,
                ],
                axis=0,
            )
        msg = "incorrect axis"
        raise ValueError(msg)

    def __str__(self) -> str:
        """Return as string object.

        x.__str__() <=> str(x)

        Returns
        -------
        str
            a string representation of VASPGrid object

        """
        outputstr: str = ""
        mesharray = self.data.reshape(self.n_frame, self.size)
        for tmp in mesharray:
            output = []
            outputstr += f"\n  {self.shape[0]}  {self.shape[1]}  {self.shape[2]}\n"
            for array in tools.each_slice(tmp, 5):
                output.append(
                    "".join(f"  {i:18.11E}" for i in array if i is not None),
                )
            outputstr += "\n".join(output)
        return outputstr + "\n"

    def average_along_axis(
        self,
        axis_name: Literal["X", "Y", "Z", "x", "y", "z"],
        frame_i: int = 0,
    ) -> NDArray[np.float64]:
        """Calculate average value of potential along 'axis'.

        Parameters
        ----------
        axis_name: str
            'X', 'Y', or 'Z'
        mode: int, optional (default is 0)
            select data by integer
        frame_i: int
            frame index

        Returns
        -------
        NDArray
            average value along the axis

        """
        data: NDArray[np.float64] = self.data[
            frame_i * self.size : (frame_i + 1) * self.size
        ].reshape((self.shape[2], self.shape[1], self.shape[0]))
        if axis_name in ("x", "X"):
            data = np.average(np.average(np.transpose(data, (2, 0, 1)), axis=2), axis=1)
        elif axis_name in ("y", "Y"):
            data = np.average(np.average(np.transpose(data, (1, 0, 2)), axis=2), axis=1)
        elif axis_name in ("z", "Z"):
            data = np.average(np.average(data, axis=2), axis=1)
        else:
            msg = "Wrong axis name set"
            raise ValueError(msg)
        return data

    def min_along_axis(
        self,
        axis_name: Literal["X", "Y", "Z", "x", "y", "z"],
        frame_i: int = 0,
    ) -> NDArray[np.float64]:
        """Calculate minimum value of potential along 'axis'.

        Parameters
        ----------
        axis_name: str
            'X', 'Y', or 'Z'
        mode: int, optional (default is 0)
            select data by integer
        frame_i: int
            frame index

        Returns
        -------
        NDArray
            minimum value along the axis

        """
        data: NDArray[np.float64] = self.data[
            frame_i * self.size : (frame_i + 1) * self.size
        ].reshape((self.shape[2], self.shape[1], self.shape[0]))
        if axis_name in ("x", "X"):
            data = np.min(np.min(np.transpose(data, (2, 0, 1)), axis=2), axis=1)
        elif axis_name in ("y", "Y"):
            data = np.min(np.min(np.transpose(data, (1, 0, 2)), axis=2), axis=1)
        elif axis_name in ("z", "Z"):
            data = np.min(np.min(data, axis=2), axis=1)
        else:
            msg = "Wrong axis name set"
            raise ValueError(msg)
        return data

    def max_along_axis(
        self,
        axis_name: Literal["X", "Y", "Z", "x", "y", "z"],
        frame_i: int = 0,
    ) -> NDArray[np.float64]:
        """Calculate maximum value of potential along 'axis'.

        Parameters
        ----------
        axis_name: str
            'X', 'Y', or 'Z'
        mode: int, optional (default is 0)
            select data by integer
        frame_i: int
            frame index

        Returns
        -------
        NDArray
            maximum value along the axis

        """
        data: NDArray[np.float64] = self.data[
            frame_i * self.size : (frame_i + 1) * self.size
        ].reshape((self.shape[2], self.shape[1], self.shape[0]))
        if axis_name in ("x", "X"):
            data = np.max(np.max(np.transpose(data, (2, 0, 1)), axis=2), axis=1)
        elif axis_name in ("y", "Y"):
            data = np.max(np.max(np.transpose(data, (1, 0, 2)), axis=2), axis=1)
        elif axis_name in ("z", "Z"):
            data = np.max(np.max(data, axis=2), axis=1)
        else:
            msg = "Wrong axis name set"
            raise ValueError(msg)
        return data

    def median_along_axis(
        self,
        axis_name: Literal["X", "Y", "Z", "x", "y", "z"],
        frame_i: int = 0,
    ) -> NDArray[np.float64]:
        """Calculate median value of potential along 'axis'.

        Parameters
        ----------
        axis_name: str
            'X', 'Y', or 'Z'
        mode: int, optional (default is 0)
            select data by integer
        frame_i: int
            frame index

        Returns
        -------
        NDArray
            median value along the axis

        """
        data: NDArray[np.float64] = self.data[
            frame_i * self.size : (frame_i + 1) * self.size
        ].reshape((self.shape[2], self.shape[1], self.shape[0]))
        if axis_name in ("x", "X"):
            data = np.median(np.median(np.transpose(data, (2, 0, 1)), axis=2), axis=1)
        elif axis_name in ("y", "Y"):
            data = np.median(np.median(np.transpose(data, (1, 0, 2)), axis=2), axis=1)
        elif axis_name in ("z", "Z"):
            data = np.median(np.median(data, axis=2), axis=1)
        else:
            msg = "Wrong axis name set"
            raise ValueError(msg)
        return data
