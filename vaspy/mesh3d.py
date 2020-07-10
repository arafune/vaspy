"""mesh3D Module to provide class VASPGRID FFT-grid NG(X,Y,Z)F

That is this is class VASPGRID is the parent class of CHGCAR,
LOCPOT, and ELFCAR.  The ELFCAR class has not yet implemented yet, though.
"""

from __future__ import division, print_function

import copy
import os
from typing import Optional, IO, Tuple, List
from nptyping import NDArray

import numpy as np

from vaspy import poscar, tools
from vaspy.tools import open_by_suffix


class VASPGrid(object):
    # Todo: Use Composite pattern!!!
    # VASPGrid should consists of POSCAR and Mesh3D object!!
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
                                            # the single blanck line
        240   240   288                     # number of gridmesh
        0.0000 0.0005 0.0002 0.0020 0.0001  # five columns in each line
        0.0030 0.0025 0.0001 0.0023 0.0003  #  ...
        ...                                 #  ...

    In CHGCAR file, information about "augmentation occupacies" is
    filled after ech charge, magnetization values. Currently,
    "augmentation occupacies" is ignored.

    In LOCPOT file, additional information (probably) about the ion
    itself by usually the value of unity is filled.  i.e. 0.1000E+1
    appears after each potential value by the times of
    self.ionnums. Currently this information is ignored.

    Attributes
    -----------------
    poscar, grid, additional

    """

    def __init__(
        self, filename: Optional[str] = None, pickleddata: Optional[str] = None
    ) -> None:
        """Initialize."""
        self.poscar = poscar.POSCAR()
        self.grid = Grid3D()
        self.additional = []
        if filename:
            self.load_file(open_by_suffix(filename), pickleddata)

    def load_file(self, thefile: IO, pickleddata: Optional[str] = None) -> None:
        """Construct the object from the file.

        Parameters
        ----------
        thefile: StringIO
            file

        pickleddata: str
            griddata stored by np.save or np.savez

        """
        separator: str = ""
        tmp = []
        griddata = ""
        # read POSCAR part
        line: str = thefile.readline()
        while not line.isspace():
            tmp.append(line.strip("\n"))
            line = thefile.readline()
        self.poscar.load_array(tmp)
        # read grid size and use it as separator
        separator = thefile.readline()
        self.grid.shape = tuple([int(string) for string in separator.split()])
        # Volumetric data
        if pickleddata:
            if os.path.splitext(pickleddata)[1] == ".npy":
                self.grid.data = np.load(pickleddata)
            elif os.path.splitext(pickleddata)[1] == ".npz":
                self.grid.data = np.load(pickleddata)["arr_0"]
            else:
                raise TypeError("Check the volmetric data type")
            thefile.close()
            return None
        griddata += next(thefile).replace("***********", "Nan")
        if self.grid.size % len(griddata.split()) == 0:
            lines_for_mesh = self.grid.size // len(griddata.split())
        else:
            lines_for_mesh = self.grid.size // len(griddata.split()) + 1
        for _ in range(lines_for_mesh - 1):  # read the first frame
            griddata += next(thefile).replace("***********", "Nan")
        section = "grid"
        for line in thefile:
            if section == "aug":
                if separator in line:
                    for _ in range(lines_for_mesh):
                        griddata += next(thefile).replace("***********", "Nan")
                    section = "grid"
                elif "augmentation occupancies " in line:
                    pass  # Used for CHGCAR, not LOCPOT. not implementd
                else:
                    pass  # Used for CHGCAR, not LOCPOT. not implementd
            elif section == "grid":
                if "augmentation occupancies " in line:
                    section = "aug"
                elif separator in line:
                    for _ in range(lines_for_mesh):
                        griddata += next(thefile).replace("***********", "Nan")
                else:
                    # for unused data stored in LOCPOT
                    self.additional.extend(line.split())
        self.grid.data = np.fromstring(griddata, dtype=float, sep=" ")
        thefile.close()

    def __str__(self) -> str:
        """String representation.

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
        -----------
        filename: str
            file name

        """
        thefile: IO
        try:  # Version safety
            thefile = open(filename, mode="w", newline="\n")
        except TypeError:
            thefile = open(filename, mode="wb")
        with thefile:
            thefile.write(str(self))

    def frame(self, frame_i: int) -> "VASPGrid":
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

    def merge(self, other: "VASPGrid") -> "VASPGrid":
        """Add density data.

        Add two VASPGrid object, but POSCAR part remains same as original.
        Use this method to calculate Bader analysis, for example.

        Parameters
        ---------------
        other: VASPGrid
            Addtion VaspGrid object

        Returns
        -----------
        VASPGrid
            Rusultant by summing two grid values

        """
        add_grid = copy.deepcopy(self)
        try:
            add_grid.grid.data = self.grid.data + other.grid.data
        except ValueError:
            raise RuntimeError("The mesh shapes are different each other")
        return add_grid

    def __add__(self, other: "VASPGrid") -> "VASPGrid":
        """Add both density and atom position.

        x.__add__(y) <=> x + y

        Parameters
        ---------------
        other: VASPGrid
            Addtion VaspGrid object

        Returns
        -----------
        Grid3D
            Rusultant by summing two grid values and poscar is also added.

        """
        add_grid = copy.deepcopy(self)
        add_grid.poscar = self.poscar + other.poscar
        try:
            add_grid.grid.data = self.grid.data + other.grid.data
        except ValueError:
            raise RuntimeError("The mesh shapes are different each other")
        return add_grid

    def __sub__(self, other: "VASPGrid") -> "VASPGrid":
        """Subtract the density.

        x.__sub__(y) <=> x - y
        Parameters
        ---------------
        other: VASPGrid
            difference VASPGrid object

        Returns
        ----------
        Grid3D
            Resultant by difference between two objects.

        Note
        --------
        The resultant grid data is the difference between two objects,
        of course. On the other hand, the atom position information
        unchange by this method.  Use the 'minuend' object.  The atom
        information in subrtrahend object is totally ignored.

        """
        diff_grid = copy.deepcopy(self)
        try:
            diff_grid.grid.data = self.grid.data - other.grid.data
        except ValueError:
            raise RuntimeError("The mesh shapes are different each other")
        return diff_grid


class Grid3D(object):
    """Class for NG(X,Y,Z)F in VASP.

    This class is used chg_array in CHGCAR, Potential in LOCPOT,
    electron localization function (ELF) in ELFCAR

    Attributes
    ----------
    size: tuple
        number of mesh in the single frame
    nframe: int
        number of frames
    shape: tuple
        shape[0], shape[1], shape[2]
    data: list or numpy.ndarray
        1D-list or 1D-numpy array.
        The length of grid is shape[0] * shape[1] * shape[2]

    """

    def __init__(
        self,
        shape: Tuple[int, int, int] = (0, 0, 0),
        data: Optional[List[float]] = None,
    ) -> None:
        """Initialize."""
        self.shape = shape
        if data is None:
            self.data: NDArray = np.array([])
        else:
            self.data = np.asarray(data)

    @property
    def size(self) -> int:
        """Return the number of meshes in the frame."""
        return self.shape[0] * self.shape[1] * self.shape[2]

    @property
    def nframe(self) -> int:
        """Return the number of grid frames."""
        return divmod(self.data.size, self.size)[0]

    def frame(self, frame_i: int) -> NDArray:
        """Return the i-th frame.

        Parameters
        -----------
        frame_i:int
            frame index

        """
        assert frame_i < self.nframe
        dest = copy.deepcopy(self)
        dest.data = self.data.reshape(self.nframe, self.size)[frame_i]
        return dest

    def slice(self, position, axis: str = "z", frame_i: int = 0) -> NDArray:
        """
        Parameters
        ----------
        axis: str
            'x', 'y', or 'z'.  Case insensitive.
        position: int
            position for slice

        Return
        ------
        numpy.ndarray
            2D numpy array that sliced from 3D mesh data.

        """
        griddata = self.data[frame_i * self.size : (frame_i + 1) * self.size]
        axis = axis.lower()
        if axis == "x":
            return griddata.reshape(self.shape[2], self.shape[1], self.shape[0])[
                :, :, position
            ]
        elif axis == "y":
            return griddata.reshape(self.shape[2], self.shape[1], self.shape[0])[
                :, position, :
            ]
        elif axis == "z":
            return griddata.reshape(self.shape[2], self.shape[1], self.shape[0])[
                position, :, :
            ]
        else:
            raise RuntimeError('axis must be "x", "y" or "z".')

    def integrate(
        self,
        axis: str,
        from_coor: Optional[int] = None,
        to_coor: Optional[int] = None,
        frame_i: int = 0,
    ) -> NDArray:
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

        Return
        ------
        numpy.ndarray
            2D numpy array that integrated from 3D mesh data

        """
        griddata = self.data[frame_i * self.size : (frame_i + 1) * self.size]
        axis = axis.lower()
        if axis == "x":
            return np.sum(
                griddata.reshape(self.shape[2], self.shape[1], self.shape[0])[
                    :, :, from_coor:to_coor
                ],
                axis=2,
            )
        elif axis == "y":
            return np.sum(
                griddata.reshape(self.shape[2], self.shape[1], self.shape[0])[
                    :, from_coor:to_coor, :
                ],
                axis=1,
            )
        elif axis == "z":
            return np.sum(
                griddata.reshape(self.shape[2], self.shape[1], self.shape[0])[
                    from_coor:to_coor, :, :
                ],
                axis=0,
            )
        else:
            raise ValueError("incorrect axis")

    def __str__(self) -> str:
        """Return as string object.

        x.__str__() <=> str(x)

        Returns
        -------
        str
            a string representation of VASPGrid object

        """
        outputstr = ""
        mesharray = self.data.reshape(self.nframe, self.size)
        for tmp in mesharray:
            output = []
            outputstr += "\n  {0}  {1}  {2}\n".format(
                self.shape[0], self.shape[1], self.shape[2]
            )
            for array in tools.each_slice(tmp, 5):
                output.append(
                    "".join("  {0:18.11E}".format(i) for i in array if i is not None)
                )
            outputstr += "\n".join(output)
        return outputstr + "\n"

    def average_along_axis(self, axis_name: str, frame_i: int = 0) -> NDArray:
        """Calculate average value of potential along 'axis'.

        Parameters
        ----------
        axis_name: str
            'X', 'Y', or 'Z'
        mode: int, optional (default is 0)
            select data by integer

        Returns
        -------
        numpy.ndarray
            average value along the axis

        """
        axis_name = axis_name.capitalize()
        data = self.data[frame_i * self.size : (frame_i + 1) * self.size].reshape(
            (self.shape[2], self.shape[1], self.shape[0])
        )
        if axis_name == "X":
            data = np.average(np.average(np.transpose(data, (2, 0, 1)), axis=2), axis=1)
        elif axis_name == "Y":
            data = np.average(np.average(np.transpose(data, (1, 0, 2)), axis=2), axis=1)
        elif axis_name == "Z":
            data = np.average(np.average(data, axis=2), axis=1)
        else:
            raise "Wrong axis name set"
        return data

    def min_along_axis(self, axis_name: str, frame_i: int = 0) -> NDArray:
        """Calculate minimum value of potential along 'axis'.

        Parameters
        -----------
        axis_name: str
            'X', 'Y', or 'Z'
        mode: int, optional (default is 0)
            select data by integer

        Returns
        -------
        numpy.ndarray
            minimum value along the axis

        """
        axis_name = axis_name.capitalize()
        data = self.data[frame_i * self.size : (frame_i + 1) * self.size].reshape(
            (self.shape[2], self.shape[1], self.shape[0])
        )
        if axis_name == "X":
            data = np.min(np.min(np.transpose(data, (2, 0, 1)), axis=2), axis=1)
        elif axis_name == "Y":
            data = np.min(np.min(np.transpose(data, (1, 0, 2)), axis=2), axis=1)
        elif axis_name == "Z":
            data = np.min(np.min(data, axis=2), axis=1)
        else:
            raise "Wrong axis name set"
        return data

    def max_along_axis(self, axis_name: str, frame_i: int = 0) -> NDArray:
        """Calculate maximum value of potential along 'axis'.

        Parameters
        -----------
        axis_name: str
            'X', 'Y', or 'Z'
        mode: int, optional (default is 0)
            select data by integer

        Returns
        -------
        numpy.ndarray
            maximum value along the axis

        """
        axis_name = axis_name.capitalize()
        data = self.data[frame_i * self.size : (frame_i + 1) * self.size].reshape(
            (self.shape[2], self.shape[1], self.shape[0])
        )
        if axis_name == "X":
            data = np.max(np.max(np.transpose(data, (2, 0, 1)), axis=2), axis=1)
        elif axis_name == "Y":
            data = np.max(np.max(np.transpose(data, (1, 0, 2)), axis=2), axis=1)
        elif axis_name == "Z":
            data = np.max(np.max(data, axis=2), axis=1)
        else:
            raise "Wrong axis name set"
        return data

    def median_along_axis(self, axis_name: str, frame_i: int = 0) -> NDArray:
        """Calculate median value of potential along 'axis'.

        Parameters
        -----------
        axis_name: str
            'X', 'Y', or 'Z'
        mode: int, optional (default is 0)
            select data by integer

        Returns
        -------
        numpy.ndarray
            median value along the axis

        """
        axis_name = axis_name.capitalize()
        data = self.data[frame_i * self.size : (frame_i + 1) * self.size].reshape(
            (self.shape[2], self.shape[1], self.shape[0])
        )
        if axis_name == "X":
            data = np.median(np.median(np.transpose(data, (2, 0, 1)), axis=2), axis=1)
        elif axis_name == "Y":
            data = np.median(np.median(np.transpose(data, (1, 0, 2)), axis=2), axis=1)
        elif axis_name == "Z":
            data = np.median(np.median(data, axis=2), axis=1)
        else:
            raise "Wrong axis name set"
        return data
