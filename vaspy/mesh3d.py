''' .. py:module:: mesh3D
mesh3D Module to provide class VASPGRID FFT-grid NG(X,Y,Z)F

That is this is class VASPGRID is the parent class of CHGCAR,
 LOCPOT, and ELFCAR.  The ELFCAR class has not yet implemented yet, though.
'''

from __future__ import division, print_function
import bz2
import copy
import os
import re
import numpy as np
from vaspy import poscar, tools

_RE_BLANK = re.compile(r'^[\s]*$')
_RE_AUG_OCC = re.compile(r'\baugmentation occupancies')


class VASPGrid(object):
    # Todo: Use Composite pattern!!!
    # VASPGrid should consists of POSCAR and Mesh3D object!!
    '''.. py:class::VaspGrid

    class for VaspGrid used in CHGCAR, LOCPOT, ELFCAR...

    General format of the file uses VaspGrid format::

           hBN-Cu                              # 1st line poscar.POSCAR[0]
           1.00000000000000                    # 2nd line poscar.POSCAR[1]
             6.762964    0.000000    0.000000  # 3rd line poscar.POSCAR[2]
             3.381482    5.856898    0.000000  # 4th line poscar.POSCAR[3]
             0.000000    0.000000   29.004836  # 5th line poscar.POSCAR[4]
           B    Cu   N    Si                   # 6th line poscar.POSCAR[5]
             7    21     7     6               # 7th line poscar.POSCAR[6]
           Direct                              # 8th line poscar.POSCAR[7]
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
    '''
    def __init__(self, filename=None):
        self.poscar = poscar.POSCAR()
        self.grid = Grid3D()
        self.additional = []
        if filename:
            self.load_from_file(filename)

    def load_from_file(self, filename):
        '''.. py:method:: load_from_file(filename)

        filename file to construct the object

        Parameters
        ---------------

        filename: str
            file name
        '''
        section = 'poscar'
        separator = None
        tmp = []
        griddata = []
        if os.path.splitext(filename)[1] == '.bz2':
            try:
                thefile = bz2.open(filename, mode='rt')
            except AttributeError:
                thefile = bz2.BZ2File(filename, mode='r')
        else:
            thefile = open(filename)
        with thefile:
            for line in thefile:
                line = line.rstrip('\n')
                if section == 'poscar':
                    if re.search(_RE_BLANK, line):
                        self.poscar.load_from_array(tmp)
                        section = 'define_separator'
                    else:
                        tmp.append(line)
                elif section == 'define_separator':
                    separator = line if separator is None else separator
                    if self.grid.shape == (0, 0, 0):
                        self.grid.shape = tuple([int(string) for string
                                                 in line.split()])
                    griddata.extend(
                        [float(i) for i in next(thefile).rstrip().replace(
                            '***********', 'Nan').split()])
                    if self.grid.size % len(griddata) == 0:
                        lines_for_mesh = self.grid.size // len(griddata)
                    else:
                        lines_for_mesh = self.grid.size // len(griddata) + 1
                    for _ in range(lines_for_mesh - 1):
                        griddata.extend([float(val) for val in
                                         next(thefile).rstrip().replace('***********',
                                                                        'Nan').split()])
                    section = 'grid'
                elif section == 'aug':
                    if separator in line:
                        for _ in range(lines_for_mesh):
                            griddata.extend([float(val) for val in
                                             next(thefile).rstrip().replace('***********',
                                                                            'Nan').split()])
                        section = 'grid'
                    elif "augmentation occupancies " in line:
                        pass  # Used for CHGCAR, not LOCPOT. not implementd
                    else:
                        pass  # Used for CHGCAR, not LOCPOT. not implementd
                elif section == 'grid':
                    if "augmentation occupancies " in line:
                        section = 'aug'
                    elif separator in line:
                        for _ in range(lines_for_mesh):
                            griddata.extend([float(val) for val in
                                             next(thefile).rstrip().replace('***********',
                                                                            'Nan').split()])
                    else:
                        # for unused data stored in LOCPOT
                        self.additional.extend(line.split())
            self.grid.data = np.array(griddata, dtype=np.float64)
            self.grid.num_frame = divmod(len(self.grid.data), self.grid.size)[0]

    def __str__(self):
        '''.. py:method:: __str__()

        x.__str__() <=> str(x)

        Returns
        -------

        str
            a string representation of VASPGrid object
        '''
        poscarstr = self.poscar.str_short()
        meshstr = self.grid.__str__()
        return poscarstr + meshstr + '\n'

    def save(self, filename):
        '''.. py:method:: save(filename)

        Save object as the same file-style

        Parameters
        -----------

        filename: str
            file name
        '''
        try:  # Version safety
            thefile = open(filename, mode='w', newline='\n')
        except TypeError:
            thefile = open(filename, mode='wb')
        with thefile:
            thefile.write(str(self))

    def merge(self, other):
        '''.. py:method:: __add__(other)

        x.__add__(y) <=> x + y

        Parameters
        ---------------

        other: VASPGrid
            Addtion VaspGrid object

        Returns
        -----------

        Grid3D
            Rusultant by summing two grid value

        '''
        add_grid = copy.deepcopy(self)
        try:
            add_grid.grid.data = self.grid.data + other.grid.data
        except ValueError:
            raise RuntimeError('The mesh shapes are different each other')
        return add_grid

    def __add__(self, other):
        '''.. py:method:: __add__(other)

        x.__add__(y) <=> x + y

        Parameters
        ---------------

        other: VASPGrid
            Addtion VaspGrid object

        Returns
        -----------

        Grid3D
            Rusultant by summing two grid value

        '''
        add_grid = copy.deepcopy(self)
        add_grid.poscar = self.poscar + other.poscar
        try:
            add_grid.grid.data = self.grid.data + other.grid.data
        except ValueError:
            raise RuntimeError('The mesh shapes are different each other')
        return add_grid

    def __sub__(self, other):
        '''.. py:method:: __sub__(other)

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
        '''

        diff_grid = copy.deepcopy(self)
        try:
            diff_grid.grid.data = self.grid.data - other.grid.data
        except ValueError:
            raise RuntimeError('The mesh shapes are different each other')
        return diff_grid


class Grid3D(object):
    '''.. py:class:: Mesh3D(size, data)
    Class for NG(X,Y,Z)F in VASP

    This class is used chg_array in CHGCAR, Potential in LOCPOT,
    electron localization function (ELF) in ELFCAR

    Parameters
    ----------

    shape: tuple
        shape[0], shape[1], shape[2]
    data: list or numpy.array
        1D-list or 1D-numpy array.
        The length of grid is shape[0] * shape[1] * shape[2]
    num_frame: int
        Number of grid frames
        for example, num_frame is 4 for CHGCAR included SOI
    '''
    def __init__(self, shape=(0, 0, 0), data=None):
        self.shape = shape
        if data is None:
            self.data = []
        else:
            self.data = np.asarray(data)
        try:
            self.num_frame = divmod(self.data.size, self.size)[0]
        except (ZeroDivisionError, AttributeError):
            self.num_frame = 0

    @property
    def size(self):
        '''Return the number of mesh in the frame'''
        return self.shape[0] * self.shape[1] * self.shape[2]

    def slice(self, axis, postition):
        '''.. py:method:: slice(axis, position)

        Parameters
        ----------

        axis: str
            'x', 'y', or 'z'.  Case insensitive.
        position: int
            position for slice

        Return
        ------

        numpy.array
            2D numpy array that sliced from 3D mesh data.
        '''
        pass

    def integrate(self, axis, from_coor, to_coor):
        '''.. py:method:: integrate(axis, from_coor, to_coor)
        Return 2D data integrated occupacy along the 'axis'
        from_coor to to_coor.

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

        numpy.array
            2D numpy array that integrated from 3D mesh data

        '''
        pass

    def __str__(self):
        '''.. py:method:: __str__()

        x.__str__() <=> str(x)

        Returns
        -------

        str
            a string representation of VASPGrid object
        '''
        outputstr = ''
        mesharray = self.data.reshape(self.data.size // self.size,
                                      self.size)
        for tmp in mesharray:
            output = []
            outputstr += '\n  {0}  {1}  {2}\n'.format(self.shape[0],
                                                      self.shape[1],
                                                      self.shape[2])
            for array in tools.each_slice(tmp, 5):
                output.append(''.join('  {0:18.11E}'.format(i)
                                      for i in array if i is not None))
            outputstr += '\n'.join(output)
        return outputstr + '\n'

    def average_along_axis(self, axis_name, mode=0):
        '''.. py:method:: average_along_axis(axis_name, mode)

        Calculate average value of potential along 'axis'

        Parameters
        ----------

        axis_name: str
             'X', 'Y', or 'Z'
        mode: int
             select data by integer

        Returns
        -------
        numpy.ndarray
            average value along the axis
        '''
        axis_name = axis_name.capitalize()
        data = self.data.reshape(self.num_frame,
                                 self.shape[2],
                                 self.shape[1],
                                 self.shape[0])[mode]
        if axis_name == 'X':
            data = np.average(np.average(
                np.transpose(data, (2, 0, 1)), axis=2), axis=1)
        elif axis_name == 'Y':
            data = np.average(np.average(
                np.transpose(data, (1, 0, 2)), axis=2), axis=1)
        elif axis_name == 'Z':
            data = np.average(np.average(data, axis=2), axis=1)
        else:
            raise 'Wrong axis name set'
        return data

    def min_along_axis(self, axis_name, mode=0):
        '''.. py:method:: min_along_axis(axis_name, potmode)

        Calculate minimum value of potential along 'axis'

        Parameters
        -----------

        axis_name: str
             'X', 'Y', or 'Z'
        mode: int
             select data by integer

        Returns
        -------
        numpy.ndarray
            minimum value along the axis
        '''
        axis_name = axis_name.capitalize()
        data = self.data.reshape(self.num_frame,
                                 self.shape[2],
                                 self.shape[1],
                                 self.shape[0])[mode]
        if axis_name == 'X':
            data = np.min(np.min(
                np.transpose(data, (2, 0, 1)), axis=2), axis=1)
        elif axis_name == 'Y':
            data = np.min(np.min(
                np.transpose(data, (1, 0, 2)), axis=2), axis=1)
        elif axis_name == 'Z':
            data = np.min(np.min(data, axis=2), axis=1)
        else:
            raise 'Wrong axis name set'
        return data

    def max_along_axis(self, axis_name, mode=0):
        '''.. py:method:: man_along_axis(axis_name, potmode)

        Calculate maximum value of potential along 'axis'

        Parameters
        -----------

        axis_name: str
             'X', 'Y', or 'Z'
        mode: int
             select data by integer

        Returns
        -------
        numpy.ndarray
            maximum value along the axis
        '''
        axis_name = axis_name.capitalize()
        data = self.data.reshape(self.num_frame,
                                 self.shape[2],
                                 self.shape[1],
                                 self.shape[0])[mode]
        if axis_name == 'X':
            data = np.max(np.max(
                np.transpose(data, (2, 0, 1)), axis=2), axis=1)
        elif axis_name == 'Y':
            data = np.max(np.max(
                np.transpose(data, (1, 0, 2)), axis=2), axis=1)
        elif axis_name == 'Z':
            data = np.max(np.max(data, axis=2), axis=1)
        else:
            raise 'Wrong axis name set'
        return data

    def median_along_axis(self, axis_name, mode=0):
        '''.. py:method:: median_along_axis(axis_name, potmode)

        Calculate median value of potential along 'axis'

        Parameters
        -----------

        axis_name: str
             'X', 'Y', or 'Z'
        mode: int
             select data by integer

        Returns
        -------
        numpy.ndarray
            median value along the axis
        '''
        axis_name = axis_name.capitalize()
        data = self.data.reshape(self.num_frame,
                                 self.shape[2],
                                 self.shape[1],
                                 self.shape[0])[mode]
        if axis_name == 'X':
            data = np.median(np.median(
                np.transpose(data, (2, 0, 1)), axis=2), axis=1)
        elif axis_name == 'Y':
            data = np.median(np.median(
                np.transpose(data, (1, 0, 2)), axis=2), axis=1)
        elif axis_name == 'Z':
            data = np.median(np.median(data, axis=2), axis=1)
        else:
            raise 'Wrong axis name set'
        return data
