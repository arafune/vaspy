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


class VASPGrid(poscar.POSCAR):
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

    mesh3d, additional, mesh_x, mesh_y, mesh_z
    '''
    def __init__(self, filename=None):
        super(VASPGrid, self).__init__()
        self.mesh_x, self.mesh_y, self.mesh_z = 0, 0, 0
        self.mesh3d = []
        self.additional = []
        if filename:
            self.load_from_file(filename)

    def load_from_file(self, filename):
        '''.. py:method:: load_from_file(Parse)

        filename file to construct the object

        Parameters
        ---------------

        filename: str
            file name
        '''
        section = 'poscar'
        separator = None
        tmp = []
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
                        self.load_from_array(tmp)
                        section = 'define_separator'
                    else:
                        tmp.append(line)
                elif section == 'define_separator':
                    separator = line if separator is None else separator
                    if self.mesh_x == self.mesh_y == self.mesh_z == 0:
                        self.mesh_x, self.mesh_y, self.mesh_z = \
                                    [int(string) for string in line.split()]
                    mesh_n = self.mesh_x * self.mesh_y * self.mesh_z
                    self.mesh3d.extend([next(thefile).rstrip().replace(
                        '***********', 'Nan').split()])
                    if mesh_n % len(self.mesh3d[0]) == 0:
                        lines_for_mesh = mesh_n // len(self.mesh3d[0])
                    else:
                        lines_for_mesh = mesh_n // len(self.mesh3d[0]) + 1
                    self.mesh3d.extend(
                        [next(thefile).rstrip().replace('***********',
                                                        'Nan').split()
                         for i in range(lines_for_mesh-1)])
                    section = 'grid'
                elif section == 'aug':
                    if separator in line:
                        self.mesh3d.extend([
                            next(thefile).rstrip().replace('***********',
                                                           'Nan').split()
                            for i in range(lines_for_mesh)])
                        section = 'grid'
                    elif "augmentation occupancies " in line:
                        pass  # Used for CHGCAR, not LOCPOT. not implementd
                    else:
                        pass  # Used for CHGCAR, not LOCPOT. not implementd
                elif section == 'grid':
                    if "augmentation occupancies " in line:
                        section = 'aug'
                    elif separator in line:
                        self.mesh3d.extend([
                            next(thefile).rstrip().replace('***********',
                                                           'Nan').split()
                            for i in range(lines_for_mesh)])
                    else:
                        # for unused data stored in LOCPOT
                        self.additional.extend(line.split())
            self.mesh3d = np.array([elem for sublist in self.mesh3d
                                    for elem in sublist], dtype=np.float64)
            self.n_frame = divmod(len(self.mesh3d),
                                  self.mesh_x * self.mesh_y * self.mesh_z)[0]
#            self.mesh3d = self.mesh3d.reshape(self.mesh3d.shape[0] *
#                                              self.mesh3d.shape[1])

    def get_mesh(self):
        '''.. py:method:: get_mesh()

        Return mesh_size

        Returns
        ------------
        tuple
           (mesh_x, mesh_y, mesh_z)
        '''
        return self.mesh_x, self.mesh_y, self.mesh_z

    def __add__(self, other):
        '''.. py:method:: __add__(other)

        x.__add__(y) <=> x + y

        Parameters
        ---------------

        other: VASPGrid
            Addtion VaspGrid object

        Returns
        -----------

        VASPGrid
            Rusultant by summing two grid value

        Notes
        ----------

        The Grid data is summed by this method, of course. And the
        atom data (POSCAR part) is also summed.
        '''
        if not isinstance(other, VASPGrid):
            return NotImplemented
        add_vaspgrid = super(VASPGrid, self).__add__(other)
        try:
            add_vaspgrid.mesh3d = self.mesh3d + other.mesh3d
        except ValueError:
            raise RuntimeError('The mesh sizes are different each other')
        return add_vaspgrid

    def __sub__(self, other):
        '''.. py:method:: __sub__(other)

        x.__sub__(y) <=> x - y

        Parameters
        ---------------

        other: VASPGrid
            difference VASPGrid object

        Returns
        ----------

        VASPGrid
            Resultant by difference between two objects.

        Notes
        --------

        The resultant grid data is the difference between two objects,
        of course. On the other hand, the atom position information
        unchange by this method.  Use the 'minuend' object.  The atom
        information in subrtrahend object is totally ignored.
        '''

        if not isinstance(other, VASPGrid):
            return NotImplemented
        diff_vaspgrid = copy.deepcopy(self)
        try:
            diff_vaspgrid.mesh3d = self.mesh3d - other.mesh3d
        except ValueError:
            raise RuntimeError('The mesh sizes are different each other')
        return diff_vaspgrid

    def merge(self, other):
        '''.. py:method:: merge(other)

        x.merge(y) -> '3D mesh data of x' + '3D mesh data of y'
        '''
        if not isinstance(other, VASPGrid):
            return NotImplemented
        merged_vaspgrid = copy.deepcopy(self)
        try:
            merged_vaspgrid.mesh3d = self.mesh3d + other.mesh3d
        except ValueError:
            raise RuntimeError('The mesh sizes are different each other')
        return merged_vaspgrid

    def __str__(self):
        '''.. py:method:: __str__()

        x.__str__() <=> str(x)

        Returns
        -------

        str
            a string representation of VASPGrid object
        '''
        mesharray = self.mesh3d.reshape(self.mesh3d.size // (self.mesh_x *
                                                             self.mesh_y *
                                                             self.mesh_z),
                                        self.mesh_x *
                                        self.mesh_y *
                                        self.mesh_z)
        outputstr = self.str_short()
        for tmp in mesharray:
            output = []
            outputstr += '\n  {0}  {1}  {2}\n'.format(self.mesh_x,
                                                      self.mesh_y,
                                                      self.mesh_z)
            for array in tools.each_slice(tmp, 5):
                output.append(''.join('  {0:18.11E}'.format(i)
                                      for i in array if i is not None))
            outputstr += '\n'.join(output)
        return outputstr + '\n'

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

    def average_along_axis(self, axis_name, mode=0):
        '''.. py:method:: average_along_axis(axis_name, mode)

        Calculate average value of potential along 'axis'

        Parameters
        ----------

        axis_name: str
             'X', 'Y', or 'Z'
        mode: int
             select meshdata by integer

        Returns
        -------
        numpy.ndarray
            average value along the axis
        '''
        axis_name = axis_name.capitalize()
        meshdata = self.mesh3d.reshape(self.n_frame,
                                       self.mesh_z,
                                       self.mesh_y,
                                       self.mesh_x)[mode]
        if axis_name == 'X':
            meshdata = np.average(np.average(
                np.transpose(meshdata, (2, 0, 1)), axis=2), axis=1)
        elif axis_name == 'Y':
            meshdata = np.average(np.average(
                np.transpose(meshdata, (1, 0, 2)), axis=2), axis=1)
        elif axis_name == 'Z':
            meshdata = np.average(np.average(meshdata, axis=2), axis=1)
        else:
            raise 'Wrong axis name set'
        return meshdata

    def min_along_axis(self, axis_name, mode=0):
        '''.. py:method:: min_along_axis(axis_name, potmode)

        Calculate minimum value of potential along 'axis'

        Parameters
        -----------

        axis_name: str
             'X', 'Y', or 'Z'
        mode: int
             select meshdata by integer

        Returns
        -------
        numpy.ndarray
            minimum value along the axis
        '''
        axis_name = axis_name.capitalize()
        meshdata = self.mesh3d.reshape(self.n_frame,
                                       self.mesh_z,
                                       self.mesh_y,
                                       self.mesh_x)[mode]
        if axis_name == 'X':
            meshdata = np.min(np.min(
                np.transpose(meshdata, (2, 0, 1)), axis=2), axis=1)
        elif axis_name == 'Y':
            meshdata = np.min(np.min(
                np.transpose(meshdata, (1, 0, 2)), axis=2), axis=1)
        elif axis_name == 'Z':
            meshdata = np.min(np.min(meshdata, axis=2), axis=1)
        else:
            raise 'Wrong axis name set'
        return meshdata

    def max_along_axis(self, axis_name, mode=0):
        '''.. py:method:: man_along_axis(axis_name, potmode)

        Calculate maximum value of potential along 'axis'

        Parameters
        -----------

        axis_name: str
             'X', 'Y', or 'Z'
        mode: int
             select meshdata by integer

        Returns
        -------
        numpy.ndarray
            maximum value along the axis
        '''
        axis_name = axis_name.capitalize()
        meshdata = self.mesh3d.reshape(self.n_frame,
                                       self.mesh_z,
                                       self.mesh_y,
                                       self.mesh_x)[mode]
        if axis_name == 'X':
            meshdata = np.max(np.max(
                np.transpose(meshdata, (2, 0, 1)), axis=2), axis=1)
        elif axis_name == 'Y':
            meshdata = np.max(np.max(
                np.transpose(meshdata, (1, 0, 2)), axis=2), axis=1)
        elif axis_name == 'Z':
            meshdata = np.max(np.max(meshdata, axis=2), axis=1)
        else:
            raise 'Wrong axis name set'
        return meshdata

    def median_along_axis(self, axis_name, mode=0):
        '''.. py:method:: median_along_axis(axis_name, potmode)

        Calculate median value of potential along 'axis'

        Parameters
        -----------

        axis_name: str
             'X', 'Y', or 'Z'
        mode: int
             select meshdata by integer

        Returns
        -------
        numpy.ndarray
            median value along the axis
        '''
        axis_name = axis_name.capitalize()
        meshdata = self.mesh3d.reshape(self.n_frame,
                                       self.mesh_z,
                                       self.mesh_y,
                                       self.mesh_x)[mode]
        if axis_name == 'X':
            meshdata = np.median(np.median(
                np.transpose(meshdata, (2, 0, 1)), axis=2), axis=1)
        elif axis_name == 'Y':
            meshdata = np.median(np.median(
                np.transpose(meshdata, (1, 0, 2)), axis=2), axis=1)
        elif axis_name == 'Z':
            meshdata = np.median(np.median(meshdata, axis=2), axis=1)
        else:
            raise 'Wrong axis name set'
        return meshdata


class FFTGRID(object):
    '''.. py:class:: FFTGRID(meshsize, meshdata)
    Class for NG(X,Y,Z)F in VASP

    This class is used chg_array in CHGCAR, Potential in LOCPOT,
    electron localization function (ELF) in ELFCAR

    Parameters
    ----------

    meshsize: tuple
        mesh_x, mesh_y, mesh_z
    meshdata: list or numpy.array
        1D-list or 1D-numpy array.  The length of grid is mesh_x *
    mesh_y * mesh_z
    '''
    def __init__(self, meshsize, meshdata):
        if len(meshdata) == meshsize(0) * meshsize(1) * meshsize(2):
            self.mesh_x = meshsize(0)
            self.mesy_y = meshsize(1)
            self.mesh_z = meshsize(2)
            self.meshdata = np.array(meshdata)
        else:
            raise RuntimeError

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
        ''' .. py:method:: __str__()

        Returns
        -------

        str
            a string representation of occupancy.
        '''
        pass
