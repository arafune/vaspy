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
                                    [int(str) for str in line.split()]
                    section = 'grid'
                elif section == 'aug':
                    if separator in line:
                        section = 'grid'
                    elif "augmentation occupancies " in line:
                        pass  # Used for CHGCAR, not LOCPOT. not implementd
                    else:
                        pass  # Used for CHGCAR, not LOCPOT. not implementd
                elif section == 'grid':
                    if "augmentation occupancies " in line:
                        section = 'aug'
                    elif separator in line:
                        pass  # Used for CHGCAR, not LOCPOT. not implementd
                    else:
                        self.mesh3d.extend(line.split())
            self.mesh3d = np.array(self.mesh3d, dtype=np.float64)
            # check the mesh size to detemine whether mesh is
            # potential or charge this algorithm is, of course, very
            # ad-hoc.
            if not len(self.mesh3d) % (self.mesh_x * self.mesh_y *
                                   self.mesh_z) == 0:
                self.mesh3d = self.mesh3d.reshape(4, self.mesh_x *
                                                  self.mesh_y *
                                                  self.mesh_z +
                                                  self.ionnums)
                self.additional = self.mesh3d[:, -self.ionnums:]
                self.mesh3d = np.delete(self.mesh3d, np.s_[-self.ionnums:])

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

    def __sub__self(self, other):
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

    def __str__(self):
        '''.. py:method:: __str__()

        x.__str__() <=> str(x)

        Returns
        -------

        str
            a string representation of VASPGrid object
        '''
        outputstring = ''
        mesharray = self.mesh3d.reshape(len(self.mesh3d) // (self.mesh_x *
                                                             self.mesh_y *
                                                             self.mesh_z),
                                        self.mesh_x *
                                        self.mesh_y *
                                        self.mesh_z)
        for tmp in mesharray:
            output = []
            outputstring += '\n  {0}  {1}  {2}\n'.format(self.mesh_x,
                                                         self.mesh_y,
                                                         self.mesh_z)
            for array in tools.each_slice(tmp, 5):
                output.append(''.join('  {0:18.11E}'.format(i)
                                      for i in array if i is not None))
            outputstring += '\n'.join(output)
        return super(VASPGrid, self).__str__() + outputstring + '\n'

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
