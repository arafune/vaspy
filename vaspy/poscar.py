#! /usr/bin/env python
# -*- conding: utf-8 -*-
''' 
This mudule provides POSCAR class

translate from poscar.rb of 2014/2/26, master branch

The below is an example of POSCAR file::

 fcc (111) surface
   3.5200000
     0.707106780  0.000000000   0.000000000
    -0.353553390  0.612372400   0.000000000
     0.000000000  0.000000000  11.547000006
   Cu
     9
 Selective dynamics
 Direct
  0.000000000  0.000000000  0.000000000   F   F   F
  0.333333333  0.666666666  0.050000000   F   F   F
  0.666666666  0.333333333  0.100000000   F   F   F
  0.000000000  0.000000000  0.150000000   T   T   T
  0.333333333  0.666666666  0.200000000   T   T   T
  0.000000000  0.000000000  0.250000000   T   T   T
  0.333333333  0.666666666  0.300000000   T   T   T
  0.666666666  0.333333333  0.350000000   T   T   T
  0.000000000  0.000000000  0.400000000   T   T   T

  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
'''

from __future__ import division, print_function  # Version safety
import bz2
import itertools as it
import copy
import re
import os
import sys

try:
    from vaspy import tools
except ImportError:
    mypath = os.readlink(__file__) if os.path.islink(__file__) else __file__
    sys.path.append(os.path.dirname(os.path.abspath(mypath)))
    import tools
import numpy as np


class POSCAR_HEAD(object):
    ''' 
    One of the parent classes of POSCAR class

    Attributes
    ----------

    system_name: str
       system name
    scaling_factor: float
       scaling factor
    iontypes : list
       list of ion name
    ionnums : list
       list of number of ions. Corresponding to `iontypes`
'''
    def __init__(self):
        self.__cell_vecs = np.array([[0., 0., 0.],
                                     [0., 0., 0.],
                                     [0., 0., 0.]])
        self.system_name = ""
        self.scaling_factor = 0.
        self.iontypes = []
        self.ionnums = []
        self.__atom_identifer = []

    @property
    def cell_vecs(self):
        '''Return the matrix of the unit cell'''
        return self.__cell_vecs

    @cell_vecs.setter
    def cell_vecs(self, vec):
        '''Setter of cell matrix'''
        if three_by_three(vec):
            self.__cell_vecs = np.array(vec)
        else:
            raise TypeError

    @property
    def atom_identifer(self):
        '''Return list style of "atom_identifer" (e.g.  "#0:Ag1")'''
        # self.__atom_identifer = []
        # ii = 1
        # for elm, n in zip(self.iontypes, self.ionnums):
        #     self.__atom_identifer.extend(
        #         '#{0}:{1}{2}'.format(ii + m, elm, m + 1) for m in range(n))
        #     ii += n
        # return self.__atom_identifer
        self.__atom_identifer = []
        atomnames = []
        for elm, ionnums in zip(self.iontypes, self.ionnums):
            for j in range(1, ionnums + 1):
                elem_num = elm + str(j)
                if elem_num not in atomnames:
                    atomnames.append(elem_num)
                else:
                    while elem_num in atomnames:
                        j = j + 1
                        elem_num = elm + str(j)
                    else:
                        atomnames.append(elem_num)
        self.__atom_identifer = [
            "#" + str(s) + ":" + a for s, a in
            zip(range(0, len(atomnames)), atomnames)]
        return self.__atom_identifer

    @atom_identifer.setter
    def atom_identifer(self, value):
        self.__atom_identifer = value


class POSCAR_POS(object):
    '''
    POSCAR_DOS Class
'''
    def __init__(self):
        self.coordinate_type = ""
        self.positions = []
        self.coordinate_changeflags = []
        self.selective = False

    def is_cartesian(self):
        '''

        Return True if Cartesian coordinate is set

        Returns
        --------

        boolean
            True if coordinate is cartesian
        '''
        return bool(re.search(r'^[ck]', self.coordinate_type, re.I))

    def is_direct(self):
        ''' 
        Return True if DIRECT coordinate is set

        Returns
        --------

        Boolean
            True if coordinate is direct (not cartesian)
        '''
        return not self.is_cartesian()


class POSCAR(POSCAR_HEAD, POSCAR_POS):
    '''

    Class for POSCAR (CONTCAR) format

    This script does *NOT* support for constructing POSCAR
    from scratch. (Use ASE for this purpose.)

    It provides a way to slightly modify POSCAR or
    CONTCAR, which has already works well.

    Attributes
    ------------
    system_name, scaling_factor, cell_vecs
    '''

    def __init__(self, arg=None):
        '''
        Parameters
        -----------

        arg: str
            POSCAR file name, or list of POSCAR text.
        '''
        super().__init__()
        POSCAR_POS.__init__(self)
        if isinstance(arg, str):
            if os.path.splitext(arg)[1] == '.bz2':
                try:
                    thefile = bz2.open(arg, mode='rt')
                except AttributeError:
                    thefile = bz2.BZ2File(arg, mode='r')
            else:
                thefile = open(arg)
            poscar = thefile.readlines()
            self.load_from_array(poscar)
        if isinstance(arg, (list, tuple)):
            self.load_from_array(arg)

    def load_from_array(self, poscar):
        '''.. :py:method:: load_from_array(poscar)

        POSCAR parser

        Parameters
        ----------

        poscar: str, list, tuple
             POSCAR data
        '''
        poscar = iter(map(str.rstrip, poscar))  # Version safety
        self.system_name = next(poscar)
        self.scaling_factor = float(next(poscar))
        self.cell_vecs[0] = [float(x) for x in next(poscar).split()]
        self.cell_vecs[1] = [float(x) for x in next(poscar).split()]
        self.cell_vecs[2] = [float(x) for x in next(poscar).split()]
        self.iontypes = next(poscar).split()
        # parse POSCAR evenif the element names are not set.
        # At present, the String representation number
        #   are used for the  dummy name.
        if self.iontypes[0].isdigit():
            self.ionnums = [int(i) for i in self.iontypes]
        else:
            self.ionnums = [int(x) for x in next(poscar).split()]
        line7 = next(poscar)
        if re.search(r'^[\s]*Selective\b', line7, re.I):
            self.selective = True
            self.coordinate_type = next(poscar)
        else:
            self.selective = False
            self.coordinate_type = line7

        for line, elem in zip(poscar, self.atom_identifer):
            # if not elem: break
            tmp = line.split()
            self.positions.append(np.float_(np.array(tmp[:3])))
            if self.selective:
                self.coordinate_changeflags.append(' '.join(tmp[3:]))

    def __iter__(self):
        for each in self.positions:
            yield each

    def sort(self, from_site=0, to_site=None, axis='z'):
        '''
        Sort positions attribute by coordinate

        Parameters
        -----------

        from_site: int
            first index # for sort

        to_site: int
            last index # for sort

        axis: str
            Axis used for sort (Default Z)


        Note
        -----

        The first site # is "0". It's the pythonic way.
        The element difference is **not** taken into account.
        '''
        if to_site is None:
            to_site = sum(self.ionnums)
        if axis == 'x' or axis == 'X' or axis == 0:
            axis = 0
        elif axis == 'y' or axis == 'Y' or axis == 1:
            axis = 1
        elif axis == 'z' or axis == 'Z' or axis == 2:
            axis = 2
        self.positions = self.positions[0:from_site] + sorted(
            self.positions[from_site:to_site],
            key=lambda sortaxis: sortaxis[axis]) + self.positions[to_site:]

    def supercell(self, n_x, n_y, n_z):
        ''' 
        Return the :math:`(n_x \\times n_y \\times n_z)` supercell

        Parameters
        -----------
        n_x: int
            repeat number along x axis
        n_y: int
            repeat number along y axis
        n_z: int
            repeat number along z axis


        Returns
        ---------

        POSCAR
            POSCAR object of the supercell
        '''
        if not isinstance(n_x, int) \
           or not isinstance(n_y, int) \
           or not isinstance(n_z, int):
            raise ValueError("arguments must be positive integer")
        if n_x <= 0 or n_y <= 0 or n_z <= 0:
            raise ValueError("arguments must be positive integer")
        #
        sposcar = copy.deepcopy(self)
        original_is_cartesian = sposcar.is_cartesian()
        if original_is_cartesian:
            sposcar.to_direct()
        sposcar.repack_in_cell()
        sposcar.cell_vecs[0] = sposcar.cell_vecs[0] * n_x
        sposcar.cell_vecs[1] = sposcar.cell_vecs[1] * n_y
        sposcar.cell_vecs[2] = sposcar.cell_vecs[2] * n_z
        sposcar.ionnums = [i * n_x * n_y * n_z for i in sposcar.ionnums]
        spositions = sposcar.positions
        sposcar.positions = []
        spositions = [np.array([x[0] / n_x, x[1] / n_y, x[2] / n_z])
                      for x in spositions]
        for spos in spositions:
            for i_z in range(0, n_z):
                for i_y in range(0, n_y):
                    for i_x in range(0, n_x):
                        sposcar.positions.append(np.array(
                            [spos[0] + i_x / n_x,
                             spos[1] + i_y / n_y,
                             spos[2] + i_z / n_z]))
        sposcar.coordinate_changeflags = []
        for flags in self.coordinate_changeflags:
            for _ in range(n_x * n_y * n_z):
                sposcar.coordinate_changeflags.append(flags)
        sposcar.atom_identifer = []
        if original_is_cartesian:
            sposcar.to_cartesian()
        return sposcar

    # class method? or independent function?
    def nearest(self, array, point):
        ''' 
        Parameters
        -----------

        array: list

        point: numpy.ndarray

        Returns
        --------

        numpy.ndarray
'''
        return min(array, key=lambda pos: np.linalg.norm(pos - point))

    # class method? or independent function?
    def make27candidate(self, position):
        '''

        Return 27 vectors set correspond the neiboring

        Parameters
        -----------

        position: numpy.ndarray, list
            atom position defined in the coordinated by
                         cell_vecs ( scaling facter is not accounted).

        Returns
        --------

        list

        '''
        position = _vectorize(position)
        candidates27 = []
        if self.is_cartesian():
            for i, j, k in it.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1]):
                candidates27.append(i * self.cell_vecs[0] +
                                    j * self.cell_vecs[1] +
                                    k * self.cell_vecs[2] +
                                    position)
        else:
            for i, j, k in it.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1]):
                candidates27.append(i * np.array([1., 0., 0.]) +
                                    j * np.array([0., 1., 0.]) +
                                    k * np.array([0., 0., 1.]) +
                                    position)
        return candidates27

    def rotate_atom(self, site, axis_name, theta_deg, center):
        '''

        Rotate the atom

        Parameters
        ----------

        site: int
            site # for rotation (The first atom is "0".).
        axis_name: str
            "X", "x", "Y", "y", "Z", or "z". Rotation axis.
        theta_deg: float
            Rotation angle (Degrees).
        center: numpy.ndarray, list, tuple
            center position for rotation.


        :Todo: check the center in the Braves lattice.
               take into account the periodic boundary.
        '''
        center = _vectorize(center)
        if len(center) != 3:
            raise ValueError
        if not point_in_box(center / self.scaling_factor, self.cell_vecs):
            raise ValueError('the center must be in the Braves lattice')
        if not isinstance(site, int):
            raise ValueError('argument error in rotate_atom method')
        if not self.is_cartesian():
            self.to_cartesian()
        position = self.positions[site]
        position -= center / self.scaling_factor
        position = globals()["rotate_" +
                             axis_name.lower()](theta_deg).dot(position)
        position += center / self.scaling_factor
        self.positions[site] = position

    def rotate_atoms(self, site_list, axis_name, theta_deg, center):
        '''
        Rotate atoms

        Parameters
        ----------

        site_list:
            list array of site for rotation.
        axis_name:
            "X", "x", "Y", "y", "Z",or "z".  Rotation axis.
        theta_deg: float
             Rotation angle (Degrees).
        center: numpy.ndarray, list, tuple
             Position of rotation center
        '''
        for site in site_list:
            self.rotate_atom(site, axis_name, theta_deg, center)

    def rotate_cell(self, theta_deg, axis_name='Z'):
        ''' 
        Rotate unit-cell (rotation angle is set by degree.)

        Parameters
        ----------

        theta_deg: float
            rotation angle
        axis_name: str
            axis name for rotation (x, y, or z)
        '''
        original_is_cartesian = False
        if self.is_cartesian():
            original_is_cartesian = True
            self.to_direct()
        axis_name = axis_name.capitalize()

        if axis_name == 'X':
            self.cell_vecs = np.dot(rotate_x(theta_deg), self.cell_vecs.T).T
        elif axis_name == 'Y':
            self.cell_vecs = np.dot(rotate_y(theta_deg), self.cell_vecs.T).T
        elif axis_name == 'Z':
            self.cell_vecs = np.dot(rotate_z(theta_deg), self.cell_vecs.T).T
        if original_is_cartesian:
            self.to_cartesian()

    def repack_in_cell(self):
        '''
        Repack all atoms in the unit cell

        No negative values in DIRECT coordinate.
        '''
        original_is_cartesian = False
        if self.is_cartesian():
            original_is_cartesian = True
            self.to_direct()
        for pos in self.positions:
            for i in (0, 1, 2):
                while pos[i] < 0.0 or pos[i] > 1.0:
                    if pos[i] < 0.0:
                        pos[i] += 1.0
                    elif pos[i] > 1.0:
                        pos[i] -= 1.0
        if original_is_cartesian:
            self.to_cartesian()

    def __add__(self, other):
        '''
        Parameters
        ----------

        other: POSCAR

        Returns
        --------
        POSCAR

        Todo
        -----

        Check the lattice vectors, coordinate_type and so on.
        '''
        if not isinstance(other, POSCAR):
            return NotImplemented
        dest_poscar = copy.deepcopy(self)
        if dest_poscar.scaling_factor != other.scaling_factor:
            raise ValueError('scaling factor is different.')
        if np.linalg.norm(dest_poscar.cell_vecs - other.cell_vecs) != 0:
            raise ValueError('lattice vectors (cell matrix) are different.')
        dest_poscar.iontypes.extend(other.iontypes)
        dest_poscar.ionnums.extend(other.ionnums)
        dest_poscar.positions.extend(other.positions)
        dest_poscar.coordinate_changeflags.extend(other.coordinate_changeflags)
        return dest_poscar

    def merge(self, other):
        '''
        lazy __add__: Return POSCAR generated from two POSCARs

        Even if the cell vectors and scaling factors are different,
        the 'merged' POSCAR is created.
        Use the cell vectors of the first POSCAR.

        Parameters
        ------------

        other: POSCAR
            POSCAR object

        Returns
        --------

        POSCAR
            added poscar object
        '''
        if not isinstance(other, POSCAR):
            return NotImplemented
        dest_poscar = copy.deepcopy(self)
        original_is_direct = False
        if dest_poscar.is_direct():
            original_is_direct = True
            dest_poscar.to_cartesian()
        other_poscar = copy.deepcopy(other)
        original_scaling_factor = dest_poscar.scaling_factor
        other_poscar.tune_scaling_factor(original_scaling_factor)
        other_poscar.to_cartesian()
        dest_poscar.iontypes.extend(other.iontypes)
        dest_poscar.ionnums.extend(other.ionnums)
        dest_poscar.positions.extend(other.positions)
        dest_poscar.coordinate_changeflags.extend(other.coordinate_changeflags)
        if original_is_direct:
            dest_poscar.to_direct()
        return dest_poscar

    def to_list(self):
        '''
        Return POSCAR object by list-style

        Returns
        ---------

        list
        a list representation of POSCAR.
        '''
        out_list = []
        out_list.append(self.system_name)
        out_list.append(self.scaling_factor)
        out_list.append(self.cell_vecs[0])
        out_list.append(self.cell_vecs[1])
        out_list.append(self.cell_vecs[2])
        if not self.iontypes[0].isdigit():
            out_list.append(self.iontypes)
        out_list.append(self.ionnums)
        if self.selective:
            out_list.append("Selective Dynamics")
        out_list.append(self.coordinate_type)
        out_list.append(self.positions)
        out_list.append(self.coordinate_changeflags)
        out_list.append(self.atom_identifer)
        return out_list

    def __str__(self):
        '''
        Returns
        --------

        str
            a string representation of POSCAR
        '''
        tmp = []
        tmp.append(self.system_name)
        tmp.append(str(self.scaling_factor))
        tmp.append(''.join('   {0:20.17f}'.format(i) for i in
                           self.cell_vecs[0]))
        tmp.append(''.join('   {0:20.17f}'.format(i) for i in
                           self.cell_vecs[1]))
        tmp.append(''.join('   {0:20.17f}'.format(i) for i in
                           self.cell_vecs[2]))
        if not self.iontypes[0].isdigit():
            tmp.append(' ' + ' '.join(self.iontypes))
        tmp.append(' ' + ' '.join(str(i) for i in self.ionnums))
        if self.selective:
            tmp.append('Selective Dynamics')
        tmp.append(self.coordinate_type)
        for pos, t_or_f, atom in tools.ZIPLONG(self.positions,
                                               self.coordinate_changeflags,
                                               self.atom_identifer,
                                               fillvalue=''):
            tmp.append(' '.join('  {0:20.17f}'.format(i) for i in pos) +
                       ' ' + t_or_f +
                       ' ' + atom)
        return '\n'.join(tmp) + '\n'

    def str_short(self):
        '''
        Returns
        ---------

        str
            a string representation of POSCAR, with short (8) digit format.
            used in CHGCAR
        '''
        tmp = []
        tmp.append(self.system_name)
        tmp.append('  {0:.14f}'.format(self.scaling_factor))
        for k in range(3):
            tmp.append(''.join('{0:12.6f}'.format(i) for i in
                               self.cell_vecs[k]))
        if not self.iontypes[0].isdigit():
            tmp.append(' ' +
                       "".join(['{0:>5}'.format(i) for i in self.iontypes]))
        tmp.append(''.join(['{0:>6}'.format(i) for i in self.ionnums]))
        tmp.append(self.coordinate_type)
        for pos in self.positions:
            tmp.append(''.join('{0:10.6f}'.format(i) for i in pos))
        return '\n'.join(tmp) + '\n'

    def tune_scaling_factor(self, new_scaling_factor=1.0):
        ''' 
        Change scaling factor to new value

        Parameters
        ------------

        new_scaling_factor: float

        Note
        -----

        **The Braves lattice are corrected (to be equal size)**

        Warning
        --------

         If you change the cell size, change scaling_factor attribute directly
        '''
        old = self.scaling_factor
        self.cell_vecs *= (old / new_scaling_factor)
        self.scaling_factor = new_scaling_factor
        if self.is_cartesian():
            self.positions = [i * old / new_scaling_factor
                              for i in self.positions]

    def to_cartesian(self):
        ''' 
        Change the coordinate to cartesian from direct
        '''
        if self.is_direct():
            self.coordinate_type = "Cartesian"
            mat = self.cell_vecs.transpose()
            self.positions = [mat.dot(v) for v in self.positions]

    def to_direct(self):
        ''' 
        Change the coordinate to direct from cartesian.
        '''
        if self.is_cartesian():
            self.coordinate_type = "Direct"
            mat = np.linalg.inv(np.transpose(self.cell_vecs))
            self.positions = [mat.dot(v) for v in self.positions]

    def guess_molecule(self, site_list, center=None):
        '''
        Arrange atom position to form a molecule.

        This method is effective to rotate a molecule.

        Parameters
        -----------

        site_list: list
            list of site number (the number begins with #1)
        center: list
            center position of "molecule" (Optional).

        Returns
        ---------

        numpy.ndarray
            Array of Vector that represents "molecule".

        Note
        -----
        When the optional argument, center, is set, the atoms are
        re-arranged as to minimize the distance from this center.
        If not, atoms are re-arranged to minimize the total bonding
        length.  As the algorithm for mimizing the total length is
        not exhaustive, the resultant atom arrangement may different
        from what you expect, in spite of time-waste.  The center
        option is highly recommended to form a molecule.
        '''
        molecule = [self.positions[j] for j in site_list]
        newposes = []
        for index, site in enumerate(site_list):
            target_atom = self.positions[site]
            atoms27 = self.make27candidate(target_atom)

            def func(pos, center):
                molecule[index] = pos
                if center is not None:  # bool([np.ndarray]) => Error
                    center = _vectorize(center)
                    return np.linalg.norm(pos - center)
                else:
                    # fixme!! when the highest symmetry point
                    # can be detemined from the position list,
                    # guess_molecule method does not require
                    # the "center" option.
                    # (molecule.
                    #     product(molecule)).inject(0.0) do | s, vectors |
                    # s+ (vectors[0]-vectors[1]).magnitude
                    s = 0.0
                    for vectors in it.product(molecule, molecule):
                        s += np.linalg.norm(vectors[0] - vectors[1])
                    return s
            newpos = min(atoms27, key=(lambda x: func(x, center)))
            newposes.append(newpos)
        for site, pos in zip(site_list, newposes):
            self.positions[site] = pos

    def translate(self, vector, atomlist):
        ''' 
        Translate the selected atom(s) by vector

        Parameters
        ----------

        vector: list, tuple, numpy.array
             translational vector (in Cartesian frame)
        atomlist: list
             list of the atom to be moved


        Note
        ------

        the first atom is "0", to follow the pythonic way.

        Returns
        --------

        numpy.ndarray
              position
        '''
        if self.is_cartesian():
            vector = _vectorize(vector)
            for i in atomlist:
                self.positions[i] = (self.positions[i] +
                                     vector / self.scaling_factor)
        else:
            vector = _vectorize(vector)
            self.to_cartesian()
            for i in atomlist:
                self.positions[i] = (self.positions[i] +
                                     vector / self.scaling_factor)
            self.to_direct()
        return self.positions

    @property
    def axes_lengthes(self):
        ''' 
        Return cell axis lengthes

        Returns
        -------

        tuple
            cell axis length of x, y, and z
        '''
        cell_x = np.linalg.norm(self.cell_vecs[0] * self.scaling_factor)
        cell_y = np.linalg.norm(self.cell_vecs[1] * self.scaling_factor)
        cell_z = np.linalg.norm(self.cell_vecs[2] * self.scaling_factor)
        return (cell_x, cell_y, cell_z)

    def translate_all(self, vector):
        ''' 
        Translate **all** atoms by vector

        Parameters
        -----------

        vector: list, numpy.array
             translational vector
        '''
        atomrange = list(range(sum(self.ionnums)))
        self.translate(vector, atomrange)

    def save(self, filename):
        ''' 
        Save POSCAR contents to the file named "filename"

        Parameters
        ----------

        filename: str
             File name for save
        '''
        try:  # Version safety
            file = open(filename, mode='w', newline='\n')
        except TypeError:
            file = open(filename, mode='wb')
        with file:
            file.write(str(self))


def point_in_box(point, cell_vecs):
    ''' 
    Return True if point is located in the box

    Parameters
    -----------

    point: numpy.ndarray, numpy.matrix, list, tuple
        vector representing the "point"
    cell_vecs: numpy.ndarray, numpy.matrix, list, tuple
        vectors defining the "box"

    Returns
    ---------

    boolean
    '''
    if three_by_three(cell_vecs):
        point = np.array(point).flatten()
        cell_vecs = np.array(cell_vecs)
        result = np.dot(np.linalg.inv(cell_vecs.T), point)
        return all((0 <= float(q) <= 1) for q in result)
    else:
        raise TypeError


def rotate_x(theta_deg):
    '''
    Rotation matrix around X-axis

    Parameters
    ----------

    theta_deg: float
        angle of rotation (Degrees)

    Returns
    -------

    numpy.ndarray
        rotation matrix

    Example
    ---------

    >>> rotate_x(60)
    array([[ 1.       ,  0.       ,  0.       ],
           [ 0.       ,  0.5      , -0.8660254],
           [ 0.       ,  0.8660254,  0.5      ]])
    '''
    degree = np.pi / 180.0
    return np.array(
        [[1.0, 0.0, 0.0],
         [0.0, np.cos(theta_deg * degree), -np.sin(theta_deg * degree)],
         [0.0, np.sin(theta_deg * degree), np.cos(theta_deg * degree)]])


def rotate_y(theta_deg):
    ''' 
    Rotation matrix around Y-axis

    Example
    --------

    >>> rotate_y(60)
    array([[ 0.5      ,  0.       ,  0.8660254],
           [ 0.       ,  1.       ,  0.       ],
           [-0.8660254,  0.       ,  0.5      ]])
    '''
    degree = np.pi / 180.0
    return np.array(
        [[np.cos(theta_deg * degree), 0.0, np.sin(theta_deg * degree)],
         [0.0, 1.0, 0.0],
         [-np.sin(theta_deg * degree), 0.0, np.cos(theta_deg * degree)]])


def rotate_z(theta_deg):
    ''' 
    Rotation matrix around Z-axis

    Example
    --------

    >>> rotate_z(60)
    array([[ 0.5      , -0.8660254,  0.       ],
           [ 0.8660254,  0.5      ,  0.       ],
           [ 0.       ,  0.       ,  1.       ]])
    '''
    degree = np.pi / 180.0
    return np.array(
        [[np.cos(theta_deg * degree), -np.sin(theta_deg * degree), 0.0],
         [np.sin(theta_deg * degree), np.cos(theta_deg * degree), 0.0],
         [0.0, 0.0, 1.0]])


def three_by_three(vec):
    ''' 
    Return True if vec can be converted into the 3x3 matrix

    Parameters
    ----------

    vec: numpy.ndarray, numpy.matrix, list, tuple
    list like object

    Returns
    --------

    boolean
    '''
    if not isinstance(vec, (np.ndarray, np.matrix, list, tuple)):
        return False
    if len(vec) != 3:
        return False
    if [3, 3, 3] == [len(i) for i in vec]:
        return True
    else:
        return False


def _vectorize(vector):
    if not isinstance(vector, (np.ndarray, np.matrix, list, tuple)):
        raise TypeError('Cannot convert into vector.')
    return np.array(vector).flatten()

# --------------------------

if __name__ == "__main__":
    import doctest
    doctest.testmod()
