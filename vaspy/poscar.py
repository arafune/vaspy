#! /usr/bin/env python
# -*- conding: utf-8 -*-
'''.. py:module:: poscar

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


class POSCAR(object):
    '''.. py:class:: POSCAR(file or array)

    Class for POSCAR (CONTCAR) format

    This script does *NOT* support for constructing POSCAR
    from scratch. (Use ASE for this purpose.)

    It provides a way to slightly modify POSCAR or
    CONTCAR, which has already works well.

    :attribute: system_name, scaling_factor, cell_vecs
    :version: 2.1
    '''

    def __init__(self, arg=None):
        '''
        :param str arg: POSCAR file name, or list of POSCAR text.
        :type arg: str
        '''
        self.system_name = ""
        self.scaling_factor = 0.
        self.__cell_vecs = np.array([[0., 0., 0.],
                                     [0., 0., 0.],
                                     [0., 0., 0.]])
        self.iontype = []
        self.ionnums = []
        self.coordinate_type = ""
        self.position = []
        self.coordinate_changeflags = []
        self.__atom_identifer = []
        self.__selective = False
        if isinstance(arg, str):
            poscar = open(arg).readlines()
            self.load_from_array(poscar)
        if isinstance(arg, (list, tuple)):
            self.load_from_array(arg)

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

    def load_from_array(self, poscar):
        '''POSCAR parser'''
        poscar = iter(map(str.rstrip, poscar))  # Version safety
        self.system_name = next(poscar)
        self.scaling_factor = float(next(poscar))
        self.cell_vecs[0] = [float(x) for x in next(poscar).split()]
        self.cell_vecs[1] = [float(x) for x in next(poscar).split()]
        self.cell_vecs[2] = [float(x) for x in next(poscar).split()]
        self.iontype = next(poscar).split()
        # parse POSCAR evenif the element names are not set.
        # At present, the String representation number
        #   are used for the  dummy name.
        if self.iontype[0].isdigit():
            self.ionnums = [int(i) for i in self.iontype]
        else:
            self.ionnums = [int(x) for x in next(poscar).split()]
        self.__atom_identifer = []
        atomnames = []
        for elm, ionnums in zip(self.iontype, self.ionnums):
            for j in range(1, ionnums + 1):
                tmp = elm + str(j)
                if tmp not in atomnames:
                    atomnames.append(tmp)
                else:
                    while tmp in atomnames:
                        j = j + 1
                        tmp = elm + str(j)
                    else:
                        atomnames.append(tmp)
        self.__atom_identifer = [
            "#" + str(s) + ":" + a for s, a in
            zip(range(1, len(atomnames) + 1), atomnames)]
        line7 = next(poscar)
        if re.search(r'^[\s]*Selective\b', line7, re.I):
            self.__selective = True
            self.coordinate_type = next(poscar)
        else:
            self.__selective = False
            self.coordinate_type = line7

        for line, elem in zip(poscar, self.atom_identifer):
            # if not elem: break
            tmp = line.split()
            self.position.append(np.float_(np.array(tmp[:3])))
            if self.is_selective:
                self.coordinate_changeflags.append(' '.join(tmp[3:]))

    @property
    def atom_identifer(self):
        '''Return list style of "atom_identifer" (e.g.  "#1:Ag1")'''
        # self.__atom_identifer = []
        # ii = 1
        # for elm, n in zip(self.iontype, self.ionnums):
        #     self.__atom_identifer.extend(
        #         '#{0}:{1}{2}'.format(ii + m, elm, m + 1) for m in range(n))
        #     ii += n
        # return self.__atom_identifer
        self.__atom_identifer = []
        atomnames = []
        for elm, ionnums in zip(self.iontype, self.ionnums):
            for j in range(1, ionnums + 1):
                tmp = elm + str(j)
                if tmp not in atomnames:
                    atomnames.append(tmp)
                else:
                    while tmp in atomnames:
                        j = j + 1
                        tmp = elm + str(j)
                    else:
                        atomnames.append(tmp)
        self.__atom_identifer = [
            "#" + str(s) + ":" + a for s, a in
            zip(range(1, len(atomnames) + 1), atomnames)]
        return self.__atom_identifer

    @atom_identifer.setter
    def atom_identifer(self, value):
        self.__atom_identifer = value

    def __iter__(self):
        for each in self.position:
            yield each

    @property
    def is_cartesian(self):
        '''.. py:method:: is_cartesian()

        Return True if Cartesian coordinate is set

        :return: True if coordinate is cartesian
        :rtype: Boolean
        '''
        return bool(re.search(r'^[ck]', self.coordinate_type, re.I))

    @property
    def is_direct(self):
        '''.. py:method:: is_direct()

        Return True if DIRECT coordinate is set

        :return: True if coordinate is direct (not cartesian)
        :rtype: Boolean
        '''
        return not self.is_cartesian

    @property
    def is_selective(self):
        '''.. py:method:: is_selective()

        Return True if "Selective Dynamcis" is set

        :return: True if "Selective Dynamics" switch on.
        :rtype: Boolean
        '''
        return self.__selective

    def pos(self, *i):
        '''.. py::method:: pos(i)

        Accessor of POSCAR.position.

        As in VASP, the atom index starts with "1", not "0".

        :param i: site indexes
        :type i: int, tuple, list, range

        :return: atom's position (When single value is set as i,
                 return just an atom position)

        :rtype: np.array or list of np.array

        .. warning:: the first site # is "1", not "0". (Follow VESTA's way.)

        .. warning:: Now, you **cannot** set the index range by using tuple.\
Use range object instead.  ex.) range(3,10) => (3, 4, 5, 6, 7, 8, 9)
        '''
        dest = []
        for thearg in i:
            if isinstance(thearg, int):
                if thearg <= 0:
                    raise ValueError
                else:
                    dest.append(self.position[thearg - 1])
            elif isinstance(thearg, (tuple, list, range)):
                for site_index in thearg:
                    if site_index <= 0:
                        raise ValueError
                    else:
                        dest.append(self.position[site_index - 1])
        if len(dest) == 1:
            dest = dest[0]
        return dest

    def average_position(self, *i):
        '''.. py:method:: average_position(*i)

        Return the average position of the sites

        :param i: site indexes
        :type i: int, tuple, list, range

        :return: atom's position

        :rtype: np.array
'''
        sitelist = []
        for thearg in i:
            if isinstance(thearg, int):
                sitelist.append(thearg)
            elif isinstance(thearg, (tuple, list, range)):
                for site_index in thearg:
                    sitelist.append(site_index)
        pos = self.pos(sitelist)
        if isinstance(pos, np.ndarray):
            return pos
        elif isinstance(pos, list):
            return sum(pos)/len(pos)

    def pos_replace(self, i, vector):
        '''.. py:method:: pos_replace(i, vector)

        :param i: site #
        :param vector: list of the i-th atom position.
        :type i: int
        :type vector: list, tuple, np.array
        :note: the first site # is "1", not "0" to follow VESTA's way.
'''
        vector = _vectorize(vector)
        if not isinstance(i, int):
            raise ValueError
        if not self.is_cartesian:
            message = 'poscar_replace method is implemented for'
            message += ' Cartesian coordinate'
            raise RuntimeError(message)
        self.position[i - 1] = vector

    def sort(self, from_index, to_index):
        '''Sort positions attribute by coordinate

        :param from_index: first index # for sort
        :param to_index: last index # for sort
        :type from_index: int
        :type to_index: last index # for sort

        :note: the first site # is "1", not "0" to follow VESTA's way.
        '''
        original_is_cartesian = False
        if self.is_cartesian:
            original_is_cartesian = True
            self.to_direct()
        self.repack_in_cell()
        poslists = self.pos(range(from_index, to_index+1))
        poslists.sort(keys=lambda x: (x[0], x[1], x[2]))
        for i in range(from_index, to_index + 1):
            self.pos_replace(i, poslists.pop(0))
        if original_is_cartesian:
            self.to_cartesian()

    def rotate_x(self, theta):
        ''' .. py:method:: rotate_x(theta)

        Rotation matrix around X-axis

        :param float theta: angle of rotation (Degrees)
        :return: rotation matrix
        :rtype: np.array

        :Example:

        >>> t = POSCAR()
        >>> t.rotate_x(60)
        array([[ 1.       ,  0.       ,  0.       ],
               [ 0.       ,  0.5      , -0.8660254],
               [ 0.       ,  0.8660254,  0.5      ]])
        '''
        degree = np.pi / 180.0
        return np.array(
            [[1.0, 0.0, 0.0],
             [0.0, np.cos(theta * degree), -np.sin(theta * degree)],
             [0.0, np.sin(theta * degree), np.cos(theta * degree)]])

    def rotate_y(self, theta):
        '''.. py:method:: rotate_y(theta)

        Rotation matrix around Y-axis

        :Example:

        >>> t = POSCAR()
        >>> t.rotate_y(60)
        array([[ 0.5      ,  0.       ,  0.8660254],
               [ 0.       ,  1.       ,  0.       ],
               [-0.8660254,  0.       ,  0.5      ]])
        '''
        degree = np.pi / 180.0
        return np.array(
            [[np.cos(theta * degree), 0.0, np.sin(theta * degree)],
             [0.0, 1.0, 0.0],
             [-np.sin(theta * degree), 0.0, np.cos(theta * degree)]])

    def rotate_z(self, theta):
        '''.. py:method:: rotate_y(theta)

        Rotation matrix around Z-axis

        :Example:

        >>> t = POSCAR()
        >>> t.rotate_z(60)
        array([[ 0.5      , -0.8660254,  0.       ],
               [ 0.8660254,  0.5      ,  0.       ],
               [ 0.       ,  0.       ,  1.       ]])
        '''
        degree = np.pi / 180.0
        return np.array(
            [[np.cos(theta * degree), -np.sin(theta * degree), 0.0],
             [np.sin(theta * degree), np.cos(theta * degree), 0.0],
             [0.0, 0.0, 1.0]])

    # class method? or independent function?
    def nearest(self, array, point):
        '''.. py:nearest::(array, point)

        :param array:
        :param point:
        :type array: list
        :type point: np.array
        :rtype: np.array
'''
        return min(array, key=lambda pos: np.linalg.norm(pos - point))

    # class method? or independent function?
    def make27candidate(self, position):
        '''.. py:make27candidate(position)

        Return 27 vectors set correspond the neiboring

        :param position: atom position defined in the coordinated by
                         cell_vecs ( scaling facter is not accounted).

        :param type: np.array, list
        :return: list of np.array
        :rtype: list

        '''
        position = _vectorize(position)
        candidates27 = []
        if self.is_cartesian:
            for l, m, n in it.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1]):
                candidates27.append(l * self.cell_vecs[0] +
                                    m * self.cell_vecs[1] +
                                    n * self.cell_vecs[2] +
                                    position)
        else:
            for l, m, n in it.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1]):
                candidates27.append(l * np.array([1., 0., 0.]) +
                                    m * np.array([0., 1., 0.]) +
                                    n * np.array([0., 0., 1.]) +
                                    position)
        return candidates27

    def atom_rotate(self, site, axis_name, theta, center):
        '''.. py:atom_rotate(site, axis_name, theta, center)

        Rotate the atom 

        :param site: site # for rotation (The first atom is "1".).
        :param axis_name: "X", "x", "Y", "y", "Z", or "z". Rotation axis.
        :param theta: Rotation angle (Degrees).
        :param center: center position for rotation.
        :type site: int
        :type axis_name: str
        :type theta: float
        :type center: np.array, list, tuple
        :todo:  check the center in the Braves lattice.
        :todo:  take into account the periodic boundary.
'''
        center = _vectorize(center)
        if len(center) != 3:
            raise ValueError
        if not self.point_in_box(center / self.scaling_factor, self.cell_vecs):
            raise ValueError('the center must be in the Braves lattice')
        if not isinstance(site, int):
            raise ValueError('argument error in atom_rotate method')
        if not self.is_cartesian:
            self.to_cartesian()
        position = self.pos(site)
        position -= center / self.scaling_factor
        position = getattr(self, 'rotate_' +
                           axis_name.lower())(theta).dot(position)
        position += center / self.scaling_factor
        self.pos_replace(site, position)

    def atoms_rotate(self, site_list, axis_name, theta, center):
        '''Rotate atoms

        :param site_list:  list array of site for rotation (The first
                               atom is "1".).

        :param axis_name: "X", "x", "Y", "y", "Z",or "z".
                          Rotation axis.

        :param theta: Rotation angle (Degrees).
        :param center: center position for rotation.
        :type site_list: list, tuple
        :type theta: float
        :type axis_name: str
        :type center: np.array, list, tuple

        '''
        for site in site_list:
            self.atom_rotate(site, axis_name, theta, center)

    def cell_rotate(self, theta, axis_name='Z'):
        '''Rotate unit-cell (rotation angle is set by degree.)

        :param theta: rotation angle
        :type theta: float
        :param axis_name: axis name for rotation (x, y, or z)
        :type axis_name: str
        :return: None'''
        original_is_cartesian = False
        if self.is_cartesian:
            original_is_cartesian = True
            self.to_direct()
        axis_name = axis_name.capitalize()

        if axis_name == 'X':
            self.cell_vecs = np.dot(self.rotate_x(theta), self.cell_vecs.T).T
        elif axis_name == 'Y':
            self.cell_vecs = np.dot(self.rotate_y(theta), self.cell_vecs.T).T
        elif axis_name == 'Z':
            self.cell_vecs = np.dot(self.rotate_z(theta), self.cell_vecs.T).T
        if original_is_cartesian:
            self.to_cartesian()

    def repack_in_cell(self):
        '''Repack all atoms in the unit cell

        No negative values in DIRECT coordinate.
        '''
        original_is_cartesian = False
        if self.is_cartesian:
            original_is_cartesian = True
            self.to_direct()
        for pos in self.position:
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
        :param other:
        :type other:  POSCAR
        :return: POSCAR
        :todo: Check the lattice vectors, coordinate_type and so on.
        '''
        if not isinstance(other, POSCAR):
            return NotImplemented
        dest_poscar = copy.deepcopy(self)
        if dest_poscar.scaling_factor != other.scaling_factor:
            raise ValueError('scaling factor is different.')
        if np.linalg.norm(dest_poscar.cell_vecs - other.cell_vecs) != 0:
            raise ValueError('lattice vectors (cell matrix) are different.')
        dest_poscar.iontype.extend(other.iontype)
        dest_poscar.ionnums.extend(other.ionnums)
        dest_poscar.position.extend(other.position)
        dest_poscar.coordinate_changeflags.extend(other.coordinate_changeflags)
#        dest_poscar.atom_identifer
        return dest_poscar

    def merge(self, other):
        '''lazy __add__: Return POSCAR generated from two POSCARs

        Even if the cell vectors and scaling factors are different,
        the 'merged' POSCAR is created.
        Use the cell vectors of the first POSCAR.

        :param other: POSCAR object
        :type other:  POSCAR
        :return: added poscar object
        :rtype: POSCAR
        '''
        if not isinstance(other, POSCAR):
            return NotImplemented
        dest_poscar = copy.deepcopy(self)
        original_is_direct = False
        if dest_poscar.is_direct:
            original_is_direct = True
            dest_poscar.to_cartesian()
        other_poscar = copy.deepcopy(other)
        original_scaling_factor = dest_poscar.scaling_factor
        other_poscar.tune_scaling_factor(original_scaling_factor)
        other_poscar.to_cartesian()
        dest_poscar.iontype.extend(other.iontype)
        dest_poscar.ionnums.extend(other.ionnums)
        dest_poscar.position.extend(other.position)
        dest_poscar.coordinate_changeflags.extend(other.coordinate_changeflags)
#        dest_poscar.atom_identifer
        if original_is_direct:
            dest_poscar.to_direct()
        return dest_poscar

    def to_list(self):
        '''Return POSCAR object by list-style

        :return: a list representation of POSCAR.
        :rtype: list
'''
        out_list = []
        out_list.append(self.system_name)
        out_list.append(self.scaling_factor)
        out_list.append(self.cell_vecs[0])
        out_list.append(self.cell_vecs[1])
        out_list.append(self.cell_vecs[2])
        if not self.iontype[0].isdigit():
            out_list.append(self.iontype)
        out_list.append(self.ionnums)
        if self.is_selective:
            out_list.append("Selective Dynamics")
        out_list.append(self.coordinate_type)
        out_list.append(self.position)
        out_list.append(self.coordinate_changeflags)
        out_list.append(self.atom_identifer)
        return out_list

    def __str__(self):
        '''
        :return: a string representation of POSCAR
        :rtype: string
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
        if not self.iontype[0].isdigit():
            tmp.append(' ' + ' '.join(self.iontype))
        tmp.append(' ' + ' '.join(str(i) for i in self.ionnums))
        if self.is_selective:
            tmp.append('Selective Dynamics')
        tmp.append(self.coordinate_type)
        for pos, t_or_f, atom in tools.ZIPLONG(self.position,  # Version safety
                                               self.coordinate_changeflags,
                                               self.atom_identifer,
                                               fillvalue=''):
            tmp.append(' '.join('  {0:20.17f}'.format(i) for i in pos) +
                       ' ' + t_or_f +
                       ' ' + atom)
        return '\n'.join(tmp) + '\n'

    def tune_scaling_factor(self, new_scaling_factor=1.0):
        '''.. py:method:: tune_scaling_factor(new_scaling_factor=1.0)

        Change scaling factor to new value

        :param new_scaling_factor:
        :type new_scaling_factor: float

        .. note::  **The Braves lattice are corrected (to be equal size)**

        .. warning:: If you change the cell size, change
                     scaling_factor attribute directly
'''
        old = self.scaling_factor
        self.cell_vecs *= (old / new_scaling_factor)
        self.scaling_factor = new_scaling_factor
        if self.is_cartesian:
            self.position = [i * old / new_scaling_factor
                             for i in self.position]

    def to_cartesian(self):
        '''.. py:method:: to_cartesian()

        Change the coordinate to cartesian from direct
'''
        if self.is_direct:
            self.coordinate_type = "Cartesian"
            mat = self.cell_vecs.transpose()
            self.position = [mat.dot(v) for v in self.position]

    def to_direct(self):
        '''.. py:method:: to_direct()

        Change the coordinate to direct from cartesian.
        '''
        if self.is_cartesian:
            self.coordinate_type = "Direct"
            mat = np.linalg.inv(np.transpose(self.cell_vecs))
            self.position = [mat.dot(v) for v in self.position]

    def guess_molecule(self, site_list, center=None):
        '''.. py:method::guess_molecule(site_list, center)

        Arrange atom position to form a molecule.

        This method is effective to rotate a molecule.

        :param site_list: list of site number
        :param center: center position of "molecule" (Optional).
        :type site_list: list
        :type center: list
        :return: Array of Vector that represents "molecule".
        :rtype: numpy.array

        .. note:: When the optional argument, center, is set, the
                  atoms are re-arranged as to minimize the distance
                  from this center.  If not, atoms are re-arranged to
                  minimize the total bonding length.  As the algorithm
                  for mimizing the total length is not exhaustive, the
                  resultant atom arrangement may different from what
                  you expect, in spite of time-waste.  The center
                  option is highly recommended to form a molecule.

        '''
        molecule = [self.pos(j) for j in site_list]
        for index, site in enumerate(site_list):
            target_atom = self.pos(site)
            atoms27 = self.make27candidate(target_atom)

            def func(pos):
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
            newpos = min(atoms27, key=func)
            molecule[index] = newpos
        return molecule

    def guess_molecule2(self, site_list):
        '''Arranges atom positions to form a molecule.

        poscar updates

        :param site_list: list of site number
        :type site_list: list
'''
        molecule = self.guess_molecule(site_list)
        for site, pos_vector in zip(site_list, molecule):
            self.pos_replace(site, pos_vector)

    def translate(self, vector, atomlist):
        '''Translate the selected atom(s) by vector

        :param vector: translational vector (in Cartesian frame)
        :param atomlist: list of the atom to be moved
        :type vector: list, tuple, numpy.array
        :type atomlist: list
        :note: the first atom is "1", not "0".
        :return: position
        :rtype: numpy.array
'''
        if self.is_cartesian:
            vector = _vectorize(vector)
            for i in atomlist:
                self.position[i - 1] = (self.position[i - 1] +
                                        vector / self.scaling_factor)
        else:
            vector = _vectorize(vector)
            self.to_cartesian()
            for i in atomlist:
                self.position[i - 1] = (self.position[i - 1] +
                                        vector / self.scaling_factor)
            self.to_direct()
        return self.position

    def translate_all(self, vector):
        '''Translate **all** atoms by vector

        :param vector: translational vector
        :type vector: list, numpy.array
'''
        atomrange = list(range(1, sum(self.ionnums) + 1))
        self.translate(vector, atomrange)

    def save(self, filename):
        '''Save POSCAR contents to the file named "filename"
'''
        try:  # Version safety
            file = open(filename, mode='w', newline='\n')
        except TypeError:
            file = open(filename, mode='wb')
        with file:
            file.write(str(self))

    def point_in_box(self, point, cell_vecs):
        '''Return True if point is located in the box

        :param point: vector representing the "point"
        :type point: numpy.ndarray, numpy.matrix, list, tuple
        :param cell_vecs: vectors defining the "box"
        :type cell_vecs: numpy.ndarray, numpy.matrix, list, tuple
        :rtype: boolean
        '''
        if three_by_three(cell_vecs):
            point = np.array(point).flatten()
            cell_vecs = np.array(cell_vecs)
            result = np.dot(np.linalg.inv(cell_vecs.T), point)
            return all((0 <= float(q) <= 1) for q in result)
        else:
            raise TypeError


def three_by_three(vec):
    '''.. py:function:three_by_three(vec)

    Return True if vec can be converted into the 3x3 matrix
    :param vec: list like object
    :type vec: numpy.ndarray, numpy.matrix, list, tuple
    :rtype: boolean
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
