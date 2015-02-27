#! /usr/bin/env python
# -*- conding: utf-8 -*-
'''This mudule provides POSCAR class

translate from poscar.rb of 2014/2/26, master branch
'''
from __future__ import division, print_function  # Version safety
import numpy as np
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

'''
 fcc (111) surface
   3.52000000000000
     0.7071067800000000    0.0000000000000000    0.0000000000000000
    -0.3535533900000000    0.6123724000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000   11.5470000000000006
   Cu
     9
Selective dynamics
Direct
  0.0000000000000000  0.0000000000000000  0.0000000000000000   F   F   F
  0.3333333300000021  0.6666666699999979  0.0499999999999972   F   F   F
  0.6666666699999979  0.3333333300000021  0.1000000000000014   F   F   F
  0.0000000000000000  0.0000000000000000  0.1499999999999986   T   T   T
  0.3333333300000021  0.6666666699999979  0.2000000000000028   T   T   T
  0.0000000000000000  0.0000000000000000  0.2500000000000000   T   T   T
  0.3333333300000021  0.6666666699999979  0.2999999999999972   T   T   T
  0.6666666699999979  0.3333333300000021  0.3500000000000014   T   T   T
  0.0000000000000000  0.0000000000000000  0.3999999999999986   T   T   T

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


class POSCAR(object):

    '''
    .. py:class:: POSCAR(object)

    Class for POSCAR (CONTCAR) format

    This script does *NOT* support for constructing POSCAR
    from scratch. (Use ASE for this purpose.)

    It provides a way to slightly modify POSCAR or
    CONTCAR, which has already works well.

    :attribute: latticeV1, latticeV1, latticeV3
    :version: 2.0
    '''

    def __init__(self, arg=None):
        '''
        :param arg: POSCAR file name, or list of POSCAR text.
        :type arg: str
        '''
        self.system_name = ""
        self.scaling_factor = 0.
        self.__latticeV1 = np.array([0., 0., 0.])
        self.__latticeV2 = np.array([0., 0., 0.])
        self.__latticeV3 = np.array([0., 0., 0.])
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
    def latticeV1(self):
        '''The vector of the first axis of the unit cell'''
        return self.__latticeV1

    @property
    def latticeV2(self):
        '''The vector of the second axis of the unit cell'''
        return self.__latticeV2

    @property
    def latticeV3(self):
        '''The vector of the third axis of the unit cell'''
        return self.__latticeV3

    def load_from_array(self, poscar):
        '''POSCAR parser'''
        poscar = iter(map(str.rstrip, poscar))  # Version safety
        self.system_name = next(poscar)
        self.scaling_factor = float(next(poscar))
        self.__latticeV1 = np.array(
            list(map(float, next(poscar).split())))  # Version safety
        self.__latticeV2 = np.array(
            list(map(float, next(poscar).split())))  # Version safety
        self.__latticeV3 = np.array(
            list(map(float, next(poscar).split())))  # Version safety
        self.iontype = next(poscar).split()
        # parse POSCAR evenif the element names are not set.
        # At present, the String representation number
        #   are used for the  dummy name.
        if self.iontype[0].isdigit():
            self.ionnums = list(map(int, self.iontype))
            #                      [int(i) for i in self.iontype]
        else:
            self.ionnums = list(
                map(int, next(poscar).split()))  # Version safety
        self.__atom_identifer = []
        atomnames = []
        for elm, n in zip(self.iontype, self.ionnums):
            for j in range(1, n + 1):
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
        for elm, n in zip(self.iontype, self.ionnums):
            for j in range(1, n + 1):
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
        '''Return True if Cartesian coordinate is set

        :return: True if coordinate is cartesian
        :rtype: Boolean
        '''
        return bool(re.search(r'^[ck]', self.coordinate_type, re.I))

    @property
    def is_direct(self):
        '''Return True if DIRECT coordinate is set

        :return: True if coordinate is direct (not cartesian)
        :rtype: Boolean
        '''
        return not self.is_cartesian

    @property
    def is_selective(self):
        '''Return True if "Selective Dynamcis" is set

        :return: True if "Selective Dynamics" switch on.
        :rtype: Boolean
        '''
        return self.__selective

    def pos(self, *i):
        '''Accessor of POSCAR.position.

        As in VASP, the atom index starts with "1", not "0".

        :param i: site index (index range can be set by tuple (from, to))
        :type i: int, or tuple of two int

        :return: atom's position (When single value is set as i,
                 return just an atom position)

        :rtype: np.array or list of np.array

        .. note:: You can set the index range by using tuple.

        :Example: (20, 25) means the index from site# 20 to site# 25.

        .. warning:: the first site # is "1", not "0". (Follow VESTA's way.)

        .. warning:: (to user who knows Python's style of list slice)
                     When set the site index range by tuple, the
                     latter index is included.  (pos(20, 25) returns
                     position[19:25], not position[19:24])

        '''
        if type(i[0]) == int:
            dest = []
            for ii in i:
                dest.append(self.position[ii - 1])
            if len(dest) == 1:
                dest = dest[0]
        elif type(i[0]) == tuple:
            if i[0] > i[1]:
                raise ValueError
            dest = self.positions[i[0] - 1: i[1]]
        return dest

    def pos_replace(self, i, vector):
        '''
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
            self.to_Direct()
        self.repack_in_cell()
        poslists = self.pos((from_index, to_index))
        poslists.sort(keys=lambda x: (x[0], x[1], x[2]))
        for i in range(from_index, to_index + 1):
            self.pos_replace(i, poslists.pop(0))
        if original_is_cartesian:
            self.to_Cartesian()

    def rotateX(self, theta):
        '''Rotation matrix around X-axis

        :param theta: angle of rotation (Degrees)
        :type theta: float
        :return: rotation matrix
        :rtype: np.array

        :Example:

        >>> t = POSCAR()
        >>> t.rotateX(60)
        array([[ 1.       ,  0.       ,  0.       ],
               [ 0.       ,  0.5      , -0.8660254],
               [ 0.       ,  0.8660254,  0.5      ]])
        '''
        degree = np.pi / 180.0
        return np.array(
            [[1.0, 0.0, 0.0],
             [0.0, np.cos(theta * degree), -np.sin(theta * degree)],
             [0.0, np.sin(theta * degree), np.cos(theta * degree)]])

    def rotateY(self, theta):
        '''Rotation matrix around Y-axis

        :Example:

        >>> t = POSCAR()
        >>> t.rotateY(60)
        array([[ 0.5      ,  0.       ,  0.8660254],
               [ 0.       ,  1.       ,  0.       ],
               [-0.8660254,  0.       ,  0.5      ]])
        '''
        degree = np.pi / 180.0
        return np.array(
            [[np.cos(theta * degree), 0.0, np.sin(theta * degree)],
             [0.0, 1.0, 0.0],
             [-np.sin(theta * degree), 0.0, np.cos(theta * degree)]])

    def rotateZ(self, theta):
        '''Rotation matrix around Z-axis

        :Example:

        >>> t = POSCAR()
        >>> t.rotateZ(60)
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
        '''

        :param array:
        :param point:
        :type array: list
        :type point: np.array
        :rtype: np.array
'''
        return min(array, key=lambda pos: np.linalg.norm(pos - point))

    # class method? or independent function?
    def make27candidate(self, position):
        '''Return 27 vector sets correspond the neiboring

        :param position: atom position defined in the coordinated by
                         latticeV1, latticeV2, latticeV3 ( scaling
                         facter is not accounted).

        :param type: np.array, list
        :return: list of np.array
        :rtype: list

        '''
        position = _vectorize(position)
        candidates27 = []
        if self.is_cartesian:
            for l, m, n in it.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1]):
                candidates27.append(l * self.latticeV1 +
                                    m * self.latticeV2 +
                                    n * self.latticeV3 +
                                    position)
        else:
            for l, m, n in it.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1]):
                candidates27.append(l * np.array([1., 0., 0.]) +
                                    m * np.array([0., 1., 0.]) +
                                    n * np.array([0., 0., 1.]) +
                                    position)
        return candidates27

    def atom_rotate(self, site, axis_name, theta, center):
        '''Rotate atom (single atom!)

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
        if not self.point_in_box(center / self.scaling_factor,
                                 self.latticeV1,
                                 self.latticeV2,
                                 self.latticeV3):
            raise ValueError('the center must be in the Braves lattice')
        if not isinstance(site, int):
            raise ValueError('argument error in atom_rotate method')
        if not self.is_cartesian:
            self.to_Cartesian()
        position = self.pos(site)
        position -= center / self.scaling_factor
        position = getattr(self, 'rotate' +
                           axis_name.capitalize())(theta).dot(position)
        position += center / self.scaling_factor
        self.pos_replace(site, position)

    def atoms_rotate(self, site_list_pack, axis_name, theta, center):
        '''Rotate atoms

        :param site_list_pack: list array of the list array (not
                               typo!) of site for rotation (The first
                               atom is "1".).

        :param axis_name: "X", "x", "Y", "y", "Z",or "z".
                          Rotation axis.

        :param theta: Rotation angle (Degrees).
        :param center: center position for rotation.
        :type site_list_pack: list, tuple
        :type theta: float
        :type axis_name: str
        :type center: np.array, list, tuple

        '''
        for site_list in site_list_pack:
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
            self.to_Direct()
        axis_name = axis_name.capitalize()
        if axis_name == 'X':
            self.__latticeV1 = np.dot(self.rotateX(theta), self.latticeV1)
            self.__latticeV2 = np.dot(self.rotateX(theta), self.latticeV2)
            self.__latticeV3 = np.dot(self.rotateX(theta), self.latticeV3)
        elif axis_name == 'Y':
            self.__latticeV1 = np.dot(self.rotateY(theta), self.latticeV1)
            self.__latticeV2 = np.dot(self.rotateY(theta), self.latticeV2)
            self.__latticeV3 = np.dot(self.rotateY(theta), self.latticeV3)
        elif axis_name == 'Z':
            self.__latticeV1 = np.dot(self.rotateZ(theta), self.latticeV1)
            self.__latticeV2 = np.dot(self.rotateZ(theta), self.latticeV2)
            self.__latticeV3 = np.dot(self.rotateZ(theta), self.latticeV3)
        if original_is_cartesian:
            self.to_Cartesian()

    def repack_in_cell(self):
        '''Repack all atoms in the unit cell

        No negative values in DIRECT coordinate.
        '''
        original_is_cartesian = False
        if self.is_cartesian:
            original_is_cartesian = True
            self.to_Direct()
        for pos in self.position:
            for i in (0, 1, 2):
                while pos[i] < 0.0 or pos[i] > 1.0:
                    if pos[i] < 0.0:
                        pos[i] += 1.0
                    elif pos[i] > 1.0:
                        pos[i] -= 1.0
        if original_is_cartesian:
            self.to_Cartesian()

    def __add__(self, other):
        '''
        :param other:
        :type other:  POSCAR
        :return: POSCAR
        :todo: Check the lattice vectors, coordinate_type and so on.
        '''
        if not isinstance(other, POSCAR):
            return NotImplemented
        destPOSCAR = copy.deepcopy(self)
        if destPOSCAR.scaling_factor != other.scaling_factor:
            raise ValueError('scaling factor is different.')
        if np.linalg.norm(destPOSCAR.latticeV1 +
                          destPOSCAR.latticeV2 +
                          destPOSCAR.latticeV3 -
                          (other.latticeV1 +
                           other.latticeV2 +
                           other.latticeV3)) != 0:
            raise ValueError('lattice vectors are different.')
        destPOSCAR.iontype.extend(other.iontype)
        destPOSCAR.ionnums.extend(other.ionnums)
        destPOSCAR.position.extend(other.position)
        destPOSCAR.coordinate_changeflags.extend(other.coordinate_changeflags)

        destPOSCAR.atom_identifer

        return destPOSCAR

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
        destPOSCAR = copy.deepcopy(self)
        original_is_direct = False
        if destPOSCAR.is_direct:
            original_is_direct = True
            destPOSCAR.to_Cartesian()
        otherPOSCAR = copy.deepcopy(other)
        original_scaling_factor = destPOSCAR.scaling_factor
        otherPOSCAR.tune_scaling_factor(original_scaling_factor)
        otherPOSCAR.to_Cartesian()
        destPOSCAR.iontype.extend(other.iontype)
        destPOSCAR.ionnums.extend(other.ionnums)
        destPOSCAR.position.extend(other.position)
        destPOSCAR.coordinate_changeflags.extend(other.coordinate_changeflags)
        destPOSCAR.atom_identifer
        if original_is_direct:
            destPOSCAR.to_Direct()
        return destPOSCAR

    def to_list(self):
        '''Return POSCAR object by list-style

        :return: a list representation of POSCAR.
        :rtype: list
'''
        out_list = []
        out_list.append(self.system_name)
        out_list.append(self.scaling_factor)
        out_list.append(self.latticeV1)
        out_list.append(self.latticeV2)
        out_list.append(self.latticeV3)
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
        tmp.append(''.join('   {0:20.17f}'.format(i) for i in self.latticeV1))
        tmp.append(''.join('   {0:20.17f}'.format(i) for i in self.latticeV2))
        tmp.append(''.join('   {0:20.17f}'.format(i) for i in self.latticeV3))
        if not self.iontype[0].isdigit():
            tmp.append(' ' + ' '.join(self.iontype))
        tmp.append(' ' + ' '.join(str(i) for i in self.ionnums))
        if self.is_selective:
            tmp.append('Selective Dynamics')
        tmp.append(self.coordinate_type)
        for pos, tf, atom in tools.ziplong(self.position,  # Version safety
                                           self.coordinate_changeflags,
                                           self.atom_identifer, fillvalue=''):
            tmp.append(' '.join('  {0:20.17f}'.format(i) for i in pos) +
                       ' ' + tf +
                       ' ' + atom)
        return '\n'.join(tmp) + '\n'

    def tune_scaling_factor(self, new_scaling_factor=1.0):
        '''Change scaling factor to new value

        :param new_scaling_factor:
        :type new_scaling_factor: float

        .. note::  **The Braves lattice are corrected (to be equal size)**

        .. warning:: If you change the cell size, change
                     scaling_factor attribute directly
'''
        old = self.scaling_factor
        self.__latticeV1 *= (old / new_scaling_factor)
        self.__latticeV2 *= (old / new_scaling_factor)
        self.__latticeV3 *= (old / new_scaling_factor)
        self.scaling_factor = new_scaling_factor
        if self.is_cartesian:
            self.position = [i * old / new_scaling_factor
                             for i in self.position]

    def to_Cartesian(self):
        '''Change the coordinate to cartesian from direct
'''
        if self.is_direct:
            self.coordinate_type = "Cartesian"
            m = np.array([self.latticeV1,
                          self.latticeV2,
                          self.latticeV3]).transpose()
            self.position = [m.dot(v) for v in self.position]

    def to_Direct(self):
        '''Change the coordinate to direct from cartesian.
        '''
        if self.is_cartesian:
            self.coordinate_type = "Direct"
            m = np.linalg.inv(np.transpose(np.array([self.latticeV1,
                                                     self.latticeV2,
                                                     self.latticeV3])))
            self.position = [m.dot(v) for v in self.position]

    def guess_molecule(self, site_list, center=None):
        '''Arrange atom position to form a molecule.

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
        for site, posVector in zip(site_list, molecule):
            self.pos_replace(site, posVector)

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
            self.to_Cartesian()
            for i in atomlist:
                self.position[i - 1] = (self.position[i - 1] +
                                        vector / self.scaling_factor)
            self.to_Direct()
        return self.position

    def translate_all(self, vector):
        '''Translate **all** atoms by vector

        :param vector: translational vector
        :type vector: list, numpy.array
'''
        atomrange = list(range(1, sum(self.ionnums) + 1))
        self.translate(vector, atomrange)

    def save(self, filename):
        try:  # Version safety
            file = open(filename, mode='w', newline='\n')
        except TypeError:
            file = open(filename, mode='wb')
        with file:
            file.write(str(self))

    def point_in_box(self, point, l1, l2, l3):
        '''Return True if point is located in the box

        .. note:: The box is defined by three vectors (l1, l2, and l3)

        :param point: vector representing the "point"
        :type point: numpy.ndarray, numpy.matrix, list, tuple
        :param l1: vector representing the "box"
        :param l2: vector representing the "box"
        :param l3: vector representing the "box"
        :type l1: numpy.ndarray, numpy.matrix, list, tuple
        :type l2: numpy.ndarray, numpy.matrix, list, tuple
        :type l3: numpy.ndarray, numpy.matrix, list, tuple
        :rtype: Boolean
        '''
        x0, y0, z0 = _vectorize(point).flatten()
        x1, y1, z1 = _vectorize(l1).flatten()
        x2, y2, z2 = _vectorize(l2).flatten()
        x3, y3, z3 = _vectorize(l3).flatten()
        l = -(-x3 * y2 * z0 + x2 * y3 * z0 + x3 * y0 * z2 - x0 * y3 * z2 -
              x2 * y0 * z3 + x0 * y2 * z3) / (x3 * y2 * z1 - x2 * y3 * z1 -
                                              x3 * y1 * z2 + x1 * y3 * z2 +
                                              x2 * y1 * z3 - x1 * y2 * z3)
        m = -(x3 * y1 * z0 - x1 * y3 * z0 - x3 * y0 * z1 + x0 * y3 * z1 +
              x1 * y0 * z3 - x0 * y1 * z3) / (x3 * y2 * z1 - x2 * y3 * z1 -
                                              x3 * y1 * z2 + x1 * y3 * z2 +
                                              x2 * y1 * z3 - x1 * y2 * z3)
        n = -(x2 * y1 * z0 - x1 * y2 * z0 - x2 * y0 * z1 + x0 * y2 * z1 +
              x1 * y0 * z2 - x0 * y1 * z2) / (-x3 * y2 * z1 + x2 * y3 * z1 +
                                              x3 * y1 * z2 - x1 * y3 * z2 -
                                              x2 * y1 * z3 + x1 * y2 * z3)
        return all((0 <= q <= 1) for q in (l, m, n))


def _vectorize(vector):
    if not isinstance(vector, _vector_acceptables):
        raise TypeError('Cannot convert into vector.')
    return np.array(vector).flatten()

_vector_acceptables = (np.ndarray, np.matrix, list, tuple)
# --------------------------

if __name__ == '__main__':
    # $-w = true
    import argparse
    import functools as ft
    arg = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""NOTE: When you running this script
on Windows Power Shell,
commas are regarded as delimiter of values.
So you must enclose values which
contains commas with quotations.
(ex.) --atom 1,2,3 -> failure / --atom "1,2,3" -> OK""")
    arg.add_argument('--atom', metavar='atoms', action='append',
                     type=tools.parse_AtomselectionNum,
                     help='''atoms specified with range using "-"
or comma-delimnated numbers.
 (ex.) --atom 1,2,7-9''')

    def split_to_float(string, n, name):
        lis = string.split(',')
        if len(lis) != n:
            message = '--{0} option requires {1} numbers'.format(name, n)
            raise argparse.ArgumentTypeError(message)
        return [float(i) for i in lis]
    arg.add_argument('--translate', metavar='x,y,z', action='append',
                     type=ft.partial(split_to_float, n=3, name='translate'),
                     help='''displacement (AA unit) by three numbers
separated by comma.''')
    arg.add_argument('--rotateX', metavar='theta,x,y,z',
                     type=ft.partial(split_to_float, n=4, name='rotateX'),
                     help='''Rotation around X-axis by "theta" at (x,y,z)
NOTE: this option is not taken into account
the periodic boundary.''')
    arg.add_argument('--rotateY', metavar='theta,x,y,z',
                     type=ft.partial(split_to_float, n=4, name='rotateY'),
                     help='''Rotation around Y-axis by "theta" at (x,y,z)
NOTE: this option is not taken into account
the periodic boundary.''')
    arg.add_argument('--rotateZ', metavar='theta,x,y,z',
                     type=ft.partial(split_to_float, n=4, name='rotateZ'),
                     help='''Rotation around Z-axis by "theta" at (x,y,z)
NOTE: this option is not taken into account
the periodic boundary.''')
    arg.add_argument('--output', metavar='file_name',
                     help='''output file name
if not specified, use standard output''')
    arg.add_argument('poscar', metavar='POSCAR_file (or CONTCAR_file)',
                     type=POSCAR)
    arguments = arg.parse_args()
    # translate option and rotate option are not set simulaneously.
    if arguments.translate and any([arguments.rotateX,
                                    arguments.rotateY,
                                    arguments.rotateZ]):
        raise RuntimeError(
            "Cannot set --translate and rotate option simultanaously.")
    # rotate options are not set multiply.
    if (arguments.rotateX,
            arguments.rotateY,
            arguments.rotateZ).count(None) < 2:
        raise RuntimeError(
            "Cannot set multiple rotate options simultanaously.")
    # print(arguments.poscar) #DEBUG
    arguments.poscar.to_Cartesian()
    # print(arguments.poscar) #DEBUG

    #
    # if "atom" option is not set, all atoms are concerned.
    #
    if not arguments.atom:
        nAtoms = sum(arguments.poscar.ionnums)
        arguments.atom = [tools.parse_AtomselectionNum('1-{0}'.format(nAtoms))]
    #
    #  Translation
    #
    if arguments.translate:
        if len(arguments.atom) != len(arguments.translate):
            raise RuntimeError
        for v, a in zip(arguments.translate, arguments.atom):
            arguments.poscar.translate(v, a)

    #
    #  Rotation
    #
    if any([arguments.rotateX, arguments.rotateY, arguments.rotateZ]):
        if len(arguments.atom) != 1:
            raise RuntimeError("--atom option set once!")
        if arguments.rotateX:
            axis_name = 'X'
            theta = arguments.rotateX[0]
            center = arguments.rotateX[1:]
        elif arguments.rotateY:
            axis_name = 'Y'
            theta = arguments.rotateY[0]
            center = arguments.rotateY[1:]
        elif arguments.rotateZ:
            axis_name = 'Z'
            theta = arguments.rotateZ[0]
            center = arguments.rotateZ[1:]
        arguments.poscar.atoms_rotate(arguments.atom, axis_name, theta, center)

    #
    #  Output result
    #
    if arguments.output is not None:
        arguments.poscar.save(arguments.output)
    else:
        print(arguments.poscar)

##
##
##

if __name__ == "__main__":
    import doctest
    doctest.testmod()
