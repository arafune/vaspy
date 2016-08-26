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
import bz2
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
        self.selective = False
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
            self.selective = True
            self.coordinate_type = next(poscar)
        else:
            self.selective = False
            self.coordinate_type = line7

        for line, elem in zip(poscar, self.atom_identifer):
            # if not elem: break
            tmp = line.split()
            self.position.append(np.float_(np.array(tmp[:3])))
            if self.selective:
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

    def is_cartesian(self):
        '''.. py:method:: is_cartesian()

        Return True if Cartesian coordinate is set

        Returns
        --------

        boolean
            True if coordinate is cartesian
        '''
        return bool(re.search(r'^[ck]', self.coordinate_type, re.I))

    def is_direct(self):
        '''.. py:method:: is_direct()

        Return True if DIRECT coordinate is set

        Returns
        --------

        Boolean
            True if coordinate is direct (not cartesian)
        '''
        return not self.is_cartesian()

    def pos(self, *i):
        '''.. py:method:: pos(i)

        Accessor of POSCAR.position.

        As in VASP, the atom index starts with "1", not "0".

        Parameters
        -----------

        i: int, tuple, list, range
             site indexes

        Returns
        --------

        numpy.ndarray, list of numpy.ndarray
             atom's position (When single value is set as i,\
                             return just an atom position)

        Warning
        --------

        the first site # is "1", not "0". (Follow VESTA's way.)

        Warning
        -------

        Now, you **cannot** set the index range by using tuple.
        Use range object instead.
        ex.) range(3,10) => (3, 4, 5, 6, 7, 8, 9)
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

    def sort(self, s_site=1, e_site=None, axis='z'):
        '''.. py:method:: sort(s_site, e_site, axis='z')

        Sort element position by the position

        Parameters
        -----------

        s_site: int
            start site index for sort

        e_site: int
            end site index for sort

        axis: str
            Axis used for sort (Default Z)

        Warning
        ------------

        The site index number starts with "1", not "0".
        The element difference is **not** taken into account.
        '''
        if e_site is None:
            e_site = sum(self.ionnums)

        if axis == 'x' or axis == 'X' or axis == 0:
            axis = 0
        elif axis == 'y' or axis == 'Y' or axis == 1:
            axis = 1
        elif axis == 'z' or axis == 'Z' or axis == 2:
            axis = 2
        self.position = self.position[0:s_site-1] + sorted(
            self.position[s_site-1:e_site],
            key=lambda sortaxis: sortaxis[axis]) + self.position[e_site:]

    def average_position(self, *i):
        '''.. py:method:: average_position(*i)

        Return the average position of the sites

        Parameters
        -----------

        i: int, tuple, list, range
            site indexes

        Returns
        --------

        numpy.ndarray
            atom's position
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

        Parameters
        -----------
        i: int
            site #
        vector: list, tuple, numpy.ndarray
            list of the i-th atom position.

        Notes
        ------

        the first site # is "1", not "0" to follow VESTA's way.
        '''
        vector = _vectorize(vector)
        if not isinstance(i, int):
            raise ValueError
        if not self.is_cartesian():
            message = 'poscar_replace method is implemented for'
            message += ' Cartesian coordinate'
            raise RuntimeError(message)
        self.position[i - 1] = vector

    def supercell(self, nx, ny, nz):
        '''.. py:method:: supercell(nx, ny, nz)

        Return the :math:`(nx \\times ny \\times nz)` supercell

        Parameters
        -----------
        nz: int
            repeat number along x axis
        ny: int
            repeat number along y axis
        nz: int
            repeat number along z axis


        Returns
        ---------

        POSCAR
            POSCAR object of the supercell
        '''
        if not isinstance(nx, int) \
           or not isinstance(ny, int) \
           or not isinstance(nz, int):
            raise ValueError("arguments must be positive integer")
        if nx <= 0 or ny <= 0 or nz <= 0:
            raise ValueError("arguments must be positive integer")
        #
        sposcar = copy.deepcopy(self)
        original_is_cartesian = sposcar.is_cartesian()
        if original_is_cartesian:
            sposcar.to_direct()
        sposcar.repack_in_cell()
        sposcar.cell_vecs[0] = sposcar.cell_vecs[0] * nx
        sposcar.cell_vecs[1] = sposcar.cell_vecs[1] * ny
        sposcar.cell_vecs[2] = sposcar.cell_vecs[2] * nz
        sposcar.ionnums = [i * nx * ny * nz for i in sposcar.ionnums]
        sposition = sposcar.position
        sposcar.position = []
        sposition = [np.array([x[0] / nx, x[1] / ny, x[2] / nz])
                     for x in sposition]
        for spos in sposition:
            for iz in range(0, nz):
                for iy in range(0, ny):
                    for ix in range(0, nx):
                        sposcar.position.append(np.array(
                            [spos[0] + ix / nx,
                             spos[1] + iy / ny,
                             spos[2] + iz / nz]))
        sposcar.coordinate_changeflags = []
        for flags in self.coordinate_changeflags:
            for i in range(nx * ny * nz):
                sposcar.coordinate_changeflags.append(flags)
        sposcar.__atom_identifer = []
        if original_is_cartesian:
            sposcar.to_cartesian()
        return sposcar

    def sort(self, from_index, to_index):
        '''.. py:method:: sort(from_index, to_index)

        Sort positions attribute by coordinate

        Parameters
        -----------

        from_index: int
            first index # for sort
        to_index: int
            last index # for sort


        Notes
        -----

        the first site # is "1", not "0" to follow VESTA's way.
        '''
        original_is_cartesian = False
        if self.is_cartesian():
            original_is_cartesian = True
            self.to_direct()
        self.repack_in_cell()
        poslists = self.pos(range(from_index, to_index+1))
        poslists.sort(keys=lambda x: (x[0], x[1], x[2]))
        for i in range(from_index, to_index + 1):
            self.pos_replace(i, poslists.pop(0))
        if original_is_cartesian:
            self.to_cartesian()

    # class method? or independent function?
    def nearest(self, array, point):
        '''.. py:method:: nearest(array, point)

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
        '''.. py:make27candidate(position)

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

        Parameters
        ----------

        site: int
            site # for rotation (The first atom is "1".).
        axis_name: str
            "X", "x", "Y", "y", "Z", or "z". Rotation axis.
        theta: float
            Rotation angle (Degrees).
        center: numpy.ndarray, list, tuple
            center position for rotation.

        Todo
        -----
        check the center in the Braves lattice.
        take into account the periodic boundary.
        '''
        center = _vectorize(center)
        if len(center) != 3:
            raise ValueError
        if not point_in_box(center / self.scaling_factor, self.cell_vecs):
            raise ValueError('the center must be in the Braves lattice')
        if not isinstance(site, int):
            raise ValueError('argument error in atom_rotate method')
        if not self.is_cartesian():
            self.to_cartesian()
        position = self.pos(site)
        position -= center / self.scaling_factor
        position = globals()["rotate_"+axis_name.lower()](theta).dot(position)
        position += center / self.scaling_factor
        self.pos_replace(site, position)

    def atoms_rotate(self, site_list, axis_name, theta, center):
        '''Rotate atoms

        Parameters
        ----------

        site_list:
            list array of site for rotation (The first atom is "1".).
        axis_name:
            "X", "x", "Y", "y", "Z",or "z".  Rotation axis.
        theta: float
             Rotation angle (Degrees).
        center: numpy.ndarray, list, tuple
             Position of rotation center
        '''
        for site in site_list:
            self.atom_rotate(site, axis_name, theta, center)

    def cell_rotate(self, theta, axis_name='Z'):
        '''.. py:method:: cell_rotate(theta, axis_name)

        Rotate unit-cell (rotation angle is set by degree.)

        Parameters
        ----------

        theta: float
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
            self.cell_vecs = np.dot(rotate_x(theta), self.cell_vecs.T).T
        elif axis_name == 'Y':
            self.cell_vecs = np.dot(rotate_y(theta), self.cell_vecs.T).T
        elif axis_name == 'Z':
            self.cell_vecs = np.dot(rotate_z(theta), self.cell_vecs.T).T
        if original_is_cartesian:
            self.to_cartesian()

    def repack_in_cell(self):
        '''.. py:method:repack_in_cell()

        Repack all atoms in the unit cell

        No negative values in DIRECT coordinate.
        '''
        original_is_cartesian = False
        if self.is_cartesian():
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
        '''.. py:method:: __add__(other)

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
        dest_poscar.iontype.extend(other.iontype)
        dest_poscar.ionnums.extend(other.ionnums)
        dest_poscar.position.extend(other.position)
        dest_poscar.coordinate_changeflags.extend(other.coordinate_changeflags)
#        dest_poscar.atom_identifer
        return dest_poscar

    def merge(self, other):
        '''.. py:method:: merge(other)

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
        dest_poscar.iontype.extend(other.iontype)
        dest_poscar.ionnums.extend(other.ionnums)
        dest_poscar.position.extend(other.position)
        dest_poscar.coordinate_changeflags.extend(other.coordinate_changeflags)
#        dest_poscar.atom_identifer
        if original_is_direct:
            dest_poscar.to_direct()
        return dest_poscar

    def to_list(self):
        '''.. py:method:: to_list()

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
        if not self.iontype[0].isdigit():
            out_list.append(self.iontype)
        out_list.append(self.ionnums)
        if self.selective:
            out_list.append("Selective Dynamics")
        out_list.append(self.coordinate_type)
        out_list.append(self.position)
        out_list.append(self.coordinate_changeflags)
        out_list.append(self.atom_identifer)
        return out_list

    def __str__(self):
        '''.. py:method:: __str__()

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
        if not self.iontype[0].isdigit():
            tmp.append(' ' + ' '.join(self.iontype))
        tmp.append(' ' + ' '.join(str(i) for i in self.ionnums))
        if self.selective:
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

        Parameters
        ------------

        new_scaling_factor: float

        Notes
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
            self.position = [i * old / new_scaling_factor
                             for i in self.position]

    def to_cartesian(self):
        '''.. py:method:: to_cartesian()

        Change the coordinate to cartesian from direct
        '''
        if self.is_direct():
            self.coordinate_type = "Cartesian"
            mat = self.cell_vecs.transpose()
            self.position = [mat.dot(v) for v in self.position]

    def to_direct(self):
        '''.. py:method:: to_direct()

        Change the coordinate to direct from cartesian.
        '''
        if self.is_cartesian():
            self.coordinate_type = "Direct"
            mat = np.linalg.inv(np.transpose(self.cell_vecs))
            self.position = [mat.dot(v) for v in self.position]

    def guess_molecule(self, site_list, center=None):
        '''.. py:method:: guess_molecule(site_list, center)

        Arrange atom position to form a molecule.

        This method is effective to rotate a molecule.

        Parameters
        -----------

        site_list: list
            list of site number
        center: list
            center position of "molecule" (Optional).

        Returns
        ---------

        numpy.ndarray
            Array of Vector that represents "molecule".

        Notes
        -----
        When the optional argument, center, is set, the atoms are
        re-arranged as to minimize the distance from this center.
        If not, atoms are re-arranged to minimize the total bonding
        length.  As the algorithm for mimizing the total length is
        not exhaustive, the resultant atom arrangement may different
        from what you expect, in spite of time-waste.  The center
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
        '''.. py:method::guess_molecule2(site_lsit)

        Arranges atom positions to form a molecule.

        poscar updates

        Parameters
        -----------

        site_list: list
            list of site number
        '''
        molecule = self.guess_molecule(site_list)
        for site, pos_vector in zip(site_list, molecule):
            self.pos_replace(site, pos_vector)

    def translate(self, vector, atomlist):
        '''.. py:method:: translate(vector, atomlist)

        Translate the selected atom(s) by vector

        Parameters
        ----------

        vector: list, tuple, numpy.array
             translational vector (in Cartesian frame)
        atomlist: list
             list of the atom to be moved


        Notes
        ------

        the first atom is "1", not "0".

        Returns
        --------

        numpy.ndarray
              position
        '''
        if self.is_cartesian():
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

    def get_axes_lengthes(self):
        '''.. py:method:: get_axes_lengthes()

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
        '''.. py:method::translate_all(vector)

        Translate **all** atoms by vector

        Parameters
        -----------

        vector: list, numpy.array
             translational vector
        '''
        atomrange = list(range(1, sum(self.ionnums) + 1))
        self.translate(vector, atomrange)

    def save(self, filename):
        '''.. py:method::samve(filename)

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
    '''.. py:function:: point_in_box(point, cell_vecs)

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


def rotate_x(theta):
    ''' .. py:function:: rotate_x(theta)

    Rotation matrix around X-axis

    Parameters
    ----------

    theta: float
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
         [0.0, np.cos(theta * degree), -np.sin(theta * degree)],
         [0.0, np.sin(theta * degree), np.cos(theta * degree)]])


def rotate_y(theta):
    '''.. py:function:: rotate_y(theta)

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
        [[np.cos(theta * degree), 0.0, np.sin(theta * degree)],
         [0.0, 1.0, 0.0],
         [-np.sin(theta * degree), 0.0, np.cos(theta * degree)]])


def rotate_z(theta):
    '''.. py:function:: rotate_y(theta)

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
        [[np.cos(theta * degree), -np.sin(theta * degree), 0.0],
         [np.sin(theta * degree), np.cos(theta * degree), 0.0],
         [0.0, 0.0, 1.0]])


def three_by_three(vec):
    '''.. py:function:three_by_three(vec)

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
