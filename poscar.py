#!/usr/bin/env python
# -*- conding: utf-8 -*-
# python 3.3.2
# translate from poscar.rb of 2014/2/26, master branch

import numpy as np
import itertools as it
import copy, re

# @example
#  l=    " 0.7071067800000000    0.0000000000000000    0.0000000000000000"
#  Vector.elements(l.split.map{|i| i.to_f}) # => Vector[0.70710678, 0.0, 0.0]

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

'''
  0.00000      0.00000      0.00000
  0.00000      1.43703      2.03227
  1.24451      0.71852      4.06454
  0.00000      0.00000      6.09682
  0.00000      1.43703      8.12909
  0.00000      0.00000     10.16136
  0.00000      1.43703     12.19363
  1.24451      0.71852     14.22590
  0.00000      0.00000     16.25818
'''

class POSCAR:
    '''
class for POSCAR (CONTCAR) format
    This script does *NOT* support for constructing POSCAR from scratch.
    (Use ASE for this purpose.)
    It provides a way to slightly modify the POSCAR(CONTCAR),
      which has already works well.
  methods list: to_Cartesian, translate, translate_all
  @version 2.0'''
    def __init__(self, arg=None):
        '''arg is POSCAR file name, or list of POSCAR text.'''
        self.system_name = ""
        self.scaling_factor = 0.
        self.__latticeV1 = np.array([[0., 0., 0.]])
        self.__latticeV2 = np.array([[0., 0., 0.]])
        self.__latticeV3 = np.array([[0., 0., 0.]])
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
        return self.__latticeV1

    @property
    def latticeV2(self):
        return self.__latticeV2

    @property
    def latticeV3(self):
        return self.__latticeV3

    def load_from_array(self, poscar):
        poscar = map(str.rstrip, poscar)
        self.system_name = next(poscar)
        self.scaling_factor = float(next(poscar))
        self.__latticeV1 = np.array([list(map(float, next(poscar).split()))])
        self.__latticeV2 = np.array([list(map(float, next(poscar).split()))])
        self.__latticeV3 = np.array([list(map(float, next(poscar).split()))])
        self.iontype = next(poscar).split()
        # parse POSCAR evenif the element names are not set.
        # At present, the String representation number 
        #   are used for the  dummy name.
        if self.iontype[0].isnum():
            self.ionnums = list(map(int, self.iontype))
            #                      [int(i) for i in self.iontype]
        else:
            self.ionnums = list(map(int, next(poscar).split()))
        ii = 1
        for elm, n in zip(self.iontype, self.ionnums):
            self.__atom_identifer.extend(
                '#{0}:{1}{2}'.format(ii+m, elm, m+1) for m in range(n))
            ii += n
        line7 = next(poscar)
        if re.search(r'^[\s]*Selective\b', line7, re.I): # re.match(r'[\s] ... )
            self.__selective = True
            self.coordinate_type = next(poscar)
        else:
            self.__selective = False
            self.coordinate_type = line7

        for line, elem in zip(poscar, self.atom_identifer):
            # if not elem: break
            tmp = line.split()
            self.position.append(np._float(np.array([tmp[:3]])))
            if self.is_selective:
                self.coordinate_changeflags.append(' '.join(tmp[3:]))

  # attribute self.atom_identifer reader
  # return [list] self.__atom_identifer 
  #   To display self.atom_identifer, the value is calculated,
  #   by using self.iontype and self.ionnums, 
  #   not just showing self.__atom_identifer
    @property
    def atom_identifer(self):
        self.__atom_identifer = []
        ii = 1
        for elm, n in zip(self.iontype, self.ionnums):
            self.__atom_identifer.extend(
                '#{0}:{1}{2}'.format(ii+m, elm, m+1) for m in range(n))
            ii += n
        return self.__atom_identifer

    @atom_identifer.setter
    def atom_identifer(self, value):
        self.__atom_identifer = value

    def __iter__(self):
        for each in self.position:
            yield each

    @property
    def is_cartesian(self):
        '''return True if coordinate is cartesian'''
        return bool(re.search(r'^[ck]', self.coordinate_type, re.I))
    
    @property
    def is_direct(self):
        '''return True if coordinate is direct (not cartesian)'''
        return not self.is_cartesian

    @property
    def is_selective(self):
        '''return True if "Selective Dynamics" switch on.'''
        return self.__selective

    def pos(self, i):
        '''
    POSCAR.pos(i): An accessor of POSCAR.position.
      As in VASP, the atom index starts with "1", not "0".  
      This method follows this manner.
    @param [Fixnum] i site #
     @note the first site # is "1", not "0". (Follow VESTA's way.)
    @return [numpy.array] Returns i-th atom's position by Vector representation.
'''
        return self.position[i - 1]

    def pos_replace(self, i, vector):
        '''
    @param [Fixnum] i site #
     @note the first site # is "1", not "0". 
        (follow VESTA's and VASP's way.)
    @param [Vector] vector  Vector representation of the i-th atom position.
'''
        vector = _vectorize(vector)
        if not isinstance(i, int): raise ValueError
        if not self.is_cartesian:
            message = 'poscar_replace method is implemented for Cartesian coordinate'
            raise RuntimeError(message)
        self.position[i -1] = vector

    def rotateX(self, theta):
        '''
    @param [float] theta angle of rotation (Degrees)
    @return [np.matrix] rotation matrix
'''
        degree = np.pi / 180.0
        return np.mat([[1.0, 0.0, 0.0],
                       [0.0, np.cos(theta * degree), -np.sin(theta.degree)],
                       [0.0, np.sin(theta * degree),  np.cos(theta.degree)]])

    def rotateY(self, theta):
        'see POSCAR.rotateX()'
        degree = np.pi / 180.0
        return np.mat([[ np.cos(theta * degree), 0.0, np.sin(theta.degree)],
                       [0.0, 1.0, 0.0],
                       [-np.sin(theta * degree), 0.0, np.cos(theta.degree)]])

    def rotateZ(self, theta):
        'see POSCAR.rotateX()'
        degree = np.pi / 180.0
        return np.mat([[np.cos(theta * degree), -np.sin(theta.degree), 0.0],
                       [np.sin(theta * degree),  np.cos(theta.degree), 0.0]
                       [0.0, 0.0, 1.0],])

    # class method? or independent function?
    def nearest(self, array, point):
        '''
    @param [list-of-np.array] array
    @param [np.array] point
    @return [np.array]
'''
        return min(array, key=lambda pos: np.linalg.norm(pos - point))

    # class method? or independent function?
    def make27candidate(self, position):
        '''
  # @param [np.array] position atom position defined
  #    in the coordinated by latticeV1, latticeV2, 
  #    latticeV3 ( scaling facter is not accounted).
  # @return [list-of-np.array]
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
                candidates27.append(l * np.array([[1., 0., 0.]]) +
                                    m * np.array([[0., 1., 0.]]) +
                                    n * np.array([[0., 0., 1.]]) +
                                    position)
        return candidates27

    def atom_rotate(self, site, axis_name, theta, center):
        '''
    Rotate atom under periodic boundary condition
    @param [Fixnum] site site # for rotation (The first atom is "1".).
    @param [String] axis_name "X", "x", "Y", "y", "Z", or "z". Rotation axis.
    @param [Float] theta Rotation angle (Degrees).
    @param [Array<Float>, Vector] center center position for rotation. 
    @todo  check the center in the Braves lattice.
    @todo  take into account the periodic boundary. 
'''
        center = _vectorize(center)
        if len(center[0]) != 3: raise ValueError
        if not self.point_in_box(center / self.scaling_factor,
                                 self.latticeV1, self.latticeV2, self.latticeV3):
            raise ValueError('the center must be in the Braves lattice')
        if not isinstance(site, int):
            raise ValueError('argument error in atom_rotate method')
        if not is_cartesian: self.to_Cartesian()
        position = self.pos(site)
        position -= center / self.scaling_factor
        position = (getattr(self, 'rotate' + axis_name.capitalize())(theta) * position.T).T
        position += center / self.scaling_factor
        self.pos_replace(site, positon)

    def atoms_rotate(self, site_list_pack, axis_name, theta, center):
        '''
  # @param [Array<Fixnum>] site_list_pack list array of the 
  #   list array (not typo!) of  site # for rotation 
  #   (The first atom is "1".).
  # @param [String] axis_name "X", "x", "Y", "y", "Z",
  #   or "z". Rotation axis.
  # @param [Float] theta Rotation angle (Degrees).
  # @param [Array<Float>, Vector] center center position for rotation. 
'''
        for site_list in site_list_pack:
            for site in site_list:
                self.atom_rotate(site, axis_name, theta, center)

    def __add__(self, other):
        '''
  # @param [POSCAR] other
  # @return [POSCAR] 
  # @todo Check the lattice vectors, coordinate_type and so on.
'''
        #if not isinstance(other, POSCAR): return NotImplemented
        destPOSCAR = copy.deepcopy(self)
        if destPOSCAR.scaling_factor != other.scaling_factor:
            raise ValueError('scaling factor is different.')
        if (destPOSCAR.latticeV1 +
            destPOSCAR.latticeV2 +
            destPOSCAR.latticeV3 -
            (other.latticeV1 +
             other.latticeV2 +
             other.latticeV3)) != np.array([[0., 0., 0.]]):
            raise ValueError('lattice vectors are different.')
        destPOSCAR.iontype.extend(other.iontype)
        destPOSCAR.ionnums.extend(other.ionnums)
        destPOSCAR.position.extend(other.position)
        destPOSCAR.coordinate_changeflags.extend(other.coordinate_changeflags)

        destPOSCAR.atom_identifer

        return destPOSCAR

    def to_list(self):
        '''
    @return [list] a list representation of POSCAR.
'''
        out_list = []
        out_list.append(self.system_name)
        out_list.append(self.scaling_factor)
        out_list.append(self.latticeV1)
        out_list.append(self.latticeV2)
        out_list.append(self.latticeV3)
        if not self.iontype[0].isnum():
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
    @return [String] a string representation of POSCAR
'''
        tmp = []
        tmp.append(self.system_name)
        tmp.append(str(self.scaling_factor))
        tmp.append(''.join('   {0:20.17f}'.format(i) for i in self.latticeV1[0]))
        tmp.append(''.join('   {0:20.17f}'.format(i) for i in self.latticeV2[0]))
        tmp.append(''.join('   {0:20.17f}'.format(i) for i in self.latticeV3[0]))
        if not self.iontype[0].isnum():
            tmp.append(' ' + ' '.join(self.iontype))
        tmp.append(' ' + ' '.join(self.ionnums))
        if self.is_selective:
            tmp.append('Selective Dynamics')
        tmp.append(self.coordinate_type)
        for pos, tf, atom in it.zip_longest(self.position,
                                            self.coordinate_changeflags,
                                            self.atom_identifer, fillvalue=''):
            tmp.append(' '.join('  {0:20.17f}'.format(i) for i in pos) +
                       ' ' + tf + ' ' + atom)
        return '\n'.join(tmp) + '\n'

    def tune_scaling_factor(self, new_scaling_factor=1.0):
        '''
    change scaling factor to new value.
        The Braves lattice are corrected (to be equal size).
    @param [Float] new_scaling_factor
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
        '''
    change the coordinate to cartesian from direct.
    @return [Boolian] true if POSCAR file is cartesian coordinate
'''
        if self.is_direct:
            self.coordinate_type = "Cartesian"
            m = np.mat([self.latticeV1, self.latticeV2, self.latticeV3]).T
            self.position = [(m * v.T).T for v in self.position]

    def to_Direct(self):
       '''
    change the coordinate to direct from cartesian.
'''
       if self.is_cartesian:
           self.coordinate_type = "Direct"
           m = (np.mat([self.latticeV1, self.latticeV2, self.latticeV3]).T).I
           self.position = [(m * v.T).T for v in self.position]

    def guess_molecule(self, site_list, center=None):
        '''
  # arrange atom position to form a molecule.  This method is effective for molecular rotation.
  # @param [Array<Fixnum>] site_list list of site #
  # @param [Vector] center center position of "molecule" (Optional).
  # @note When the optional argument, center, is set, the 
  #   atoms are re-arranged as to
  #   minimize the distance from this center.  If not the
  #   center is not set, atoms are re-arranged to minimize 
  #   the total bonding length.  As the algorithm for mimizing the
  #   total length is not exhaustive, the resultant atom arrangement
  #   may different from what you expect, in spite of time-waste.
  #   So, the center option is highly recommended to form a molecule.
  # @return [Array<Vector>] Array of Vector that represents "molecule".
'''
        #list for atom positions for "molecule"
        molecule = [self.pos(j) for j in site_list]
        for index, site in enumerate(site_list):
            target_atom = self.pos(site)
            atoms27 = self.make27candidate(target_atom)
            def func(pos):
                molecule[index] = pos
                if center is not None: # bool([np.ndarray]) => Error
                    center = _vectorize(center)
                    return np.linalg.norm(pos - center)
                else:
                    # fixme!! when the highest symmetry point can be detemined from the position list,
                    #              guess_molecule method does not require the "center" option.
                    #(molecule.product(molecule)).inject(0.0) do | s, vectors |
                    # s+ (vectors[0]-vectors[1]).magnitude
                    s = 0.0
                    for vectors in it.product(molecule, molecule):
                        s += np.linalg.norm(vectors[0] - vectors[1])
                    return s
            newpos = min(atoms27, key=func)
            molecule[index] = newpos
        return molecule

    def guess_molecule2(self, site_list):
        '''
  # arrange atom position to form a molecule.
  # poscar updates
  # @param [Array<Fixnum>] site_list list of site #
'''
        molecule = guess_molecule(site_list)
        for site, posVector in zip(site_list, molecule):
            self.pos_replace(site, posVector)

    def translate(self, vector, atomlist):
        '''  # Translate the selected atom(s) by vector
  # @param [Vector] vector translational vector
  # @param [Array] atomlist list of the atome for moving 
  #  @note the first atom is "1", not "0".
  # @return [Array]
'''
        if self.is_cartesian:
            vector = _vectorize(vector)
            for i in atomlist:
                self.position[i - 1] = self.position[i - 1] + vector / self.scaling_factor
        else:
            ###
            print('line is not correct')
        return self.postion

    def translate_all(self, vector):
        '''  # Translate all atoms by vector
  # @param [Vector] vector translational vector
  # @return [Array]
'''
        atomrange = list(range(1, sum(self.ionnums) + 1))
        translate(vector, atomrange)

    def save(self, filename):
        with open(filename, 'w') as f:
            f.write(str(self))

    def point_in_box(self, point, l1, l2, l3):
        x0, y0, z0 = _vectorize(point).flatten()
        x1, y1, z1 = _vectorize(l1).flatten()
        x2, y2, z2 = _vectorize(l2).flatten()
        x3, y3, z3 = _vectorize(l3).flatten()
        l = -(-x3*y2*z0+x2*y3*z0+x3*y0*z2-x0*y3*z2-x2*y0*z3+x0*y2*z3)/(x3*y2*z1-x2*y3*z1-x3*y1*z2+x1*y3*z2+x2*y1*z3-x1*y2*z3)
        m = -(x3*y1*z0-x1*y3*z0-x3*y0*z1+x0*y3*z1+x1*y0*z3-x0*y1*z3)/(x3*y2*z1-x2*y3*z1-x3*y1*z2+x1*y3*z2+x2*y1*z3-x1*y2*z3)
        n = -(x2*y1*z0-x1*y2*z0-x2*y0*z1+x0*y2*z1+x1*y0*z2-x0*y1*z2)/(-x3*y2*z1+x2*y3*z1+x3*y1*z2-x1*y3*z2-x2*y1*z3+x1*y2*z3)
        return all((0 <= q <= 1) for q in (l, m, n))

#--------------------------
def _vectorize(vector):
    if not isinstance(vector, _vector_acceptables):
        raise TypeError('Cannot convert into vector.')
    return np.array([np.array(vector).flatten()])

_vector_acceptables = (np.ndarray, np.matrix, list, tuple)
#--------------------------
#

if __name__ == '__main__':
    # $-w = true
    import argparse, os, sys
    import functools as ft
    mypath = os.readlink(__file__) if os.path.islink(__file__) else __file__
    sys.path.append(os.path.dirname(os.path.abspath(mypath)))
    import argParsePostVASP as appv
    arg = appv.APPV()
    arg.add_argument('--atom', metavar='atoms', action='append',
                     type=arg.parse_AtomselectionNum,
                     help='''atoms specified with range using "-"
or comma-delimnated numbers.
 (ex.) --atom 1,2,7-9''')
    def split_to_float(string, n, name):
        lis = string.split(',')
        if len(lis) != n:
            message = '--{0} option requires {1} numbers'.format(name, n)
            raise argparse.ArgumenTypeError(message)
        return [float(i) for i in lis]
    arg.add_argument('--translate', metavar='x,y,z', action='append',
                     type=ft.partial(split_to_float, n=3, name='translate'),
                     help='displacement (AA unit) by three numbers separated by comma.')
    arg.add_argument('--rotateX', metavar='theta,x,y,z',
                     type=ft.partial(split_to_float, n=4, name='rotateX'),
                     help='''Rotation around X-axis by "theta" at (x,y,z)
NOTE: this option is not taken into account the periodic boundary.''')
    arg.add_argument('--rotateY', metavar='theta,x,y,z',
                     type=ft.partial(split_to_float, n=4, name='rotateY'),
                     help='''Rotation around Y-axis by "theta" at (x,y,z)
NOTE: this option is not taken into account the periodic boundary.''')
    arg.add_argument('--rotateZ', metavar='theta,x,y,z',
                     type=ft.partial(split_to_float, n=4, name='rotateZ'),
                     help='''Rotation around Z-axis by "theta" at (x,y,z)
NOTE: this option is not taken into account the periodic boundary.''')
    arg.add_argument('--output', metavar='file_name',
                     help='''output file name
if not specified, use standard output''')
    arg.add_argument('poscar', metavar='POSCAR_file (or CONTCAR_file)',
                     type=POSCAR)
    arguments = arg.parse_args()
    # translate option and rotate option are not set simulaneously.
    if arguments.translate and any([arguments.rotateX, arguments.rotateY, arguments.rotateZ]):
        raise RuntimeError("Cannot set --translate and rotate option simultanaously.")
    # rotate options are not set multiply.
    if (arguments.rotateX, arguments.rotateY, arguments.rotateZ).count(None) < 2:
        raise RuntimeError("Cannot set multiple rotate options simultanaously.")
    #
    
    ############

    print(arguments.poscar)
    arguments.poscar.to_Cartesian()
    print(arguments.poscar)

    #
    #  if "atom" option is not set, all atoms are concerned. 
    #
    if not arguments.atom:
        nAtoms = sum(arguments.poscar.ionnums)
        arguments.atom = [arg.parse_AtomselectionNum('1-{0}'.format(nAtoms))]
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
    if arguments.output:
        arguments.poscar.save(arguments.output)
    else:
        print(arguments.poscar)
