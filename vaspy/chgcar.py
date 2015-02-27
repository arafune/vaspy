#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Module for CHGCAR class

translate from chgcar.rb in scRipt4VASP, 2014/2/26 master branch
'''
from __future__ import division, print_function  # Version safety
import re
import copy
import os
import sys

try:
    from vaspy import poscar, tools
except ImportError:
    mypath = os.readlink(__file__) if os.path.islink(__file__) else __file__
    sys.path.append(os.path.dirname(os.path.abspath(mypath)))
    import poscar
    import tools


_re_blank = re.compile(r'^[\s]*$')
_re_aug_occ = re.compile(r'\baugmentation occupancies')


class CHGCAR(poscar.POSCAR):

    '''Class for CHGCAR

     An example of the first few lines of CHGCAR. ::

           hBN-Cu                                  #1st line   poscar.POSCAR[0]
           1.00000000000000                        #2nd line   poscar.POSCAR[1]
             6.762964    0.000000    0.000000      #3rd line   poscar.POSCAR[2]
             3.381482    5.856898    0.000000      #4th line   poscar.POSCAR[3]
             0.000000    0.000000   29.004836      #5th line   poscar.POSCAR[4]
           B    Cu   N    Si                       #6th line   poscar.POSCAR[5]
             7    21     7     6                   #7th line   poscar.POSCAR[6]
           Direct                                  #8th line   poscar.POSCAR[7]
             0.047680  0.261795  0.361962          #9th line   poscar.POSCAR[8]
             ....

    :attribute: chgArray, meshX, meshY, meshZ, spininfo
    :version: 1.0.0

    .. note:: the current verstion ignores
              "augmentation occupacies".
'''
    # accessor: chgArray, meshX-Y-Z

    def __init__(self, arg=None):
        super(CHGCAR, self).__init__(None)
        self.__meshX = 0
        self.__meshY = 0
        self.__meshZ = 0
        self.__spininfo = 0
        self.__chgArray = []
        if arg:
            self.load_from_file(arg)

    def load_from_file(self, chgcarfile):
        '''Parse CHGCAR file to construct CHGCAR object

        :param chgcarfile: CHGCAR file name
        :type chgcarfile: str
        '''
        section = 'poscar'
        separator = None
        tmp = []
        with open(chgcarfile) as f:
            for line in f:
                line = line.rstrip('\n')
                if section == 'poscar':
                    if re.search(_re_blank, line):
                        self.load_from_array(tmp)
                        section = 'define_separator'
                    else:
                        tmp.append(line)
                elif section == 'define_separator':
                    separator = line if separator is None else separator
                    if self.meshX == self.meshY == self.meshZ == 0:
                        self.__meshX, self.__meshY, self.__meshZ = \
                            list(map(int, line.split()))
                    section = 'grid'
                elif section == 'aug':
                    if separator in line:
                        section = 'grid'
                    elif "augmentation occupancies " in line:
                        pass  # not implemented
                    else:
                        pass  # not implemented
                elif section == 'grid':
                    if "augmentation occupancies " in line:
                        section = 'aug'
                    elif separator in line:
                        pass
                    else:
                        self.__chgArray.extend(map(float, line.split()))
        if len(self.chgArray) % (self.meshX * self.meshY * self.meshZ) != 0:
            print(len(self.chgArray), ": Should be ",
                  self.meshX * self.meshY * self.meshZ, "x 1, 2, 4")
            print(len(self.chgArray) % (self.meshX * self.meshY * self.meshZ),
                  ": Should be zero")
            raise RuntimeError('Failed: Construction')
        self.__spininfo = len(self.chgArray) // (self.meshX *
                                                 self.meshY *
                                                 self.meshZ)
        if self.spininfo == 1:
            self.__spininfo = [""]
        elif self.spininfo == 2:
            self.__spininfo = ["up+down", "up-down"]
        elif self.spininfo == 4:
            self.__spininfo = ["mT", "mX", "mY", "mZ"]
        else:
            raise RuntimeError("CHGCAR is correct?")

    @property
    def meshX(self):
        '''Number of mesh along the first axis of the cell'''
        return self.__meshX

    @property
    def meshY(self):
        '''Number of mesh along the second axis of the cell'''
        return self.__meshY

    @property
    def meshZ(self):
        '''Number of mesh along the third axis of the cell'''
        return self.__meshZ

    @property
    def spininfo(self):
        '''Spin property of CHGCAR

        *  for ``ISPIN = 1``, [""]

        *  for ``ISPIN = 2`` (but ``LSORBIT=.FALSE.``),
           ["up+down", "up-down"]

        *  for ``ISPIN = 2`` (but ``LSORBIT=.TRUE.``),
           ["mT", "mX", "mY", "mZ"]

        '''
        return self.__spininfo

    @property
    def chgArray(self):
        '''charge data'''
        return self.__chgArray

    def magnetization(self, direction=None):
        '''
        Return CHGCAR for magnetization

        For collinear spin-polarized calculations
        (``ISPIN=2`` but ``LSORBIT=.FALSE.``),
        two sets of data are found in CHGCAR file. The first set
        is the total charge density (spin-up plus spin-down),
        the second one the magnetization density (spin-up minus spin-down).
        For non-collinear spin-polarized calculations
        (``ISPIN=2`` and ``LSORBIT=.TRUE.``),
        CHGCAR file stores the total charge density and the
        magnetisation density in the x, y and z direction in this order.

        For collinear spinpolarized calculation the argument does
        not make a sense.  For non-collinear CHGCAR, direction
        should be one of 'x', 'y', and 'z'

        :param direction: specify x, y, or z in noncollinear calculation
        :type direction: str
        :return: CHGCAR of the spin-distribution
        :rtype: CHGCAR
'''
        if len(self.spininfo) == 1:
            raise RuntimeError("This CHGCAR is not spinresolved version")
        destCHGCAR = copy.deepcopy(self)
        if len(self.spininfo) == 2:
            s1, s2 = tools.each_slice(self.chgArray,
                                      self.meshX * self.meshY * self.meshZ)
            destCHGCAR.__chgArray = list(s2)
            destCHGCAR.__spininfo = ["up-down"]
        elif len(self.spininfo) == 4:
            s1, s2, s3, s4 = tools.each_slice(self.chgArray,
                                              self.meshX * self.meshY * self.meshZ)
            if direction is None:
                direction = 'x'
            if direction == 'x':
                destCHGCAR.__chgArray = list(s2)
                destCHGCAR.__spininfo = ["mX"]
            elif direction == 'y':
                destCHGCAR.__chgArray = list(s3)
                destCHGCAR.__spininfo = ["mY"]
            elif direction == 'z':
                destCHGCAR.__chgArray = list(s4)
                destCHGCAR.__spininfo = ["mZ"]
        return destCHGCAR

    def majorityspin(self):
        '''Return CHGCAR for majority spin

        This method is for CHGCAR given by ``ISPIN=2`` but not-SOI
        calculations.

        :return: CHGCAR for the majority spin charge
        :rtype:  CHGCAR

        '''
        if len(self.spininfo) != 2:
            raise RuntimeError('This CHGCAR is not spinresolved version')
        destCHGCAR = copy.deepcopy(self)
        s1, s2 = tools.each_slice(self.chgArray,
                                  self.meshX * self.meshY * self.meshZ)
        destCHGCAR.__chgArray = [(up + down) / 2 for up, down in zip(s1, s2)]
        destCHGCAR.__spininfo = ["up"]
        return destCHGCAR

    def minorityspin(self):
        '''Return CHGCAR for minority spin

        This method is for CHGCAR given by ``ISPIN=2`` but not-SOI
        calculations.

        :return: CHGCAR for the minority  spin charge
        :rtype: CHGCAR

        '''
        if len(self.spininfo) != 2:
            raise RuntimeError('This CHGCAR is not spinresolved version')
        destCHGCAR = copy.deepcopy(self)
        s1, s2 = tools.each_slice(self.chgArray,
                                  self.meshX * self.meshY * self.meshZ)
        destCHGCAR.__chgArray = [(up - down) / 2 for up, down in zip(s1, s2)]
        destCHGCAR.__spininfo = ["down"]
        return destCHGCAR

    def __add__(self, other):
        '''
        x.__add__(y) <=> x + y

        :param CHGCAR: other addition CHGCAR object
        :return: CHGCAR of the result by summing two CHGCARs:
        :rtype: CHGCAR
        :note:
        in the returned CHGCAR :
        the charge distribution is just summantion of two CHGCARs,
        and the atoms are also summantion of two CHGCARs.
'''
        # augend + aggend
        if not isinstance(other, CHGCAR):
            return NotImplemented
        addCHGCAR = super(CHGCAR, self).__add__(other)
        if any([self.meshX != other.meshX,
                self.meshY != other.meshY,
                self.meshZ != other.meshZ]):
            raise RuntimeError('Mesh sizes are inconsistent')
        augend = self.chgArray
        addend = other.chgArray
        if len(augend) == len(addend):
            addCHGCAR.__chgArray = [x + y for x, y in zip(augend, addend)]
        else:
            raise RuntimeError('the mesh sies are different.')
        return addCHGCAR

    def __sub__(self, other):
        '''x.__sub__y <=> x - y

        :param other: difference CHGCAR object
        :type other: CHGCAR
        :return: CHGCAR of the result of difference between two CHGCARs:
        :rtype: CHGCAR
        :note:  in the returned CHGCAR :
        the charge distribution is just difference
        of two CHGCARs, and the atoms are used
        for "munuend" CHGCAR, not difference.
        The atoms in subtrahend CHGCAR are totally ignored.
'''
        # minuend - subtrahend
        if not isinstance(other, CHGCAR):
            return NotImplemented
        diffCHGCAR = copy.deepcopy(self)
        if any([self.meshX != other.meshX,
                self.meshY != other.meshY,
                self.meshZ != other.meshZ]):
            raise RuntimeError('Mesh sizes are incinsistent')
        minuend = self.chgArray
        subtrahend = other.chgArray
        if len(minuend) == len(subtrahend):
            diffCHGCAR.__chgArray = [x - y for x, y in
                                     zip(minuend, subtrahend)]
        else:
            raise RuntimeError('the mesh sies are different.')
        return diffCHGCAR

    def __str__(self):
        '''x.__str__() <=> str(x)

        :return: a string representation of CHGCAR.
        :rtype: str
'''
        outputstring = ''
        tmp = self.chgArray
        for tmp in tools.each_slice(self.chgArray,
                                    self.meshX * self.meshY * self.meshZ):
            output = []
            outputstring += '\n  {0}  {1}  {2}\n'.format(self.meshX,
                                                         self.meshY,
                                                         self.meshZ)
            for array in tools.each_slice(tmp, 5):
                output.append(''.join('  {0:18.11E}'.format(i)
                                      for i in array if i is not None))
            outputstring += '\n'.join(output)
        return super(CHGCAR, self).__str__() + outputstring + '\n'

    def save(self, filename):
        '''Save CHGCAR object as CHGCAR file style

        :param filename: file name
        :type filename: str
'''
        try:  # Version safety
            file = open(filename, mode='w', newline='\n')
        except TypeError:
            file = open(filename, mode='wb')
        with file:
            file.write(str(self))

# ------------------------- Main
if __name__ == '__main__':
    import argparse
    arg = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter)
    group = arg.add_mutually_exclusive_group(required=True)
    group.add_argument('--add', action='store_true', default=False,
                       help="Add two CHGCAR files")
    group.add_argument('--diff', action='store_true', default=False,
                       help="Get difference of two CHGCAR files")
    group.add_argument('--spin', metavar='spin_operation',
                       help="""spin-relatated operation.
when this option is set --add, -diff are ignored,
and CHGCAR_file_2 must not be set.
spin operation is one of the followings:
mag : show the magnetisation
      density (for spin resolved calculations)
magX : show the magnetisation density in
       the X direction (for non collinear calc.)
magY : show the magnetisation density in
       the Y direction (for non collinear calc.)
magZ : show the magnetisation density in
       the Z direction (for non collinear calc.)
majority : extract the part for the
           majority spin (for spin resolved calc.)
minority : extract the part for the
           inority spin (for spin resolved calc.)""")
    arg.add_argument('--output', metavar='file_name',
                     help="""output file name
if not specified, use standard output""")
    arg.add_argument('CHGCAR_file_1', type=CHGCAR)
    arg.add_argument('CHGCAR_file_2', type=CHGCAR, nargs='?')
    # if CHGCAR_file_2 is not specified,
    # *None* is stored in arguments.CHGCAR_file_2, not CHGCAR(None)
    arguments = arg.parse_args()
    #
    if arguments.spin is not None:
        if arguments.CHGCAR_file_2 is not None:
            raise RuntimeError("Only one CHGCAR file for --spin operations")
        if arguments.spin == "mag":
            c = arguments.CHGCAR_file_1.magnetization()
        elif arguments.spin == "magX":
            c = arguments.CHGCAR_file_1.magnetization('x')
        elif arguments.spin == "magY":
            c = arguments.CHGCAR_file_1.magnetization('y')
        elif arguments.spin == "magZ":
            c = arguments.CHGCAR_file_1.magnetization('z')
        elif arguments.spin == "majority":
            c = arguments.CHGCAR_file_1.majorityspin()
        elif arguments.spin == "minority":
            c = arguments.CHGCAR_file_1.minorityspin()
        else:
            raise RuntimeError("Such spin operation parameter is not defined.")
    #
    if arguments.add or arguments.diff:
        if arguments.CHGCAR_file_2 is None:
            raise RuntimeError('Two CHGCAR files are required.')
        if arguments.add:
            c = arguments.CHGCAR_file_1 + arguments.CHGCAR_file_2
        else:
            c = arguments.CHGCAR_file_1 - arguments.CHGCAR_file_2
    #
    if arguments.output is not None:
        c.save(arguments.output)
    else:
        print(c)
