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

    def __init__(self, chgcar_file=None):
        super(CHGCAR, self).__init__(None)
        self.__meshX = 0
        self.__meshY = 0
        self.__meshZ = 0
        self.__spininfo = 0
        self.__chgArray = []
        if chgcar_file:
            self.load_from_file(chgcar_file)

    def load_from_file(self, chgcarfile):
        '''Parse CHGCAR file to construct CHGCAR object

        :param chgcarfile: CHGCAR file name
        :type chgcarfile: str
        '''
        section = 'poscar'
        separator = None
        tmp = []
        with open(chgcarfile) as thefile:
            for line in thefile:
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
        dest_chgcar = copy.deepcopy(self)
        if len(self.spininfo) == 2:
            total, sd1 = tools.each_slice(self.chgArray,
                                          self.meshX *
                                          self.meshY *
                                          self.meshZ)
            dest_chgcar.__chgArray = list(sd1)
            dest_chgcar.__spininfo = ["up-down"]
        elif len(self.spininfo) == 4:
            total, sd1, sd2, sd3 = tools.each_slice(self.chgArray,
                                                    self.meshX *
                                                    self.meshY *
                                                    self.meshZ)
            if direction is None:
                direction = 'x'
            if direction == 'x':
                dest_chgcar.__chgArray = list(sd1)
                dest_chgcar.__spininfo = ["mX"]
            elif direction == 'y':
                dest_chgcar.__chgArray = list(sd2)
                dest_chgcar.__spininfo = ["mY"]
            elif direction == 'z':
                dest_chgcar.__chgArray = list(sd3)
                dest_chgcar.__spininfo = ["mZ"]
        return dest_chgcar

    def majorityspin(self):
        '''Return CHGCAR for majority spin

        This method is for CHGCAR given by ``ISPIN=2`` but not-SOI
        calculations.

        :return: CHGCAR for the majority spin charge
        :rtype:  CHGCAR

        '''
        if len(self.spininfo) != 2:
            raise RuntimeError('This CHGCAR is not spinresolved version')
        dest_chgcar = copy.deepcopy(self)
        total, magnetization = tools.each_slice(self.chgArray,
                                                self.meshX *
                                                self.meshY *
                                                self.meshZ)
        dest_chgcar.__chgArray = [
            (up + down) / 2 for up, down in zip(total, magnetization)]
        dest_chgcar.__spininfo = ["up"]
        return dest_chgcar

    def minorityspin(self):
        '''Return CHGCAR for minority spin

        This method is for CHGCAR given by ``ISPIN=2`` but not-SOI
        calculations.

        :return: CHGCAR for the minority  spin charge
        :rtype: CHGCAR

        '''
        if len(self.spininfo) != 2:
            raise RuntimeError('This CHGCAR is not spinresolved version')
        dest_chgcar = copy.deepcopy(self)
        total, magnetization = tools.each_slice(self.chgArray,
                                                self.meshX *
                                                self.meshY *
                                                self.meshZ)
        dest_chgcar.__chgArray = [
            (up - down) / 2 for up, down in zip(total, magnetization)]
        dest_chgcar.__spininfo = ["down"]
        return dest_chgcar

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
        add_chgcar = super(CHGCAR, self).__add__(other)
        if any([self.meshX != other.meshX,
                self.meshY != other.meshY,
                self.meshZ != other.meshZ]):
            raise RuntimeError('Mesh sizes are inconsistent')
        augend = self.chgArray
        addend = other.chgArray
        if len(augend) == len(addend):
            add_chgcar.__chgArray = [x + y for x, y in zip(augend, addend)]
        else:
            raise RuntimeError('the mesh sies are different.')
        return add_chgcar

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
        diff_chgcar = copy.deepcopy(self)
        if any([self.meshX != other.meshX,
                self.meshY != other.meshY,
                self.meshZ != other.meshZ]):
            raise RuntimeError('Mesh sizes are incinsistent')
        minuend = self.chgArray
        subtrahend = other.chgArray
        if len(minuend) == len(subtrahend):
            diff_chgcar.__chgArray = [x - y for x, y in
                                      zip(minuend, subtrahend)]
        else:
            raise RuntimeError('the mesh sizes are different.')
        return diff_chgcar

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
            thefile = open(filename, mode='w', newline='\n')
        except TypeError:
            thefile = open(filename, mode='wb')
        with thefile:
            thefile.write(str(self))
