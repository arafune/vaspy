#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
# translate from chgcar.rb in scRipt4VASP, 2014/2/26 master branch

import re, copy, os, sys
import itertools as it
mypath = os.readlink(__file__) if os.path.islink(__file__) else __file__
sys.path.append(os.path.dirname(os.path.abspath(mypath)))
import poscar


_re_blank = re.compile(r'^[\s]*$')
_re_aug_occ = re.compile(r'\baugmentation occupancies')

def _each_slice(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return it.zip_longest(*args, fillvalue=fillvalue)

def _removeall(L, value):
    'remove all *value* in [list] L'
    while L.count(value):
        L.remove(value)
    return L

class CHGCAR(poscar.POSCAR):
    '''
    class for CHGCAR format
    @version 1.0.0
    @note the current verstion does not take account "augmentation occupacies".

     An example of the first few line of the CHGCAR. :
      hBN-Cu                                  #1st line   @poscar[0]
      1.00000000000000                        #2nd line   @poscar[1]
        6.762964    0.000000    0.000000      #3rd line   @poscar[2]
        3.381482    5.856898    0.000000      #4th line   @poscar[3]
        0.000000    0.000000   29.004836      #5th line   @poscar[4]
      B    Cu   N    Si                       #6th line   @poscar[5]
        7    21     7     6                   #7th line   @poscar[6]
      Direct                                  #8th line   @poscar[7]
        0.047680  0.261795  0.361962          #9th line   @poscar[8]
        ....
'''
    # accessor: chgArray, meshX-Y-Z
    def __init__(self, arg=None):
        poscar.POSCAR.__init__(self, None)
        self.__meshX = 0
        self.__meshY = 0
        self.__meshZ = 0
        self.__spininfo = 0
        self.__chgArray = []
        if arg:
            self.load_from_file(arg)

    def load_from_file(chgcarfile):
        '''
    @param [String] chgcarfile CHGCAR file name
    @return [CHGCAR]
'''
        section = 'poscar'
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
                    separator = separator if 'separator' in locals() else line
                    if self.meshX == self.meshY == self.meshZ == 0:
                        self.__meshX, self.__meshY, self.__meshZ = list(map(int, line.split))
                    section = 'grid'
                elif section == 'aug':
                    if re.search(r'#{0}'.format(separator), line):
                        section = 'grid'
                    elif re.search(_re_aug_occ, line):
                        pass # not implemented
                    else:
                        pass # not implemented
                elif section == 'grid':
                    if re.search(_re_aug_occ, line):
                        section = 'aug'
                    else:
                        self.__chgArray.extend(map(float, line.split()))
        if len(self.chgArray) % (self.meshX * self.meshY * self.meshZ) != 0:
            raise RuntimeError('Failed: Construction')
        self.__spininfo = int(len(self.chgArray) / (self.meshX * self.meshY * self.meshZ))
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
        return self.__meshX

    @property
    def meshY(self):
        return self.__meshY

    @property
    def meshZ(self):
        return self.__meshZ

    @property
    def spininfo(self):
        return self.__spininfo

    @property
    def chgArray(self):
        return self.__chgeArray

    def magnetization(self, direction=None):
        '''
  # CHGCAR#magnetization(direction=nil)
  # @return [CHGCAR] returns CHGCAR of the spin-distribution 
  #
  #   For spinpolarized calculations, two sets of data can be found in
  #   the CHGCAR file. The first set contains the total charge density
  #   (spin up plus spin down), the second one the magnetization
  #   density (spin up minus spin down). For non collinear
  #   calculations the CHGCAR file contains the total charge density
  #   and the magnetisation density in the x, y and z direction in
  #   this order.

  # for spinpolarized calculation the argument does not make a sense.
  # For non collinear CHGCAR direction should be one of 'x', 'y' 'z'
'''
        if len(self.spininfo) == 1:
            raise RuntimeError("This CHGCAR is not spinresolved version")
        destCHGCAR = copy.deepcopy(self)
        s1, s2, s3, s4 = _each_slice(self.chgArray, self.meshX * self.meshY * self.meshZ)
        if len(self.spininfo) == 2:
            desfCHGCAR.__chgArray = list(s2)
            destCHGCAR.__spininfo = ["up-down"]
        elif len(self.spininfo) == 4:
            if direction is None: direction = 'x'
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
        '''
  # CHGCAR#majorityspin
  # @return [CHGCAR] returns CHGCAR for the majority spin charge
  #   from CHGCAR given by ISPIN=2 but not-SOI calculations.
  #   According to Dr. Minamitani, the former part of charge 
  #   distribution corresponds for majority spin + minority spin,
  #   the latter part for  majority - minority
'''
        if len(self.spininfo) != 2:
            raise RuntimeError('This CHGCAR is not spinresolved version')
        destCHGCAR = copy.deepcopy(self)
        s1, s2 = _each_slice(self.meshX * self.meshY * self.meshZ)
        destCHGCAR.__chgArray = [(up + down) / 2. for up, down in zip(s1, s2)]
        destCHGCAR.__spininfo = ["up"]
        return destCHGCAR

    def minorityspin(self):
        '''  # CHGCAR#minorityspin
  # @return [CHGCAR] returns CHGCAR for the minority  spin charge
  #   from CHGCAR given by ISPIN=2 but not-SOI calculations.
  #   According to Dr. Minamitani, the former part of charge 
  #   distribution corresponds for majority spin + minority spin,
  #   the latter part for  majority - minority
'''
        if len(self.spininfo) != 2:
            raise RuntimeError('This CHGCAR is not spinresolved version')
        destCHGCAR = copy.deepcopy(self)
        s1, s2 = _each_slice(self.meshX * self.meshY * self.meshZ)
        destCHGCAR.__chgArray = [(up - down) / 2. for up, down in zip(s1, s2)]
        destCHGCAR.__spininfo = ["down"]
        return destCHGCAR

    def __add__(self, other):
        '''x.__add__(y) <=> x + y
    @param [CHGCAR] other addition CHGCAR object
    @return [CHGCAR] returns CHGCAR of the result by summing two CHGCARs:
    @note
      in the returned CHGCAR :
      the charge distribution is just summantion of two CHGCARs,
      and the atoms are also summantion of two CHGCARs.
'''
        # augend + aggend
        if not isinstance(other, CHGCAR):
            return NotImplemented
        addCHGCAR = poscar.POSCAR.__add__(self, other)
        if any([self.meshX != other.meshX, self.meshY != other.meshY, self.meshZ != other.meshZ]):
            raise RuntimeError('Mesh sizes are incinsistent')
        augend = self.chgArray
        addend = other.chgArray
        if len(augend) == len(addend):
            addCHGCAR.__chgArray = [x + y for x, y in zip(augend, addend)]
        else:
            raise RuntimeError('the mesh sies are different.')
        return addCHGCAR

    def __sub__(self, other):
        '''x.__sub__y <=> x - y
    @param [CHGCAR] other difference CHGCAR object
    @return [CHGCAR] returns CHGCAR of the result 
      of difference between two CHGCARs:
    @note
      in the returned CHGCAR :
      the charge distribution is just difference 
      of two CHGCARs, and the atoms are used 
      for "munuend" CHGCAR, not difference.  
      The atoms in subtrahend CHGCAR are totally ignored.
'''
        # minuend - subtrahend
        if not isinstance(other, CHGCAR):
            return NotImplemented
        diffCHGCAR = copy.deepcopy(self)
        if any([self.meshX != other.meshX, self.meshY != other.meshY, self.meshZ != other.meshZ]):
            raise RuntimeError('Mesh sizes are incinsistent')
        minuend = self.chgArray
        subtrahend = other.chgArray
        if len(minuend) == len(subtrahend):
            diffCHGCAR.__chgArray = [x - y for x, y in zip(minuend, subtrahend)]
        else:
            raise RuntimeError('the mesh sies are different.')
        return diffCHGCAR

    def __str__(self):
        '''x.__str__() <=> str(x)
  # @return [String] return a string representation of CHGCAR.
'''
        outputstring = '\n'
        outputstring += '  {0}  {1}  {2}\n'.format(self.meshX, self.meshY, self.meshZ)
       #outputstring += '  ' + str(self.meshX) + '  ' + self.meshY + '  ' + self.meshZ + '\n'
        tmp = self.chgArray
        output = []
        for array in _each_slice(tmp, 5):
            array = _removeall(array, None)
            output.append(''.join('  {0:E}'.format(i) for i in array))
        outputstring += '\n'.join(output)
        return poscar.POSCAR.__str__(self) + outputstring + '\n'

    def save(self, filename):
        '''
    @param [String] filename
    save operated CHGCAR to the file.
'''
        with open(filename, mode='w') as file:
            file.write(str(self))

# ------------------------- Main
if __name__ == '__main__':
    import argparse
    arg = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
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
    # if CHGCAR_file_2 is not specified, *None* is stored in arguments.CHGCAR_file_2, not CHGCAR(None)
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
            c = arguments.CHGCAR_file_1 + argumetns.CHGCAR_file_2
        else:
            c = arguments.CHGCAR_file_1 - argumetns.CHGCAR_file_2
    # 
    if arguments.output is not None:
        c.save(arguments.output)
    else:
        print(c)

