#!/usr/bin/env/ python
# -*- coding: utf-8 -*-

import re, copy
import itertools as it
import poscar

re_blank = re.compile(r'^[\s]*$')

def _each_slice(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

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
    # accessor: chargeArray, meshX-Y-Z
    def __init__(self, arg=None):
        poscar.POSCAR.__init__(self, None)
        self.__meshX = 0
        self.__meshY = 0
        self.__meshZ = 0
        self.__chargeArray = []
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
                if section == 'poscar'
                    if re.search(re_blank, line):
                        self.load_from_array(tmp)
                        section = 'grid'
                    else:
                        tmp.append(line)
                elif section is 'grid':
                    if self.meshX == 0 and self.meshY == 0 and self.meshZ == 0:
                        self.__meshX, self.__meshY, self.__meshZ = list(
                            map(int, line.split))
                    elif len(self.chargeArray) < self.meshX * self.meshY * self.meshZ:
                        self.__chargeArray.extend(map(float, line.split()))
                    else:
                        section = 'aocc'
        if len(self.chargeArray) != self.meshX * self.meshY * self.meshZ:
            raise RuntimeError('Failed: Construction')

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
    def chargeArray(self):
        return self.__chargeArray

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
        augend = self.chargeArray
        addend = other.chargeArray
        if len(augend) == len(addend):
            addCHGCAR.__chargeArray = [x + y for x, y in zip(augend, addend)]
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
        minuend = self.chargeArray
        subtrahend = other.chargeArray
        if len(minuend) == len(subtrahend):
            diffCHGCAR.__chargeArray = [x - y for x, y in zip(minuend, subtrahend)]
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
        tmp = self.chargeArray
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
    arg = argparse.ArgumentParser()
    arg.add_argument('--add', action='store_true', default=False,
                     help="Add two CHGCAR files")
    arg.add_argument('--diff', action='store_true', default=False,
                     help="Get difference of two CHGCAR files")
    arg.add_argument('--output', metavar='file_name',
                     help="""output file name
if not specified, use standard output""")
    arg.add_argument('CHGCAR_file_1')
    arg.add_argument('CHGCAR_file_2')
    arguments = arg.parse_args()
    if arguments.add == arguments.diff:
        raise RuntimeError('\n--add or --diff ? Choose which one.\n')
    if (arguments.CHGCAR_file_1 is None) or (arguments.CHGCAR_file_2 is None):
        # parhaps in such situation, parser raise error automatically.
        raie RuntimeError('Two CHGCAR files are required.')
    a = CHGCAR(arguments.CHGCAR_file_1)
    b = CHGCAR(argumetns.CHGCAR_file_2)
    if arguments.add:
        c = a + b
    else:
        c = a - b
    if arguments.output is not None:
        c.save(arguments.output)
    else:
        print(c)

