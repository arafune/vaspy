#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''.. py:module:: chgcar

Module for CHGCAR class

translate from chgcar.rb in scRipt4VASP, 2014/2/26 master branch
'''
from __future__ import division, print_function  # Version safety
import re
import copy
import os
import sys
import bz2
import numpy as np

try:
    from vaspy import poscar, tools
except ImportError:
    MYPATH = os.readlink(__file__) if os.path.islink(__file__) else __file__
    sys.path.append(os.path.dirname(os.path.abspath(MYPATH)))
    import poscar
    import tools


_RE_BLANK = re.compile(r'^[\s]*$')
_RE_AUG_OCC = re.compile(r'\baugmentation occupancies')


class CHGCAR(poscar.POSCAR):
    '''.. py:class:: CHGCAR

    Class for CHGCAR

     An example of the first few lines of CHGCAR. ::

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


    Attributes
    ----------

    chg_array, mesh_x, mesh_y, mesh_z, spininfo

    Notes
    -----

    the current verstion ignores "augmentation occupacies".
    '''
    # accessor: chg_array, mesh_x-y-z

    def __init__(self, chgcar_file=None):
        super(CHGCAR, self).__init__(None)
        self.mesh_x = 0
        self.mesh_y = 0
        self.mesh_z = 0
        self.spininfo = 0
        self.chg_array = []
        if chgcar_file:
            self.load_from_file(chgcar_file)

    def load_from_file(self, chgcarfile):
        '''.. py:method:: load_from_file(chgcarfile)

        Parse CHGCAR file to construct CHGCAR object

        Parameters
        ----------

        chgcarfile: str
            CHGCAR file name
        '''
        section = 'poscar'
        separator = None
        tmp = []
        if os.path.splitext(chgcarfile)[1] == '.bz2':
            try:
                thefile = bz2.open(chgcarfile, mode='rt')
            except AttributeError:
                thefile = bz2.BZ2File(chgcarfile, mode='r')
        else:
            thefile = open(chgcarfile)
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
                        self.mesh_x, self.mesh_y, self.mesh_z = [
                            int(str) for str in line.split()]
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
                        self.chg_array.extend(line.split())
        #
        self.chg_array = np.array(self.chg_array, dtype=np.float64)
        #
        if divmod(len(self.chg_array),
                  self.mesh_x * self.mesh_y * self.mesh_z) == (1, 0):
            self.spininfo = [""]
        elif divmod(len(self.chg_array),
                    self.mesh_x * self.mesh_y * self.mesh_z) == (2, 0):
            self.spininfo = ["up+down", "up-down"]
        elif divmod(len(self.chg_array),
                    self.mesh_x * self.mesh_y * self.mesh_z) == (4, 0):
            self.spininfo = ["mT", "mX", "mY", "mZ"]
        else:
            raise RuntimeError("CHGCAR is correct?")

    def magnetization(self, direction=None):
        '''.. py:method:: magnetization(direction)

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
        should be one of 'x', 'y', 'z' and 't'

        Parameters
        ----------

        direction: str
            specify x, y, z or t in noncollinear calculation.
            't' means the total.

        Returns
        -------

        CHGCAR
             CHGCAR of the spin-distribution
        '''
        if len(self.spininfo) == 1:
            raise RuntimeError("This CHGCAR is not spinresolved version")
        dest_chgcar = copy.deepcopy(self)
        if len(self.spininfo) == 2:
            dest_chgcar.chg_array = dest_chgcar.chg_array.reshape(2,
                                                                  self.mesh_x,
                                                                  self.mesh_y,
                                                                  self.mesh_z)[1]
            dest_chgcar.spininfo = ["up-down"]
        elif len(self.spininfo) == 4:
            dest_chgcar.chg_array.reshape(4,
                                          self.mesh_x,
                                          self.mesh_y,
                                          self.mesh_z)
            if direction is None or direction == 't':
                dest_chgcar.chg_array = dest_chgcar.chg_array.reshape(4,
                                                                      self.mesh_x,
                                                                      self.mesh_y,
                                                                      self.mesh_z)[0]
                dest_chgcar.spininfo = ["mT"]
            if direction == 'x':
                dest_chgcar.chg_array = dest_chgcar.chg_array.reshape(4,
                                                                      self.mesh_x,
                                                                      self.mesh_y,
                                                                      self.mesh_z)[1]
                dest_chgcar.spininfo = ["mX"]
            elif direction == 'y':
                dest_chgcar.chg_array = dest_chgcar.chg_array.reshape(4,
                                                                      self.mesh_x,
                                                                      self.mesh_y,
                                                                      self.mesh_z)[2]
                dest_chgcar.spininfo = ["mY"]
            elif direction == 'z':
                dest_chgcar.chg_array = dest_chgcar.chg_array.reshape(4,
                                                                      self.mesh_x,
                                                                      self.mesh_y,
                                                                      self.mesh_z)[3]
                dest_chgcar.spininfo = ["mZ"]
        return dest_chgcar

    def majorityspin(self):
        '''.. py:method:: majorityspin()

        Return CHGCAR for majority spin

        This method is for CHGCAR given by ``ISPIN=2`` but not-SOI
        calculations.

        Returns
        -------

        CHGCAR
            CHGCAR for the majority spin charge
        '''
        if len(self.spininfo) != 2:
            raise RuntimeError('This CHGCAR is not spinresolved version')
        dest_chgcar = copy.deepcopy(self)
        tmp = dest_chgcar.chg_array.reshape(2,
                                            self.mesh_x,
                                            self.mesh_y,
                                            self.mesh_z)
        dest_chgcar.chg_array = (tmp[0] + tmp[1]) / 2
        dest_chgcar.spininfo = ["up"]
        return dest_chgcar

    def minorityspin(self):
        '''.. py:method:: majorityspin()

        Return CHGCAR for minority spin

        This method is for CHGCAR given by ``ISPIN=2`` but not-SOI
        calculations.

        Returns
        -------

        CHGCAR
            CHGCAR for the minority  spin charge
        '''
        if len(self.spininfo) != 2:
            raise RuntimeError('This CHGCAR is not spinresolved version')
        dest_chgcar = copy.deepcopy(self)
        tmp = dest_chgcar.chg_array.reshape(2,
                                            self.mesh_x,
                                            self.mesh_y,
                                            self.mesh_z)
        dest_chgcar.chg_array = (tmp[0] - tmp[1]) / 2
        dest_chgcar.spininfo = ["up"]
        return dest_chgcar

    def __add__(self, other):
        '''.. py:method:: __add__(other)

        x.__add__(y) <=> x + y

        Parameters
        ----------

        other : CHGCAR
            addition CHGCAR object

        Returns
        -------

        CHGCAR
            CHGCAR of the result by summing two CHGCARs:

        Notes
        -----

        In the returned CHGCAR, \
        the charge distribution is just summantion of two CHGCARs, \
        and the atoms are also summantion of two CHGCARs.
        '''
        # augend + aggend
        if not isinstance(other, CHGCAR):
            return NotImplemented
        add_chgcar = super(CHGCAR, self).__add__(other)
        if any([self.mesh_x != other.mesh_x,
                self.mesh_y != other.mesh_y,
                self.mesh_z != other.mesh_z]):
            raise RuntimeError('Mesh sizes are inconsistent')
        try:
            add_chgcar.chg_array = self.chg_array + other.chg_array
        except ValueError:
            raise RuntimeError('the mesh sies are different.')
        return add_chgcar

    def __sub__(self, other):
        '''.. py:method:: __sub__(other)

        x.__sub__y <=> x - y

        Parameters
        -----------

        other: CHGCAR
            difference CHGCAR object

        Returns
        -------

        CHGCAR
            CHGCAR of the result of difference between two CHGCARs:

        Notes
        -----

        In the returned CHGCAR, the charge distribution is \
        just difference of two CHGCARs, and the atoms are used \
        for "minuend" CHGCAR, not difference. \
        The atoms in subtrahend CHGCAR are totally ignored.
        '''
        # minuend - subtrahend
        if not isinstance(other, CHGCAR):
            return NotImplemented
        diff_chgcar = copy.deepcopy(self)
        if any([self.mesh_x != other.mesh_x,
                self.mesh_y != other.mesh_y,
                self.mesh_z != other.mesh_z]):
            raise RuntimeError('Mesh sizes are incinsistent')
        try:
            diff_chgcar.chg_array = self.chg_array - other.chg_array
        except ValueError:
            raise RuntimeError('the mesh sizes are different.')
        return diff_chgcar

    def __str__(self):
        '''.. py:method:: __str__()

        x.__str__() <=> str(x)

        Returns
        -------

        str
            a string representation of CHGCAR.
        '''
        outputstring = ''
        chgarray = self.chg_array.reshape(len(self.spininfo),
                                          self.mesh_x *
                                          self.mesh_y *
                                          self.mesh_z)
        for tmp in chgarray:
            output = []
            outputstring += '\n  {0}  {1}  {2}\n'.format(self.mesh_x,
                                                         self.mesh_y,
                                                         self.mesh_z)
            for array in tools.each_slice(tmp, 5):
                output.append(''.join('  {0:18.11E}'.format(i)
                                      for i in array if i is not None))
            outputstring += '\n'.join(output)
        return super(CHGCAR, self).__str__() + outputstring + '\n'

    def save(self, filename):
        '''.. py:method:: save(filename)

        Save CHGCAR object as CHGCAR file style

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
