#!/usr/bin/env python 
# -*- coding: utf-8 -*-

from __future__ import print_function, division
import re as _re
import number.Number as _Number


class DOSCAR(object):
    '''# class for DOSCAR file
       # @author Ryuichi Arafune
    '''
    # attr_accessor :nAtom, :dos_container

    def __init__(self, arg=None):
        self.nAtom = 0
        self.__nEnergypoints = 0
        self.dos_container = list()
        
        if arg is not None:
            self.load_from_file(arg)

    
    def load_from_file(self, doscar_file):
        """# @param [String] doscar_file file name of "DOSCAR"."""
        with open(doscar_file) as f:
            separate_text = ""
            aDOS = list()
            for idx, line in enumerate(f, 1):
                if idx == 1:
                    self.nAtom = int(line.split()[0])
                elif 2<= idx <= 5:
                    continue
                elif idx == 6:
                    separate_text = line.rstrip('\n')
                    self.__nEnergypoints = int(separate_text.split()[2])
                    continue
                else: # idx >= 7
                    if _re.search(separate_text, line):
                        self.dos_container.append(aDOS)
                        aDOS = list()
                    else:
                        aDOS.append([float(i) for i in line.split()])
            else:
                self.dos_container.append(aDOS)



class DOS(object):
    '''# Class for DOS
       #
       #  attribute : dos
       #       Array object.  that consists arrays of two elements.
       #       the first element is the float representing the energy.
       #       the other element is still array that contains the density.
       #
       # @author Ryuichi Arafune
    '''
    # attr_accessor :dos
    def __init__(self, array=None):
        self.__dos = list()
        if array is not None:
            for a in array:
                a_float = [float(i) for i in a]
                self.__dos.append((a_float[0], a_float[1:]))
    

    @property
    def dos(self):
        return self.__dos

    
    def __copy__(self, orig):
        self.__dos = self.dos.copy()

    
    def __iter__(self):
        for ith_point in self.dos:
            yield ith_point
    

    def append(self, dos_data):
        '''# @param [Array] dos_data
        '''
        self.__dos.append(_filter_dos_data(dos_data))
    

    def pop(self, i=-1):
        '''# @return [Array] return and remove the last element 
           #    (the highest energy data) of the DOS object.
        '''
        self.__dos.pop(i)
    

    # self.unshift(X) => self.[0:0] = X
    # self.shift => self.pop(0)
    

    def __len__(self):
        """x.__len__() <=> len(x)"""
        return len(self.dos)

    
    def __setitem__(self, idx, value):
        """x.__setitem__(i, y) <=> x[i]=y"""
        if isinstance(idx, int):
            value = _filter_dos_data(value)
        if isinstance(idx, slice):
            value = [_filter_dos_data(each) for each in value]
        self.__dos.__setitem__(idx, value)

    
    def __getitem__(self, idx):
        """x.__getitem__(i) <=> x[i]"""
        self.__dos.__getitem__(idx)

    
    def fermilevel_correction(self, fermi):
        '# @param [Float] fermi fermi level'
        self.__dos = [(each[0]-fermi, each[1]) for each in self.dos]

    
    def energies(self, i=None):
        '''# @param [Fixnum] i 
           # @return [Float, Array] if i is set, return the energy value of the i-th point.
           #   if not, return the all energies in DOS object by Array representation.
           # @param [Fixnum] i 
           # @return [Float, Array] if i is set, return the energy value of the i-th point.
           #   if not, return the all energies in DOS object by Array representation.
        '''
        if i is None:
            return [each[0] for each in self.dos]
        else:
            return self.dos[i][0]


    def densities(self, i=None):
        '''# @param [Fixnum] i 
           # @return [Array]
        '''
        if i is None:
            return [each[1] for each in self.dos]
        else:
            return self.dos[i][1]

    
    def __str__(self):
        """x.__str__() <=> str(x)
           # @return [String] returns string representation of DOS object.
           # csv-like (tab-deliminated) format.
        """
        return '\n'.join(str(line[0]) + '\t' + '\t'.join(line[1])
                            for line in self.dos)


def _filter_dos_data(data):
    listlikes = (list, tuple)
    if not isinstance(data, listlikes):
        raise RuntimeError("Invalid argument")
    if not isinstance(data[0], _Number):
        raise RuntimeError("non Numeric instance in header.")
    if not isinstance(data[1], listlikes):
        raise RuntimeError("non list/tuple instance in data.")
    return data

"""
From VASP webpage:

The file DOSCAR contains the DOS and integrated DOS The units are
"number of states/unit cell". For dynamic simulations and relaxations,
an averaged DOS and an averaged integrated DOS is written to the file.
For a description of how the averaging is done see 7.18, 7.32).
The first few lines of the DOSCAR file are made up by a header,
which is followed by NDOS lines holding three data

 energy     dos     integrated_dos

For spin-polarized calculations each line holds five data

 energy     dos(up) dos(dwn)  integrated_dos(up) integrated_dos(dwn)

If RWIGS (Wigner Seitz radii, see section 7.29) is set in the INCAR file,
a l- and site-projected DOS is calculated and also written to the file
DOSCAR. One set of data is written for each ion,
each set of data holds NDOS lines with the following data

 energy s-dos p-dos d-dos

and

 energy s-dos(up) p-dos(up) d-dos(up) s-dos(dwn) p-dos(dwn) d-dos(dwn)

for the non spin-polarized and spin polarized case respectively.
The units of the l- and site projected DOS are states/atom.
Please mind, that the site projected DOS is not evaluated
in the parallel version if NPAR tex2html_wrap_inline5201 1.

Mind: For relaxations the DOSCAR is usually useless.
If you want to get an accurate DOS for the final configuration
copy CONTCAR to POSCAR and make another static (ISTART=1; NSW=0)
calculation. 
"""