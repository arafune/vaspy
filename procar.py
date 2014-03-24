#!/usr/bin/env python
# -*- coding: utf-8 -*-
# translate from procar.rb of scRipt4VASP 2014/2/26 master branch

from __future__ import division, print_function
import re
import numpy as np

class PROCAR:
    '''Class for PROCAR file
 
   PROCAR consists of these lines.  Appeer once per file.
   1 # the first line
    ex.)   PROCAR lm decomposed + phase  
   2 set number of k-points, bands and ions.  
      Appear once when spin-integrated, twice when spin-resolved.
    ex.)   # of k-points:   50         # of bands: 576         # of ions:  98
   3  k-point character
    ex.)  k-point    1 :    0.00000000 0.00000000 0.00000000     weight = 0.02000000
    note that the first character is "blank".
   4 band character
    ex.)  band   1 # energy  -11.87868466 # occ.  2.00000000
   5 orbital contribution.
    ex.)  1  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000
  @author Ryuichi Arafune
'''
    #attr_reader :orbital, :phase, :spininfo, :numk,
    #:nBands, :nAtoms, :kvectors, :energies, :oritalname
    def __init__(self, arg=None, phase_read=None):
        self.__orbital = list()
        self.__phase = list()
        self.__spininfo = 0 # nospin: 1, spinresolved: 2, soi: 4
        self.__numk = 0
        self.__nBands = 0
        self.__nAtoms = 0
        self.__kvectors = list()
        self.__energies = list()

        if isinstance(arg, str):
            self.load_from_file(arg)
        elif isinstance(arg, (list, tuple)):
            self.load_from_array(arg)

    @property
    def orbital(self):
        return self.__orbital

    @property
    def phase(self):
        return self.__phase

    @property
    def spininfo(self):
        return self.__spininfo

    @property
    def numk(self):
        return self.__numk

    @property
    def nBands(self):
        return self.__nBands

    @property
    def nAtoms(self):
        return self.__nAtoms

    @property
    def kvectors(self):
        return self.__kvectors

    @property
    def energies(self):
        return self.__energies

    @property
    def orbitalname(self):
        return self.__orbitalname

    def load_from_file(self, file, phase_read=None):
        '''PROCAR#load_from_file
    A virtual parser of PROCAR
    @param [String] file filename of *PROCAR* file
    @param [Boolian] phase_read Switch for loading phase characters
'''
        section = list()
        separator_to_orbital = None
        with open(file) as f:
            for line in f:
                if re.findall(r'^[\s]*$', line): continue
                elif re.findall(r'^#', line):
                    self.__numk, self.__nBands, self.__nAtoms = [int(i) for i in line.split() if i.isdigit()]
                elif re.findall(r'\bk-points\b', line):
                    self.__kvectors.append([float(i) for i in line.split()[3:6]]) # check data
                    section.pop()
                elif re.findall(r'^band\b', line):
                    self.__energies.append(float(line.split()[4]))
                    section.pop()
                elif re.findall(r'^ion\b', line):
                    separator_to_orbital = separator_to_orbital if 'separator_to_orbital' in locals() else line.rstrip('\n')
                    separator_to_phase = separator_to_phase if 'separator_to_phase' in locals() else separator_to_orbital[0:-7]
                    self.__orbitalname = self.__orbitalname if hasattr(self, 'orbitalname') else separator_to_orbital.split()
                    if re.findall(separator_to_orbital, line):
                        section = ['orbital']
                    elif re.findall(separator_to_phase, line):
                        section = ['phase']
                else:
                    if section == ['orbital']:
                        if re.findall(r'\btot\b', line): continue
                        tmp = [float(i) for i in line.split()]
                        tmp[0] = int(tmp[0])
                        self.__orbital.append(tmp)
                    elif section == ['phase']:
                        if re.findall(r'\btot\b', line): continue
                        if not phase_read: continue
                        tmp = [float(i) for i in line.split()]
                        tmp[0] = int(tmp[0])
                        self.__phase.append(tmp)
        self.__spininfo = len(self.orbital) // (self.numk * self.nBands * self.nAtoms)
        if len(self.orbital) % (self.numk * self.nBands * self.nAtoms) != 0:
            raise RuntimeError("PROCAR file may be broken")
        if self.spininfo == 1:
            self.__spininfo = ['']
        elif self.spininfo == 2:
            self.__spininfo = ['_up', '_down']
        elif self.spininfo == 4:
            self.__spininfo = ['_mT', '_mX', '_mY', '_mZ']
        # orbitalname
        tmpOrb = Orbital()
        tmpOrb.redefine_orbital_list(self.orbitalname)

    def __str__(self):
        '''x.__str__() <=> str(x)

    show the PROCAR character, not contents.
'''
        template = '''The properties of this procar:
  # of k-points: {0.numk}
  # of bands: {0.nBands}
  # of ions: {0.nAtoms}
  # of kvectors: {1}
  # of energies: {2}
     ((# of k-points) * (# of bands) = {0.numk}*{0.nBands}={3}
  # of orbital component: {4}
     ((# of k-points) * (# of bands) * (# of ions) = {0.numk}*{0.nBands}*{0.nAtoms}={5}
  # of phase component: {6}
  Orbitals are: {0.orbitalname}
  spininfo: {0.spininfo}
'''
        return template.format(self,
                               len(self.kvectors),
                               len(self.energies),
                               self.numk * self.nBands,
                               len(self.orbital),
                               self.numk * self.nBansd * self.nAtoms,
                               len(phase))

    def __iter__(self):
        for line in self.orbital:
            yield line

    def to_band(self):
        band = Band.(self.kvectors[0:self.numk]) # for no-spin and soi,
                                                 # [0:self.numk] does not affect
                                                 # the results.  However for
                                                 # spin-case, remove the
                                                 # identical latter kvectors.
        for i, orbital in enumerate(self.orbital):
            if len(spininfo) == 1: # spin integrated
                kindex, x = divmod(i, self.nBands * self.nAtoms)
                bindex = x // self.nAtoms
                eigenvalue = self.energies[i // self.nAtoms]
                spininfo = self.spininfo[0]
                aState = State(kindex, bindex, eigenvalue, spininfo, orital)
                band.push(aState)
            elif len(spininfo) == 2: # spin resolved
                kindex, x = divmod(i, self.nBands * self.nAtoms)
                if i >= self.numk * self.nBands * self.nAtoms: kindex -= self.numk
                bindex = x // self.nAtoms
                eigenvalue = self.energies[i // self.nAtoms]
                if i >= self.numk * self.nBands * self.nAtoms:
                    spininfo = self.spininfo[1]
                else:
                    spininfo = self.spininfo[0]
                aState = State(kindex, bindex, eigenvalue, spininfo, orbital)
            elif len(spininfo) == 4: # SOI
                kindex, x = divmod(i, self.nBands * self.nAtoms * 4) ## 4 is mT, mX, mY, mZ
                bindex = x // (self.nAtoms * 4)                      ## 4 is mT, mX, mY, mZ
                eigenvalue = self.energies[i // (self.nAtoms * 4)]   ## 4 is mT, mX, mY, mZ
                spininfo = self.spininfo[i % 4]                      ## 4 is mT, mX, mY, mZ
                aState = State(kindex, bindex, eigenvalue, spininfo, orbital)
                band.append(aState)
        return band

#
