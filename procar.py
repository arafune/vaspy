#!/usr/bin/env python
# -*- coding: utf-8 -*-
# translate from procar.rb of scRipt4VASP 2014/2/26 master branch

from __future__ import division, print_function
import re, copy, os, sys
import itertools as it
import functools as ft
import numpy as np
mypath = os.readlink(__file__) if os.path.islink(__file__) else __file__
sys.path.append(os.path.dirname(os.path.abspath(mypath)))
import tools

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
    def __init__(self, arg=None, phase_read=False):
        self.__orbital = list()
        self.__orbitalname = list()
        self.__phase = list()
        self.__spininfo = 0 # nospin: 1, spinresolved: 2, soi: 4
        self.__numk = 0
        self.__nBands = 0
        self.__nAtoms = 0
        self.__kvectors = list()
        self.__energies = list()

        if isinstance(arg, str):
            self.load_from_file(arg, phase_read)
        elif isinstance(arg, (list, tuple)):
            self.load_from_array(arg, phase_read)

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

    def load_from_file(self, file, phase_read=False):
        '''PROCAR#load_from_file
    A virtual parser of PROCAR
    @param [String] file filename of *PROCAR* file
    @param [Boolian] phase_read Switch for loading phase characters
'''
        section = list()
        separator_to_orbital = None
        separator_to_phase = None
        with open(file) as f:
            for line in f:
                if re.findall(r'^[\s]*$', line): continue
                elif re.findall(r'^#', line):
                    self.__numk, self.__nBands, self.__nAtoms = [int(i) for i in line.split() if i.isdigit()]
                elif re.findall(r'\bk-points\b', line):
                    self.__kvectors.append(np.asarray([list(float(i) for i in line.split()[3:6])])) # check data
                    section.pop()
                elif re.findall(r'^band\b', line):
                    self.__energies.append(float(line.split()[4]))
                    section.pop()
                elif re.findall(r'^ion\b', line):
                    separator_to_orbital = separator_to_orbital or line.rstrip('\n')
                    separator_to_phase = separator_to_phase or separator_to_orbital[0:-7]
                    self.__orbitalname = self.orbitalname or separator_to_orbital.split()
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
            self.__spininfo = ('',)
        elif self.spininfo == 2:
            self.__spininfo = ('_up', '_down')
        elif self.spininfo == 4:
            self.__spininfo = ('_mT', '_mX', '_mY', '_mZ')
        # orbitalname
        tmpOrb = Orbital()
        tmpOrb.redefine_orbital_list(self.orbitalname)

    def load_from_array(procar, phase_read=False):
        '''This method effectively acts as a parser of PROCAR.
    @param [Array] procar Array of PROCAR file (IO.readlines(procar file))
    @param [Boolian] phase_read Switch for loading phase characters
    @deprecated Too sloooow!! (about 10 times slower than load_from_file)
'''
        section = list()
        #
        self.__numk, self.__nBands, self.__nAtoms = [int(i) for i in procar[1].split() if i.isdigit()]
        separator_to_orbital = procar[7].rstrip('\n')
        self.__orbitalname = separator_to_orbital.split()
        separator_to_phase = separator_to_orbital[0:-7]
        for line in procar:
            if re.findall(r'^[\s]*$', line): continue
            elif re.findall(r'^#', line): continue
            elif re.findall(r'\bk-points\b', line):
                self.__kvectors.append([float(i) for i in line.split()[3:6]]) # check data
                section.pop()
            elif re.findall(r'^band\b', line):
                self.__energies.append(float(line.split()[4]))
                section.pop()
            elif re.findall(separator_to_orbital, line):
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
            self.__spininfo = ('',)
        elif self.spininfo == 2:
            self.__spininfo = ('_up', '_down')
        elif self.spininfo == 4:
            self.__spininfo = ('_mT', '_mX', '_mY', '_mZ')
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
        return iter(self.__orbital)

    def to_band(self):
        band = Band(self.kvectors[0:self.numk]) # for no-spin and soi,
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

class Band:
    '''Band class
  container of State objects
  class variable : @@kvectors
  @author Ryuichi Arafune
'''
    # attr_accessor :band
    __kvectors = list()
    __distance = list()

    def __init__(self, arg=None):
        '''Band(list-of-kvectors) => Band object
    Band#initialize: Sets @@kvectors and then calculate @@distance
    @param[Array] arg kvectors array
'''
        self.__band = list()
        if arg is not None:
            self.__kvectors = arg
            self.__distance.extend([None] * len((self.kvectors) - len(self.distance)))
            for i, k in enumerate(self.kvectors):
                if i == 0:
                    self.__distance[i] = 0.0
                else:
                    self.__distance[i] = self.distance[i - 1] + np.linalg.norm(self.kvectors[i - 1] - k)
    
    @property
    def kvectors(self):
        '''Band#kvectors
  # @return [Array<Float, Float, Float>] Returns @@kvectors
'''
        return self.__kvectors

    @property
    def distance(self):
        '''@return [Array<Float>] Returns @@distance
'''
        return self.__distance

    def initialize_copy(self, org):
        '''for Band.dup()'''
        self.band = self.band.copy()

    def __getitem__(self, i):
        '''x.__getitem__(i) <==> x[i]
  # same as list[i] or tuple[i]
  # @return[State] i-th State in *band* object
'''
        return self.band[i]

    def __setitem__(self, index, value):
        '''x.__setitem__(i, y) <==> x[i]=y
    value must be a procar.State (for single index)
    or iterable of procar.State objects (ranged indices)'''
        if isinstance(index, slice):
            value = tuple(value)
            if all(isinstance(item, State) for item in value):
                self.__band[index] = value
            else:
                raise TypeError('argument must be an iterable of procar.State objects.')
        else:
            if isinstance(value, State):
                self.__band[index] = value
            else:
                raise TypeError('argument must be a procar.State object.')

    def __iter__(self):
        '''x.__iter__() <==> iter(x)
to work as iterable
return iterator object'''
        return iter(self.__band)

    def fermilevel_correction(self, ef): # check whether it works correctly
        '''
  # eigenvalues are corrected by fermi level.
  # Usually the value are obtained from OUTCAR.
  #
  #  @param[Float] ef
'''
        for aState in self:
            aState.fermilevel_correction(ef)

    def append(self, aState):
        '''same as list.append(), but argument must be a State object.
    @param [State] aState
    '''
        if not isinstance(aState, State):
            raise TypeError('argument must be a procar.State object.')
        self.__band.append(aState)

    def extend(self, iterable_of_States):
        '''same as list.extend(), but argument must be an iterable of State objects.
'''
        statetuple = tuple(iterable_of_States)
        if not all(isinstance(item, State) for item in statetuple):
            raise TypeError("All items in iterable must be a State object.")
        self.__band.extend(statetuple)

    def pop(self, n=-1):
        return self.__band.pop(n)
    
    # Band.shift() [Ruby] <==> Band.pop(0) [Python]
    # Band.shift(n) [Ruby] <==> [Band.pop(0) for i in range(n)] # n != 1 [Python]
    # Array.unshift(a, b, c...) [Ruby] <==> list[:0] = (a, b, c...) [Python]
    
    def insert(self, index, item):
        self.__band.insert(index, item)

    def __len__(self):
        '''x.__len__() <==> len(x)'''
        return len(self.band)

    def select_by_site(self, *sites):
        sites = flatten(sites)
        sites = sorted(set(sites))
        dest = Band()
        #dest.band = self.find_all{|s| sites.include?(s[:ion])}
        return dest

    def __add__(self, other):
        '''x.__add__(y) <==> x+y
'''
        other = _convert_other_band(other)
        if other is NotImplemented: return other
        dest = Band()
        dest.__band = self.band + other.band
        return dest

    def sort(self, key=None, reverse=False):
        '''In-place sorting of self.__band.
Use sorted(self) for not In-place sorting.'''
        self.__band.sort(key, reverse)

    def site_integrate(self, *sites):
        '''  # @param [Fixnum, Range, Array] sites sites are
  #     specified by comma-separated number, array,
  #     or Range object.
  # @return [Band] site-integrated band object
'''
        tmpBand = self.select_by_site(*sites)
        siteGroup = dict()
        for aState in tmpBand:
            x = (aState.kindex, aState.bandindex, aState.spininfo)
            if siteGroup.get(x) is None:
                siteGroup[x] = aState
            else:
                siteGroup[x] += aState
        dest = Band()
        dest.__band = list(siteGroup.values())
        return dest

    def extract_orbitals_in_place(self, *orbital_symbols):
        '''  # extract certain orbital component in place
  # @param [Array] orbital_symbols
  # @return [self]
'''
        for aState in self.band:
            aState.extract_orbitals(*orbital_symbols)
    
    def extract_orbitals(self, *orbital_symbols):
        '''  # extract certain orbital component in place
  # Receiver itself does not change
  # @param [Array] orbital_symbols
  # @return [Band] orbital_extracted Band
'''
        dest = copy.deepcopy(self)
        for aState in dest.band:
            aState.extract_orbitals(*orbital_symbols)
        return dest

    def header(self):
        '''  # @param[Array] orbital_symbols
  # @return [Array] Returns arrray represents header for 
  #   Band#to_a 
'''
        nSites = self.number_of_sites()
        nSpintype = self.number_of_spintype()
        tmp = sorted(self)[0:(nSites * nSpintype)]
        tmp = [str(aState.site) + str(orbital) + aState.spininfo
               for aState in tmp
               for orbital in aState.orbital_keys()]
        tmp[0:0] = ['k', 'energy']
        return tmp

    def tolist(self): # orbital_symbols?
        '''x.tolist() => list

  # @param [Array] orbital_symbols
  # @return [Array] Returns *2D array* (Array of array) for output
'''
        nBands = self.number_of_bands()
        nSites = self.number_of_sites()
        nSpintype = self.number_of_spintype()
        orbital_symbols = self[0].defined_orbital_list()
        dest = list()
        distance = self.distance * nBands
        self.sort()

        for array in tools.each_slice(self, nSites * nSpintype):
            dest.append([aState.send(orbital) for aState in array
                        for orbital in aState.orbital_keys ].insert(0, array[0].eigenvalue))
        return list(tools.flatten(i) for i in zip(distance, dest))

    def __str__(self):
        '''x.__str__() <==> str(x)

  # @return[String] Returns cvs formatted text'''
        array = self.tolist()
        header = self.header()
        nColumn = len(header)
        text = '\t'.join(header) + '\n'
        for i, e in enumerate(array):
            if i % len(self.distance) == 0 and i > 0:
                text += '\t' * (nColumn - 1) + '\n'
                text += '\t'.join(e) + '\n'
            else:
                text += '\t'.join(e) + '\n'
        return text

    def save(self, filename):
        '''Store bandstructure with orbital
  #   contribution to the file.
  # @param [String] filename
'''
        with open(filename, mode='w') as file:
            dum = file.write(str(self))

    def rename_site(*new_names):
        '''
  # @param [Array] new_names
'''
        sites = self.sites
        if len(sites) != len(new_names): raise RuntimeError("Number of sites and new names are must be the same.")
        rule = dict(zip(sites, new_names))
        for aState in self:
            aState['ion'] = rule[aState['ion']]

    def group_by_spintype(self):
        '''  # Returns Array consists of spin selected band.  
  # @return[Array<Band>] 
'''
        tmp = dict()
        for aState in self:
            tmp[aState.spininfo] = tmp.get(aState.spininfo, []).append(aState)
        dest = dict()
        for key, value in tmp.items():
            band = Band()
            band.__band = value
            dest[key] = band
        if '_mT' in dest:
            destArray = [dest['_mT'], dest['_mX'], dest['_mY'], dest['_mZ']]
        elif '_up' in dest:
            destArray = [dest['_up'], dest['_down']]
        else:
            destArray = list(dest.values())
        return destArray

    def number_of_bands(self):
        '''  # @return [Fixnum] Returns the number of bands in *band*'''
        return len(set(aState.bandindex for aState in self))

    def number_of_sites(self):
        '''  # @return [Fixnum] Returns the number of sites in *band*'''
        return len(set(aState['ion'] for aState in self))

    def number_of_spintype(self):
        '''  # @return [Fixnum] Returns the number of spintype in *band*'''
        return len(set(aState.spininfo for aState in self))

    def sites(self):
        '''  # @return [Array] Return *array* that consists of site names'''
        return set(aState['ion'] for aState in self)

    def select_by_band(self, *bandindices):
        '''
  # @param [Fixnum, Range, Array] bandindexes bandindexes
  #   are specified by comma-separated number, array,
  #   or Range object. 
  # @example Band#select_by_band(1,2,3,4) : same as 
  #   Band#select_by_band(1..4)/Band#select_by_band([1,2,3,4]) 
  # @return [Band] Returns a new band object that consists 
  #    of States specified by bandindex.
'''
        bandindices = sorted(set(tools.flatten(bandindices)))
        dest = Band()
        dest.__band = [s for s in self if s.bandindex in bandindices]
        return dest

    def dump(self, filename):
        pass
        

class State:

    def extract_orbitals(self, *orbital_symbols):
        pass
    
    def orbital_keys(self):
        pass

def _convert_other_band(other):
    if isinstance(other, Band):
        return other
    return NotImplemented
