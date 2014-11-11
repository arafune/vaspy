#!/usr/bin/env python
# -*- coding: utf-8 -*-
# translate from procar.rb of scRipt4VASP 2014/2/26 master branch

from __future__ import division, print_function # Version safety
import re, copy, os, sys, csv
import functools as ft
if sys.version_info[0] >= 3: # Version safety
    from io import StringIO
else:
    from cStringIO import StringIO
import numpy as np
mypath = os.readlink(__file__) if os.path.islink(__file__) else __file__
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(mypath))))
from vaspy import tools


class PROCAR(object): # Version safety
    '''
    Class for PROCAR file

    .. py:class:: PROCAR(object)

    PROCAR consists of these lines.  Appeer once per file.

    1. # the first line
     ex.)   PROCAR lm decomposed + phase
    2. set number of k-points, bands and ions.
      Appear once when spin-integrated, twice when spin-resolved.
     ex.)   # of k-points:   50         # of bands: 576         # of ions:  98
    3.  k-point character
     ex.)  k-point    1 :    0.00000000 0.00000000 0.00000000     weight = 0.02000000
     note that the first character is "blank".
    4. band character
     ex.)  band   1 # energy  -11.87868466 # occ.  2.00000000
    5. orbital contribution.
     ex.)  1  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000

    :author: Ryuichi Arafune
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
                    self.__numk, self.__nBands, self.__nAtoms = [
                        int(i) for i in line.split() if i.isdigit()]
                elif re.findall(r'\bk-point\b', line):
                    self.__kvectors.append(np.asarray(
                        [list(float(i) for i in line.split()[3:6])]))
                    section = []
                elif re.findall(r'^band\b', line):
                    self.__energies.append(float(line.split()[4]))
                    section = []
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

    def load_from_array(self, procar, phase_read=False):
        '''This method effectively acts as a parser of PROCAR.
        @param [Array] procar Array of PROCAR file (IO.readlines(procar file))
        @param [Boolian] phase_read Switch for loading phase characters
        @deprecated Too sloooow!! (about 10 times slower than load_from_file)
        '''
        section = list()
        #
        self.__numk, self.__nBands, self.__nAtoms = [
            int(i) for i in procar[1].split() if i.isdigit()]
        separator_to_orbital = procar[7].rstrip('\n')
        self.__orbitalname = separator_to_orbital.split()
        separator_to_phase = separator_to_orbital[0:-7]
        for line in procar:
            if re.findall(r'^[\s]*$', line): continue
            elif re.findall(r'^#', line): continue
            elif re.findall(r'\bk-points\b', line):
                self.__kvectors.append(np.asarray(
                        [list(float(i) for i in line.split()[3:6])]))
                section = []
            elif re.findall(r'^band\b', line):
                self.__energies.append(float(line.split()[4]))
                section = []
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
     ((# of k-points) * (# of bands) = {0.numk}*{0.nBands}={3})
  # of orbital component: {4}
     ((# of k-points) * (# of bands) * (# of ions) = {0.numk}*{0.nBands}*{0.nAtoms}={5})
  # of phase component: {6}
  Orbitals are: {0.orbitalname}
  spininfo: {0.spininfo}
'''
        return template.format(self,
                               len(self.kvectors),
                               len(self.energies),
                               self.numk * self.nBands,
                               len(self.orbital),
                               self.numk * self.nBands * self.nAtoms,
                               len(self.phase))

    def __iter__(self):
        return iter(self.__orbital)

    def to_band(self):
        band = Band(self.kvectors[0:self.numk]) # for no-spin and soi,
                                                # [0:self.numk] does not affect
                                                # the results.  However for
                                                # spin-case, remove the
                                                # identical latter kvectors.
        for i, orbital in enumerate(self.orbital):
            if len(self.spininfo) == 1: # spin integrated
                kindex, x = divmod(i, self.nBands * self.nAtoms)
                bindex = x // self.nAtoms
                eigenvalue = self.energies[i // self.nAtoms]
                spininfo = self.spininfo[0]
                aState = State(kindex, bindex, eigenvalue, spininfo, orbital)
                band.append(aState)
            elif len(self.spininfo) == 2: # spin resolved
                kindex, x = divmod(i, self.nBands * self.nAtoms)
                if i >= self.numk * self.nBands * self.nAtoms: kindex -= self.numk
                bindex = x // self.nAtoms
                eigenvalue = self.energies[i // self.nAtoms]
                if i >= self.numk * self.nBands * self.nAtoms:
                    spininfo = self.spininfo[1]
                else:
                    spininfo = self.spininfo[0]
                aState = State(kindex, bindex, eigenvalue, spininfo, orbital)
                band.append(aState)
            elif len(self.spininfo) == 4: # SOI
                kindex, x = divmod(i, self.nBands * self.nAtoms * 4) ## 4 is mT, mX, mY, mZ
                bindex = x // (self.nAtoms * 4)                      ## 4 is mT, mX, mY, mZ
                eigenvalue = self.energies[i // (self.nAtoms * 4)]   ## 4 is mT, mX, mY, mZ
                spininfo = self.spininfo[i % 4]                      ## 4 is mT, mX, mY, mZ
                aState = State(kindex, bindex, eigenvalue, spininfo, orbital)
                band.append(aState)
        return band

class Band(object): # Version safety
    '''
    .. py:class:: Band class

    container of State objects

    :class variable: kvectors
    :author: Ryuichi Arafune
    '''
    # attr_accessor :band
    __kvectors = list()
    __distance = list()

    def __init__(self, arg=None):
        '''
        Band(list-of-kvectors) => Band object
        Band#initialize: Sets @@kvectors and then calculate @@distance

        :param Array arg: kvectors array
        '''
        self.__band = list()
        if arg is not None:
            self.__kvectors = arg
            self.__distance.extend([None] * (len(self.kvectors) - len(self.distance)))
            for i, k in enumerate(self.kvectors):
                if i == 0:
                    self.__distance[i] = 0.0
                else:
                    self.__distance[i] = self.distance[i - 1] + np.linalg.norm(self.kvectors[i - 1] - k)
    
    @property
    def kvectors(self):
        '''
        Band#kvectors

        :return: kvectors
        :rtype: Array<Float, Float, Float>
        '''
        return self.__kvectors

    @property
    def distance(self):
        '''
        :return: distance
        :rtype: Array<Float>
        '''
        return self.__distance

    @property
    def band(self):
        return self.__band

    def __copy__(self):
        'x.__copy__() <==> copy.copy(x)'
        dest = Band()
        dest.__band = self.band[:]
        return dest

#    def __deepcopy__(self, memo):
#        'x.__deepcopy__() <==> copy.deepcopy(x)'
#        dest = Band()
#        dest.__band == copy.deepcopy(self.band, memo)

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
        eigenvalues are corrected by fermi level.  Usually the value are obtained from OUTCAR.

        :param float ef:
        '''
        for aState in self:
            aState.fermilevel_correction(ef)

    def append(self, aState):
        '''
        same as list.append(), but argument must be a State object.
        :param State aState:
        '''
        if not isinstance(aState, State):
            raise TypeError('argument must be a procar.State object.')
        self.__band.append(aState)

    def extend(self, iterable_of_States):
        '''
        same as list.extend().
        but argument must be a finite iterable of State objects.
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
        sites = tools.flatten(sites)
        sites = sorted(set(sites))
        dest = Band()
        dest.__band = [aState for aState in self if aState['ion'] in sites]
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
        self.__band.sort(key=key, reverse=reverse)

    def site_integrate(self, *sites):
        '''
        :param [Fixnum, Range, Array] sites: sites are specified by comma-separated number, array, or Range object.
        :return:  site-integrated band object
        :rtype: Band
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
        '''
        Extract certain orbital component in place.  Receiver itself does not change
        :param Array orbital_symbols:
        :return: orbital_extracted Band
        :rtype: Band
        '''
        dest = copy.deepcopy(self) ### check it is really deep copy ###
        for aState in dest.band:
            aState.extract_orbitals(*orbital_symbols)
        return dest

    def header(self):
        '''
        :param Array orbital_symbols:
        :return: Arrray represents header for Band#to_a 
        :rtype: Array 
        '''
        nSites = self.number_of_sites()
        nSpintype = self.number_of_spintype()
        tmp = sorted(self)[0:(nSites * nSpintype)]
        tmp = [str(aState.site()) + str(orbital) + aState.spininfo
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
        #orbital_symbols = self[0].defined_orbital_list()
        dest = list()
        distance = self.distance * nBands
        self.sort()
        
        for array in tools.each_slice(self, nSites * nSpintype):
            dest.append([getattr(aState, orbital) for aState in array
                        for orbital in aState.orbital_keys()])
            dest[-1].insert(0, array[0].eigenvalue)
        return list(tools.flatten(i) for i in zip(distance, dest))

    def __str__(self):
        '''x.__str__() <==> str(x)

        # @return[String] Returns csv formatted text'''
        with StringIO() as stream:
            self.export_csv(stream, delimiter='\t', lineterminator='\n')
            return stream.getvalue()


    def export_csv(self, file, **kwargs):
        """Export data to file object (or file-like object) as csv format.
        kwargs are keyword options of csv.writer().
        see help(csv.writer) for detail.
        """
        csvwriter = csv.writer(file, **kwargs)
        array = self.tolist()
        header = self.header()
        nColumn = len(header)
        csvwriter.writerow(header)
        for i, each in enumerate(array):
            if i % len(self.distance) == 0 and i > 0:
                csvwriter.writerow([None] * (nColumn-1))
                csvwriter.writerow(str(x) for x in each)
            else:
                csvwriter.writerow(str(x) for x in each)


    def save(self, filename):
        '''Store bandstructure with orbital
        #   contribution to the file.
        # @param [String] filename
        '''
        try: # Version safety
            file = open(filename, mode='w', newline='')
        except TypeError:
            file = open(filename, mode='wb')
        with file:
            self.export_csv(file, delimiter='\t', lineterminator='\n')


    def rename_site(self, *new_names):
        '''
        # @param [Array] new_names
        '''
        sites = self.sites()
        if len(sites) != len(new_names): 
            raise RuntimeError(
            "Number of sites and new names are must be the same.")
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
        return sorted(set(aState['ion'] for aState in self))

    def select_by_band(self, *bandindices):
        '''
        :param [Fixnum, Range, Array] bandindexes: bandindexes are specified by comma-separated number, array, or Range object. 
        :example: Band#select_by_band(1,2,3,4) : same as Band#select_by_band(1..4)/Band#select_by_band([1,2,3,4]) 
        :return: a new band object that consists of States specified by bandindex.
        :rtype: Band
        '''
        bandindices = sorted(set(tools.flatten(bandindices)))
        dest = Band()
        dest.__band = [s for s in self if s.bandindex in bandindices]
        return dest

    def dump(self, filename):
        pass

class Orbital(object): # Version safety
    '''
    .. py:class:: Orbital(object)
    
    Electronic orbital contribution.
    
    A "functionalized" Hash

    :author: Ryuichi Arafune
    :version: 2.0
    '''
    __orbital_list = ('ion', 's', 'py', 'pz', 'px', 'dxy',
                      'dyz', 'dz2', 'dxz', 'dx2', 'tot')

    def __init__(self, arg=None):
        self.orbital = dict()
        if arg is not None: self.load_from_line(arg)

    def __copy__(self):
        'x.__copy__() <==> copy.copy(x)'
        dest = Orbital()
        dest.orbital = copy.copy(self.orbital)
        return dest

    def orbital_keys(self):
        keys = list(self.orbital.keys())
        keys.remove('ion')
        return sorted(keys)

    def load_from_line(self, arg):
        '''
        Read values of orbital contribution from the line.  (Essentially, this method acts as a parser.) The order of the orbitals is determined with the orbitalname array.  In PROCAR file, the table that corresponds the orbital contribution is like this::

          ion    s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot
          1  0.001  0.000  0.001  0.000  0.000  0.000  0.000  0.000  0.000  0.002
          2  0.003  0.000  0.001  0.000  0.000  0.000  0.000  0.000  0.000  0.005
          3  0.002  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.004
          4  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.001
          5  0.005  0.000  0.001  0.000  0.000  0.000  0.000  0.000  0.000  0.007
          6  0.003  0.000  0.001  0.000  0.000  0.000  0.001  0.000  0.000  0.005
          ...


          An array called orbitalname is created from the first line of the table 
          that starts with "ion"


        :param [String,Array,Hash,Orbital] arg:

        the format is like this::

           1  0.001  0.000  0.001  0.000  0.000  0.000  0.000  0.000  0.000  0.002
                #   (the first column indicates the atom #.)

           (Each element is separated by the space.)

        :return: atom #.
        :rtype:  string
        '''
        if isinstance(arg, str):
            self.orbital = dict(zip(self.defined_orbital_list(),
                                    (float(i) for i in arg.split())))
        elif isinstance(arg, (list, tuple)):
            self.orbital = dict(zip(self.defined_orbital_list(), arg))
        elif isinstance(arg, dict):
            if 'ion' in arg:
                self.orbital = arg
            else:
                raise ValueError(
                "dict() of {0} cannot be Orbital object.".format(arg))
        elif isinstance(arg, Orbital):
            self.orbital = arg.todict()
        else:
            raise ValueError(
            "Cannot create Orbital object from {0}() object.".format(type(arg)))
        self.orbital['ion'] = int(self.orbital['ion'])

    def __iter__(self):
        return iter(self.orbital)

    def names(self):
        return self.orbital.keys()

    def contributions(self):
        return self.orbital.values()

    def items(self):
        return self.orbital.items()

    def defined_orbital_list(self):
        return self.__orbital_list[:]

    def __getitem__(self, orbital_symbol):
        '''x.__getitem__(i) <==> x[i]

        # Returns contribution of specified orbital.
        # @param [Symbol] orbital_symbol Orbital name by symbol format
        #     (ex.):s, :sp.  if orbitalsymol is not exist in the object, it
        #     tries to generate the value of contribution.  
        # @return [Float] contribution of the specified orbital.
        '''
        if orbital_symbol in self.orbital:
            return self.orbital[orbital_symbol]
        else:
            return getattr(self, orbital_symbol)

    def __setitem__(self, orbital_symbol, value):
        '''x.__setitem__(i, y) <==> x[i]=y

        # Associate the *value* of contribution with the specified *orbital*.
        # @param [Symbol] orbitalsymbol Orbital name by symbol format
        #  (ex.):s, :sp.
        # @param [Float] value of contribution for the specified orbital.
        '''
        self.orbital[orbital_symbol] = value

    def get(self, symbol, default=None):
        '''Orbital.get(s[, D]) -> Orbital[s] if exist, else D'''
        return self.orbital.get(symbol, default)

    def setdefault(self, symbol, default=None):
        '''X.setdefault(s[, D]) -> X.get(s, D), also set X[k]=D if k not in X
        '''
        return self.orbital.setdefault(symbol, default)

    def site(self): return self.orbital['ion']

    #def orbital(self): return self

    def todict(self): return self.orbital

    @property
    def s(self): return self.orbital['s']

    @property
    def px(self): return self.orbital['px']

    @property
    def py(self): return self.orbital['py']

    @property
    def pz(self): return self.orbital['pz']

    @property
    def dxy(self): return self.orbital['dxy']

    @property
    def dyz(self): return self.orbital['dyz']

    @property
    def dz2(self): return self.orbital['dz2']

    @property
    def dxz(self): return self.orbital['dxz']

    @property
    def dx2(self): return self.orbital['dx2']

    @property
    def p(self):
        return self.setdefault('p', self.px + self.py + self.pz)

    @property
    def pxpy(self):
        return self.setdefault('pxpy', self.px + self.py)

    pypx = pxpy
    
    @property
    def pypz(self):
        return self.setdefault('pypz', self.py + self.pz)

    pzpy = pypz

    @property
    def pzpx(self):
        return self.setdefault('pzpx', self.pz + self.px)

    pxpz = pzpx

    @property
    def sp(self):
        return self.setdefault('sp', self.s + self.p)

    @property
    def d(self):
        return self.setdefault('d', self.dxy + self.dyz + 
                                    self.dz2 + self.dxz + self.dx2)

    @property
    def tot(self): return self.orbital['tot']

    spd = total = tot

    def is_empty(self):
        '''  # @return [Boolean] true if #Orbital is empty.
        '''
        return len(self.orbital) == 0

    def __add__(self, other):
        other = _convert_other_orbital(other)
        if other is NotImplemented: return other
        dest = Orbital()
        if self.is_empty():
            dest = copy.copy(other)
        else:
            for orbital_symbol, value in self.orbital.items():
                if orbital_symbol == 'ion':
                    dest[orbital_symbol] = (str(value) + ',' +
                                            str(other['ion']))
                else:
                    dest[orbital_symbol] = value + other[orbital_symbol]
        return dest

    def redefine_orbital_list(self, arg):
        '''Redefine orbital list by arg
        #   The machine default is set as above.
        # @param [Array] arg Symbol list consisting orbital name.
        #   First element must be :ion
        #   Final element must be :tot
        '''
        if not isinstance(arg, (list, tuple)):
            raise RuntimeError("redefine_orbital_list fail")
        if arg[0] != 'ion': 
            raise ValueError("First element of argument mst be 'ion'.")
        if arg[-1] != 'tot':
            raise ValueError("Last element of argument mst be 'tot'.")
        self.__orbital_list = arg

class State(Orbital):
    '''# Class for electronic state
    # @author Ryuichi Arafune
    # @version 2.0
    '''
    #  attr_accessor :orbital
    #  attr_reader :kindex, :bandindex, :spininfo, :eigenvalue
    def __init__(self, kindex=None, bandindex=None,
                 eigenvalue=None, spininfo=None, orbital=None):
        super(State, self).__init__(orbital)
        self.__kindex = kindex
        self.__bandindex = bandindex
        self.__eigenvalue = eigenvalue
        self.__spininfo = spininfo

    @property
    def kindex(self):
        return self.__kindex

    @property
    def bandindex(self):
        return self.__bandindex

    @property
    def eigenvalue(self):
        return self.__eigenvalue

    @property
    def spininfo(self):
        return self.__spininfo

    def __copy__(self):
        'x.__copy__() <==> copy.copy(x)'
        dest = State()
        dest.orbital = copy.copy(self.orbital)
        dest.__kindex = self.kindex
        dest.__bandindex = self.bandindex
        dest.__eigenvalue = self.eigenvalue
        dest.__spininfo = copy.copy(self.spininfo)
        return dest

    # python 3.x don't have cmp(),
    # python 2.6 or older don't have functools.total_odering
    # => define all rich comparision method.

    def __cmp__(self, other): # Version safety
        '''x.__cmp__(y) <==> cmp(x, y) # python 2.x
        '''
        if not isinstance(other, State): return NotImplemented
        cmp = np.sign(self.bandindex - other.bandindex)
        if cmp != 0: return cmp
        cmp = np.sign(self.kindex - other.kindex)
        if cmp != 0: return cmp
        cmp = np.sign(self.site - other.site)
        if cmp != 0: return cmp
        return np.sign(self.spininfo - other.spininfo)

    def __eq__(self, other): # Version safety
        'x.__eq__(y) <==> x==y'
        cmp = self.__cmp__(other)
        if cmp is NotImplemented: return cmp
        else: return cmp == 0

    def __ne__(self, other): # Version safety
        'x.__ne__(y) <==> x!=y'
        return self.__cmp__(other) != 0

    def __lt__(self, other): # Version safety
        'x.__lt__(y) <==> x<y'
        cmp = self.__cmp__(other)
        if cmp is NotImplemented: return cmp
        return cmp == -1

    def __le__(self, other): # Version safety
        'x.__le__(y) <==> x<=y'
        cmp = self.__cmp__(other)
        if cmp is NotImplemented: return cmp
        return cmp != 1

    def __gt__(self, other): # Version safety
        'x.__gt__(y) <==> x>y'
        cmp = self.__cmp__(other)
        if cmp is NotImplemented: return cmp
        return cmp == 1

    def __ge__(self, other): # Version safety
        'x.__ge__(y) <==> x>=y'
        cmp = self.__cmp__(other)
        if cmp is NotImplemented: return cmp
        return cmp != -1

    def fermilevel_correction(self, ef):
        self.__eigenvalue -= ef

    def __add__(self, other):
        '''x.__add__(y) <==> x+y

        #  eigenvalue must be identical with each other.
        #  In many cases, kvector and bandindex should be same, but not checked. 
        #  You may calculate k-integrated density of states from *band*...
        # @param [State] other'''
        if self.is_empty():
            dest = copy.copy(other)
        else:
            if self.eigenvalue != other.eigenvalue:
                raise ValueError("Eigenvalues are different.")
            dest = copy.copy(self)
            orbital = super().__add__(other)
            dest.orbital = orbital.todict()
        return dest

    def extract_orbitals(self, *orbital_symbols):
        '''
        #   Extract specified orbitals 
        #   (remove other orbital contributions)
        # @param [Array] orbital_symbols
        # @return [self] 
        '''
        orbital_symbols = tools.flatten(orbital_symbols)
        orbital_symbols.append('ion')
        self.orbital = dict((key, value) for key, value
                            in self.items() if key in orbital_symbols)
    
  # State#<<() :  
  #  eigenvalue must be identical with each other.
  #  In many cases, kvector and bandindex should be same, but not checked. 
  #  You may calculate k-integrated density of states from *band*...
  # @param [State] other
  #
  #  def <<(other)
  #    @@orbital_list[1..-1].each do orbital
  #      self[orbital]+= other[orbital]
  #    end
  #    self[:ion] = self[:ion].to_s
  #    self[:ion] << "_"<< other[:ion].to_s
  #  end
  #  self
  #end 
    
def _convert_other_band(other):
    if isinstance(other, Band):
        return other
    return NotImplemented

def _convert_other_orbital(other):
    if isinstance(other, Orbital):
        return other
    return NotImplemented

# -------------------------------------------------------------
#
# Main (test) routine
#
if __name__ == '__main__':
    import argparse
    from outcar import OUTCAR
    
    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--outcar', metavar='outcar_file',
                        help='''Use "OUTCAR" for the Fermi level correction.
outcar_file must be specified.
NOTE:  E-fermi of OUTCAR generated in
Band-calculation may NOT be reliable.''')
    parser.add_argument('--fermi', metavar='value', type=float,
                        help='''Fermi level correction
Energy shifts by this value
if --outcar is set, this option is ignored''')
    parser.add_argument('--site', metavar='atom_indices', dest='atomname',
                        action='append',
                        type=tools.parse_AtomselectionNum,
                        help='''atom index specifed with range.
Use "-" or ","
 (ex.) --site 1,2,7-9''')
    parser.add_argument('--as', metavar='name', nargs='+', dest='atomsetname',
                        action='append',
                        help='''the name of the sites identified
         by --site option
         the name is used in the title of the column''')
    parser.add_argument('--orbital', metavar='orbitals', action='append',
                        type=ft.partial(re.split, r'[,:]'),
                        help='''orbital name
         deliminated by ":" or ",".
         orbital names are:
         s, p, pxpy, pz, d, dxy, dyz, dz2, dxz, dx2
         (ex.) --orbital s:pxpy:d''')
    parser.add_argument('procar', metavar='PROCAR_file',
                        help='''PROCAR file''')
    
    args = parser.parse_args()

    # ---
    if not (len(args.atomname) == len(args.orbital) == len(args.atomsetname)):
        raise parser.error("--atom, --as and --orbital are mismatched.")
    # ---

    if args.outcar is not None:
        outcar = OUTCAR(args.outcar)
        fermi = outcar.fermi
    elif args.fermi is not None:
        fermi = args.fermi
    else:
        fermi = 0.0

    procar = PROCAR(args.procar)
    band = procar.to_band()
    
    if fermi != 0.0:
        band.fermilevel_correction(fermi)
    
    output_band = Band()
    for site, name, orbital in zip(args.atomname,
                                   args.atomsetname,
                                   args.orbital):
        tmpBand = band.site_integrate(site)
        tmpBand.rename_site(*name)
        tmpBand.extract_orbitals_in_place(orbital)
        output_band += tmpBand
    output_band.sort()
    print(output_band)
        
