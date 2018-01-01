#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module provides OUTCAR class
"""

from __future__ import unicode_literals  # Version safety
from __future__ import print_function  # Version safety
import os.path
import bz2


class OUTCAR(object):  # Version safety
    '''
    Class for OUTCAR file that stores calculation details
    and/or calculation progress

    Attributes
    -------------

    ntom: int
        number of atom
    iontypes: list
        list of ion name
    ionnums: list
        list of number of ions
    kvecs: list
        kvector
    weights: list
        weight list
    '''

    def __init__(self, arg=None):
        self.natom = 0
        self.iontypes = []
        self.ionnums = []
        self.posforce = []
        self.posforce_title = []
        self.atom_names = []
        self.fermi = 0.0
        self.atom_identifer = []
        self.numk = 0
        self.nkdim = 0
        self.nbands = 0
        self.magnetization = []
        self.tot_chage = []
        self.kvecs = []
        self.weights = []
        if arg is not None:
            self.load_file(arg)

    def set_atom_names(self):
        '''
        build atom_names (the list of atomname_with_index)
        '''
        self.atom_names = []
        for elm, ionnum in zip(self.iontypes, self.ionnums):
            for j in range(1, ionnum + 1):
                tmp = elm + str(j)
                if tmp not in self.atom_names:
                    self.atom_names.append(tmp)
                else:
                    #                    jj = j
                    while tmp in self.atom_names:
                        j = j + 1
                        tmp = elm + str(j)
                    else:
                        self.atom_names.append(tmp)
        return self.atom_names

    def set_posforce_title(self):
        """
        build posforce_title
        """
        self.set_atom_names()
        self.posforce_title = [[i + "_x",
                                i + "_y",
                                i + "_z",
                                i + "_fx",
                                i + "_fy",
                                i + "_fz", ]
                               for i in self.atom_names]

    def load_file(self, arg):
        '''
        Effectively, this is a constructor of OUTCAR object.

        Parameters
        ----------

        arg: str
            File name of "OUTCAR"
        '''
        # local variables
        section = []
        posforce = []
        magnetization = []
        kvec_weight = []
        # parse
        if os.path.splitext(arg)[1] == '.bz2':
            try:
                thefile = bz2.open(arg, mode='rt')
            except AttributeError:
                thefile = bz2.BZ2File(arg, mode='r')
        else:
            thefile = open(arg)
        for line in thefile:
            if section == ["force"]:
                if "total drift" in line:
                    section.pop()
                elif "---------------" in line:
                    pass
                elif "total drift:" in line:
                    section.pop()
                else:
                    posforce.append([float(x) for x in line.split()])
            elif section == ["magnetization"]:
                if "---------------------------------" in line:
                    pass
                elif "# of ion" in line:
                    pass
                elif "tot    " in line:
                    self.magnetization.append(magnetization)
                    section.pop()
                elif len(line) == 2:
                    pass
                else:
                    magnetization.append([
                        float(x) for x in line.split()[1:4]])
            elif section == ['kvec_weight']:
                if len(line) > 3:
                    kvec_weight.append(
                        [float(x) for x in line.split()])
                else:
                    section.pop()
            else:
                if "number of dos" in line:
                    self.natom = int(line.split()[-1])
                elif "TITEL  =" in line:
                    self.iontypes.append(line.split()[3])
                elif "ions per type " in line:
                    self.ionnums = [int(x) for x in line.split()[4:]]
                elif "POSITION" in line and "TOTAL-FORCE" in line:
                    section.append("force")
                elif "E-fermi" in line:
                    self.fermi = float(line.split()[2])
                elif "NBANDS" in line:
                    self.numk = int(line.strip().split()[3])
                    self.nkdim = int(line.strip().split()[9])
                    self.nbands = int(line.strip().split()[14])
                elif "reciprocal lattice vectors" in line:
                    self.recvec = [[float(i) for i in
                                    next(thefile).strip().split()[3:]]
                                   for i in range(3)]
                elif " magnetization (x)" in line:
                    magnetization = []
                    section.append("magnetization")
                elif " Following reciprocal coordinates:" in line:
                    next(thefile)
                    kvec_weight = []
                    section.append('kvec_weight')
                else:
                    pass
        self.atom_identifer = [name + ":#" + str(index + 1)
                               for (index, name)
                               in enumerate(
                                   [elm + str(j)
                                    for (elm, n) in zip(self.iontypes,
                                                        self.ionnums)
                                    for j in range(1, int(n) + 1)])]
        self.posforce = [posforce[i:i + self.natom]
                         for i in range(0, len(posforce), self.natom)]
        self.set_atom_names()
        self.set_posforce_title()
        for i in kvec_weight:
            self.kvecs.append([i[0], i[1], i[2]])
            self.weights.append(i[3])

    def select_posforce_header(self, posforce_flag, *sites):
        '''
        '''
        if sites == () or sites[0] == []:
            sites = range(1, self.natom + 1)
        if isinstance(sites[0], (list, tuple)):
            sites = [n for n in sites[0]]
        return [posforce for (index, site)
                in enumerate(self.posforce_title)
                for (posforce, boolian) in zip(site, posforce_flag)
                if boolian and (index + 1 in sites)]
# return [posforce for (posforce, boolian) in zip(ithAtom, poforce_flag)
# if boolian==True for ithAtom in self.posforce_title  ] #which is
# correct?

    def select_posforce(self, posforce_flag, *sites):
        '''
        Return the posforce corresponding the posforce_flag

        Note
        -------
        posforce_flag: An 6-element True/False list that indicates \
                  the output (ex.) [True, True, False, True, True, False]
        '''
        if sites == () or sites[0] == []:
            sites = range(1, self.natom + 1)
        if isinstance(sites[0], (list, tuple)):
            sites = [n for n in sites[0]]
        return [[posforce for (index, site) in enumerate(one_cycle)
                 for (posforce, boolian) in zip(site, posforce_flag)
                 if boolian
                 if index + 1 in sites]
                for one_cycle in self.posforce]
