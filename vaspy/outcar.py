#! /usr/bin/env python
# -*- coding: utf-8 -*-
""".. py:module:: outcar

This module provides OUTCAR class
"""

from __future__ import unicode_literals  # Version safety
from __future__ import print_function  # Version safety


class OUTCAR(object):  # Version safety
    '''.. py:class:: OUTCAR(outcar_file)

    Class for OUTCAR file that stores calculation details
    and/or calculation progress
    '''

    def __init__(self, arg=None):
        self.nions = 0
        self.iontype = []
        self.ionnums = []
        self.posforce = []
        self.posforce_title = []
        self.atom_names = []
        self.fermi = 0.0
        self.atom_identifer = []
        if arg is not None:
            self.load_from_file(arg)

    def set_atom_names(self):
        '''.. py:method:: set_atom_names()

        build atom_names (the list of atomname_with_index)
        '''
        self.atom_names = []
        for elm, ionnum in zip(self.iontype, self.ionnums):
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
        """.. py:method:: set_posforce_title()

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

    def load_from_file(self, arg):
        ''' .. py:method:: load_from_file(outcar_filename)

        Effectively, this is a constructor of OUTCAR object.

        :param arg: File name of "OUTCAR"
        '''
        # local variables
        section = []
        posforce = []
        # parse
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
            else:
                if "number of dos" in line:
                    self.nions = int(line.split()[-1])
                elif "TITEL  =" in line:
                    self.iontype.append(line.split()[3])
                elif "ions per type " in line:
                    self.ionnums = [int(x) for x in line.split()[4:]]
                elif "POSITION" in line and "TOTAL-FORCE" in line:
                    section.append("force")
                elif "E-fermi" in line:
                    self.fermi = float(line.split()[2])
                else:
                    pass
        self.atom_identifer = [name + ":#" + str(index + 1)
                               for (index, name)
                               in enumerate(
                                   [elm + str(j)
                                    for (elm, n) in zip(self.iontype,
                                                        self.ionnums)
                                    for j in range(1, int(n) + 1)])]
        self.posforce = [posforce[i:i + self.nions]
                         for i in range(0, len(posforce), self.nions)]
        self.set_atom_names()
        self.set_posforce_title()

    def select_posforce_header(self, posforce_flag, *sites):
        '''.. py:method:: select_posforce_header(posforce_flag, *sites)        
'''
        if sites == () or sites[0] == []:
            sites = range(1, self.nions + 1)
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
        '''.. py:method:: select_posforce(posforce_flag, *sites)

        Return the posforce corresponding the posforce_flag

        .. note:: posforce_flag: An 6-element True/False list that indicates
                  the output (ex.) [True, True, False, True, True, False]
        '''
        if sites == () or sites[0] == []:
            sites = range(1, self.nions + 1)
        if isinstance(sites[0], (list, tuple)):
            sites = [n for n in sites[0]]
        return [[posforce for (index, site) in enumerate(one_cycle)
                 for (posforce, boolian) in zip(site, posforce_flag)
                 if boolian
                 if index + 1 in sites]
                for one_cycle in self.posforce]
