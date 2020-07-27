#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module provides OUTCAR class
"""


from __future__ import unicode_literals  # Version safety
from __future__ import annotations

from vaspy.tools import open_by_suffix
from typing import List, Optional, Sequence, Union, IO


class OUTCAR(object):  # Version safety
    """Class for OUTCAR file.

    Attributes
    -------------
    natom: int
        number of atom
    atomtypes: list
        list of ion name
    atomnums: list
        list of number of ions
    kvecs: list
        kvector
    weights: list
        weight list
    fermi: float
        fermi energy
    numk: int
        number of k points
    nbands: int
        number of bands
    posforce : list
        evolution of position and force
    magnetization: list
        evolution of magnetization
    total_charge: list
        evolution of total charge


    Todo
    ------
    posforce should be devided to positions and forces.

    """

    def __init__(self, filename: Optional[str] = None) -> None:
        """Initialize."""
        self.natom = 0
        self.atomtypes: List[str] = []
        self.atomnums: List[int] = []
        self.posforce: List[List[List[float]]] = []
        self.posforce_title: List[List[str]] = []
        self.atom_names: List[str] = []
        self.fermi = 0.0
        self.site_label: List[str] = []
        self.numk = 0
        self.nkdim = 0
        self.nbands = 0
        self.magnetizations: List[List[List[float]]] = []
        self.total_charges: List[List[List[float]]] = []
        self.kvecs: List[List[float]] = []
        self.weights: List[float] = []
        if filename:
            self.load_file(open_by_suffix(filename))

    def set_atom_names(self) -> List[str]:
        """Build atom_names (the list of atomname_with_index)."""
        self.atom_names = []
        for elm, ionnum in zip(self.atomtypes, self.atomnums):
            for j in range(1, ionnum + 1):
                tmp = elm + str(j)
                if tmp not in self.atom_names:
                    self.atom_names.append(tmp)
                else:
                    while tmp in self.atom_names:
                        j = j + 1
                        tmp = elm + str(j)
                    else:
                        self.atom_names.append(tmp)
        return self.atom_names

    def set_posforce_title(self) -> None:
        """Build posforce_title."""
        self.set_atom_names()
        self.posforce_title = [
            [i + "_x", i + "_y", i + "_z", i + "_fx", i + "_fy", i + "_fz"]
            for i in self.atom_names
        ]

    def load_file(self, thefile: Union[IO[str]]) -> None:
        """Parse OUTCAR file.

        Parameters
        ----------
        thefile: StringIO
            "OUTCAR" file

        """
        # local variables
        section: List[str] = []
        posforce = []
        magnetizations: List[List[float]] = []
        total_charges: List[List[float]] = []
        kvec_weight = []
        # parse
        for line in thefile:
            if section == ["force"]:
                if "total drift" in line:
                    section.pop()
                elif "---------------" in line:
                    pass
                else:
                    posforce.append([float(x) for x in line.split()])
            elif section == ["magnetization"]:
                if "---------------------------------" in line:
                    pass
                elif "# of ion" in line:
                    pass
                elif "tot    " in line:
                    self.magnetizations.append(magnetizations)
                    section.pop()
                elif len(line) == 2:
                    pass
                else:
                    magnetizations.append([float(x) for x in line.split()[1:4]])
                    if self.natom == 1:
                        section.pop()
            elif section == ["total_charge"]:
                if "---------------------------------" in line:
                    pass
                elif "# of ion" in line:
                    pass
                elif "tot    " in line:
                    self.total_charges.append(total_charges)
                    section.pop()
                elif len(line) == 2:
                    pass
                else:
                    total_charges.append([float(x) for x in line.split()[1:4]])
                    if self.natom == 1:
                        section.pop()
            elif section == ["kvec_weight"]:
                if len(line) > 3:
                    kvec_weight.append([float(x) for x in line.split()])
                else:
                    section.pop()
            else:
                if "number of dos" in line:
                    self.natom = int(line.split()[-1])
                elif "TITEL  =" in line:
                    self.atomtypes.append(line.split()[3])
                elif "ions per type " in line:
                    self.atomnums = [int(x) for x in line.split()[4:]]
                elif "POSITION" in line and "TOTAL-FORCE" in line:
                    section.append("force")
                elif "E-fermi" in line:
                    self.fermi = float(line.split()[2])
                elif "NBANDS" in line:
                    self.numk = int(line.split()[3])
                    self.nkdim = int(line.split()[9])
                    self.nbands = int(line.split()[14])
                elif "reciprocal lattice vectors" in line:
                    self.recvec = [
                        [float(i) for i in next(thefile)[43:].split()] for i in range(3)
                    ]
                elif " magnetization (x)" in line:
                    magnetizations = []
                    section.append("magnetization")
                elif " total charge     " in line:
                    total_charges = []
                    section.append("total_charge")
                elif " Following reciprocal coordinates:" in line:
                    next(thefile)
                    kvec_weight = []
                    section.append("kvec_weight")
                else:
                    pass
        self.site_label = [
            name + ":#" + str(index + 1)
            for (index, name) in enumerate(
                [
                    elm + str(j)
                    for (elm, n) in zip(self.atomtypes, self.atomnums)
                    for j in range(1, int(n) + 1)
                ]
            )
        ]
        self.posforce = [
            posforce[i : i + self.natom] for i in range(0, len(posforce), self.natom)
        ]
        self.set_atom_names()
        self.set_posforce_title()
        for i in kvec_weight:
            self.kvecs.append([i[0], i[1], i[2]])
            self.weights.append(i[3])
        thefile.close()

    def select_posforce_header(
        self, posforce_flag: Sequence[bool], *sites: int
    ) -> List[str]:
        """Return the position and force header selected."""
        selected_sites: Sequence[int]
        if sites == () or sites[0] == []:
            selected_sites = range(1, self.natom + 1)
        if isinstance(sites[0], (list, tuple)):
            selected_sites = [n for n in sites[0]]
        return [
            posforce
            for (index, site) in enumerate(self.posforce_title)
            for (posforce, boolian) in zip(site, posforce_flag)
            if boolian and (index + 1 in selected_sites)
        ]

    # return [posforce for (posforce, boolian) in zip(ithAtom, poforce_flag)
    # if boolian==True for ithAtom in self.posforce_title  ] #which is
    # correct?

    def select_posforce(
        self, posforce_flag: Sequence[bool], *sites: int
    ) -> List[List[float]]:
        """Return the position and force selected by posforce_flag.

        Notes
        -------
        posforce_flag: An 6-element True/False list that indicates
                the output (ex.) [True, True, False, True, True, False]

        """
        if sites == () or sites[0] == []:
            selected_sites = range(1, self.natom + 1)
        if isinstance(sites[0], (list, tuple)):
            selected_sites = [n for n in sites[0]]
        return [
            [
                posforce
                for (index, site) in enumerate(one_cycle)
                for (posforce, boolian) in zip(site, posforce_flag)
                if boolian
                if index + 1 in selected_sites
            ]
            for one_cycle in self.posforce
        ]
