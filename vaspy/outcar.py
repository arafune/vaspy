"""OUTCAR class."""


from __future__ import annotations

from typing import IO, TYPE_CHECKING

from vaspy.tools import open_by_suffix

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path


class OUTCAR:  # Version safety
    """Class for OUTCAR file.

    Attributes:
    ----------
    n_atom: int
        number of atom
    atom_types: list
        list of ion name
    atomnums: list
        list of number of ions
    k_vectors: list
        kvector
    weights: list
        weight list
    fermi: float
        fermi energy
    num_k: int
        number of k points
    n_bands: int
        number of bands
    posforce : list
        evolution of position and force
    magnetization: list
        evolution of magnetization
    total_charge: list
        evolution of total charge


    Todo:
    ----
    posforce should be divided to positions and forces.

    """

    def __init__(self, filename: str | Path = "") -> None:
        """Initialize."""
        self.n_atom = 0
        self.atom_types: list[str] = []
        self.atomnums: list[int] = []
        self.posforce: list[list[list[float]]] = []
        self.posforce_title: list[list[str]] = []
        self.atom_names: list[str] = []
        self.fermi = 0.0
        self.site_label: list[str] = []
        self.num_k = 0
        self.nkdim = 0
        self.n_bands = 0
        self.magnetizations: list[list[list[float]]] = []
        self.total_charges: list[list[list[float]]] = []
        self.k_vectors: list[list[float]] = []
        self.weights: list[float] = []
        self.totens: list[float] = []
        if filename:
            self.load_file(open_by_suffix(str(filename)))

    def set_atom_names(self) -> list[str]:
        """Build atom_names (the list of atomname_with_index)."""
        self.atom_names = []
        for elm, ionnum in zip(self.atom_types, self.atomnums, strict=True):
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

    def load_file(self, the_file: IO[str]) -> None:
        """Parse OUTCAR file.

        Parameters
        ----------
        the_file: StringIO
            "OUTCAR" file

        """
        # local variables
        section: list[str] = []
        posforce: list[list[float]] = []
        magnetizations: list[list[float]] = []
        total_charges: list[list[float]] = []
        kvec_weight = []
        # parse
        for line in the_file:
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
                    if self.n_atom == 1:
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
                    if self.n_atom == 1:
                        section.pop()
            elif section == ["kvec_weight"]:
                if len(line) > 3:
                    kvec_weight.append([float(x) for x in line.split()])
                else:
                    section.pop()
            else:
                if "number of dos" in line:
                    self.n_atom = int(line.split()[-1])
                elif "TITEL  =" in line:
                    self.atom_types.append(line.split()[3])
                elif "ions per type " in line:
                    self.atomnums = [int(x) for x in line.split()[4:]]
                elif "POSITION" in line and "TOTAL-FORCE" in line:
                    section.append("force")
                elif "E-fermi" in line:
                    self.fermi = float(line.split()[2])
                elif "NBANDS=" in line:
                    self.num_k = int(line.split()[3])
                    self.nkdim = int(line.split()[9])
                    self.n_bands = int(line.split()[14])
                elif "free energy    TOTEN  =" in line:
                    self.totens.append(float(line.split("=")[-1].split()[0]))
                elif "reciprocal lattice vectors" in line:
                    self.recvec = [
                        [float(_) for _ in next(the_file)[43:].split()]
                        for _ in range(3)
                    ]
                elif " magnetization (x)" in line:
                    magnetizations = []
                    section.append("magnetization")
                elif " total charge     " in line:
                    total_charges = []
                    section.append("total_charge")
                elif " Following reciprocal coordinates:" in line:
                    next(the_file)
                    kvec_weight = []
                    section.append("kvec_weight")
                else:
                    pass
        self.site_label = [
            name + ":#" + str(index + 1)
            for (index, name) in enumerate(
                [
                    elm + str(j)
                    for (elm, n) in zip(self.atom_types, self.atomnums)
                    for j in range(1, int(n) + 1)
                ],
            )
        ]
        self.posforce = [
            posforce[i : i + self.n_atom] for i in range(0, len(posforce), self.n_atom)
        ]
        self.set_atom_names()
        self.set_posforce_title()
        for i in kvec_weight:
            self.k_vectors.append([i[0], i[1], i[2]])
            self.weights.append(i[3])
        the_file.close()

    def select_posforce_header(
        self,
        posforce_flag: Sequence[bool],
        *sites: tuple[int | list[int]],
    ) -> list[str]:
        """Return the position and force header selected."""
        selected_sites: Sequence[int]
        if sites == () or sites[0] == []:
            selected_sites = range(1, self.n_atom + 1)
        if isinstance(sites[0], list | tuple):
            selected_sites = list(sites[0])
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
        self,
        posforce_flag: Sequence[bool],
        *sites: tuple[int | list[int]],
    ) -> list[list[float]]:
        """Return the position and force selected by posforce_flag.

        Notes
        -----
        posforce_flag: An 6-element True/False list that indicates
                the output (ex.) [True, True, False, True, True, False]

        """
        if sites == () or sites[0] == []:
            selected_sites: Sequence[int] = range(1, self.n_atom + 1)
        if isinstance(sites[0], list | tuple):
            selected_sites = list(sites[0])
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
