# -*- coding: utf-8 -*-
# translate from procar.rb of scRipt4VASP 2014/2/26 master branch
"""
This module provides PROCAR, ProjectionBand classes.

* PROCAR class is used to stored the PROCAR information in the memory.
* Band_with_Projection class is used in vaspy-procar.py script.
    *  This class is essentially used on ipython, but not so easy for use.
* EnergyBand class is used for drawing the energy band.
* Projection class is used for storing the orbital projection data
"""

import csv
import re
from logging import DEBUG, INFO, Formatter, StreamHandler, getLogger
from typing import IO, List, Optional, Sequence, Tuple, Union

import numpy as np

from vaspy.eigenval import EnergyBand
from vaspy.tools import open_by_suffix

# logger
LOGLEVEL = INFO
logger = getLogger(__name__)
fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
formatter = Formatter(fmt)
handler = StreamHandler()
handler.setLevel(LOGLEVEL)
logger.setLevel(LOGLEVEL)
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False


class ProjectionBand(EnergyBand):
    """Class for band structure including orbital projection.

    This is the child class of eigenval.EnergyBand which represents
    the band structure, to represent the orbital.

    Attributes
    -----------
    natom: int
        number of atoms
    proj: numpy.ndarray
        Orbital projection. proj[spin_i, k_i, band_i, site_i, orbital_i]
    phase: numpy.ndarray
        Phase data.  phase[spin_i, k_i, band_i, site_i, orbital_i]

    """

    def __init__(self, kvecs=(), energies=(), proj=(), phase=(), nspin=1) -> None:
        """Initialize."""
        super(ProjectionBand, self).__init__()
        self.natom = 0
        self.proj = proj
        self.phase = phase
        self.kvecs = []

    def append_sumsite(self, sites, site_name):
        """Append site-sum results.

        After this method, shape changes as following
        self[nspin, numk, nbands, natom + 1, norbital]

        Parameters
        ------------
        sites: tuple
            site index for summention
        site_name: str
            label name for summed site, such as 'silicene', and 'SurfaceAu'

        """
        # As the original label['site'] is just a number beginnig from zero
        logger.debug("self.label: {}".format(self.label))
        logger.debug("self.natom: {}".format(self.natom))
        logger.debug("sites: {}".format(sites))

        if len(sites) == 1 and sites[0] in self.label["site"][: self.natom]:
            self.label["site"][sites[0]] = site_name
            logger.debug("self.label: {}".format(self.label))
        if site_name in self.label["site"]:
            return None
        self.label["site"].append(site_name)
        logger.debug("self.label: {}".format(self.label))
        #    spin, k, band, atom
        sumsite = self.proj[:, :, :, sites, :].sum(axis=-2, keepdims=True)
        self.proj = np.concatenate((self.proj, sumsite), axis=-2)
        return sumsite

    def append_sumorbital(self, orbitals, orbital_name):
        """Append orbital-sum results.

        After this method, shape changes as following
        self[nspin, numk, nbands, natom, norbital + 1]

        Parameters
        -----------
        orbitals: tuple
            orbital index for summention
        oribtal_name : str
            label name for summed orbital, such as 'p' and 'sp'

        """
        if orbital_name in self.label["orbital"]:
            return None
        self.label["orbital"].append(orbital_name)
        #    spin, k, band, atom
        sumorbital = self.proj[:, :, :, :, orbitals].sum(axis=-1, keepdims=True)
        self.proj = np.concatenate((self.proj, sumorbital), axis=-1)
        return sumorbital

    def orbital_index(self, arg):
        """Return the indexes corresponding orbital names.

        This method returns the tuple of orbital number in
        self.label['orbital'].
        (i.e. self.label['orbital'].index(orbitalname).  If the
        orbital name has not been in self.orbital_names (i.e. if the
        orbital name is not used in PROCAR file) but the orbital name
        is appropriate as the composed orbital ((ex.) sp, pxpy),
        returns the indexes of the orbitals to be composed as the
        tuple.


        Parameters
        ----------
        arg: str
            name of (composed) orbital


        Returns
        --------
        tuple
            number corresponding to (composed) orbital name.

        """
        orbital_names = self.label["orbital"]
        orbital_name = check_orbital_name(arg)
        if orbital_name in orbital_names:
            return (orbital_names.index(orbital_name),)
        elif orbital_name == "p":
            return (
                orbital_names.index("px"),
                orbital_names.index("py"),
                orbital_names.index("pz"),
            )
        elif orbital_name == "d":
            return (
                orbital_names.index("dxy"),
                orbital_names.index("dyz"),
                orbital_names.index("dx2"),
                orbital_names.index("dxz"),
                orbital_names.index("dz2"),
            )
        elif orbital_name == "sp":
            return (
                orbital_names.index("s"),
                orbital_names.index("px"),
                orbital_names.index("py"),
                orbital_names.index("pz"),
            )
        elif orbital_name == "pxpy":
            return (orbital_names.index("px"), orbital_names.index("py"))
        elif orbital_name == "pypz":
            return (orbital_names.index("py"), orbital_names.index("pz"))
        elif orbital_name == "pxpz":
            return (orbital_names.index("px"), orbital_names.index("pz"))
        else:
            err = str(orbital_name) + " is not a proper (composed) orbital name."
        raise RuntimeError(err)

    def make_label(self, site_indexes=None, orbital_indexes_sets=None):
        """Return array the used for **label** for CSV-like data.

        Parameters
        ----------
        site_indexes: tuple
            key tuple used for label
        orbital_indexes_sets: tuple
            index tuple for output

        """
        label_list = super(ProjectionBand, self).make_label("k", "energy")
        if not (site_indexes and orbital_indexes_sets):
            return label_list
        for site_i, orbitals in zip(site_indexes, orbital_indexes_sets):
            for orbital_i in orbitals:
                for spin_i in self.label["spin"]:
                    try:
                        label = (
                            str(self.label["site"][site_i])
                            + spin_i
                            + "_"
                            + self.label["orbital"][orbital_i]
                        )
                    except TypeError:
                        label = (
                            str(self.label["site"][site_i])
                            + spin_i
                            + "_"
                            + self.label["orbital"][orbital_i[0]]
                        )
                    label_list.append(label)
        return label_list

    def to_3dlist(
        self,
        site_indexes: Sequence[int] = (),
        orbital_indexes_sets: Sequence[Sequence[int]] = (),
    ):
        """Return 3D list data that are easily converted to txt data for csv.

        Parameters
        ------------
        site_indexes: list or tuple that contains int
            site name for output  (the elements must be in self.label['site'])
            e.g., (3, 5)
        orbitals: list or tuple that contains list or tuple of int
            tuple (list)  of tuple (list) for output
            e.g., ((1 ,5 , 11), (0, 3))

        Returns
        --------
        list:
            3D list, the first dimenstion corresponds to the data for band_i
            and each band_i contains kdistances, energy (or energies), and
            orbital data.

        """
        assert len(site_indexes) == len(
            orbital_indexes_sets
        ), "must len(sites)==len(orbitals)"
        if not (site_indexes and orbital_indexes_sets):
            return super(ProjectionBand, self).to_3dlist()
        # make np.array named projband from proj by using
        # site_indexes, orbital_indexes_sets values:
        # projband[spin_i, k_i, band_i, site_orbital]
        # Note that projband.ndim is 4
        array_list = []
        for site_i, orbitals in zip(site_indexes, orbital_indexes_sets):
            for orbital_i in orbitals:
                array_list.append(
                    self.proj[:, :, :, site_i, orbital_i][:, :, :, np.newaxis]
                )
        # このtranspose で、
        # A_mT_orb1, A_mX_orb1, A_mY_orb1, A_mZ_orb1, サイトB_mT_orb2...
        # というフォーマットになる。
        output_proj = (
            np.concatenate(tuple(array_list), axis=-1)
            .transpose(2, 1, 3, 0)
            .reshape(self.nbands, self.numk, -1)
        )
        kvalues = np.asarray((self.kdistances.tolist()) * self.nbands)[
            :, np.newaxis
        ].reshape(self.nbands, self.numk, 1)
        projband = np.concatenate((kvalues, self.energies.T, output_proj), axis=-1)
        return projband.tolist()

    def to_csv(
        self, csv_file, site_indexes=(), orbital_indexes_sets=(), blankline: bool = True
    ):
        """Write data to csv file.

        Parameters
        ----------
        csv_file: str
            filename for output
        label_str: str
            string for label (put it on the first line)
        blankline: boolean
            It True (default), the blank line is inserted between band data

        """
        label_str = (
            "\t".join(self.make_label(site_indexes, orbital_indexes_sets)) + "\n"
        )
        with open(csv_file, "w") as fhandle:
            fhandle.writelines(label_str)
            writer = csv.writer(fhandle, delimiter="\t")
            for band_i in self.to_3dlist(site_indexes, orbital_indexes_sets):
                writer.writerows(band_i)
                if blankline:
                    fhandle.writelines("\n")

    def text_sheet(
        self,
        site_indexes: Sequence[int] = (),
        orbital_indexes_sets: Sequence[Sequence[int]] = (),
    ):
        """Return csv-like text data.

        Parameters
        -----------
        site_indexes: list or tuple that contains int
            site name for output  (the elements must be in self.label['site'])
            e.g., (3, 5)
        orbital_indexes_sets: list or tuple that contains list or tuple of int
            tuple (list)  of tuple (list) for output
            e.g., ((1 ,5 , 11), (0, 3))

        Returns
        --------
        str

        """
        assert len(site_indexes) == len(
            orbital_indexes_sets
        ), "must len(site_indexes)==len(orbitals_indexes_sets)"
        output = "\t".join(self.make_label(site_indexes, orbital_indexes_sets))
        output += "\n"
        list3d = self.to_3dlist(site_indexes, orbital_indexes_sets)
        for band_i in list3d:
            for line in band_i:
                output += "{0: .8e}".format(line[0])
                for element in line[1:]:
                    output += "\t{0:.8e}".format(element)
                output += "\n"
            output += "\n"
        return output


class PROCAR(ProjectionBand):  # Version safety
    """Class for storing the data saved in PROCAR file.

    PROCAR consists of the following lines.  Appear once per file.

    1. The first line is used just as a comment.

        :Example:   PROCAR lm decomposed + phase

    2. Number of k-points, bands and ions.
        (Appear once when spin-integrated, twice when spin-resolved.)

        :Example:
            # of k-points:   50         # of bands: 576         # of ions:  98

    3.  k-point character

        :Example:
            k-point    1 :    0.00000 0.00000 0.00000 weight = 0.02000000

        .. Note::  That the first character must be "blank".

    4. Band character

        :Example:

        band   1 # energy  -11.87868466 # occ.  2.00000000

    5. orbital contribution.

        :Example:
        1 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000

    Attributes
    ----------
    orbital: list
        Store orbital character

    phase: list
        list (Not used at present)

    Parameters
    ----------
    filename: str
        Filename of *PROCAR* file
    phase_read: boolean
        True if read phase data (default is false)

    """

    def __init__(
        self, filename: Optional[str] = None, phase_read: bool = False
    ) -> None:
        """Initialize."""
        super(PROCAR, self).__init__()
        if filename:
            self.load_file(open_by_suffix(filename), phase_read)

    def load_file(
        self, thefile: Union[IO[str], IO[bytes]], phase_read: bool = False
    ) -> None:
        """Parse PROCAR.

        Parameters
        ----------
        thefile: StringIO
            'PROCAR' file
        phase_read: boolean
            Switch for loading phase characters

        """
        first_line = next(thefile)
        assert (
            "PROCAR lm decomposed + phase" in first_line
        ), "This PROCAR is not a proper format\n \
            Check your INCAR the calculations.\n"

        self.kvecs = []
        energies = []
        orbitals = ""
        phase_r = ""
        phase_i = ""
        self.label["orbital"] = []
        for line in thefile:
            if line.isspace():
                continue
            elif "k-points: " in line:
                self.numk, self.nbands, self.natom = [
                    int(i) for i in re.split(r"\s|:", line) if i.isdigit()
                ]
            elif "k-point " in line:
                try:
                    self.kvecs.append([float(i) for i in line.split()[3:6]])
                except ValueError:
                    self.kvecs.append(
                        [
                            np.float_(line[18:29]),
                            np.float_(line[29:40]),
                            np.float_(line[40:51]),
                        ]
                    )
            elif "band " in line:
                energies.append(float(line.split()[4]))
            elif "ion" in line:
                if "tot" in line and not self.label["orbital"]:
                    self.label["orbital"] = line.split()[1:]
                line = next(thefile)
                while "ion " not in line:
                    if "tot " not in line:
                        orbitals += line[3:]
                    line = next(thefile)
                if phase_read:
                    line = next(thefile)
                    for _ in range(self.natom):
                        try:
                            phase_r += line[3:]
                            line = next(thefile)
                            phase_i += line[3:]
                            line = next(thefile)
                        except StopIteration:
                            continue
        #
        self.kvecs = np.asarray(self.kvecs[: self.numk])
        self.proj = np.fromstring(orbitals, dtype=float, sep=" ")
        del orbitals
        norbital = len(self.label["orbital"])
        self.label["site"] = list(range(self.natom))
        self.nspin = self.proj.size // (self.numk * self.nbands * self.natom * norbital)
        if phase_read:
            self.phase = np.fromstring(phase_r, dtype=float, sep=" ") + (
                0 + 1.0j
            ) * np.fromstring(phase_i, dtype=float, sep=" ")
            del phase_r, phase_i
        #
        if self.nspin == 1:  # standard
            self.label["spin"] = [""]
            self.label["energy"] = ["Energy"]
            self.energies = np.asarray(energies).reshape(1, self.numk, self.nbands)
            self.proj = self.proj.reshape(
                (self.nspin, self.numk, self.nbands, self.natom, norbital)
            )
            if phase_read:
                self.phase = self.phase.reshape(
                    (self.nspin, self.numk, self.nbands, self.natom, norbital - 1)
                )
        elif self.nspin == 2:  # collinear
            self.label["spin"] = ["_up", "_down"]
            self.label["energy"] = ["Energy_up", "Energy_down"]
            self.energies = np.asarray(energies).reshape(2, self.numk, self.nbands)
            self.proj = self.proj.reshape(
                (self.nspin, self.numk, self.nbands, self.natom, norbital)
            )
            if phase_read:
                self.phase = self.phase.reshape(
                    (self.nspin, self.numk, self.nbands, self.natom, norbital - 1)
                )
        elif self.nspin == 4:  # non-collinear
            self.label["spin"] = ["_mT", "_mX", "_mY", "_mZ"]
            self.label["energy"] = ["Energy"]
            self.energies = np.asarray(energies).reshape(1, self.numk, self.nbands)
            self.proj = self.proj.reshape(
                self.numk, self.nbands, self.nspin, self.natom, norbital
            ).transpose((2, 0, 1, 3, 4))
            if phase_read:
                self.phase = self.phase.reshape(
                    (1, self.numk, self.nbands, self.natom, norbital - 1)
                )
        else:
            raise ValueError
        del energies
        thefile.close()

    def set_spin_character(self, phase_read: bool = False):
        """Set label of Energy and Spin character.

        Parameters
        ----------
        phase_read: boolean
            If True, treat the shape of self.phase

        """
        pass

    def __repr__(self) -> str:
        """__str__() <=> str(x).

        Show the PROCAR character, not contents.
        """
        template1 = """The properties of this procar:
    # of k-points: {0.numk}
    # of bands: {0.nbands}
    # of atoms: {0.natom}
    # of spin: {0.nspin}
    # of kvecs: {1}
    # of energies: {2}
    ((# of k-points) * (# of bands) = {0.numk}*{0.nbands}={3})
    # of orbital component: {4}
    ((# of k-points) * (# of bands) * (# of ions) =
        {0.numk}*{0.nbands}*{0.natom}={5})
    # of phase component: {6}"""
        string = ""
        for orb in self.label["orbital"]:
            string += "{0}  ".format(orb)
        template2 = """
    # Orbitals are: {0}
        """.format(
            str
        )
        return (
            template1.format(
                self,
                len(self.kvecs),
                len(self.energies),
                self.numk * self.nbands,
                len(self.proj),
                self.numk * self.nbands * self.natom,
                len(self.phase),
            )
            + template2
        )


def check_orbital_name(arg):
    """Return arg without change if arg is a member of the 'orbital name'.

    If arg is an alias of the (more appropriate) orbital
    name, return it as is.  If arg is neither the appropriate
    orbital name nor the alias, raise ValueError.

    Parameters
    ----------
    arg: str
        the string to be checked as the orbital name

    Returns
    --------
    str

    """
    translate_dict = {
        "pypx": "pxpy",
        "pzpx": "pxpz",
        "pzpy": "pypz",
        "pxpypz": "p",
        "pxpzpy": "p",
        "pypxpz": "p",
        "pypzpx": "p",
        "pzpxpy": "p",
        "pzpypx": "p",
        "spd": "tot",
    }
    proper_orbital_name_list = [
        "s",
        "py",
        "pz",
        "px",
        "dxy",
        "dyz",
        "dz2",
        "dxz",
        "dx2",
        "tot",
    ] + ["sp", "p", "pxpy", "pxpz", "pypz", "spd", "d"]
    if arg in translate_dict.keys():
        arg = translate_dict[arg]
    if arg in proper_orbital_name_list:
        return arg
    errmsg = arg + ": (composed) orbital name was not defined."
    raise ValueError(errmsg)


def shortcheck(procar) -> Tuple[int, int, int, List[str], bool]:
    """Check whether PROCAR file is good.

    Return numk, nbands, nion, orbital_names and
    True/False if collienar calculation

    """
    if "PROCAR lm decomposed + phase" not in next(procar):
        procar.close()
        raise RuntimeError(
            "This PROCAR is not a proper format\n \
                            Check your INCAR the calculations.\n"
        )
    tmp = next(procar)
    numk, nbands, natom = [int(i) for i in (tmp[14:20], tmp[39:43], tmp[62:-1])]
    _ = [next(procar) for i in range(5)]
    section = []
    orbitals = []
    phases = []
    for line in procar:
        if line.isspace():
            break
        elif "ion" in line and "tot" in line:
            orbitalnames = line.split()[1:]
            section = ["orbital"]
        elif "ion" in line and "tot" not in line:
            section.pop()
            section = ["phase"]
        elif "tot" in line and "ion" not in line:
            continue
        elif section == ["orbital"]:
            orbitals.append([float(i) for i in line.split()[1:]])
        elif section == ["phase"]:
            phases.append([float(i) for i in line.split()[1:]])
    if len(orbitals) == natom:
        collinear = True
    elif len(orbitals) == natom * 4:
        collinear = False
    else:
        raise RuntimeError("PROCAR is not proper format")
    procar.seek(0)
    return numk, nbands, natom, orbitalnames, collinear
