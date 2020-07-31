# -*- conding: utf-8 -*-
"""VASPY class for INCAR file"""

from __future__ import annotations

from logging import DEBUG, INFO, Formatter, StreamHandler, getLogger
from collections.abc import Mapping
import pprint
import re
from pathlib import Path
from typing import Any, Dict, Generator, IO, List, Tuple, Union
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


tags: Dict[str, List[str]] = {}
tags["generic"] = [
    "SYSTEM",
    "ISTART",
    "ICHARG",
    "ENCUT",
    "ISMEAR",
    "SIGMA",
    "PREC",
    "ADDGRID",
    "NBANDS",
    "LREAL",
    "ALGO",
]
tags["electron"] = ["NELMIN", "NELM", "EDIFF", "NELMDL", "RWIGS", "LASPH", "LEPSILON"]
tags["spin"] = ["ISPIN", "MAGMOM"]
tags["core"] = ["IBRION", "NSW", "EDIFFG", "POTIM", "NFREE"]
tags["dipole"] = ["IDIPOL", "LDIPOL"]
tags["calc"] = ["NCORE", "NPAR"]
tags["tuning"] = ["AMIX", "BMIX", "AMIX_MAG", "BMIX_MAG", "MAXMIX"]
tags["grid"] = ["NGXF", "NGYF", "NGZF", "NGX", "NGY", "NGZ"]
tags["lda+u"] = ["LDAU", "LDATYPE", "LDAUL", "LDAUU", "LDAUJ", "LDAUPRINT"]
tags["output"] = [
    "LWAVE",
    "LCHARG",
    "LVHAR",
    "LELF",
    "LORBIT",
    "NEDOS",
    "EMIN",
    "EMAX",
    "LAECHG",
    "LVTOT",
]
tags["soi"] = ["LSORBIT", "LMAXMIX", "SAXIS", "NBANDS", "ISYM"]
tags["vdw"] = ["LUSE_VDW", "GGA", "AGGAC", "PARAM1", "PARAM2", "ZAB_VDW"]
tags["stm"] = ["LPARD", "EINT", "NBMOD", "IBAND"]

float_ = [
    "ENCUT",
    "SIGMA",
    "EDIFFG",
    "POTIM",
    "EDIFF",
    "AMIX",
    "BMIX",
    "AMIX_MAG",
    "BMIX_MAG",
    "AGGAC",
    "PARAM1",
    "PARAM2",
    "ZAB_VDW",
    "NMIN",
    "NMAX",
]
bool_ = [
    "LDIPOL",
    "LSORBIT",
    "LVHAR",
    "LWAVE",
    "LCHARG",
    "ADDGRID",
    "LUSE_VDW",
    "LPARD",
    "LDAU",
    "LVTOT",
    "LAECHG",
    "LASPH",
    "LEPSILON",
]
str_ = ["SYSTEM", "PREC", "LREAL", "GGA", "ALGO"]
# RWIGS and EINT are Tuple of float, SAXIS is Tuple of int
str_2 = ["RWIGS", "SAXIS", "EINT", "LDATYPE", "LDAUL", "LDAUU", "LDAUJ", "MAGMOM"]


class Incar(Mapping):
    """General class for INCAR file
    
    """

    def __init__(self, filename: Union[Path, str, None] = None) -> None:
        """
        Parameters
        ----------
        filename: str, pathlib.Path
            filename of INCAR
        """
        self._incar: Dict[str, Tuple[Union[str, int, float], bool]] = {}
        self.additional_comments: Dict[str, str] = {}

        if filename:
            self.load_file(open_by_suffix(str(filename)))

    def load_file(self, thefile: IO[str]) -> None:
        """Load INCAR file.

        Parameters
        -----------
        thefile: StringIO
            "INCAR" file
        """
        incar: Dict[str, Tuple[str, bool]] = {}
        for line in thefile:
            line, active = remove_sharp(line)
            if "=" in line:
                tag, value_str = line.split("=")
                incar[tag.upper().strip()] = value_str.strip(), active
        for tag, value in incar.items():
            if tag in str_:
                self[tag] = value
            elif tag in bool_:
                self[tag] = value
            elif tag in str_2:
                self[tag] = value
            elif tag in float_:
                try:
                    self[tag] = (float(value[0]), value[1])
                except ValueError:
                    if "#" in value[0]:
                        val_com: str = re.sub("#", "!", value[0])
                    else:
                        val_com = value[0]
                    real_value, com = val_com.split("!", maxsplit=1)
                    self[tag] = (float(real_value), value[1])
                    self.additional_comments[tag] = com
            else:  # Default value is int
                try:
                    self[tag] = (int(value[0]), value[1])
                except ValueError:
                    if "#" in value[0]:
                        val_com = re.sub("#", "!", value[0])
                    else:
                        val_com = value[0]
                    logger.debug(tag, val_com)
                    real_value, com = val_com.split("!", maxsplit=1)
                    self[tag] = (int(real_value), value[1])
                    self.additional_comments[tag] = com

    def __iter__(self) -> Generator:
        for key in self._incar:
            yield key

    def __getitem__(self, key_item: str) -> Any:
        return self._incar.__getitem__(key_item)

    def __setitem__(
        self, key_item: str, value_item: Tuple[Union[str, int, float], bool]
    ) -> None:
        self._incar.__setitem__(key_item, value_item)

    def __len__(self) -> int:
        return len(self._incar)

    def __repr__(self) -> str:
        return pprint.pformat(self._incar)

    def __str__(self) -> str:
        output: str = ""
        incar = self._incar.copy()
        for tag_type in tags.keys():
            if len(set(tags[tag_type]) & set(incar.keys())) > 0:
                output += " {}:\n".format(tag_type)
                for tag in tags[tag_type]:
                    if tag in incar:
                        if tag in ["EDIFF"]:
                            if incar[tag][1]:
                                output += "    {} = {:.2E}".format(tag, incar[tag][0])
                            else:
                                output += "#   {} = {:.2E}".format(tag, incar[tag][0])
                        else:
                            if incar[tag][1]:
                                output += "    {} = {}".format(tag, incar[tag][0])
                            else:
                                output += "#   {} = {}".format(tag, incar[tag][0])
                        if tag in self.additional_comments:
                            output += "   ! {}\n".format(self.additional_comments[tag])
                        else:
                            output += "\n"
                        del incar[tag]
        if len(incar) != 0:
            print(incar)
            raise RuntimeError(
                "Unkonwn tags are used!!! Check your INCAR, or the script"
            )
        return output

    def active(self, keyword: str) -> bool:
        """Return True if keyword is active.  False if keyword is not set or comment out

        Parameters
        ----------
        keyword: str
            INCAR keyword
        """
        if (keyword in self) and self[keyword][1]:
            return True
        return False

    def lint_all(self) -> str:
        msg_lint = ""
        msg_lint += self.lint0()
        msg_lint += self.lint1()
        msg_lint += self.lint2()
        msg_lint += self.lint3()
        msg_lint += self.lint4()
        return msg_lint

    def lint0(self) -> str:
        """Check WAVECAR/CHGCAR should be written or not"""
        msg: str = ""
        if self["ICHARG"] == (11, True) and (
            self["LWAVE"] == (".TRUE.", True) or self["LCHARG"] == (".TRUE.", True)
        ):
            msg = "Want to draw the band ?.\n"
            msg += 'Recommend "LWAVE = .FALSE, LCHARG = .FALSE"\n'
        return msg

    def lint1(self) -> str:
        """if IBRION = 7, 8, switch off NCORE / NPAR"""
        msg: str = ""
        if (self["IBRION"][0] > 4 and self["IBRION"][0]) and (
            self.active("NPAR") or self.active("NCORE")
        ):
            msg = 'if IBRION > 4 (to DFPT), Remove "NPAR/NCORE" keyword\n'
        return msg

    def lint2(self) -> str:
        """dipole correction should be start with ISTART=2"""
        msg = ""
        if self["ISTART"] != (2, True) and self.active("LDIPOL"):
            msg = "DIPOLE correction is hard for not ISTART=2."
        return msg

    def lint3(self) -> str:
        """when LDIPOL is set, IDIPOL must be set, too"""
        msg = ""
        if self.active("LDIPOL") ^ self.active("IDIPOL"):
            msg = 'For dipoloe correction, need both "LDIPOL" and "IDIPOL"'
        return msg

    def lint4(self) -> str:
        """To draw the simulated stm image, ISTART must be 2"""
        msg = ""
        if self.active("LPARD") and self["ISTART"] != (2, True):
            msg = 'For LPARD to get partial charge densities, ISTART must be "2"'
        return msg


def remove_sharp(str_: str) -> Tuple[str, bool]:
    """ Remvoe # / ! from the string head.
    
    Parameters
    ----------
        str_: str
    """
    active: bool = True
    str_ = str_.strip()
    if re.search("!|#", str_):
        str_ = re.sub(r"^[!|#|\s]+", "", str_).strip()
        active = False
    logger.debug(str_)
    return str_, active
