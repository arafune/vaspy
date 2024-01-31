# -*- conding: utf-8 -*-
"""VASPY class for INCAR file."""

from __future__ import annotations

import pprint
import re
from collections.abc import Generator, Mapping
from logging import INFO, Formatter, StreamHandler, getLogger
from typing import IO, TYPE_CHECKING

from vaspy.tools import open_by_suffix

if TYPE_CHECKING:
    from pathlib import Path

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


tags: dict[str, list[str]] = {}
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
# RWIGS and EINT are tuple of float, SAXIS is tuple of int
str_2 = ["RWIGS", "SAXIS", "EINT", "LDATYPE", "LDAUL", "LDAUU", "LDAUJ", "MAGMOM"]


def remove_sharp(str_: str) -> tuple[str, bool]:
    """Remove # / ! from the string head.

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


class Incar(Mapping):
    """General class for INCAR file."""

    def __init__(self, filename: Path | str = "") -> None:
        """Parameters
        ----------
        filename: str, pathlib.Path
            filename of INCAR
        """
        self._incar: dict[str, tuple[str | int | float, bool]] = {}
        self.additional_comments: dict[str, str] = {}
        #
        if filename:
            self.load_file(open_by_suffix(str(filename)))

    def load_file(self, the_file: IO[str]) -> None:
        """Load INCAR file.

        Parameters
        ----------
        the_file: StringIO
            "INCAR" file

        """
        incar: dict[str, tuple[str, bool]] = {}
        for line in the_file:
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
        yield from self._incar

    def __getitem__(self, key_item: str) -> tuple[str | float, bool]:
        return self._incar.__getitem__(key_item)

    def __setitem__(self, key_item: str, value_item: tuple[str | float, bool]) -> None:
        self._incar.__setitem__(key_item, value_item)

    def __len__(self) -> int:
        return len(self._incar)

    def __repr__(self) -> str:
        return pprint.pformat(self._incar)

    def __str__(self) -> str:
        output: str = ""
        incar = self._incar.copy()
        for tag_type in tags:
            if len(set(tags[tag_type]) & set(incar.keys())) > 0:
                output += f" {tag_type}:\n"
                for tag in tags[tag_type]:
                    if tag in incar:
                        if tag in ["EDIFF"]:
                            if incar[tag][1]:
                                output += f"    {tag} = {incar[tag][0]:.2E}"
                            else:
                                output += f"#   {tag} = {incar[tag][0]:.2E}"
                        else:
                            if incar[tag][1]:
                                output += f"    {tag} = {incar[tag][0]}"
                            else:
                                output += f"#   {tag} = {incar[tag][0]}"
                        if tag in self.additional_comments:
                            output += f"   ! {self.additional_comments[tag]}\n"
                        else:
                            output += "\n"
                        del incar[tag]
        if len(incar) != 0:
            print(incar)
            msg = "Unknown tags are used!!! Check your INCAR, or the script"
            raise RuntimeError(
                msg,
            )
        return output

    def active(self, keyword: str) -> str | float | bool:
        """Return True if keyword is active.  False if keyword is not set or comment out.

        Parameters
        ----------
        keyword: str
            INCAR keyword

        """
        if (keyword in self) and self[keyword][1]:
            return self[keyword][0]
        return False

    def lint_all(self) -> str:
        """Tiny lint for vasp.

        Returns
        -------
        str
            Check messages

        """
        checks: dict[str, bool | int | float | str] = {
            'When ICHARG = 11, Recommend "LWAVE = .FALSE, LCHARG = .FALSE"\n': (
                self.active("ICHARG") == 11
                and (
                    self.active("LWAVE") == ".TRUE."
                    or self.active("LCHARG") == ".TRUE."
                )
            ),
        }
        checks["DIPOLE correction is hard for not ISTART=2."] = (
            self.active("ISTART") != 2 and self.active("LDIPOL") == ".TRUE."
        )
        checks['For dipoloe correction, need both "LDIPOL" and "IDIPOL"\n'] = (
            self.active("LDIPOL") == ".TRUE."
        ) ^ self.active("IDIPOL")
        checks['For LPARD to get partial charge densities, ISTART must be "2"'] = (
            self.active("LPARD") == ".TRUE." and self.active("ISTART") != 2
        )
        checks["For SOI calculation, ISYM = -1 is recommended."] = (
            self.active("LSORBIT") == ".TRUE." and self.active("ISYM") != -1
        )
        checks['"NPAR" is not recommend. Consider to use "NCORE"\n'] = self.active(
            "NPAR",
        )
        checks['if IBRION > 4 (to DFPT), Remove "NPAR/NCORE" keyword\n'] = (
            self.active("IBRION") > 4
        ) and (self.active("NPAR") or self.active("NCORE"))

        msg_lint = ""
        for mesg, check_point in checks.items():
            if check_point:
                msg_lint += mesg
        return msg_lint
