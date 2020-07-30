# -*- conding: utf-8 -*-
"""VASPY class for INCAR file"""

from __future__ import annotations

from logging import DEBUG, INFO, Formatter, StreamHandler, getLogger
import re
from pathlib import Path
from typing import Dict, IO, List, Tuple, Union
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
    "SIGMA",
    "PREC",
    "ADDGRID",
    "NBANDS",
    "LREAL",
    "ENCUT",
]
tags["electron"] = ["EDIFF", "NELMIN", "NELM"]
tags["core"] = ["IBRION", "POTIM", "NSW", "EDIFFG"]
tags["dipole"] = ["IDIPOL", "LDIPOL"]
tags["calc"] = ["NCORE", "NPAR"]
tags["output"] = ["LWAVE", "LCHARG", "LVHAR", "LELF", "LORBIT", "NEDOS"]
tags["soi"] = ["LSORBIT", "LMAXMIX", "SAXIS", "NBANDS", "ISYM"]
tags["vdw"] = ["LUSE_VDW", "GGA", "AGGAC", "PARAM1", "PARAM2", "ZAB_VDW"]
tags["stm"] = ["LPARD", "EINT", "NBMOD"]

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
]
bool_ = ["LDIPOL", "LSORBIT", "LVHAR", "LWAVE", "LCHARG", "ADDGRID", "LUSE_VDW"]
str_ = ["SYSTEM", "PREC", "LREAL"]
# RWIGS and EINT are Tuple of float, SAXIS is Tuple of int
str_2 = ["RWIGS", "SAXIS", "EINT", "MAGMON"]


class Incar:
    """General class for INCAR file
    
    """

    def __init__(self, filename: Union[Path, str, None] = None) -> None:
        """
        Parameters
        ----------
        filename: str, pathlib.Path
            filename of INCAR
        """
        self._incar = {}
        self.section = ["generic:", "core:", "spin:", "tuning:", "SOI:", "VDW:"]

        if filename:
            self.load_file(open_by_suffix(str(filename)))

    def load_file(self, thefile: IO[str]) -> None:
        """Load INCAR file.

        Parameters
        -----------
        thefile: StringIO
            "INCAR" file
        """
        incar: Dict[str, Tuple[str, str]] = {}
        for line in thefile:
            line, active = remove_sharp(line)
            if "=" in line:
                tag, value = line.split("=")
                incar[tag.upper().strip()] = value.strip(), active
        for tag, value in incar.items():
            if tag in str_:
                self._incar[tag] = value
            elif tag in bool_:
                self._incar[tag] = value
            elif tag in str_2:
                self._incar[tag] = value
            elif tag in float_:
                self._incar[tag] = (float(value[0]), value[1])
            else:  # Default value is int
                self._incar[tag] = (int(value[0]), value[1])

    def repr(self) -> str:
        pass

    def __str__(self) -> str:
        output: str = ""
        incar = self._incar
        for tag_type in tags.keys():
            if len(set(tags[tag_type]) & set(incar.keys())) > 0:
                output += tag_type + ":\n"
                for tag in tags[tag_type]:
                    if tag in incar:
                        if tag in ["EDIFF"]:
                            output += "{} = {:.5E}\n".format(tag, incar[tag])
                        else:
                            output += "{} = {}\n".format(tag, incar[tag])
                        del incar[tag]
        if len(incar) != 0:
            print(incar)
            raise RuntimeError(
                "Unkonwn tags are used!!! Check your INCAR, or the script"
            )
        return output


def remove_sharp(str_: str) -> Tuple[str, bool]:
    """ Remvoe # / ! from the string.
    
    Parameters
    str_: str
    """
    active: bool = True
    str_ = str_.strip()
    if re.search("!|#", str_):
        str_ = re.sub("!|#", "", str_).strip()
        active = False
    return str_, active
