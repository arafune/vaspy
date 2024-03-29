#!/usr/bin/env python
"""Script to use(demonstrate) vaspy.doscar."""

import argparse
import os
import sys
from pathlib import Path

from vaspy import tools
from vaspy.doscar import DOSCAR, PDOS, TDOS
from vaspy.outcar import OUTCAR

mypath = os.readlink(__file__) if Path(__file__).is_symlink() else __file__
sys.path.append(os.path.dirname(os.path.abspath(mypath)))

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("doscar", metavar="DOSCAR_file")
parser.add_argument(
    "-o",
    "--outcar",
    metavar="OUTCAR_file",
    nargs="?",
    default=NotImplemented,
    help="""Use OUTCAR file.
If file path is not given,
try to open OUTCAR file
in the same directory to DOSCAR.""",
)
# TIPS
# -o data/OUTCAR => args.outcar == "data/OUTCAR"
# -o             => args.outcar == None
# (not given)    => args.outcar == NotImplemented (default)
parser.add_argument(
    "-f",
    "--fermi",
    metavar="value",
    type=float,
    default=0.0,
    help="""Fermi level correction.
Energy shifts by this value.
if --outcar is set, this option is ignored.""",
)
parser.add_argument(
    "-s",
    "--site",
    metavar="atoms",
    action="append",
    dest="atomset",
    type=tools.atom_selection_to_list,
    help="""atom # specified with range.
Use "-" or ","
(ex.) --site 1,2,7-9""",
)
parser.add_argument(
    "-a",
    "--as",
    metavar="name",
    action="append",
    dest="atomsetname",
    help="""the name of the range identified by --site.
(ex.) --as layer1
the name is used in the output filename.""",
)

args = parser.parse_args()

doscar = DOSCAR(args.doscar)
atomlist: list[str] = []
if args.outcar is not NotImplemented:
    if args.outcar is None:
        args.outcar = str(Path(args.doscar).parent / "OUTCAR")
    outcar = OUTCAR(args.outcar)
    atomlist = outcar.atom_names
    args.fermi = outcar.fermi
#
if atomlist == []:
    atomlist.extend("atom" + str(i) for i in range(1, doscar.n_atom + 1))
#
# construct TDOS & PDOS objects
#
doses = [TDOS(doscar.dos_container.pop(0))]
#
doses.extend(
    PDOS(*each) for each in zip(doscar.dos_container, atomlist, strict=True)
)  # tmp[1:] ?
#
# Fermi level correction
#
[thedos.fermi_correction(args.fermi) for thedos in doses]
#
if args.atomset is not None:  # atomset and atomsetname are given by
    # the command line argument.
    if len(args.atomset) == len(args.atomsetname):
        pdoses = []
        for site, name in zip(args.atomset, args.atomsetname, strict=True):
            each = PDOS()
            for atom_number in site:
                each += doses[atom_number]
            each.site = name
            pdoses.append(each)
        for a_pdos in pdoses:
            filename = a_pdos.site + ".dat"
            a_pdos.export_csv(filename)
    else:
        msg = "Check your command."
        raise ValueError(msg)
#
doses[0].export_csv("total.dat")
for thedos, atom_index in zip(doses[1:], atomlist, strict=True):
    if isinstance(thedos, PDOS) and thedos.site == "":
        thedos.site = atom_index
    filename = atom_index + "_dos.dat"
    thedos.export_csv(filename)
