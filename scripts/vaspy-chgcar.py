#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script to use(demonstrate) vaspy.chgcar module."""

import argparse

import numpy as np

from vaspy.chgcar import CHGCAR

# HINT list methods for --spin option below.
# they are called with just one argument, CHGCAR instance.
spin_methods = {
    "mag": CHGCAR.magnetization,
    "magX": (lambda chg: chg.magnetization("x")),  # < Error !! FixMe
    "magY": (lambda chg: chg.magnetization("y")),  # < Error !! FixMe
    "magZ": (lambda chg: chg.magnetization("z")),  # < Error !! FixMe
    "majority": CHGCAR.majority_spin,
    "minority": CHGCAR.minority_spin,
}

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
    "--add", action="store_true", default=False, help="Add two CHGCAR files"
)
group.add_argument(
    "--merge",
    action="store_true",
    default=False,
    help="Add two CHGCAR files (ion position is not modified)",
)
group.add_argument(
    "--diff",
    action="store_true",
    default=False,
    help="Get difference of two CHGCAR files",
)
group.add_argument(
    "--spin",
    metavar="spin_operation",
    help="""spin-relatated operation.
when this option is set --add, -diff are ignored,
and CHGCAR_file_2 must not be set.
spin operation is one of the followings:
mag : show the magnetization
        density (for spin resolved calculations)
magX : show the magnetization density in
        the X direction (for non collinear calc.)
magY : show the magnetization density in
        the Y direction (for non collinear calc.)
magZ : show the magnetization density in
        the Z direction (for non collinear calc.)
majority : extract the part for the
            majority spin (for spin resolved calc.)
minority : extract the part for the
            minority spin (for spin resolved calc.)""",
)

parser.add_argument(
    "--output",
    metavar="file_name",
    help="""output file name
if not specified, use standard output""",
)
#
parser.add_argument("CHGCAR_file_1", type=CHGCAR)
parser.add_argument("CHGCAR_file_2", type=CHGCAR, nargs="?")
# TIP: if CHGCAR_file_2 is not specified,
# *None* is stored in arguments.CHGCAR_file_2, not CHGCAR(None)

args = parser.parse_args()


def check_position_axes(chgcar1: CHGCAR, chgcar2: CHGCAR) -> bool:
    """Check the cell vectors and atom positions are same in two CHGCAR.

    Parameters
    -----------
    chgcar1, chgcar2: vaspy.CHGCAR

    Returns
    -------
    bool

    """
    cell1 = chgcar1.poscar.cell_vecs
    cell2 = chgcar2.poscar.cell_vecs
    pos1 = np.array(chgcar1.poscar.positions)
    pos2 = np.array(chgcar2.poscar.positions)
    assert np.allclose(cell1, cell2), "UnitCells are inconsistent.  Abort."
    assert np.allclose(pos1, pos2), "Atom positions are inconsistent!!! Abort."
    return True


#
if args.spin is not None:
    if args.CHGCAR_file_2 is not None:
        parser.error("Only one CHGCAR file for --spin operations")
    if args.spin in spin_methods:
        if args.spin == "magX":
            dest_chgcar = args.CHGCAR_file_1.magnetization("x")
        elif args.spin == "magY":
            dest_chgcar = args.CHGCAR_file_1.magnetization("y")
        elif args.spin == "magZ":
            dest_chgcar = args.CHGCAR_file_1.magnetization("z")
        elif args.spin == "mag":
            dest_chgcar = args.CHGCAR_file_1.magnetization()
        elif args.spin == "majority":
            dest_chgcar = args.CHGCAR_file_1.majority_spin()
        elif args.spin == "minority":
            dest_chgcar = args.CHGCAR_file_1.minority_spin()
    else:
        parser.error("Such spin operation parameter is not defined.")
#
if args.add or args.diff or args.merge:
    if args.CHGCAR_file_2 is None:
        raise parser.error("Two CHGCAR files are required.")
    if args.add:
        dest_chgcar = args.CHGCAR_file_1 + args.CHGCAR_file_2
    elif args.diff:
        check_position_axes(args.CHGCAR_file_1, args.CHGCAR_file_2)
        dest_chgcar = args.CHGCAR_file_1 - args.CHGCAR_file_2
    elif args.merge:
        check_position_axes(args.CHGCAR_file_1, args.CHGCAR_file_2)
        dest_chgcar = args.CHGCAR_file_1.merge(args.CHGCAR_file_2)
#
if args.output is not None:
    dest_chgcar.save(args.output)
else:
    print(dest_chgcar)
