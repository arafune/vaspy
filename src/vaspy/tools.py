#!/usr/bin/env python
"""Module for tools used in vaspy."""
from __future__ import annotations

import bz2
import re
from itertools import zip_longest
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from numpy.typing import NDArray


def open_by_suffix(filename: str) -> IO[str]:
    """Open file."""
    if Path(filename).suffix == ".bz2":
        the_file = bz2.open(filename, mode="rt")
    else:
        the_file = Path(filename).open()
    return the_file


def each_slice(
    iterable: Iterable,
    n: int,
    fillvalue: float | None = None,
) -> Iterable[Any]:
    """each_slice(iterable, n[, fillvalue]) => iterator.

    make new iterator object which get n item from [iterable] at once.
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


_RERANGE = re.compile(r"(\d+)-(\d+)")
_RESINGLE = re.compile(r"\d+")


def atom_selection_to_list(input_str: str, number: bool = True) -> list[int | str]:
    """Return list of ordered "String" represents the number.

    Parameters
    ----------
    input_str: str
        range of the atoms. the numbers deliminated by "-" or ","

    Returns
    -------
    list
        ordered "String" represents the number.

    Example
    --------

    >>> atom_selection_to_list("1-5,8,8,9-15,10", False)
    ['1', '10', '11', '12', '13', '14', '15', '2', '3', '4', '5', '8', '9']
    >>> atom_selection_to_list("1-5,8,8,9-15,10")
    [1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15]

    """
    array = input_str.split(",")
    output = set()
    for each in array:
        if re.search(_RERANGE, each):
            start, stop = re.findall(_RERANGE, each)[0]
            # Version safety
            output |= {str(i) for i in range(int(start), int(stop) + 1)}
        elif re.search(_RESINGLE, each):
            output.add(each)
    if number:
        return sorted(int(i) for i in output)
    return sorted(output)


def atom_types_atomnums_to_atoms(
    atom_types: Iterable[str],
    atomnums: Iterable[int],
) -> tuple[str, ...]:
    """Return list representation for atom in use.

    Parameters
    ----------
    atom_types: list
        atom names
    atomnums: list
        atom numbers

    Examples
    --------
    >>> test_nums = [2, 3, 2, 1]
    >>> test_elements = ['Si', 'Ag', 'H', 'Si']
    >>> atom_types_atomnums_to_atoms(test_elements, test_nums)
    ('Si', 'Si', 'Ag', 'Ag', 'Ag', 'H', 'H', 'Si')

    """
    atoms = []
    for elem, nums in zip(atom_types, atomnums, strict=True):
        for _ in range(nums):
            atoms.append(elem)
    return tuple(atoms)


def atoms_to_atom_types_atomnums(atoms: list[str]) -> tuple[list[str], list[int]]:
    r"""Return atomnums and atom_types list.

    Returns
    -------
    atomnums
        list of number of atoms
    atom_types
        list of atomnames


    Examples
    --------
    >>> test = ['Si', 'Si', 'Ag', 'Ag', 'Ag', 'H', 'H', 'Si']
    >>> atoms_to_atom_types_atomnums(test)
    (['Si', 'Ag', 'H', 'Si'], [2, 3, 2, 1])

    """
    thelast = ""
    atomnums: list[int] = []
    atom_types: list[str] = []
    while atoms:
        atom = atoms.pop(0)
        if thelast == atom:
            atomnums[-1] = atomnums[-1] + 1
        else:
            atomnums.append(1)
            atom_types.append(atom)
        thelast = atom
    return atom_types, atomnums


def cuboid(crystal_axes: Sequence[float]) -> NDArray[np.float64]:
    """Return the coordinates for cuboid.

    Return the coordinates for the cuboid that includes tetrahedron
    represented by vectors.

    Parameters
    ----------
    vectors: array-like.
        Three vectors for tetrahedron.  (Crystal axis a,b,c)

    Return

    """
    a: NDArray[np.float64] = np.array(crystal_axes[0])
    b: NDArray[np.float64] = np.array(crystal_axes[1], dtype=np.float64)
    c: NDArray[np.float64] = np.array(crystal_axes[2])
    o: NDArray[np.float64] = np.array((0, 0, 0))
    points: NDArray[np.float64] = np.array((o, a, b, c, a + b, a + c, b + c, a + b + c))
    box: NDArray[np.float64] = np.array(
        (np.min(points.T, axis=1), np.max(points.T, axis=1)),
    ).T
    return box


if __name__ == "__main__":
    import argparse

    def EACH_SLICE_DEMO(L: Sequence[Any], n: int) -> list[Any]:
        return list(each_slice(L, n))

    demo = {
        "EACH_SLICE_DEMO": (range(10), 3),
        "removeall": ([1, 0, 0, 1, 0, 1, 0, 0], 0),
        "flatten": (
            (1, [range(2), 3, {4, 5}, [6]], frozenset([7, 8])),  # Version safety
        ),
        "parse_Atomselection": ("1-5,8,9,11-15",),
    }
    argcounts = {
        "EACH_SLICE_DEMO": 2,
        "removeall": 2,
        "flatten": 1,
        "parse_Atomselection": 1,
    }
    available = ["all", *list(demo.keys())]
    parser = argparse.ArgumentParser(
        description="""collection of tools used in this package.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""-a/--args option arguments are interpreted by eval().
so strings must be given with quotations("" or '').
Because command line regards spaces as break,
list argument must be written without any space.
(i.e. [1,2,3,4,5] is valid, while [1, 2, 3, 4, 5] is invalid.)""",
    )
    parser.add_argument(
        "choice",
        metavar="funcname",
        nargs="+",
        choices=available,
        help="""Demonstrate the function.
*all* shows all function in the choice.
If -a/--args option is given, get argument(s) from command line.
Otherwise use prepared argument(s).""",
    )
    parser.add_argument(
        "-a",
        "--args",
        metavar="values",
        nargs="+",
        action="append",
        dest="values",
        help="""Use given argument(s) for demonstration.
You have to use this option for each function.
See epilog for notices for argument notation.""",
    )
    args = parser.parse_args()

    choices = demo.keys() if "all" in args.choice else args.choice
    index = 0
    for func in choices:
        if args.values is None or len(args.values) <= index:
            values = demo[func]
        else:
            if argcounts[func] != len(args.values[index]):
                print(
                    f"""argument number not match (require {argcounts[func]}, given {len(args.values[index])})
use default values.""",
                )
                values = demo[func]
            else:
                values = tuple([eval(s) for s in args.values[index]])
            index += 1
        print(f"Demonstrate function {func}()\ninput:  ", end="")
        print(*values, sep=" : ")
        print("output: ", end="")
        print(vars()[func](*values))
