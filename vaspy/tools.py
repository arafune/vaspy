#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module for tools used in vaspy
"""

from __future__ import division, print_function  # Version safety

import bz2
import itertools as it
import os
import re
from collections import Iterable

# Version safety
ZIPLONG = it.izip_longest if hasattr(it, 'izip_longest') else it.zip_longest
FLATTEN_IGNORE = (dict, str, bytes, bytearray)  # or (dict, basestring)


def open_by_suffix(filename):
    """Open file."""
    if os.path.splitext(filename)[1] == '.bz2':
        try:
            thefile = bz2.open(filename, mode='rt')
        except AttributeError:
            thefile = bz2.BZ2File(filename, mode='r')
    else:
        thefile = open(filename)
    return thefile


def each_slice(iterable, n, fillvalue=None):
    """each_slice(iterable, n[, fillvalue]) => iterator

    make new iterator object which get n item from [iterable] at once.
    """
    args = [iter(iterable)] * n
    return ZIPLONG(*args, fillvalue=fillvalue)  # Version safety


def flatten(nested, target=Iterable, ignore=FLATTEN_IGNORE):
    """Flatten iterabable object.

    flatten(iterable) => list

    flatten nested iterable.
    expand object which is included in *target*
    AND NOT included in *ignore*.
    default *include* is all iterable object.
    defualt *ignore* is dict() and str-like objects.
    (i.e. str, bytes, unicode(2.x), bytearray(3.x))

    >>> flatten((1, [range(2), 3, set([4, 5]), [6]], frozenset([7, 8])))
    [1, 0, 1, 3, 4, 5, 6, 8, 7]
    """
    if (isinstance(nested, target) and not isinstance(nested, ignore)):
        nested = list(nested)
    i = 0
    while i < len(nested):
        while (isinstance(nested[i], target)
               and not isinstance(nested[i], ignore)):
            if not nested[i]:
                nested.pop(i)
                i -= 1
                break
            else:
                nested[i:i + 1] = nested[i]
        i += 1
    return nested


_RERANGE = re.compile(r'(\d+)-(\d+)')
_RESINGLE = re.compile(r'\d+')


def atom_selection_to_list(input_str, number=True):
    """Return list of ordered "String" represents the number.

    Parameters
    ----------
    input_str: str
        range of the atoms. the numbers deliminated by "-" or ","

    Returns
    --------
    list
       ordered "String" represents the number.

    Example
    --------

    >>> atom_selection_to_list("1-5,8,8,9-15,10", False)
    ['1', '10', '11', '12', '13', '14', '15', '2', '3', '4', '5', '8', '9']
    >>> atom_selection_to_list("1-5,8,8,9-15,10")
    [1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15]

    """
    array = input_str.split(',')
    output = set()
    for each in array:
        if re.search(_RERANGE, each):
            start, stop = re.findall(_RERANGE, each)[0]
            # Version safety
            output |= set(str(i) for i in range(int(start), int(stop) + 1))
        elif re.search(_RESINGLE, each):
            output.add(each)
    if number:
        return sorted(int(i) for i in output)
    return sorted(output)


def atomtypes_atomnums_to_atoms(atomtypes, atomnums):
    """Return list representation for atom in use.

    Parameters
    ------------
    atomtypes: list
        atom names
    atomnums: list
        atom numbers

    Examples
    --------
    >>> test_nums = [2, 3, 2, 1]
    >>> test_elements = ['Si', 'Ag', 'H', 'Si']
    >>> atomtypes_atomnums_to_atoms(test_elements, test_nums)
    ('Si', 'Si', 'Ag', 'Ag', 'Ag', 'H', 'H', 'Si')

    """
    atoms = []
    for elem, nums in zip(atomtypes, atomnums):
        for _ in range(nums):
            atoms.append(elem)
    return tuple(atoms)


def atoms_to_atomtypes_atomnums(atoms):
    r"""Return atomnums and atomtypes list.

    Returns
    --------
    atomnums
        list of number of atoms
    atomtypes
        list of atomnames


    Examples
    --------
    >>> test = ['Si', 'Si', 'Ag', 'Ag', 'Ag', 'H', 'H', 'Si']
    >>> atoms_to_atomtypes_atomnums(test)
    (['Si', 'Ag', 'H', 'Si'], [2, 3, 2, 1])

    """
    thelast = ''
    atomnums = []
    atomtypes = []
    while atoms:
        atom = atoms.pop(0)
        if thelast == atom:
            atomnums[-1] = atomnums[-1] + 1
        else:
            atomnums.append(1)
            atomtypes.append(atom)
        thelast = atom
    return atomtypes, atomnums


if __name__ == '__main__':
    import argparse

    def EACH_SLICE_DEMO(L, n):
        return list(each_slice(L, n))

    demo = {
        'EACH_SLICE_DEMO': (range(10), 3),
        'removeall': ([1, 0, 0, 1, 0, 1, 0, 0], 0),
        'flatten': (
            (
                1,
                [range(2), 3, set([4, 5]), [6]],  # Version safety
                frozenset([7, 8])), ),
        'parse_Atomselection': ('1-5,8,9,11-15', ),
    }
    argcounts = {
        'EACH_SLICE_DEMO': 2,
        'removeall': 2,
        'flatten': 1,
        'parse_Atomselection': 1
    }
    available = ['all'] + list(demo.keys())
    parser = argparse.ArgumentParser(
        description="""collection of tools used in this package.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""-a/--args option arguments are interpleted by eval().
so strings must be given with quotations("" or '').
Because command line regards spaces as break,
list argument must be written without any space.
(i.e. [1,2,3,4,5] is valid, while [1, 2, 3, 4, 5] is invalid.)""")
    parser.add_argument('choice',
                        metavar='funcname',
                        nargs='+',
                        choices=available,
                        help="""Demonstrate choosen function.
*all* shows all function in the choice.
If -a/--args option is given, get argument(s) from command line.
Otherwise use prepared argument(s).""")
    parser.add_argument('-a',
                        '--args',
                        metavar='values',
                        nargs='+',
                        action='append',
                        dest='values',
                        help="""Use given argument(s) for demonstration.
You have to use this option for each function.
See epilog for notices for argument notation.""")
    args = parser.parse_args()

    if 'all' in args.choice:
        choices = demo.keys()
    else:
        choices = args.choice
    index = 0
    for func in choices:
        if args.values is None or len(args.values) <= index:
            values = demo[func]
        else:
            if argcounts[func] != len(args.values[index]):
                print("""argument number not match (require {0}, given {1})
use default values.""".format(argcounts[func], len(args.values[index])))
                values = demo[func]
            else:
                values = [eval(s) for s in args.values[index]]
            index += 1
        print('Demonstrate function {0}()\ninput:  '.format(func), end='')
        print(*values, sep=' : ')
        print('output: ', end='')
        print(vars()[func](*values))
