#!/usr/env python
# -*- coding: utf-8 -*-
''''''

from __future__ import print_function, division
import collections, re
import itertools as it

ziplong = it.izip_longest if hasattr(it, 'izip_longest') else it.zip_longest

def each_slice(iterable, n, fillvalue=None):
    '''each_slice(iterable, n[, fillvalue]) => iterator

make new iterator object which get n item from [iterable] at once.
'''
    args = [iter(iterable)] * n
    return ziplong(*args, fillvalue=fillvalue)

def removeall(L, value):
    'remove all *value* in [list] L'
    while L.count(value):
        L.remove(value)
    return L

def flatten(nested):
    '''flatten(iterable) => list

flatten nested iterables.
str and bytes are regarded as non-iterable'''
    basestring = basestring if hasattr(__builtins__, 'basestring') else (str, bytes)
    i = 0
    while i < len(nested):
        while isinstance(nested[i], collections.Iterable) and not isinstance(nested[i], basestring):
            if not nested[i]:
                nested.pop(i)
                i -= 1
                break
            else:
                nested[i:i+1] = nested[i]
        i += 1
    return nested

_rerange = re.compile(r'(\d+)-(\d+)')
_resingle = re.compile(r'\d+')

def parse_Atomselection(L):
    array = L.split(',')
    output = set()
    for each in array:
        if re.search(_rerange, each):
            start, stop = re.findall(_rerange, each)[0]
            output |= {str(i) for i in range(int(start), int(stop)+1)}
        elif re.search(_resingle, each):
            output.add(each)
    return sorted(output)

def parse_AtomselectionNum(L):
    return sorted(int(i) for i in parse_Atomselection(L))

