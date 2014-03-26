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
    '''
  # @param [String] l String to represent the range
  #   of the numbers deliminated by "-" or ",". 
  #   (ex.) "1-5,8,9,11-15" 
  #   @example l="1-5,8,8,9-15,10" #=> >["1", "10", "11", "12", "13", "14", "15", "2", "3", "4", "5", "8", "9"]}
  # @return [Array<String>] Returns the array that consists 
  #   of the ordered "String" represents the number.
'''
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

if __name__ == '__main__':
    import argparse
    each_slice_demo = lambda L, n: list(each_slice(L, n))
    demo = {'each_slice_demo': (range(10), 3),
            'removeall': ([1, 0, 0, 1, 0, 1, 0, 0], 0),
            'flatten': ([1, [2, 3, [4, 5], [6]], 7, [8, 9]],),
            'parse_Atomselection': ('1-5, 8, 9, 11-15',),
            }
    argcounts = {'each_slice_demo': 2,
                 'removeall': 2,
                 'flatten': 1,
                 'parse_Atomselection': 1}
    available = ['all'] + list(demo.keys())
    parser = argparse.ArgumentParser(description='''collection of tools used in this package.''',
                 formatter_class=argparse.RawTextHelpFormatter,
                 epilog='''-a/--args option arguments are interpleted by eval().
so strings must be given with quotations("" or '').
Because command line regards spaces as break,
list argument must be written without any space.
(i.e. [1,2,3,4,5] is valid, while [1, 2, 3, 4, 5] is invalid.)''')
    parser.add_argument('choice', metavar='funcname', nargs='+',
                        choices=available,
                        help='''Demonstrate choosen function.
*all* shows all function in the choice.
If -a/--args option is given, get argument(s) from command line.
Otherwise use prepared argument(s).''')
    parser.add_argument('-a', '--args', metavar='values', nargs='+',
                        action='append', dest='values',
                        help='''Use given argument(s) for demonstration.
You have to use this option for each function.
See epilog for notices for argument notation.''')
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
                print('argument number not match (require {}, given {})\nuse default values.'.format(
                      argcounts[func], len(args.values[index])))
                values = demo[func]
            else:
                values = [eval(s) for s in args.values[index]]
            index += 1
        print('Demonstrate function {}()\ninput:  '.format(func), end='')
        print(*values, sep=' : ')
        print('output: ', end='')
        print(vars()[func](*values))
