#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re


class APPV(argparse.ArgumentParser):
    __doc__ = '''Add methods parse_Atomselection, parse_AtonselectionNum to
ArgumentParser class ''' + argparse.ArgumentParser.__doc__

    def parse_Atomselection(self, l):
        '''Return list of ordered "String" represents the number

        :param l: String for the range of the ions (deliminator: "-" or ",")
        :type l: str
        :return: list of the ordered "String" represents the number.
        :rtype: list

        >>> a = APPV()
        >>> a.parse_Atomselection("1-5,8,9,11-15")
        ['1', '11', '12', '13', '14', '15', '2', '3', '4', '5', '8', '9']
        '''
        a = l.split(',')
        output = set()
        reRange = re.compile(r'(\d+)-(\d+)')
        reSingle = re.compile(r'\d+')
        for i in a:
            if re.search(reRange, i):
                start, stop = re.findall(reRange, i)[0]
                output |= {str(i) for i in range(int(start), int(stop) + 1)}
            elif re.search(reSingle, i):
                output.add(i)
        return sorted(output)

    def parse_AtomselectionNum(self, l):
        '''Very similar with parse_Atomselection but returns the array
        of the number not array of the string.

        :param l: String for the range of the ions (deliminator: "-" or ",")
        :type l: str
        :return: list of the ordered 'integer' represents the number.
        :rtype: list

        >>> a = APPV()
        >>> a.parse_AtomselectionNum("1-5,8,9,11-15")
        [1, 2, 3, 4, 5, 8, 9, 11, 12, 13, 14, 15]
        '''
        tmp = self.parse_Atomselection(l)
        return sorted(int(i) for i in tmp)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    opt = APPV()
    opt.add_argument('--atom', dest='atomrange', help='atom#')
    args = opt.parse_args()
    print(opt.parse_Atomselection(args.atomrange))
