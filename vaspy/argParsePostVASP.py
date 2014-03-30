#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse, re

class APPV(argparse.ArgumentParser):
    __doc__ = '''
# Add methods parse_Atomselection, parse_AtonselectionNum
#   to ArgumentParser class
''' + argparse.ArgumentParser.__doc__

    def parse_Atomselection(self, l):
        '''
  # @param [String] l String to represent the range
  #   of the numbers deliminated by "-" or ",". 
  #   (ex.) "1-5,8,9,11-15" 
  #   @example l="1-5,8,8,9-15,10" #=> >["1", "10", "11", "12", "13", "14", "15", "2", "3", "4", "5", "8", "9"]}
  # @return [Array<String>] Returns the array that consists 
  #   of the ordered "String" represents the number.
'''
        a = l.split(',')
        output = set()
        reRange = re.compile(r'(\d+)-(\d+)')
        reSingle = re.compile(r'\d+')
        for i in a:
            if re.search(reRange, i):
                start, stop = re.findall(reRange, i)[0]
                output |= {str(i) for i in range(int(start), int(stop)+1)}
            elif re.search(reSingle, i):
                output.add(i)
        return sorted(output)

    def parse_AtomselectionNum(self, l):
        tmp = self.parse_Atomselection(l)
        return sorted(int(i) for i in tmp)

if __name__ == '__main__':
    opt = APPV()
    opt.add_argument('--atom', dest='atomrange', help='atom#')
    args = opt.parse_args()
    print(opt.parse_Atomselection(args.atomrange))
        
                
