#! /usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

from __future__ import unicode_literals  # Version safety
from __future__ import print_function  # Version safety
import re


class OUTCAR(object):  # Version safety

    """Class for OUTCAR file that stores calculation details and/or calculation progress

    :Author: Ryuichi Arafune
    """

    def __init__(self, arg=None):
        self.nions = 0
        self.iontype = []
        self.ionnums = []
        self.posforce = []
        self.posforce_title = []
        self.atom_names = []
        self.fermi = 0.0
        self.atom_identifer = []
        if arg is not None:
            self.load_from_file(arg)

    def set_atom_names(self):
        """ build atom_names (the list of atomname_with_index)
        """
        self.atom_names = []
        for elm, n in zip(self.iontype, self.ionnums):
            for j in range(1, n + 1):
                tmp = elm + str(j)
                if tmp not in self.atom_names:
                    self.atom_names.append(tmp)
                else:
                    #                    jj = j
                    while tmp in self.atom_names:
                        j = j + 1
                        tmp = elm + str(j)
                    else:
                        self.atom_names.append(tmp)
        return self.atom_names

    def set_posforce_title(self):
        """build posforce_title
        """
        self.set_atom_names()
        self.posforce_title = [[i + "_x",
                                i + "_y",
                                i + "_z",
                                i + "_fx",
                                i + "_fy",
                                i + "_fz", ]
                               for i in self.atom_names]

    def load_from_file(self, arg):
        """Effectively, this is a constructor of OUTCAR object.
        arg is the file object of "OUTCAR"
        """
        # local variables
        section = []
        posforce = []
        # parse
        f = open(arg)
        for line in f:
            if section == ["force"]:
                if "total drift" in line:
                    section.pop()
                elif "---------------" in line:
                    pass
                elif "total drift:" in line:
                    section.pop()
                else:
                    posforce.append([float(x) for x in line.split()])
            else:
                if "number of dos" in line:
                    self.nions = int(line.split()[-1])
                elif "TITEL  =" in line:
                    self.iontype.append(line.split()[3])
                elif "ions per type " in line:
                    self.ionnums = [int(x) for x in line.split()[4:]]
                elif "POSITION" in line and "TOTAL-FORCE" in line:
                    section.append("force")
                elif "E-fermi" in line:
                    self.fermi = float(line.split()[2])
                else:
                    pass

        self.atom_identifer = [name + ":#" + str(index + 1)
                               for (index, name)
                               in enumerate(
            [elm + str(j)
             for (elm, n) in zip(self.iontype,
                                 self.ionnums)
             for j in range(1, int(n) + 1)])]
        #            i=1
        #            for elm, n in zip(self.iontype, self.ionnums):
        #                for j in range(1, int(n)+1):
        #                    self.atom_identifer.append(elem+str(j)+":#"+str(i))
        #                    i = i+1
        self.posforce = [posforce[i:i + self.nions]
                         for i in range(0, len(posforce), self.nions)]
        self.set_atom_names()
        self.set_posforce_title()

    def select_posforce_header(self, posforce_flag, *sites):
        if sites == () or sites[0] == []:
            sites = range(1, self.nions + 1)
        # if type(sites[0])==list or type(sites[0])==tuple:
        if isinstance(sites[0], (list, tuple)):
            sites = [n for n in sites[0]]
        return [posforce for (index, site) in enumerate(self.posforce_title)
                for (posforce, boolian) in zip(site, posforce_flag)
                if boolian and (index + 1 in sites)]
# return [posforce for (posforce, boolian) in zip(ithAtom, poforce_flag)
# if boolian==True for ithAtom in self.posforce_title  ] #which is
# correct?

    def select_posforce(self, posforce_flag, *sites):
        """Return the posforce corresponding the posforce_flag

        .. note:: posforce_flag: An 6-element True/False list that indicates the output (ex.) [True, True, False, True, True, False]
        """

        if sites == () or sites[0] == []:
            sites = range(1, self.nions + 1)
        if type(sites[0]) == list or type(sites[0]) == tuple:
            sites = [n for n in sites[0]]
        return [[posforce for (index, site) in enumerate(one_cycle)
                 for (posforce, boolian) in zip(site, posforce_flag)
                 if boolian
                 if index + 1 in sites]
                for one_cycle in self.posforce]


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Show position/force evolution in VASP calculation')

    parser.add_argument("-x", "--posx", action="store_true", dest="posx",
                        help="X-axis position", default=False)
    parser.add_argument("-y", "--posy", action="store_true", dest="posy",
                        help="Y-axis position", default=False)
    parser.add_argument("-z", "--posz", action="store_true", dest="posz",
                        help="Z-axis position", default=False)
    parser.add_argument("-X", "--forcex", action="store_true", dest="forcex",
                        help="Force along X-axis", default=False)
    parser.add_argument("-Y", "--forcey", action="store_true", dest="forcey",
                        help="Force along Y-axis", default=False)
    parser.add_argument("-Z", "--forcez", action="store_true", dest="forcez",
                        help="Force along Z-axis", default=False)
    parser.add_argument("--site", action="store", dest="site",
                        help="numbers deliminated by comma and hyphen for SITE you want to see", type=str)
    parser.add_argument("outcarfiles", type=str, nargs="+", metavar="OUTCAR",
                        help="OUTCAR_file")
    parser.add_argument("--plot", action="store_true",
                        help="show plot")
    args = parser.parse_args()
    pff = [args.posx, args.posy, args.posz,
           args.forcex, args.forcey, args.forcez]  # pff is "posforce_flag"
    if pff == [False, False, False, False, False, False]:
        pff = [True, True, True, True, True, True, True]

    def parse_SiteSelection(string):
        if string is None:
            return []
        else:
            array = string.split(",")
            output = []
            for i in array:
                if re.findall(r"(\d+)-(\d+)", i):
                    tmp = re.findall(r"(\d+)-(\d+)", i)
                    a_range = range(int(tmp[0][0]), int(tmp[0][1]) + 1)
                    for j in a_range:
                        output.append(j)
                elif re.findall(r"\d+", i):
                    tmp = re.findall(r"\d+", i)
                    output.append(int(tmp[0]))
            output = list(set(output))
            output.sort()
            return output

    sites = parse_SiteSelection(args.site)

    o = OUTCAR(args.outcarfiles.pop(0))
    headers = o.select_posforce_header(pff, sites)
    posforce = o.select_posforce(pff, sites)
    for outcarfile in args.outcarfiles:
        o = OUTCAR(outcarfile)
        posforce.extend(o.select_posforce(pff, sites))

    if args.plot:
        import matplotlib.pyplot as plt
        data = map(list, zip(*posforce))
        for i in zip(headers, data):
            plt.plot(i[1], label=i[0])
        plt.legend()
        plt.show()
    else:
        print("\t".join(headers))
        for item in posforce:
            print("\t".join([str(i) for i in item]))
