# -*- codinng: utf-8 -*-
'''This module provides vsim_asc class

This module might not be a member of vaspy. But the input file .*ascii
is generated usually by phonopy, and the input files needed to phonopy
calculation is provided from the VASP calculations.

Purpose of this module
This module generates

1) POSCAR files for phonon displacement
2) POVRAY scene file for animation

The first is absolutely required.
 '''

import os.path
import bz2
import numpy as np


class VSIM_ASC(object):
    def __init__(self, filename=None):
        self.ionnums = []
        self.iontypes = []
        self.system_name = ""
        #
        self.qpts = []
        self.freqs = []
        #
        if filename:
            if os.path.splitext(filename)[1] == '.bz2':
                try:
                    thefile = bz2.open(filename, mode='rt')
                except AttributeError:
                    thefile = bz2.BZ2File(filename, mode='r')
            else:
                thefile = open(filename)
            self.load_file(thefile)

    def load_file(self, thefile):
        '''vsim.ascii parser

        Parameters
        ----------

        thefile: StringIO
            "VSIM.ascii" file
'''
        phonon_lines = []
        # the first line is system name
        self.system_name = next(thefile)[1:].strip()
        # the 2nd line represents dxx, dyx, dyy
        dxx, dyx, dyy = [float(x) for x in next(thefile).split()]
        # the 3rd line represents dzx, dzy, dzz
        dzx, dzy, dzz = [float(x) for x in next(thefile).split()]
        self.lattice_vector = np.array([[dxx, 0, 0],
                                        [dyx, dyy, 0],
                                        [dzx, dzy, dzz]])
        ions = []
        positions = []
        for line in thefile:
            line = line.strip()
            if line[0] == "#" or line[0] == "!":
                phonon_lines.append(line[1:].strip())
            else:
                x, y, z, ion = line.split()
                ions.append(ion)
                positions.append(np.array([float(x), float(y), float(z)]))
        self.ionnums, self.iontypes = ions_to_iontypes_ionnums(ions)
        #
        dis_vectors = []
        #
        for line in phonon_lines:
            if 'metaData' in line:
                aline = line[15:]
                if aline[-1] == '\\':
                    aline = aline[:-1]
                modedata = aline.split(';')
                qpt = np.array([float(x) for x in modedata[0:3]])
                freq = float(modedata[3])
                self.qpts.append(qpt)
                self.freqs.append(freq)
            elif ']' in line:
                pass
            else:  # displacement vector
                vectors = [float(x) for x in line[1:-1].split(';')]
                dis_vectors.append([vectors[0] + vectors[3]*1j,
                                    vectors[1] + vectors[4]*1j,
                                    vectors[2] + vectors[5]*1j])
        n_phonons = len(self.freqs)
        self.d_vectors = np.array(dis_vectors).reshape(n_phonons,
                                                       sum(self.ionnums),
                                                       3)
        self.freqs = np.array(self.freqs)


def ions_to_iontypes_ionnums(ions):
    '''Return ionnums and iontypes list

    Returns
    --------

    ionnums
        list of number of ions

    iontypes
        list of ionnames


    Example
    -----------

    >>> ions_to_iontypes_ionnums(['Si', 'Si', 'Ag', 'Ag', 'Ag', \
                                  'Ag', 'H', 'H', 'Si'])
    ([2, 4, 2, 1], ['Si', 'Ag', 'H', 'Si'])
'''
    thelast = ''
    ionnums = []
    iontypes = []
    while ions:
        ion = ions.pop(0)
        if thelast == ion:
            ionnums[-1] = ionnums[-1] + 1
        else:
            ionnums.append(1)
            iontypes.append(ion)
        thelast = ion
    return ionnums, iontypes


if __name__ == '__main__':
    import doctest
    doctest.testmod()
