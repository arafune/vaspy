# -*- coding:  utf-8 -*-
'''This module provide constant data such as atom mass etc'''
radii = {
    'Ac': 2.15,
    'Ag': 1.45,
    'Al': 1.21,
    'Am': 1.8,
    'Ar': 1.06,
    'As': 1.19,
    'At': 1.5,
    'Au': 1.36,
    'B': 0.84,
    'Ba': 2.15,
    'Be': 1.96,
    'Bi': 1.48,
    'Br': 1.2,
    'C': 0.76,
    'Ca': 1.76,
    'Cd': 1.44,
    'Ce': 2.04,
    'Cl': 1.02,
    'Cm': 1.69,
    'Co': 1.26,
    'Cr': 1.39,
    'Cs': 2.44,
    'Cu': 1.32,
    'Dy': 1.92,
    'Er': 1.89,
    'Eu': 1.98,
    'F': 0.57,
    'Fe': 1.32,
    'Fr': 2.6,
    'Ga': 1.22,
    'Gd': 1.96,
    'Ge': 1.2,
    'H': 0.31,
    'He': 0.28,
    'Hf': 1.75,
    'Hg': 1.32,
    'Ho': 1.92,
    'I': 1.39,
    'In': 1.42,
    'Ir': 1.41,
    'K': 2.03,
    'Kr': 1.16,
    'La': 2.07,
    'Li': 1.28,
    'Lu': 1.87,
    'Mg': 1.41,
    'Mn': 1.39,
    'Mo': 1.54,
    'N': 0.71,
    'Na': 1.66,
    'Nb': 1.64,
    'Nd': 2.01,
    'Ne': 0.58,
    'Ni': 1.24,
    'Np': 1.9,
    'O': 0.66,
    'Os': 1.44,
    'P': 1.07,
    'Pa': 2.0,
    'Pb': 1.46,
    'Pd': 1.39,
    'Pm': 1.99,
    'Po': 1.49,
    'Pr': 2.03,
    'Pt': 1.36,
    'Pu': 1.87,
    'Ra': 2.21,
    'Rb': 2.2,
    'Re': 1.51,
    'Rh': 1.42,
    'Rn': 1.5,
    'Ru': 1.46,
    'S': 1.05,
    'Sb': 1.39,
    'Sc': 1.7,
    'Se': 1.2,
    'Si': 1.11,
    'Sm': 1.98,
    'Sn': 1.39,
    'Sr': 1.95,
    'Ta': 1.7,
    'Tb': 1.94,
    'Tc': 1.47,
    'Te': 1.38,
    'Th': 2.06,
    'Ti': 1.6,
    'Tl': 1.45,
    'Tm': 1.9,
    'U': 1.96,
    'V': 1.53,
    'W': 1.62,
    'Xe': 1.4,
    'Y': 1.9,
    'Yb': 1.87,
    'Zn': 1.22,
    'Zr': 1.75,
}

# Relative atomic masses from J. S. Coursey, D. J. Schwab, J. J. Tsai,
# and R. A. Dragoset, NIST Physical Measurement Laboratory
# www.nist.gov/pml/data/comp
#
# Where this reference gives a range and/or the relative abundance of
# isotopes is unknown, a simple mean was taken

masses = {
    'Ac': 227.0,
    'Ag': 107.8682,
    'Al': 26.9815385,
    'Am': 242.0591053,
    'Ar': 39.948,
    'As': 74.921595,
    'At': 210.0,
    'Au': 196.966569,
    'B': 10.8135,
    'Ba': 137.327,
    'Be': 9.0121831,
    'Bh': 272.13826,
    'Bi': 208.9804,
    'Bk': 248.0726475,
    'Br': 79.904,
    'C': 12.0106,
    'Ca': 40.078,
    'Cd': 112.414,
    'Ce': 140.116,
    'Cf': 250.578118975,
    'Cl': 35.451499999999996,
    'Cm': 245.0654423,
    'Cn': 285.17712,
    'Co': 58.933194,
    'Cr': 51.9961,
    'Cs': 132.90545196,
    'Cu': 63.546,
    'D': 2.01410177812,
    'Db': 268.12567,
    'Ds': 281.16451,
    'Dy': 162.5,
    'Er': 167.259,
    'Es': 252.08298,
    'Eu': 151.964,
    'F': 18.998403163,
    'Fe': 55.845,
    'Fl': 289.19042,
    'Fm': 257.0951061,
    'Fr': 223.0,
    'Ga': 69.723,
    'Gd': 157.25,
    'Ge': 72.63,
    'H': 1.00782503223,
    'He': 4.002,
    'Hf': 178.4,
    'Hg': 200.5,
    'Ho': 164.93,
    'Hs': 270.19,
    'I': 126.90,
    'In': 114.8,
    'Ir': 192.2,
    'K': 39.098,
    'Kr': 83.79,
    'La': 138.97,
    'Li': 6.967,
    'Lr': 262.11,
    'Lu': 174.9,
    'Lv': 293.29,
    'Md': 259.14,
    'Mg': 24.30,
    'Mn': 54.934,
    'Mo': 95.95,
    'Mt': 276.19,
    'N': 14.006,
    'Na': 22.98928,
    'Nb': 92.90,
    'Nd': 144.2,
    'Ne': 20.17,
    'Ni': 58.69,
    'No': 259.13,
    'Np': 237.0,
    'O': 15.999,
    'Os': 190.2,
    'P': 30.973998,
    'Pa': 231.08,
    'Pb': 207.2,
    'Pd': 106.4,
    'Pm': 145.0,
    'Po': 209.0,
    'Pr': 140.96,
    'Pt': 195.0,
    'Pu': 244.0,
    'Ra': 226.0,
    'Rb': 85.4678,
    'Re': 186.207,
    'Rf': 267.121,
    'Rg': 280.165,
    'Rh': 102.905,
    'Rn': 222.0,
    'Ru': 101.07,
    'S': 32.0674999999995,
    'Sb': 121.76,
    'Sc': 44.9559,
    'Se': 78.971,
    'Sg': 271.13393,
    'Si': 28.085,
    'Sm': 150.36,
    'Sn': 118.71,
    'Sr': 87.62,
    'T': 3.0160492779,
    'Ta': 180.94788,
    'Tb': 158.92535,
    'Tc': 98.0,
    'Te': 127.6,
    'Th': 232.0377,
    'Ti': 47.867,
    'Tl': 204.3835,
    'Tm': 168.93422,
    'U': 238.02891,
    'Uuo': 294.21392,
    'Uup': 288.19274,
    'Uus': 292.20746,
    'Uut': 284.17873,
    'V': 50.9415,
    'W': 183.84,
    'Xe': 131.293,
    'pY': 88.90584,
    'Yb': 173.054,
    'Zn': 65.38,
    'Zr': 91.224,
}
