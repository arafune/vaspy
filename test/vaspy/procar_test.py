#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import os
import tempfile
import numpy as np
import vaspy.procar as procar


class TestPROCAR(unittest.TestCase):
    def setUp(self):
        # single state
        global procar_single
        filePROCAR = tempfile.mkstemp()
        f = open(filePROCAR[1], 'w')
        f.write(procar_single)
        f.close()
        self.procar_single = procar.PROCAR(filePROCAR[1],phase_read = True)
        os.remove(filePROCAR[1])
        # without Spi/Spin
        global procar_woSpi
        filePROCAR = tempfile.mkstemp()
        f = open(filePROCAR[1], 'w')
        f.write(procar_woSpi)
        f.close()
        self.procar_std = procar.PROCAR(filePROCAR[1], phase_read = True)
        os.remove(filePROCAR[1])
        # with Spin
        global procar_Spin
        filePROCAR = tempfile.mkstemp()
        f = open(filePROCAR[1], 'w')
        f.write(procar_Spin)
        f.close()
        self.procar_spin = procar.PROCAR(filePROCAR[1], phase_read = True)
        os.remove(filePROCAR[1])

    def procar_load_test(self):
        pass

    def procar_std_print_test(self):
        self.assertEqual(output_print_procar_std, self.procar_std.__str__())   

    def procar_spin_print_test(self):
        self.assertEqual(output_print_procar_spin, self.procar_spin.__str__())

output_print_procar_std="""The properties of this procar:
  # of k-points: 3
  # of bands: 2
  # of ions: 3
  # of kvectors: 3
  # of energies: 6
    ((# of k-points) * (# of bands) = 3*2=6)
  # of orbital component: 18
    ((# of k-points) * (# of bands) * (# of ions) =
        3*2*3=18)
  # of phase component: 0
  Orbitals are: ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']
  spininfo: ('',)
"""

output_print_procar_spin="""The properties of this procar:
  # of k-points: 3
  # of bands: 2
  # of ions: 3
  # of kvectors: 6
  # of energies: 12
    ((# of k-points) * (# of bands) = 3*2=6)
  # of orbital component: 36
    ((# of k-points) * (# of bands) * (# of ions) =
        3*2*3=18)
  # of phase component: 0
  Orbitals are: ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']
  spininfo: ('_up', '_down')
"""

if __name__ == '__main__':
    unittest.main()

procar_woSpi = """PROCAR lm decomposed + phase
# of k-points:    3         # of bands: 2         # of ions: 3

 k-point    1 :    0.00000000 0.00000000 0.00000000     weight = 0.33333333

band   1 # energy  -10.000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  0.0000  0.0001  0.0002  0.0003  0.0004  0.0005  0.0006  0.0007  0.0008  0.0036
  2  0.0010  0.0011  0.0012  0.0013  0.0014  0.0015  0.0016  0.0017  0.0018  0.0126
  3  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

band   2 # energy  -5.000000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  0.0030  0.0031  0.0032  0.0033  0.0034  0.0035  0.0036  0.0037  0.0038  0.0306
  2  0.0040  0.0041  0.0042  0.0043  0.0044  0.0045  0.0046  0.0047  0.0048  0.0396
  3  0.0050  0.0051  0.0052  0.0053  0.0054  0.0055  0.0056  0.0057  0.0058  0.0486
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

 k-point    2 :    0.2500000 0.2500000 0.00000000     weight = 0.33333333

band   1 # energy  -7.00000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  0.0060  0.0061  0.0062  0.0063  0.0064  0.0065  0.0066  0.0067  0.0068  0.0576
  2  0.0070  0.0071  0.0072  0.0073  0.0074  0.0075  0.0076  0.0077  0.0078  0.0666
  3  0.0080  0.0081  0.0082  0.0083  0.0084  0.0085  0.0086  0.0087  0.0088  0.0756
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

band   2 # energy  -4.000000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  0.0090  0.0091  0.0092  0.0093  0.0094  0.0095  0.0096  0.0097  0.0098  0.0846
  2  0.0100  0.0101  0.0102  0.0103  0.0104  0.0105  0.0106  0.0107  0.0108  0.0936
  3  0.0110  0.0111  0.0112  0.0113  0.0114  0.0115  0.0116  0.0117  0.0118  0.1026
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

 k-point    3 :    0.50000000 0.50000000 0.00000000     weight = 0.33333333

band   1 # energy  -6.0000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  0.0120  0.0121  0.0122  0.0123  0.0124  0.0125  0.0126  0.0127  0.0128  0.1116
  2  0.0130  0.0131  0.0132  0.0133  0.0134  0.0135  0.0136  0.0137  0.0138  0.1206
  3  0.0140  0.0141  0.0142  0.0143  0.0144  0.0145  0.0146  0.0147  0.0148  0.1296
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

band   2 # energy  -1.00000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  0.0150  0.0151  0.0152  0.0153  0.0154  0.0155  0.0156  0.0157  0.0158  0.1386
  2  0.0160  0.0161  0.0162  0.0163  0.0164  0.0165  0.0166  0.0167  0.0168  0.1476
  3  0.0170  0.0171  0.0172  0.0173  0.0174  0.0175  0.0176  0.0177  0.0178  0.1566
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
"""

procar_Spin = """PROCAR lm decomposed + phase
# of k-points:    3         # of bands: 2         # of ions: 3

 k-point    1 :    0.00000000 0.00000000 0.00000000     weight = 0.33333333

band   1 # energy  -10.000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  0.0000  0.0001  0.0002  0.0003  0.0004  0.0005  0.0006  0.0007  0.0008  0.0036 
  2  0.0010  0.0011  0.0012  0.0013  0.0014  0.0015  0.0016  0.0017  0.0018  0.0126 
  3  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

band   2 # energy  -5.000000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  0.0030  0.0031  0.0032  0.0033  0.0034  0.0035  0.0036  0.0037  0.0038  0.0306
  2  0.0040  0.0041  0.0042  0.0043  0.0044  0.0045  0.0046  0.0047  0.0048  0.0396
  3  0.0050  0.0051  0.0052  0.0053  0.0054  0.0055  0.0056  0.0057  0.0058  0.0486
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

 k-point    2 :    0.2500000 0.2500000 0.00000000     weight = 0.33333333

band   1 # energy  -7.00000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  0.0060  0.0061  0.0062  0.0063  0.0064  0.0065  0.0066  0.0067  0.0068  0.0576 
  2  0.0070  0.0071  0.0072  0.0073  0.0074  0.0075  0.0076  0.0077  0.0078  0.0666 
  3  0.0080  0.0081  0.0082  0.0083  0.0084  0.0085  0.0086  0.0087  0.0088  0.0756
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

band   2 # energy  -4.000000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  0.0090  0.0091  0.0092  0.0093  0.0094  0.0095  0.0096  0.0097  0.0098  0.0846 
  2  0.0100  0.0101  0.0102  0.0103  0.0104  0.0105  0.0106  0.0107  0.0108  0.0936 
  3  0.0110  0.0111  0.0112  0.0113  0.0114  0.0115  0.0116  0.0117  0.0118  0.1026 
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

 k-point    3 :    0.50000000 0.50000000 0.00000000     weight = 0.33333333

band   1 # energy  -6.0000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  0.0120  0.0121  0.0122  0.0123  0.0124  0.0125  0.0126  0.0127  0.0128  0.1116 
  2  0.0130  0.0131  0.0132  0.0133  0.0134  0.0135  0.0136  0.0137  0.0138  0.1206 
  3  0.0140  0.0141  0.0142  0.0143  0.0144  0.0145  0.0146  0.0147  0.0148  0.1296 
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

band   2 # energy  -1.00000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  0.0150  0.0151  0.0152  0.0153  0.0154  0.0155  0.0156  0.0157  0.0158  0.1386 
  2  0.0160  0.0161  0.0162  0.0163  0.0164  0.0165  0.0166  0.0167  0.0168  0.1476 
  3  0.0170  0.0171  0.0172  0.0173  0.0174  0.0175  0.0176  0.0177  0.0178  0.1566 
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  1  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  2  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
  3  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

# of k-points:    3         # of bands: 2         # of ions: 3

 k-point    1 :    0.00000000 0.00000000 0.00000000     weight = 0.33333333

band   1 # energy  -10.5000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  1.0000  1.0001  1.0002  1.0003  1.0004  1.0005  1.0006  1.0007  1.0008  9.0036 
  2  1.0010  1.0011  1.0012  1.0013  1.0014  1.0015  1.0016  1.0017  1.0018  9.0126 
  3  1.0020  1.0021  1.0022  1.0023  1.0024  1.0025  1.0026  1.0027  1.0028  9.0216 
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  1  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  2  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  2  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  3  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  3  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000

band   2 # energy  -5.500000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  1.0030  1.0031  1.0032  1.0033  1.0034  1.0035  1.0036  1.0037  1.0038  9.0306 
  2  1.0040  1.0041  1.0042  1.0043  1.0044  1.0045  1.0046  1.0047  1.0048  9.0396 
  3  1.0050  1.0051  1.0052  1.0053  1.0054  1.0055  1.0056  1.0057  1.0058  9.0486 
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  1  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  2  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  2  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  3  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  3  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000

 k-point    2 :    1.2500000 1.2500000 1.00000000     weight = 1.33333333

band   1 # energy  -7.500000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  1.0060  1.0061  1.0062  1.0063  1.0064  1.0065  1.0066  1.0067  1.0068  9.0576 
  2  1.0070  1.0071  1.0072  1.0073  1.0074  1.0075  1.0076  1.0077  1.0078  9.0666 
  3  1.0080  1.0081  1.0082  1.0083  1.0084  1.0085  1.0086  1.0087  1.0088  9.0756 
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  1  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  2  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  2  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  3  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  3  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000

band   2 # energy  -4.5000000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  1.0090  1.0091  1.0092  1.0093  1.0094  1.0095  1.0096  1.0097  1.0098  9.0846 
  2  1.0100  1.0101  1.0102  1.0103  1.0104  1.0105  1.0106  1.0107  1.0108  9.0936 
  3  1.0110  1.0111  1.0112  1.0113  1.0114  1.0115  1.0116  1.0117  1.0118  9.1026 
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2  
  1  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000   
  1  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000   
  2  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000   
  2  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000   
  3  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000   
  3  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000   

 k-point    3 :    1.50000000 1.50000000 1.00000000     weight = 1.33333333

band   1 # energy  -6.50000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  1.0120  1.0121  1.0122  1.0123  1.0124  1.0125  1.0126  1.0127  1.0128  9.1116 
  2  1.0130  1.0131  1.0132  1.0133  1.0134  1.0135  1.0136  1.0137  1.0138  9.1206 
  3  1.0140  1.0141  1.0142  1.0143  1.0144  1.0145  1.0146  1.0147  1.0148  9.1296 
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  1  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  2  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  2  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  3  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  3  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000

band   2 # energy  -1.500000000 # occ.  1.00000000

ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2     tot
  1  1.0150  1.0151  1.0152  1.0153  1.0154  1.0155  1.0156  1.0157  1.0158  9.1386
  2  1.0160  1.0161  1.0162  1.0163  1.0164  1.0165  1.0166  1.0167  1.0168  9.1476
  3  1.0170  1.0171  1.0172  1.0173  1.0174  1.0175  1.0176  1.0177  1.0178  9.1566
tot  0.0020  0.0021  0.0022  0.0023  0.0024  0.0025  0.0026  0.0027  0.0028  0.0216
ion       s      py      pz      px     dxy     dyz     dz2     dxz     dx2
  1  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  1  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000   
  2  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  2  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  3  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
  3  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
"""
