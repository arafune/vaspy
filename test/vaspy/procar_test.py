#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Test for PROCAR class'''
# import unittest
import os
# import tempfile
from nose.tools import eq_, ok_
from nose.tools import with_setup, assert_equal, raises
from nose.tools import assert_false
import numpy as np
import vaspy.procar as procar


class TestSinglePROCAR(object):
    '''Class for Test of PROCAR module

    Use PROCAR_single for test data (dummy PROCAR)'''
    def setup(self):
        '''PROCAR object by using file of "PROCAR_single" in data directory'''
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + "PROCAR_single"
        self.singleprocar = procar.PROCAR(data_file)

    @with_setup(setup=setup)
    def test_sigleprocar_firstcheck(self):
        '''Load test for PROCAR_single'''
        eq_([''], self.singleprocar.label['spin'])
        eq_(3, self.singleprocar.natom)
        eq_(1, self.singleprocar.nbands)
        eq_(1, self.singleprocar.numk)
        np.testing.assert_array_equal([[[-15.0]]], self.singleprocar.energies)
        np.testing.assert_array_equal(
            ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot'],
            self.singleprocar.label['orbital'])
        np.testing.assert_array_equal([[0.0, 0.0, 0.0]],
                                      self.singleprocar.kvecs)

    @with_setup(setup=setup)
    def test_projection1(self):
        '''Test Projection class'''
        # testarray is taken from PROCAR_single
        eq_(5, self.singleprocar.proj.ndim)
        np.testing.assert_array_equal((1, 1, 1, 3, 10),
                                      self.singleprocar.proj.shape)
        np.testing.assert_array_equal([0.0000, 0.0001, 0.0002, 0.0003,
                                       0.0004, 0.0005, 0.0006, 0.0007,
                                       0.0008, 0.0036],
                                      self.singleprocar.proj[0, 0, 0, 0])
        np.testing.assert_array_equal([0.0010, 0.0011, 0.0012, 0.0013,
                                       0.0014, 0.0015, 0.0016, 0.0017,
                                       0.0018, 0.0126],
                                      self.singleprocar.proj[0, 0, 0, 1])
        np.testing.assert_array_equal([0.0020, 0.0021, 0.0022, 0.0023,
                                       0.0024, 0.0025, 0.0026, 0.0027,
                                       0.0028, 0.0216],
                                      self.singleprocar.proj[0, 0, 0, 2])

    def test_singleprocar_fermi_correction(self):
        '''test for fermi_correction
        '''
        self.singleprocar.fermi_correction(1.0)
        np.testing.assert_array_equal(self.singleprocar.energies,
                                      [[[-16.0]]])

    @with_setup(setup=setup)
    def test_singleprocar_sum_site(self):
        '''test for append_sumsite  (single)'''
        self.singleprocar.append_sumsite((0, 2), 'zero-two')
        eq_(self.singleprocar.label['site'][-1], 'zero-two')
        np.testing.assert_array_almost_equal(
            self.singleprocar.proj[0, 0, 0, 3],
            [0.0020, 0.0022, 0.0024,
             0.0026, 0.0028, 0.0030,
             0.0032, 0.0034, 0.0036,
             0.0252])

    @with_setup(setup=setup)
    def test_singleprocar_sum_orbital(self):
        '''test for append_sumorbital  (single)
        '''
        self.singleprocar.append_sumsite((0, 2), 'zero-two')
        self.singleprocar.append_sumorbital((1, 3), 'pxpy')
        eq_(self.singleprocar.label['orbital'][-1], 'pxpy')
        np.testing.assert_array_almost_equal(
            self.singleprocar.proj[0, 0, 0, 3],
            [0.0020, 0.0022, 0.0024,
             0.0026, 0.0028, 0.0030,
             0.0032, 0.0034, 0.0036,
             0.0252, 0.0048])

    @with_setup(setup=setup)
    def test_make_label(self):
        '''test for make label from PROCAR_single w/o append_sum*
        '''
        label = self.singleprocar.make_label((0, 2), ((0, 3, 1), (4, 7)))
        eq_(label,
            ['#k', 'Energy',
             '0_s', '0_px', '0_py', '2_dxy', '2_dxz'])
#            ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot'],

#     @with_setup(setup=setup)
#     def test_singleprocar_band_sum_orbital1(self):
#         '''test for Band_with_projection.sum_orbital ()

#         raise RuntimeError when no item in sitecomposed
#         '''
#         self.singleband.sum_site((0, 2))
#         self.singleband.sum_orbital(('p', 'pxpy', 'd'))
#         np.testing.assert_allclose(self.singleband.sitecomposed,
#                                    [[[[[0.0020, 0.0022, 0.0024,
#                                         0.0026, 0.0028, 0.0030,
#                                         0.0032, 0.0034, 0.0036,
#                                         0.0252,
#                                         0.0072, 0.0048, 0.0160]]]]])
#         eq_(self.singleband.orb_names,
#             ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2',
#              'dxz', 'dx2', 'tot', 'p', 'pxpy', 'd'])

#     @with_setup(setup=setup)
#     def test_singleprocar_band_orbnums(self):
#         '''test for Band_with_projection.get_orbnums
#         '''
#         self.singleband.sum_site((0, 2))
#         self.singleband.sum_orbital(('p', 'pxpy', 'd'))
#         result = self.singleband.orbnums((('p', 'pxpy', 'd'),))
#         eq_(((10, 11, 12), ), result)
#         result = self.singleband.orbnums((('d', 'pxpy', 'p'),))
#         eq_(((12, 11, 10), ), result)
#         result = self.singleband.orbnums((['p', 'pxpy', 'd'],))
#         eq_(((10, 11, 12), ), result)
#         result = self.singleband.orbnums([('d', 'pxpy', 'p')])
#         eq_(((12, 11, 10), ), result)

#     @with_setup(setup=setup)
#     def test_singleprocar_band_orb_index(self):
#         '''test for Band_with_projection.get_orbnums
#         '''
#         self.singleband.sum_site((0, 2))
#         self.singleband.sum_orbital(('p', 'pxpy', 'd'))
#         eq_(self.singleband.orb_index('p'), (3, 1, 2))
#         eq_(self.singleband.orb_index('d'), (4, 5, 8, 7, 6))
#         eq_(self.singleband.orb_index('pxpy'), (3, 1))

#     @raises(ValueError)
#     @with_setup(setup=setup)
#     def test_singleprocar_band_check_orb_name(self):
#         '''test for Band_with_projection.check_orb_name
#         '''
#         self.singleband.sum_site((0, 2))
#         eq_(procar.check_orb_name('p'), 'p')
#         procar.check_orb_name('k')

#     @with_setup(setup=setup)
#     def test_singleprocar_setheader(self):
#         '''test for Band_with_projection.set_header
#         '''
#         self.singleband.sum_site((0, 2))
#         self.singleband.sum_orbital(('p', 'pxpy', 'd'))
#         eq_(self.singleband.set_header((('test'), ), (('p', 'pxpy', 'd'), )),
#             ['#k', 'energy', 'test_p', 'test_pxpy', 'test_d'])

# #     @with_setup(setup=setup)
# #     def test_to_band(self):
# #         '''test for simple band output'''
# #         single_onlyband = self.singleprocar.to_band()
# #         eq_(single_onlyband.__str__(),
# #             '''#k	Energy
# # 0.000000000e+00	-1.500000000e+01

# # ''')

# # ------------------------------


class TestSpinPolarizedPROCAR(object):
    '''Class for Test of PROCAR module

    Use PROCAR_spin_dummy for test data (dummy PROCAR)'''
    def setup(self):
        '''PROCAR object by using file of "PROCAR_single" in data directory'''
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + "PROCAR_spin_dummy"
        self.spinprocar = procar.PROCAR(data_file)

    @with_setup(setup=setup)
    def test_spinprocar_firstcheck(self):
        '''Load test for PROCAR_spin_dummy'''
        eq_(['_up', '_down'], self.spinprocar.label['spin'])
        eq_(3, self.spinprocar.natom)
        eq_(4, self.spinprocar.nbands)
        eq_(3, self.spinprocar.numk)
        np.testing.assert_array_almost_equal(
            np.array([[[-10.,  -5.,   0.,   5.],
                       [ -7.,  -4.,  -1.,   4.],
                       [ -6.,  -1.,  -3.,   0.]],
                      [[-10.5, -5.5, -0.5, -5.5],
                       [ -7.5, -4.5, -1.5, -4.5],
                       [ -6.5, -1.5, -3.5, -0.5]]]),
            self.spinprocar.energies)
        np.testing.assert_array_equal(
            ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot'],
            self.spinprocar.label['orbital'])
        np.testing.assert_array_equal([[0.0, 0.0, 0.0],
                                       [0.25, 0.25, 0.00],
                                       [0.50, 0.50, 0.00],
                                       [0.0, 0.0, 0.0],
                                       [1.25, 1.25, 1.00],
                                       [1.50, 1.50, 1.00]],
                                      self.spinprocar.kvecs)
        np.testing.assert_array_equal(720, self.spinprocar.proj.size)
        np.testing.assert_array_equal([0.0000, 0.0001, 0.0002, 0.0003,
                                       0.0004, 0.0005, 0.0006, 0.0007,
                                       0.0008, 0.0036],
                                      self.spinprocar.proj[0, 0, 0, 0])
        np.testing.assert_array_equal([0.0010, 0.0011, 0.0012, 0.0013,
                                       0.0014, 0.0015, 0.0016, 0.0017,
                                       0.0018, 0.0126],
                                      self.spinprocar.proj[0, 0, 0, 1])
        np.testing.assert_array_equal([0.0020, 0.0021, 0.0022, 0.0023,
                                       0.0024, 0.0025, 0.0026, 0.0027,
                                       0.0028, 0.0216],
                                      self.spinprocar.proj[0, 0, 0, 2])

    @with_setup(setup=setup)
    def test_make_label(self):
        '''test for make label from PROCAR_spin w/o append_sum*
        '''
        label = self.spinprocar.make_label((0, 2), ((0, 3, 1), (4, 7)))
        eq_(label,
            ['#k', 'Energy_up','Energy_down',
             '0_up_s', '0_down_s',
             '0_up_px', '0_down_px',
             '0_up_py','0_down_py',
             '2_up_dxy', '2_down_dxy',
             '2_up_dxz', '2_down_dxz'])

        
#            ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot'],


#     @with_setup(setup=setup)
#     def test_spinprocar_band(self):
#         '''Band_with_projection object test generated from PROCAR_spin_dummy'''
#         spinband = self.spinprocar.band()
#         eq_(spinband.nbands, self.spinprocar.nbands)
#         np.testing.assert_array_almost_equal(spinband.kvecs,
#                                              self.spinprocar.kvecs[0:3])
#         eq_(spinband.spininfo, self.spinprocar.spininfo)
#         ok_(spinband.isready())
#         np.testing.assert_array_almost_equal(spinband.kdistance,
#                                              [0., 0.353553, 0.707107])
#         '''test for Band_with_projection.energies setter (SPIN)
#         '''
#         spinband = self.spinprocar.band()
#         eq_(spinband.energies.shape, (2, 3, 4))
#         np.testing.assert_array_equal(spinband.energies,
#                                       np.array([[[-10., -5., 0., 5.],
#                                                  [-7., -4., -1., 4.],
#                                                  [-6., -1., -3., 0.]],
#                                                 [[-10.5, -5.5, -0.5, -5.5],
#                                                  [-7.5, -4.5, -1.5, -4.5],
#                                                  [-6.5, -1.5, -3.5, -0.5]]]))

    def test_spinband_fermi_correction(self):
        '''test for fermi_correction (SPIN)
        '''
        self.spinprocar.fermi_correction(1.0)
        np.testing.assert_array_equal(self.spinprocar.energies,
                                      np.array([[[-11.,  -6., -1., 4.],
                                                 [ -8.,  -5., -2., 3.],
                                                 [ -7.,  -2., -4., -1.]],
                                                [[-11.5, -6.5, -1.5, -6.5],
                                                 [ -8.5, -5.5, -2.5, -5.5],
                                                 [ -7.5, -2.5, -4.5, -1.5]]]))
#     @with_setup(setup=setup)
#     def test_spinprocar_band_orbitalread(self):
#         '''test for Band_with_projection.orbitals setter (SPIN)'''
#         spinband = self.spinprocar.band()
#         ok_(isinstance(spinband.orbitals, list))
#         eq_(spinband.orbitals[0].shape, (3, 4, 3, 10))
#         eq_(spinband.orbitals[1].shape, (3, 4, 3, 10))
#         # for up spin, ik = 1, ib = 0, iatom = 2
#         #            ->  (In PROCAR, k# =2, band# = 1, atom# = 3)
#         #
#         np.testing.assert_array_almost_equal(
#             spinband.orbitals[0][1][0][2],
#             [0.0080, 0.0081, 0.0082, 0.0083,
#              0.0084, 0.0085, 0.0086, 0.0087,
#              0.0088, 0.0756])

#     @with_setup(setup=setup)
#     def test_spinprocar_band_sum_orbital1(self):
#         '''test for Band_with_projection.sum_orbital (SOI)

#         raise RuntimeError when no item in sitecomposed
#         '''
#         spinband = self.spinprocar.band()
#         spinband.sum_site((0, 2))
#         spinband.sum_orbital(('p', 'pxpy', 'd'))
#         np.testing.assert_allclose(spinband.sitecomposed[0][0][0][0],
#                                    [0.0020, 0.0022, 0.0024,
#                                     0.0026, 0.0028, 0.0030,
#                                     0.0032, 0.0034, 0.0036, 0.0252,
#                                     0.0072, 0.0048, 0.0160])
#         np.testing.assert_allclose(spinband.sitecomposed[1][0][0][0],
#                                    [2.0020, 2.0022, 2.0024,
#                                     2.0026, 2.0028, 2.0030,
#                                     2.0032, 2.0034, 2.0036, 18.0252,
#                                     6.0072, 4.0048, 10.0160])
#         eq_(spinband.orb_names,
#             ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2',
#              'dxz', 'dx2', 'tot', 'p', 'pxpy', 'd'])

#     @with_setup(setup=setup)
#     def test_spinprocar_setheader(self):
#         '''test for Band_with_projection.set_header  (SPIN)
#         '''
#         spinband = self.spinprocar.band()
#         spinband.sum_site((0, 2))
#         spinband.sum_orbital(('p', 'pxpy', 'd'))
#         eq_(spinband.set_header((('test'), ), (('p', 'pxpy', 'd'), )),
#             ['#k', 'energy_up', 'test_p_up', 'test_pxpy_up', 'test_d_up',
#              'energy_down', 'test_p_down', 'test_pxpy_down', 'test_d_down'])

# #     @with_setup(setup=setup)
# #     def test_toband(self):
# #         '''test for simple band data output (Spin polarized)'''
# #         onlyband = self.spinprocar.to_band()
# #         eq_(onlyband.__str__(),
# #             '''#k	Energy_up	Energy_down
# # 0.000000000e+00	-1.000000000e+01	-1.050000000e+01
# # 3.535533906e-01	-7.000000000e+00	-7.500000000e+00
# # 7.071067812e-01	-6.000000000e+00	-6.500000000e+00

# # 0.000000000e+00	-5.000000000e+00	-5.500000000e+00
# # 3.535533906e-01	-4.000000000e+00	-4.500000000e+00
# # 7.071067812e-01	-1.000000000e+00	-1.500000000e+00

# # 0.000000000e+00	0.000000000e+00	-1.000000000e+01
# # 3.535533906e-01	-1.000000000e+00	-1.000000000e+01
# # 7.071067812e-01	-3.000000000e+00	-1.000000000e+01

# # 0.000000000e+00	5.000000000e+00	-5.000000000e+00
# # 3.535533906e-01	4.000000000e+00	-5.000000000e+00
# # 7.071067812e-01	0.000000000e+00	-5.000000000e+00

# # ''')

# # -------------------------


class TestSOIPROCAR(object):
    '''Class for Test of PROCAR module

    Use PROCAR_SOI_dummy for test data (dummy PROCAR)'''
    def setup(self):
        '''PROCAR object by using file of "PROCAR_single" in data directory'''
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + "PROCAR_SOI_dummy"
        self.soiprocar = procar.PROCAR(data_file)

    def test_projection1(self):
        '''Test Projection class  (SOI)'''
        testarray = [
            0.0000,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0036,
            0.0010,0.0011,0.0012,0.0013,0.0014,0.0015,0.0016,0.0017,0.0018,0.0126,
            0.0020,0.0021,0.0022,0.0023,0.0024,0.0025,0.0026,0.0027,0.0028,0.0216,
            1.0000,1.0001,1.0002,1.0003,1.0004,1.0005,1.0006,1.0007,1.0008,1.0036,
            1.0010,1.0011,1.0012,1.0013,1.0014,1.0015,1.0016,1.0017,1.0018,1.0126,
            1.0020,1.0021,1.0022,1.0023,1.0024,1.0025,1.0026,1.0027,1.0028,1.0216,
            2.0000,2.0001,2.0002,2.0003,2.0004,2.0005,2.0006,2.0007,2.0008,2.0036,
            2.0010,2.0011,2.0012,2.0013,2.0014,2.0015,2.0016,2.0017,2.0018,2.0126,
            2.0020,2.0021,2.0022,2.0023,2.0024,2.0025,2.0026,2.0027,2.0028,2.0216,
            3.0000,3.0001,3.0002,3.0003,3.0004,3.0005,3.0006,3.0007,3.0008,3.0036,
            3.0010,3.0011,3.0012,3.0013,3.0014,3.0015,3.0016,3.0017,3.0018,3.0126,
            3.0020,3.0021,3.0022,3.0023,3.0024,3.0025,3.0026,3.0027,3.0028,3.0216,
            0.0030,0.0031,0.0032,0.0033,0.0034,0.0035,0.0036,0.0037,0.0038,0.0306,
            0.0040,0.0041,0.0042,0.0043,0.0044,0.0045,0.0046,0.0047,0.0048,0.0396,
            0.0050,0.0051,0.0052,0.0053,0.0054,0.0055,0.0056,0.0057,0.0058,0.0486,
            1.0030,1.0031,1.0032,1.0033,1.0034,1.0035,1.0036,1.0037,1.0038,1.0306,
            1.0040,1.0041,1.0042,1.0043,1.0044,1.0045,1.0046,1.0047,1.0048,1.0396,
            1.0050,1.0051,1.0052,1.0053,1.0054,1.0055,1.0056,1.0057,1.0058,1.0486,
            2.0030,2.0031,2.0032,2.0033,2.0034,2.0035,2.0036,2.0037,2.0038,2.0306,
            2.0040,2.0041,2.0042,2.0043,2.0044,2.0045,2.0046,2.0047,2.0048,2.0396,
            2.0050,2.0051,2.0052,2.0053,2.0054,2.0055,2.0056,2.0057,2.0058,2.0486,
            3.0030,3.0031,3.0032,3.0033,3.0034,3.0035,3.0036,3.0037,3.0038,3.0306,
            3.0040,3.0041,3.0042,3.0043,3.0044,3.0045,3.0046,3.0047,3.0048,3.0396,
            3.0050,3.0051,3.0052,3.0053,3.0054,3.0055,3.0056,3.0057,3.0058,3.0486,
            0.0060,0.0061,0.0062,0.0063,0.0064,0.0065,0.0066,0.0067,0.0068,0.0576,
            0.0070,0.0071,0.0072,0.0073,0.0074,0.0075,0.0076,0.0077,0.0078,0.0666,
            0.0080,0.0081,0.0082,0.0083,0.0084,0.0085,0.0086,0.0087,0.0088,0.0756,
            1.0060,1.0061,1.0062,1.0063,1.0064,1.0065,1.0066,1.0067,1.0068,1.0576,
            1.0070,1.0071,1.0072,1.0073,1.0074,1.0075,1.0076,1.0077,1.0078,1.0666,
            1.0080,1.0081,1.0082,1.0083,1.0084,1.0085,1.0086,1.0087,1.0088,1.0756,
            2.0060,2.0061,2.0062,2.0063,2.0064,2.0065,2.0066,2.0067,2.0068,2.0576,
            2.0070,2.0071,2.0072,2.0073,2.0074,2.0075,2.0076,2.0077,2.0078,2.0666,
            2.0080,2.0081,2.0082,2.0083,2.0084,2.0085,2.0086,2.0087,2.0088,2.0756,
            3.0060,3.0061,3.0062,3.0063,3.0064,3.0065,3.0066,3.0067,3.0068,3.0576,
            3.0070,3.0071,3.0072,3.0073,3.0074,3.0075,3.0076,3.0077,3.0078,3.0666,
            3.0080,3.0081,3.0082,3.0083,3.0084,3.0085,3.0086,3.0087,3.0088,3.0756,
            0.0090,0.0091,0.0092,0.0093,0.0094,0.0095,0.0096,0.0097,0.0098,0.0846,
            0.0100,0.0101,0.0102,0.0103,0.0104,0.0105,0.0106,0.0107,0.0108,0.0936,
            0.0110,0.0111,0.0112,0.0113,0.0114,0.0115,0.0116,0.0117,0.0118,0.1026,
            1.0090,1.0091,1.0092,1.0093,1.0094,1.0095,1.0096,1.0097,1.0098,1.0846,
            1.0100,1.0101,1.0102,1.0103,1.0104,1.0105,1.0106,1.0107,1.0108,1.0936,
            1.0110,1.0111,1.0112,1.0113,1.0114,1.0115,1.0116,1.0117,1.0118,1.1026,
            2.0090,2.0091,2.0092,2.0093,2.0094,2.0095,2.0096,2.0097,2.0098,2.0846,
            2.0100,2.0101,2.0102,2.0103,2.0104,2.0105,2.0106,2.0107,2.0108,2.0936,
            2.0110,2.0111,2.0112,2.0113,2.0114,2.0115,2.0116,2.0117,2.0118,2.1026,
            3.0090,3.0091,3.0092,3.0093,3.0094,3.0095,3.0096,3.0097,3.0098,3.0846,
            3.0100,3.0101,3.0102,3.0103,3.0104,3.0105,3.0106,3.0107,3.0108,3.0936,
            3.0110,3.0111,3.0112,3.0113,3.0114,3.0115,3.0116,3.0117,3.0118,3.1026,
            0.0120,0.0121,0.0122,0.0123,0.0124,0.0125,0.0126,0.0127,0.0128,0.1116,
            0.0130,0.0131,0.0132,0.0133,0.0134,0.0135,0.0136,0.0137,0.0138,0.1206,
            0.0140,0.0141,0.0142,0.0143,0.0144,0.0145,0.0146,0.0147,0.0148,0.1296,
            1.0120,1.0121,1.0122,1.0123,1.0124,1.0125,1.0126,1.0127,1.0128,1.1116,
            1.0130,1.0131,1.0132,1.0133,1.0134,1.0135,1.0136,1.0137,1.0138,1.1206,
            1.0140,1.0141,1.0142,1.0143,1.0144,1.0145,1.0146,1.0147,1.0148,1.1296,
            2.0120,2.0121,2.0122,2.0123,2.0124,2.0125,2.0126,2.0127,2.0128,2.1116,
            2.0130,2.0131,2.0132,2.0133,2.0134,2.0135,2.0136,2.0137,2.0138,2.1206,
            2.0140,2.0141,2.0142,2.0143,2.0144,2.0145,2.0146,2.0147,2.0148,2.1296,
            0.0120,0.0121,0.0122,0.0123,0.0124,0.0125,0.0126,0.0127,0.0128,0.1116,
            0.0130,0.0131,0.0132,0.0133,0.0134,0.0135,0.0136,0.0137,0.0138,0.1206,
            0.0140,0.0141,0.0142,0.0143,0.0144,0.0145,0.0146,0.0147,0.0148,0.1296,
            0.0150,0.0151,0.0152,0.0153,0.0154,0.0155,0.0156,0.0157,0.0158,0.1386,
            0.0160,0.0161,0.0162,0.0163,0.0164,0.0165,0.0166,0.0167,0.0168,0.1476,
            0.0170,0.0171,0.0172,0.0173,0.0174,0.0175,0.0176,0.0177,0.0178,0.1566,
            1.0150,1.0151,1.0152,1.0153,1.0154,1.0155,1.0156,1.0157,1.0158,1.1386,
            1.0160,1.0161,1.0162,1.0163,1.0164,1.0165,1.0166,1.0167,1.0168,1.1476,
            1.0170,1.0171,1.0172,1.0173,1.0174,1.0175,1.0176,1.0177,1.0178,1.1566,
            2.0150,2.0151,2.0152,2.0153,2.0154,2.0155,2.0156,2.0157,2.0158,2.1386,
            2.0160,2.0161,2.0162,2.0163,2.0164,2.0165,2.0166,2.0167,2.0168,2.1476,
            2.0170,2.0171,2.0172,2.0173,2.0174,2.0175,2.0176,2.0177,2.0178,2.1566,
            0.0150,0.0151,0.0152,0.0153,0.0154,0.0155,0.0156,0.0157,0.0158,0.1386,
            0.0160,0.0161,0.0162,0.0163,0.0164,0.0165,0.0166,0.0167,0.0168,0.1476,
            0.0170,0.0171,0.0172,0.0173,0.0174,0.0175,0.0176,0.0177,0.0178,0.1566]
        self.testProjection=procar.Projection(testarray,
                                              natom=3, nbands=2, numk=3, nspin=4)

    @with_setup(setup=setup)
    def test_soiprocar_firstcheck(self):
        '''Load test for PROCAR_SOI_dummy'''
        eq_(['_mT', '_mX', '_mY', '_mZ'], self.soiprocar.label['spin'])
        eq_(3, self.soiprocar.natom)
        eq_(2, self.soiprocar.nbands)
        eq_(3, self.soiprocar.numk)
        np.testing.assert_array_equal(np.array([[[-10.0, -5.0],
                                                 [-7.0, -4.0],
                                                 [-6.0, -1.0]]]),
                                      self.soiprocar.energies)
        np.testing.assert_array_equal(
            ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot'],
            self.soiprocar.label['orbital'])
        np.testing.assert_array_equal([[0.0, 0.0, 0.0],
                                       [0.25, 0.25, 0.00],
                                       [0.50, 0.50, 0.00]],
                                      self.soiprocar.kvecs)
        # 72 = natom * nbands * numk * 4
        np.testing.assert_array_equal(720, self.soiprocar.proj.size)
        np.testing.assert_array_equal([0.0000, 0.0001, 0.0002, 0.0003,
                                       0.0004, 0.0005, 0.0006, 0.0007,
                                       0.0008, 0.0036],
                                      self.soiprocar.proj[0, 0, 0, 0])
        np.testing.assert_array_equal([0.0010, 0.0011, 0.0012, 0.0013,
                                       0.0014, 0.0015, 0.0016, 0.0017,
                                       0.0018, 0.0126],
                                      self.soiprocar.proj[0, 0, 0, 1])
        np.testing.assert_array_equal([0.0020, 0.0021, 0.0022, 0.0023,
                                       0.0024, 0.0025, 0.0026, 0.0027,
                                       0.0028, 0.0216],
                                      self.soiprocar.proj[0, 0, 0, 2])

    @with_setup(setup=setup)
    def test_make_label(self):
        '''test for make label from PROCAR_soi w/o append_sum*
        '''
        label = self.soiprocar.make_label((0, 2), ((0, 3, 1), (4, 7)))
        eq_(label,
            ['#k', 'Energy',
             '0_mT_s', '0_mX_s','0_mY_s', '0_mZ_s',
             '0_mT_px', '0_mX_px', '0_mY_px', '0_mZ_px',
             '0_mT_py','0_mX_py','0_mY_py','0_mZ_py',
             '2_mT_dxy', '2_mX_dxy', '2_mY_dxy', '2_mZ_dxy',
             '2_mT_dxz', '2_mX_dxz', '2_mY_dxz', '2_mZ_dxz'])
#            ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot'],

#     @with_setup(setup=setup)
#     def test_soiprocar_band(self):
#         '''Band_with_projection object test generated from PROCAR_soi_dummy'''
#         eq_(self.soiband.nbands, self.soiprocar.nbands)
#         np.testing.assert_array_almost_equal(self.soiband.kvecs,
#                                              self.soiprocar.kvecs[0:3])
#         eq_(self.soiband.spininfo, self.soiprocar.spininfo)
#         ok_(self.soiband.isready())
#         np.testing.assert_array_almost_equal(self.soiband.kdistance,
#                                              [0., 0.353553, 0.707107])
#         '''test for Band_with_projection.energies setter (SOI)'''
#         eq_(self.soiband.energies.shape, (self.soiband.numk,
#                                           self.soiband.nbands))
#         np.testing.assert_array_equal(self.soiband.energies,
#                                       [[-10, -5], [-7, -4], [-6, -1]])
#         '''test for Band_with_projection.fermi_correction
#         '''
#         self.soiband.fermi_correction(1.0)
#         np.testing.assert_array_equal(self.soiband.energies,
#                                       [[-11, -6], [-8, -5], [-7, -2]])

#     @with_setup(setup=setup)
#     def test_soiprocar_band_orbitalread(self):
#         '''test for Band_with_projection.orbitals setter (SOI)'''
#         eq_(self.soiband.orbitals.shape, (self.soiband.numk,
#                                           self.soiband.nbands,
#                                           self.soiband.natom * 4,
#                                           10))
#         # for ik = 0, ib = 1, atom=2, spin=mY,
#         #                      (k#=1,  band#=2, atom#=3, spin=mY)
#         np.testing.assert_array_equal(
#             self.soiband.orbitals[0][1][2 + 2 * self.soiband.natom],
#             #                               ^ This two means "mY"
#             [2.0050, 2.0051, 2.0052, 2.0053,
#              2.0054, 2.0055, 2.0056, 2.0057,
#              2.0058, 2.0486])
#         # for ik = 1, ib = 0, atom=1, spin=mZ,
#         #                      (k#=2,  band#=1, atom#=2, spin=mZ)
#         np.testing.assert_array_equal(
#             self.soiband.orbitals[1][0][1 + 3 * self.soiband.natom],
#             #                               ^ This three means "mZ"
#             [3.0070, 3.0071, 3.0072, 3.0073,
#              3.0074, 3.0075, 3.0076, 3.0077,
#              3.0078, 3.0666])

#     @with_setup(setup=setup)
#     def test_soiprocar_band_sum_site(self):
#         '''test for Band_with_projection.sum_site (SOI)'''
#         self.soiband.sum_site((0, 2))
#         np.testing.assert_allclose(self.soiband.sitecomposed[0][0][0][0],
#                                    [0.0020, 0.0022, 0.0024,
#                                     0.0026, 0.0028, 0.0030,
#                                     0.0032, 0.0034, 0.0036, 0.0252])

#     @raises(RuntimeError)
#     @with_setup(setup=setup)
#     def test_soiprocar_band_sum_orbital0(self):
#         '''test for Band_with_projection.sum_orbital (0) (SOI)

#         raise RuntimeError when no item in sitecomposed
#         '''
#         self.soiband.sum_orbital(('p', 'pxpy', 'd'))

#     @with_setup(setup=setup)
#     def test_soiprocar_band_sum_orbital1(self):
#         '''test for Band_with_projection.sum_orbital (SOI)

#         raise RuntimeError when no item in sitecomposed
#         '''
#         self.soiband.sum_site((0, 2))
#         self.soiband.sum_orbital(('p', 'pxpy', 'd'))
#         np.testing.assert_allclose(self.soiband.sitecomposed[0][0][0][0],  # mT
#                                    [0.0020, 0.0022, 0.0024,
#                                     0.0026, 0.0028, 0.0030,
#                                     0.0032, 0.0034, 0.0036, 0.0252,
#                                     0.0072, 0.0048, 0.0160])
#         np.testing.assert_allclose(self.soiband.sitecomposed[2][0][0][0],  # mY
#                                    [4.002, 4.0022, 4.0024, 4.0026, 4.0028,
#                                     4.0030, 4.0032, 4.0034, 4.0036, 4.0252,
#                                     12.0072, 8.0048, 20.016])
#         eq_(self.soiband.orb_names,
#             ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2',
#              'dxz', 'dx2', 'tot', 'p', 'pxpy', 'd'])

#     @with_setup(setup=setup)
#     def test_soiprocar_setheader(self):
#         '''test for Band_with_projection.set_header  (SOI)
#         '''
#         self.soiband.sum_site((0, 2))
#         self.soiband.sum_orbital(('p', 'pxpy', 'd'))
#         eq_(self.soiband.set_header((('test'), ), (('p', 'pxpy', 'd'), )),
#             ['#k', 'energy',
#              'test_p_mT', 'test_pxpy_mT', 'test_d_mT',
#              'test_p_mX', 'test_pxpy_mX', 'test_d_mX',
#              'test_p_mY', 'test_pxpy_mY', 'test_d_mY',
#              'test_p_mZ', 'test_pxpy_mZ', 'test_d_mZ'])

#     @with_setup(setup=setup)
#     def test_to_band_soi(self):
#         '''test for soi energies (SOI)'''
#         np.testing.assert_array_almost_equal(
#             np.array([[[-10.0, -5.0],
#                        [-7.0, -4.0],
#                        [-6.0, -1.0]]]),
#             self.soiprocar.energies)
#         # soiband = self.soiprocar.to_band()
#         # np.testing.assert_array_almost_equal(
#         #     [[-10., -5.],
#         #      [ -7., -4.],
#         #      [ -6., -1.]],
#         #     soiband.energies)

# #     @with_setup(setup=setup)
# #     def test_onlyband(self):
# #         onlyband = self.soiprocar.to_band()
# #         eq_(onlyband.__str__(),
# #             '''#k	Energy
# # 0.000000000e+00	-1.000000000e+01
# # 3.535533906e-01	-7.000000000e+00
# # 7.071067812e-01	-6.000000000e+00

# # 0.000000000e+00	-5.000000000e+00
# # 3.535533906e-01	-4.000000000e+00
# # 7.071067812e-01	-1.000000000e+00

# # ''')


class test_functions_in_procarpy(object):
    def test_shortcheck_function(self):
        '''Test for shortfunction'''
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_single = open(datadir + "PROCAR_single")
        data_spin = open(datadir + "PROCAR_spin_dummy")

        result_single = procar.shortcheck(data_single)
        eq_(1, result_single[0])  #  numk
        eq_(1, result_single[1])  #  nbands
        eq_(3, result_single[2])  #  natom
        eq_(['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot'],
            result_single[3])
        ok_(result_single[4])     #  collinear
        result_spin = procar.shortcheck(data_spin)
        eq_(3, result_spin[0])  #  numk
        eq_(4, result_spin[1])  #  nbands
        eq_(3, result_spin[2])  #  natom
        eq_(['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot'],
            result_spin[3])
        ok_(result_spin[4])     #  collinear
        data_soi = open(datadir + "PROCAR_SOI_dummy")
        result_soi = procar.shortcheck(data_soi)
        eq_(3, result_soi[0])  #  numk
        eq_(2, result_soi[1])  #  nbands
        eq_(3, result_soi[2])  #  natom
        eq_(['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot'],
            result_soi[3])
        assert_false(result_soi[4])     #  collinear

