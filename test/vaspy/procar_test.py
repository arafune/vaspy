#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Test for PROCAR class'''
# import unittest
import os
# import tempfile
from nose.tools import eq_, ok_
from nose.tools import with_setup, assert_equal, raises
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
        self.singleband = self.singleprocar.band()

    @with_setup(setup=setup)
    def test_sigleprocar_firstcheck(self):
        '''Load test for PROCAR_single'''
        eq_(('',), self.singleprocar.spininfo)
        eq_(3, self.singleprocar.n_atoms)
        eq_(1, self.singleprocar.n_bands)
        eq_(1, self.singleprocar.numk)
        eq_([-15.0], self.singleprocar.energies)
        eq_(('s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot'),
            self.singleprocar.orb_names)
        np.testing.assert_array_equal([[0.0, 0.0, 0.0]],
                                      self.singleprocar.kvectors)
        np.testing.assert_array_equal(3, len(self.singleprocar.orbital))
        np.testing.assert_array_equal([0.0000, 0.0001, 0.0002, 0.0003,
                                       0.0004, 0.0005, 0.0006, 0.0007,
                                       0.0008, 0.0036],
                                      self.singleprocar.orbital[0])
        np.testing.assert_array_equal([0.0010, 0.0011, 0.0012, 0.0013,
                                       0.0014, 0.0015, 0.0016, 0.0017,
                                       0.0018, 0.0126],
                                      self.singleprocar.orbital[1])
        np.testing.assert_array_equal([0.0020, 0.0021, 0.0022, 0.0023,
                                       0.0024, 0.0025, 0.0026, 0.0027,
                                       0.0028, 0.0216],
                                      self.singleprocar.orbital[2])

    @with_setup(setup=setup)
    def test_singleprocar_band(self):
        '''Band_with_projection object test generated from PROCAR_single'''
        eq_(self.singleband.n_bands, self.singleprocar.n_bands)
        np.testing.assert_array_equal(self.singleband.kvectors,
                                      self.singleprocar.kvectors)
        eq_(self.singleband.spininfo, self.singleprocar.spininfo)
        ok_(self.singleband.isready())
        np.testing.assert_array_equal(self.singleband.kdistance,
                                      [0.])

    @with_setup(setup=setup)
    def test_singleprocar_band_energies(self):
        '''test for Band_with_projection.energies setter'''
        np.testing.assert_array_equal(self.singleband.energies,
                                      [[-15.0]])

    def test_singleprocar_fermi_correction(self):
        '''test for Band_with_projection.fermi_correction
        '''
        self.singleband.fermi_correction(1.0)
        np.testing.assert_array_equal(self.singleband.energies,
                                      [[-16.0]])

    @with_setup(setup=setup)
    def test_singleprocar_band_orbitalread(self):
        '''test for Band_with_projection.orbitals setter'''
        np.testing.assert_array_equal(self.singleband.orbitals[0][0],
                                      self.singleprocar.orbital)

    @with_setup(setup=setup)
    def test_singleprocar_band_sum_site(self):
        '''test for Band_with_projection.sum_site'''
        self.singleband.sum_site((0, 2))
        np.testing.assert_array_almost_equal(self.singleband.sitecomposed,
                                             [[[[[0.0020, 0.0022, 0.0024,
                                                  0.0026, 0.0028, 0.0030,
                                                  0.0032, 0.0034, 0.0036,
                                                  0.0252]]]]])

        np.testing.assert_allclose(self.singleband.sitecomposed,
                                   [[[[[0.0020, 0.0022, 0.0024,
                                        0.0026, 0.0028, 0.0030,
                                        0.0032, 0.0034, 0.0036,
                                        0.0252]]]]])

    @raises(RuntimeError)
    @with_setup(setup=setup)
    def test_singleprocar_band_sum_orbital0(self):
        '''test for Band_with_projection.sum_orbital (0)

        raise RuntimeError when no item in sitecomposed
        '''
        self.singleband.sum_orbital(('p', 'pxpy', 'd'))

    @with_setup(setup=setup)
    def test_singleprocar_band_sum_orbital1(self):
        '''test for Band_with_projection.sum_orbital ()

        raise RuntimeError when no item in sitecomposed
        '''
        self.singleband.sum_site((0, 2))
        self.singleband.sum_orbital(('p', 'pxpy', 'd'))
        np.testing.assert_allclose(self.singleband.sitecomposed,
                                   [[[[[0.0020, 0.0022, 0.0024,
                                        0.0026, 0.0028, 0.0030,
                                        0.0032, 0.0034, 0.0036,
                                        0.0252,
                                        0.0072, 0.0048, 0.0160]]]]])
        eq_(self.singleband.orb_names,
            ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2',
             'dxz', 'dx2', 'tot', 'p', 'pxpy', 'd'])

    @with_setup(setup=setup)
    def test_singleprocar_band_get_orbnums(self):
        '''test for Band_with_projection.get_orbnums
        '''
        self.singleband.sum_site((0, 2))
        self.singleband.sum_orbital(('p', 'pxpy', 'd'))
        result = self.singleband.get_orbnums((('p', 'pxpy', 'd'),))
        eq_(((10, 11, 12), ), result)
        result = self.singleband.get_orbnums((('d', 'pxpy', 'p'),))
        eq_(((12, 11, 10), ), result)
        result = self.singleband.get_orbnums((['p', 'pxpy', 'd'],))
        eq_(((10, 11, 12), ), result)
        result = self.singleband.get_orbnums([('d', 'pxpy', 'p')])
        eq_(((12, 11, 10), ), result)

    @with_setup(setup=setup)
    def test_singleprocar_band_get_orb_index(self):
        '''test for Band_with_projection.get_orbnums
        '''
        self.singleband.sum_site((0, 2))
        self.singleband.sum_orbital(('p', 'pxpy', 'd'))
        eq_(self.singleband.get_orb_index('p'), (3, 1, 2))
        eq_(self.singleband.get_orb_index('d'), (4, 5, 8, 7, 6))
        eq_(self.singleband.get_orb_index('pxpy'), (3, 1))

    @raises(ValueError)
    @with_setup(setup=setup)
    def test_singleprocar_band_check_orb_name(self):
        '''test for Band_with_projection.check_orb_name
        '''
        self.singleband.sum_site((0, 2))
        eq_(procar.check_orb_name('p'), 'p')
        procar.check_orb_name('k')

    @with_setup(setup=setup)
    def test_singleprocar_setheader(self):
        '''test for Band_with_projection.set_header
        '''
        self.singleband.sum_site((0, 2))
        self.singleband.sum_orbital(('p', 'pxpy', 'd'))
        eq_(self.singleband.set_header((('test'), ), (('p', 'pxpy', 'd'), )),
            ['#k', 'energy', 'test_p', 'test_pxpy', 'test_d'])

    @with_setup(setup=setup)
    def test_to_band(self):
        '''test for simple band output'''
        single_onlyband = self.singleprocar.to_band()
        eq_(single_onlyband.__str__(),
            '''#k	Energy
0.000000000e+00	-1.500000000e+01

''')

# ------------------------------


class TestSpinPolarizedPROCAR(object):
    '''Class for Test of PROCAR module

    Use PROCAR_spin_dummy for test data (dummy PROCAR)'''
    def setup(self):
        '''PROCAR object by using file of "PROCAR_single" in data directory'''
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + "PROCAR_spin_dummy"
        self.spinprocar = procar.PROCAR(data_file)
        self.spinband = self.spinprocar.band()

    @with_setup(setup=setup)
    def test_spinprocar_firstcheck(self):
        '''Load test for PROCAR_spin_dummy'''
        eq_(('_up', '_down'), self.spinprocar.spininfo)
        eq_(3, self.spinprocar.n_atoms)
        eq_(2, self.spinprocar.n_bands)
        eq_(3, self.spinprocar.numk)
        eq_([-10.0, -5.0, -7.0, -4.0, -6.0, -1.0,
             -10.5, -5.5, -7.5, -4.5, -6.5, -1.5],
            self.spinprocar.energies)
        eq_(('s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot'),
            self.spinprocar.orb_names)
        np.testing.assert_array_equal([[0.0, 0.0, 0.0],
                                       [0.25, 0.25, 0.00],
                                       [0.50, 0.50, 0.00],
                                       [0.0, 0.0, 0.0],
                                       [1.25, 1.25, 1.00],
                                       [1.50, 1.50, 1.00]],
                                      self.spinprocar.kvectors)
        np.testing.assert_array_equal(36, len(self.spinprocar.orbital))
        np.testing.assert_array_equal([0.0000, 0.0001, 0.0002, 0.0003,
                                       0.0004, 0.0005, 0.0006, 0.0007,
                                       0.0008, 0.0036],
                                      self.spinprocar.orbital[0])
        np.testing.assert_array_equal([0.0010, 0.0011, 0.0012, 0.0013,
                                       0.0014, 0.0015, 0.0016, 0.0017,
                                       0.0018, 0.0126],
                                      self.spinprocar.orbital[1])
        np.testing.assert_array_equal([0.0020, 0.0021, 0.0022, 0.0023,
                                       0.0024, 0.0025, 0.0026, 0.0027,
                                       0.0028, 0.0216],
                                      self.spinprocar.orbital[2])

    @with_setup(setup=setup)
    def test_spinprocar_band(self):
        '''Band_with_projection object test generated from PROCAR_spin_dummy'''
        eq_(self.spinband.n_bands, self.spinprocar.n_bands)
        np.testing.assert_array_almost_equal(self.spinband.kvectors,
                                             self.spinprocar.kvectors[0:3])
        eq_(self.spinband.spininfo, self.spinprocar.spininfo)
        ok_(self.spinband.isready())
        np.testing.assert_array_almost_equal(self.spinband.kdistance,
                                             [0., 0.353553, 0.707107])

    @with_setup(setup=setup)
    def test_spinprocar_band_energies(self):
        '''test for Band_with_projection.energies setter (SPIN)
        '''
        eq_(self.spinband.energies.shape, (2, 3, 2))
        np.testing.assert_array_equal(self.spinband.energies,
                                      [[[-10, -5], [-7, -4], [-6, -1]],
                                       [[-10.5, -5.5],
                                        [-7.5, -4.5],
                                        [-6.5, -1.5]]])

    def test_spinprocar_fermi_correction(self):
        '''test for Band_with_projection.fermi_correction (SPIN)
        '''
        self.spinband.fermi_correction(1.0)
        np.testing.assert_array_equal(self.spinband.energies,
                                      [[[-11, -6], [-8, -5], [-7, -2]],
                                       [[-11.5, -6.5],
                                        [-8.5, -5.5],
                                        [-7.5, -2.5]]])

    @with_setup(setup=setup)
    def test_spinprocar_band_orbitalread(self):
        '''test for Band_with_projection.orbitals setter (SPIN)'''
        ok_(isinstance(self.spinband.orbitals, list))
        eq_(self.spinband.orbitals[0].shape, (3, 2, 3, 10))
        eq_(self.spinband.orbitals[1].shape, (3, 2, 3, 10))
        # for up spin, ik = 1, ib = 0, iatom = 2
        #            ->  (In PROCAR, k# =2, band# = 1, atom# = 3)
        #
        np.testing.assert_array_almost_equal(
            self.spinband.orbitals[0][1][0][2],
            [0.0080, 0.0081, 0.0082, 0.0083,
             0.0084, 0.0085, 0.0086, 0.0087,
             0.0088, 0.0756])

    @with_setup(setup=setup)
    def test_spinprocar_band_sum_orbital1(self):
        '''test for Band_with_projection.sum_orbital (SOI)

        raise RuntimeError when no item in sitecomposed
        '''
        self.spinband.sum_site((0, 2))
        self.spinband.sum_orbital(('p', 'pxpy', 'd'))
        np.testing.assert_allclose(self.spinband.sitecomposed[0][0][0][0],
                                   [0.0020, 0.0022, 0.0024,
                                    0.0026, 0.0028, 0.0030,
                                    0.0032, 0.0034, 0.0036, 0.0252,
                                    0.0072, 0.0048, 0.0160])
        np.testing.assert_allclose(self.spinband.sitecomposed[1][0][0][0],
                                   [2.0020, 2.0022, 2.0024,
                                    2.0026, 2.0028, 2.0030,
                                    2.0032, 2.0034, 2.0036, 18.0252,
                                    6.0072, 4.0048, 10.0160])
        eq_(self.spinband.orb_names,
            ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2',
             'dxz', 'dx2', 'tot', 'p', 'pxpy', 'd'])

    @with_setup(setup=setup)
    def test_spinprocar_setheader(self):
        '''test for Band_with_projection.set_header  (SPIN)
        '''
        self.spinband.sum_site((0, 2))
        self.spinband.sum_orbital(('p', 'pxpy', 'd'))
        eq_(self.spinband.set_header((('test'), ), (('p', 'pxpy', 'd'), )),
            ['#k', 'energy_up', 'test_p_up', 'test_pxpy_up', 'test_d_up',
             'energy_down', 'test_p_down', 'test_pxpy_down', 'test_d_down'])

    @with_setup(setup=setup)
    def test_to_band_soi(self):
        '''test for soi energies (SPIN)'''
        np.testing.assert_array_almost_equal(
            [[-10., -10.5],
             [-5., -5.5],
             [-7., -7.5],
             [-4., -4.5],
             [-6., -6.5],
             [-1., -1.5]],
            self.spinprocar.energies)
        band = self.spinprocar.to_band()
        np.testing.assert_array_almost_equal(
            [[[-10., -5.],
              [ -7., -4.],
              [ -6., -1.]],
             [[-10.5, -5.5],
              [ -7.5, -4.5],
              [ -6.5, -1.5]]],
            band.energies)

    @with_setup(setup=setup)
    def test_toband(self):
        '''test for simple band data output (Spin polarized)'''
        onlyband = self.spinprocar.to_band()
        eq_(onlyband.__str__(),
            '''#k	Energy_up	Energy_down
0.000000000e+00	-1.000000000e+01	-1.050000000e+01
3.535533906e-01	-7.000000000e+00	-7.500000000e+00
7.071067812e-01	-6.000000000e+00	-6.500000000e+00

0.000000000e+00	-5.000000000e+01	-5.500000000e+01
3.535533906e-01	-4.000000000e+00	-4.500000000e+00
7.071067812e-01	-1.000000000e+00	-1.500000000e+00

''')

# -------------------------


class TestSOIPROCAR(object):
    '''Class for Test of PROCAR module

    Use PROCAR_SOI_dummy for test data (dummy PROCAR)'''
    def setup(self):
        '''PROCAR object by using file of "PROCAR_single" in data directory'''
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        data_file = datadir + "PROCAR_SOI_dummy"
        self.soiprocar = procar.PROCAR(data_file)
        self.soiband = self.soiprocar.band()

    @with_setup(setup=setup)
    def test_soiprocar_firstcheck(self):
        '''Load test for PROCAR_SOI_dummy'''
        eq_(('_mT', '_mX', '_mY', '_mZ'), self.soiprocar.spininfo)
        eq_(3, self.soiprocar.n_atoms)
        eq_(2, self.soiprocar.n_bands)
        eq_(3, self.soiprocar.numk)
        eq_([-10.0, -5.0, -7.0, -4.0, -6.0, -1.0],
            self.soiprocar.energies)
        eq_(('s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot'),
            self.soiprocar.orb_names)
        np.testing.assert_array_equal([[0.0, 0.0, 0.0],
                                       [0.25, 0.25, 0.00],
                                       [0.50, 0.50, 0.00]],
                                      self.soiprocar.kvectors)
        # 72 = n_atom * n_bands * numk * 4
        np.testing.assert_array_equal(72, len(self.soiprocar.orbital))
        np.testing.assert_array_equal([0.0000, 0.0001, 0.0002, 0.0003,
                                       0.0004, 0.0005, 0.0006, 0.0007,
                                       0.0008, 0.0036],
                                      self.soiprocar.orbital[0])
        np.testing.assert_array_equal([0.0010, 0.0011, 0.0012, 0.0013,
                                       0.0014, 0.0015, 0.0016, 0.0017,
                                       0.0018, 0.0126],
                                      self.soiprocar.orbital[1])
        np.testing.assert_array_equal([0.0020, 0.0021, 0.0022, 0.0023,
                                       0.0024, 0.0025, 0.0026, 0.0027,
                                       0.0028, 0.0216],
                                      self.soiprocar.orbital[2])

    @with_setup(setup=setup)
    def test_soiprocar_band(self):
        '''Band_with_projection object test generated from PROCAR_soi_dummy'''
        eq_(self.soiband.n_bands, self.soiprocar.n_bands)
        np.testing.assert_array_almost_equal(self.soiband.kvectors,
                                             self.soiprocar.kvectors[0:3])
        eq_(self.soiband.spininfo, self.soiprocar.spininfo)
        ok_(self.soiband.isready())
        np.testing.assert_array_almost_equal(self.soiband.kdistance,
                                             [0., 0.353553, 0.707107])

    @with_setup(setup=setup)
    def test_soiprocar_band_energies(self):
        '''test for Band_with_projection.energies setter (SOI)'''
        eq_(self.soiband.energies.shape, (self.soiband.numk,
                                          self.soiband.n_bands))
        np.testing.assert_array_equal(self.soiband.energies,
                                      [[-10, -5], [-7, -4], [-6, -1]])

    def test_soiprocar_fermi_correction(self):
        '''test for Band_with_projection.fermi_correction
        '''
        self.soiband.fermi_correction(1.0)
        np.testing.assert_array_equal(self.soiband.energies,
                                      [[-11, -6], [-8, -5], [-7, -2]])

    @with_setup(setup=setup)
    def test_soiprocar_band_orbitalread(self):
        '''test for Band_with_projection.orbitals setter (SOI)'''
        eq_(self.soiband.orbitals.shape, (self.soiband.numk,
                                          self.soiband.n_bands,
                                          self.soiband.n_atoms * 4,
                                          10))
        # for ik = 0, ib = 1, atom=2, spin=mY,
        #                      (k#=1,  band#=2, atom#=3, spin=mY)
        np.testing.assert_array_equal(
            self.soiband.orbitals[0][1][2 + 2 * self.soiband.n_atoms],
            #                               ^ This two means "mY"
            [2.0050, 2.0051, 2.0052, 2.0053,
             2.0054, 2.0055, 2.0056, 2.0057,
             2.0058, 2.0486])
        # for ik = 1, ib = 0, atom=1, spin=mZ,
        #                      (k#=2,  band#=1, atom#=2, spin=mZ)
        np.testing.assert_array_equal(
            self.soiband.orbitals[1][0][1 + 3 * self.soiband.n_atoms],
            #                               ^ This three means "mZ"
            [3.0070, 3.0071, 3.0072, 3.0073,
             3.0074, 3.0075, 3.0076, 3.0077,
             3.0078, 3.0666])

    @with_setup(setup=setup)
    def test_soiprocar_band_sum_site(self):
        '''test for Band_with_projection.sum_site (SOI)'''
        self.soiband.sum_site((0, 2))
        np.testing.assert_allclose(self.soiband.sitecomposed[0][0][0][0],
                                   [0.0020, 0.0022, 0.0024,
                                    0.0026, 0.0028, 0.0030,
                                    0.0032, 0.0034, 0.0036, 0.0252])

    @raises(RuntimeError)
    @with_setup(setup=setup)
    def test_soiprocar_band_sum_orbital0(self):
        '''test for Band_with_projection.sum_orbital (0) (SOI)

        raise RuntimeError when no item in sitecomposed
        '''
        self.soiband.sum_orbital(('p', 'pxpy', 'd'))

    @with_setup(setup=setup)
    def test_soiprocar_band_sum_orbital1(self):
        '''test for Band_with_projection.sum_orbital (SOI)

        raise RuntimeError when no item in sitecomposed
        '''
        self.soiband.sum_site((0, 2))
        self.soiband.sum_orbital(('p', 'pxpy', 'd'))
        np.testing.assert_allclose(self.soiband.sitecomposed[0][0][0][0],  # mT
                                   [0.0020, 0.0022, 0.0024,
                                    0.0026, 0.0028, 0.0030,
                                    0.0032, 0.0034, 0.0036, 0.0252,
                                    0.0072, 0.0048, 0.0160])
        np.testing.assert_allclose(self.soiband.sitecomposed[2][0][0][0],  # mY
                                   [4.002, 4.0022, 4.0024, 4.0026, 4.0028,
                                    4.0030, 4.0032, 4.0034, 4.0036, 4.0252,
                                    12.0072, 8.0048, 20.016])
        eq_(self.soiband.orb_names,
            ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2',
             'dxz', 'dx2', 'tot', 'p', 'pxpy', 'd'])

    @with_setup(setup=setup)
    def test_soiprocar_setheader(self):
        '''test for Band_with_projection.set_header  (SOI)
        '''
        self.soiband.sum_site((0, 2))
        self.soiband.sum_orbital(('p', 'pxpy', 'd'))
        eq_(self.soiband.set_header((('test'), ), (('p', 'pxpy', 'd'), )),
            ['#k', 'energy',
             'test_p_mT', 'test_pxpy_mT', 'test_d_mT',
             'test_p_mX', 'test_pxpy_mX', 'test_d_mX',
             'test_p_mY', 'test_pxpy_mY', 'test_d_mY',
             'test_p_mZ', 'test_pxpy_mZ', 'test_d_mZ'])

    @with_setup(setup=setup)
    def test_to_band_soi(self):
        '''test for soi energies (SOI)'''
        np.testing.assert_array_almost_equal(
            [-10.0, -5.0, -7.0, -4.0, -6.0, -1.0],
            self.soiprocar.energies)
        soiband = self.soiprocar.to_band()
        np.testing.assert_array_almost_equal(
            [[-10., -5.],
             [ -7., -4.],
             [ -6., -1.]],
            soiband.energies)

    @with_setup(setup=setup)
    def test_onlyband(self):
        onlyband = self.soiprocar.to_band()
        eq_(onlyband.__str__(),
            '''#k	Energy
0.000000000e+00	-1.000000000e+01
3.535533906e-01	-7.000000000e+00
7.071067812e-01	-6.000000000e+00

0.000000000e+00	-5.000000000e+00
3.535533906e-01	-4.000000000e+00
7.071067812e-01	-1.000000000e+00

''')
