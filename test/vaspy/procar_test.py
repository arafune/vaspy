#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
        eq_(self.singleband.kvectors, self.singleprocar.kvectors)
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
    def test_singleprocar_band_compose_sites(self):
        '''test for Band_with_projection.compose_sites'''
        self.singleband.compose_sites((0, 2))
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
    def test_singleprocar_band_compose_orbital0(self):
        '''test for Band_with_projection.compose_orbital (0)

        raise RuntimeError when no item in sitecomposed
        '''
        self.singleband.compose_orbital(('p', 'pxpy', 'd'))

    @with_setup(setup=setup)
    def test_singleprocar_band_compose_orbital1(self):
        '''test for Band_with_projection.compose_orbital ()

        raise RuntimeError when no item in sitecomposed
        '''
        self.singleband.compose_sites((0, 2))
        self.singleband.compose_orbital(('p', 'pxpy', 'd'))
        np.testing.assert_allclose(self.singleband.sitecomposed,
                                   [[[[[0.0020, 0.0022, 0.0024,
                                        0.0026, 0.0028, 0.0030,
                                        0.0032, 0.0034, 0.0036,
                                        0.0252,
                                        0.0072, 0.0048, 0.0160]]]]])
        eq_(self.singleband.orb_names,
            ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2',
             'dxz', 'dx2', 'tot', 'p', 'pxpy', 'd'])


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
        eq_(self.spinband.kvectors, self.spinprocar.kvectors[0:3])
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
    def test_spinprocar_band_compose_orbital1(self):
        '''test for Band_with_projection.compose_orbital (SOI)

        raise RuntimeError when no item in sitecomposed
        '''
        self.spinband.compose_sites((0, 2))
        self.spinband.compose_orbital(('p', 'pxpy', 'd'))
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
        eq_(self.soiband.kvectors, self.soiprocar.kvectors[0:3])
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
    def test_soiprocar_band_compose_sites(self):
        '''test for Band_with_projection.compose_sites (SOI)'''
        self.soiband.compose_sites((0, 2))
        np.testing.assert_allclose(self.soiband.sitecomposed[0][0][0][0],
                                   [0.0020, 0.0022, 0.0024,
                                    0.0026, 0.0028, 0.0030,
                                    0.0032, 0.0034, 0.0036, 0.0252])

    @raises(RuntimeError)
    @with_setup(setup=setup)
    def test_soiprocar_band_compose_orbital0(self):
        '''test for Band_with_projection.compose_orbital (0) (SOI)

        raise RuntimeError when no item in sitecomposed
        '''
        self.soiband.compose_orbital(('p', 'pxpy', 'd'))

    @with_setup(setup=setup)
    def test_soiprocar_band_compose_orbital1(self):
        '''test for Band_with_projection.compose_orbital (SOI)

        raise RuntimeError when no item in sitecomposed
        '''
        self.soiband.compose_sites((0, 2))
        self.soiband.compose_orbital(('p', 'pxpy', 'd'))
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

'''

    def test_procar_std_print(self):
        self.assertEqual(output_print_procar_std, self.procar_std.__str__())

    def test_procar_spin_print(self):
        self.assertEqual(output_print_procar_spin, self.procar_spin.__str__())

    def test_band_check_orb_name(self):
        self.assertEqual('sp', self.nullband.check_orb_name('sp'))
        self.assertEqual('pxpy', self.nullband.check_orb_name('pypx'))
        self.assertRaises(ValueError, self.nullband.check_orb_name, 'nonexist')

    def test_BScreation_fromnull(self):
        self.nullband.kvectors = self.procar_std.kvectors
        np.testing.assert_array_equal([np.array([0.0, 0.0, 0.0]),
                                       np.array([0.25, 0.25, 0.]),
                                       np.array([0.5, 0.5, 0.])],
                                      self.nullband.kvectors)
        self.assertEqual([0.0,
                          np.sqrt(0.25 ** 2 + 0.25 ** 2),
                          np.sqrt(0.5 ** 2 + 0.5 ** 2)],
                         self.nullband.kdistance)
        self.assertEqual(3, self.nullband.numk)

    def test_BScreation_from_procar_std(self):
        testBand = self.procar_std.band()
        self.assertEqual(3, testBand.n_atoms)
        self.assertEqual(2, testBand.n_bands)
        self.assertEqual(("",), testBand.spininfo)
        np.testing.assert_array_equal([np.array([0.0, 0.0, 0.0]),
                                       np.array([0.25, 0.25, 0.]),
                                       np.array([0.5, 0.5, 0.])],
                                      testBand.kvectors)
        self.assertEqual((0,), testBand.get_orb_index('s',))
        self.assertEqual((3, 1), testBand.get_orb_index('pxpy'))
        self.assertEqual((3, 2), testBand.get_orb_index('pzpx'))
        self.assertEqual((9,), testBand.get_orb_index('spd'))
        self.assertEqual((0, 3, 1, 2),
                         testBand.get_orb_index('sp'))

    def test_BScreation_from_procar_single(self):
        testBand = self.procar_single.band()
        self.assertEqual(3, testBand.n_atoms)
        self.assertEqual(1, testBand.n_bands)
        np.testing.assert_array_equal([np.array([0.0, 0.0, 0.0])],
                                      testBand.kvectors)
        np.testing.assert_array_equal(
            np.array([[[[0.0000, 0.0001, 0.0002, 0.0003, 0.0004,
                         0.0005, 0.0006, 0.0007, 0.0008, 0.0036],
                        [0.0010, 0.0011, 0.0012, 0.0013, 0.0014,
                         0.0015, 0.0016, 0.0017, 0.0018, 0.0126],
                        [0.0020, 0.0021, 0.0022, 0.0023, 0.0024,
                         0.0025, 0.0026, 0.0027, 0.0028, 0.0216]]]]),
            testBand.orbitals)
        np.testing.assert_array_equal(
            np.array([[[[0.0000 + 0.0010J, 0.0001 + 0.0011J,
                         0.0002 + 0.0012J, 0.0003 + 0.0013J,
                         0.0004 + 0.0014J, 0.0005 + 0.0015J,
                         0.0006 + 0.0016J, 0.0007 + 0.0017J,
                         0.0008 + 0.0018J],
                        [0.0020 + 0.0030J, 0.0021 + 0.0031J,
                         0.0022 + 0.0032J, 0.0023 + 0.0033J,
                         0.0024 + 0.0034J, 0.0025 + 0.0035J,
                         0.0026 + 0.0036J, 0.0027 + 0.0037J,
                         0.0028 + 0.0038J],
                        [0.0040 + 0.0050J, 0.0041 + 0.0051J,
                         0.0042 + 0.0052J, 0.0043 + 0.0053J,
                         0.0044 + 0.0054J, 0.0045 + 0.0055J,
                         0.0046 + 0.0056J, 0.0047 + 0.0057J,
                         0.0048 + 0.0058J]]]]),
            testBand.phases)

    def test_BScreation_from_procar_spin(self):
        testBand = self.procar_spin.band()
        self.assertEqual(3, testBand.n_atoms)
        self.assertEqual(2, testBand.n_bands)
        self.assertEqual(("_up", "_down"), testBand.spininfo)
        np.testing.assert_array_equal([np.array([0.0, 0.0, 0.0]),
                                       np.array([0.25, 0.25, 0.]),
                                       np.array([0.5, 0.5, 0.])],
                                      testBand.kvectors)
        self.assertEqual((3, 2, 3, 9), testBand.phases[0].shape)
        self.assertEqual((3, 2, 3, 9), testBand.phases[1].shape)
        self.assertEqual((3, 2, 3, 10), testBand.orbitals[0].shape)
        self.assertEqual((3, 2, 3, 10), testBand.orbitals[1].shape)

    def test_BScreation_from_procar_soi(self):
        testBand = self.procar_soi.band()
        self.assertEqual(3, testBand.n_atoms)
        self.assertEqual(2, testBand.n_bands)
        self.assertEqual(("_mT", "_mX", "_mY", "_mZ"), testBand.spininfo)
        np.testing.assert_array_equal([np.array([0.0, 0.0, 0.0]),
                                       np.array([0.25, 0.25, 0.]),
                                       np.array([0.5, 0.5, 0.])],
                                      testBand.kvectors)

    def test_siteCompose_band_from_procar_std(self):
        testBand = self.procar_std.band()
        testBand.compose_sites((0, 2))
        self.assertEqual((3, 2, 1, 10), testBand.sitecomposed[0].shape)
        np.testing.assert_array_almost_equal(
            np.array([[0.0020, 0.0022, 0.0024, 0.0026, 0.0028,
                       0.0030, 0.0032, 0.0034, 0.0036, 0.0252]]),
            testBand.sitecomposed[0][0, 0])
        np.testing.assert_array_almost_equal(
            np.array([[0.020, 0.0202,  0.0204, 0.0206, 0.0208,
                       0.021, 0.0212, 0.0214, 0.0216, 0.1872]]),
            testBand.sitecomposed[0][1, 1])

    def test_composeOrbital_from_procar_std(self):
        testband = self.procar_std.band()
        testband.compose_sites((0, 1))
        testband.compose_sites((1, 2))
        testband.compose_orbital(('s', 'pxpy', 'd'))
        self.assertEqual((3, 2, 2, 12), testband.sitecomposed[0].shape)
        np.testing.assert_array_almost_equal(np.array([0.021, 0.0212, 0.0214,
                                                       0.0216, 0.0218, 0.022,
                                                       0.0222, 0.0224, 0.0226,
                                                       0.1962, 0.0428, 0.111]),
                                             testband.sitecomposed[0][1, 1, 1])
        testheader = testband.set_header(('Ag', '1st'),
                                         (('s', 'pxpy', 'tot'), ('d', 'pxpy')))
        self.assertEqual(testheader, ['#k', 'energy',
                                      'Ag_s', 'Ag_pxpy', 'Ag_tot',
                                      '1st_d', '1st_pxpy'])
        testlist = testband.list_sitecomposed_data((('s', 'pxpy', 'tot'),
                                                    ('d', 'pxpy')))
        self.assertEqual(len(testlist[0]), 7)

    def test_composeOrbital_from_procar_spin(self):
        testband = self.procar_spin.band()
        testband.compose_sites((0, 1))
        testband.compose_sites((1, 2))
        testband.compose_orbital(('s', 'pxpy', 'd'))
        self.assertEqual((3, 2, 2, 12), testband.sitecomposed[0].shape)
        self.assertEqual((3, 2, 2, 12), testband.sitecomposed[1].shape)
        np.testing.assert_array_almost_equal(np.array([0.021, 0.0212, 0.0214,
                                                       0.0216, 0.0218, 0.022,
                                                       0.0222, 0.0224, 0.0226,
                                                       0.1962, 0.0428, 0.111]),
                                             testband.sitecomposed[0][1, 1, 1])
        testheader = testband.set_header(('Ag', '1st'),
                                         (('s', 'pxpy', 'tot'), ('d', 'pxpy')))
        self.assertEqual(testheader, ['#k', 'energy_up',
                                      'Ag_s_up', 'Ag_pxpy_up', 'Ag_tot_up',
                                      '1st_d_up', '1st_pxpy_up',
                                      'energy_down',
                                      'Ag_s_down', 'Ag_pxpy_down',
                                      'Ag_tot_down',
                                      '1st_d_down', '1st_pxpy_down'])
        testlist = testband.list_sitecomposed_data((('s', 'pxpy', 'tot'),
                                                    ('d', 'pxpy')))
        self.assertEqual(len(testlist[0]), 13)

    def test_composeOrbital_from_procar_soi(self):
        testband = self.procar_soi.band()
        testband.compose_sites((0, 1))
        testband.compose_sites((1, 2))
        testband.compose_orbital(('s', 'pxpy', 'd'))
        self.assertEqual((3, 2, 2, 12), testband.sitecomposed[0].shape)
        self.assertEqual((3, 2, 2, 12), testband.sitecomposed[1].shape)
        self.assertEqual((3, 2, 2, 12), testband.sitecomposed[2].shape)
        self.assertEqual((3, 2, 2, 12), testband.sitecomposed[3].shape)
        np.testing.assert_array_almost_equal(np.array([0.021, 0.0212, 0.0214,
                                                       0.0216, 0.0218, 0.022,
                                                       0.0222, 0.0224, 0.0226,
                                                       0.1962, 0.0428, 0.111]),
                                             testband.sitecomposed[0][1, 1, 1])
        testheader = testband.set_header(('Ag', '1st'),
                                         (('s', 'pxpy', 'tot'), ('d', 'pxpy')))
        self.assertEqual(testheader, ['#k', 'energy',
                                      'Ag_s_mT', 'Ag_pxpy_mT', 'Ag_tot_mT',
                                      '1st_d_mT', '1st_pxpy_mT',
                                      'Ag_s_mX', 'Ag_pxpy_mX', 'Ag_tot_mX',
                                      '1st_d_mX', '1st_pxpy_mX',
                                      'Ag_s_mY', 'Ag_pxpy_mY', 'Ag_tot_mY',
                                      '1st_d_mY', '1st_pxpy_mY',
                                      'Ag_s_mZ', 'Ag_pxpy_mZ', 'Ag_tot_mZ',
                                      '1st_d_mZ', '1st_pxpy_mZ'])
        testlist = testband.list_sitecomposed_data((('s', 'pxpy', 'tot'),
                                                    ('d', 'pxpy')))
        self.assertEqual(len(testlist[0]), 22)
'''
