#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from pathlib import Path

import pytest
import vaspy
from vaspy.doscar import DOSCAR, PDOS


class TestDOSCAR(object):
    def setUp(self):
        currentpath = os.path.abspath(os.path.dirname(__file__))
        testDOS = currentpath + "/DOSCAR_dummy"
        self.testDOSCAR = DOSCAR(testDOS)


datadir = Path(__file__).parent / "data"


@pytest.fixture
def ispin1() -> vaspy.doscar.DOSCAR:
    return vaspy.load(str(datadir / "DOSCAR_dummy"))


@pytest.fixture
def ispin2() -> vaspy.doscar.DOSCAR:
    return vaspy.load(str(datadir / "DOSCAR_spin"))


class Test_ispin1_TDOS:
    def test_basic_property(self, ispin1: DOSCAR):
        assert ispin1.n_atom == 39
        assert len(ispin1.energies) == 301
        assert len(ispin1.tdos) == 301
        assert ispin1.tdos[0][0] == 0.0
        assert ispin1.tdos[0] == (0.0,)
        assert ispin1.tdos[-1][0] == 0.0
        assert ispin1.tdos[-13][0] == -0.2781e-16
        assert ispin1.tdos.header == ["Energy", "TDOS"]


class Test_ispin2_TDOS:
    def test_basic_property(self, ispin2: DOSCAR):
        assert ispin2.n_atom == 57
        assert len(ispin2.energies) == 31
        assert len(ispin2.tdos) == 31
        assert ispin2.tdos[0][0] == 0.0
        assert ispin2.tdos[0][0] == 0.0
        assert ispin2.tdos[-3] == (0.6508e01, -0.6492e01)
        assert ispin2.tdos.header == ["Energy", "TDOS_up", "TDOS_down"]


class Test_ispin2_PDOS(object):
    def test_pdos_sign(self, ispin2: DOSCAR):
        assert ispin2.pdoses[0][-2] == [
            0.3502e-04,
            -0.3458e-04,
            0.2131e-04,
            -0.2032e-04,
            0.3049e-03,
            -0.3362e-03,
            0.5880e-04,
            -0.5700e-04,
            0.0000e00,
            -0.0000e00,
            0.0000e00,
            -0.0000e00,
            0.0000e00,
            -0.0000e00,
            0.0000e00,
            -0.0000e00,
            0.0000e00,
            -0.0000e00,
        ]

    def test_pdos_add(self, ispin2: DOSCAR):
        twicedos: PDOS = ispin2.pdoses[0] + ispin2.pdoses[0]
        assert twicedos[-2] == [i * 2 for i in ispin2.pdoses[0][-2]]

    def test_pdos_projected(self, ispin2: DOSCAR):
        assert ispin2.pdoses[0].projected("s_down") == ispin2.pdoses[0].projected(1)
        assert ispin2.pdoses[0].projected("s_down") == [
            -0.0,
            -0.0,
            -0.0698,
            -0.02912,
            -1.173e-13,
            -4.605e-06,
            -0.00121,
            -0.0,
            -0.02841,
            -0.0,
            -0.0581,
            -8.137e-12,
            -0.01212,
            -0.004887,
            -0.01728,
            -0.0006637,
            -0.00955,
            -0.00152,
            -4.348e-07,
            -0.008619,
            -0.0009346,
            -0.0005141,
            -0.0,
            -5.045e-22,
            -0.001406,
            -0.0,
            -2.951e-18,
            -0.0,
            -0.0002834,
            -3.458e-05,
            -0.0,
        ]

    def test_pdos_total(self, ispin2: DOSCAR):
        assert ispin2.pdoses[0].total == [
            (0.0, -0.0),
            (0.0, -0.0),
            (0.07552, -0.07585),
            (0.00015729400000000002, -0.0358072),
            (1.78758e-11, -2.36926e-13),
            (0.0009279900000000001, -0.00017783),
            (0.014800899999999999, -0.015593),
            (0.0, -0.0),
            (0.0299106, -0.029691500000000003),
            (0.0, -0.0),
            (0.0770382, -0.07541179999999999),
            (1.3100700000000002e-10, -1.6231599999999998e-11),
            (0.049339999999999995, -0.04812),
            (0.013352799999999998, -0.01311858),
            (0.07481900000000001, -0.07617),
            (0.04378570000000001, -0.04423870000000001),
            (0.10165099999999999, -0.095474),
            (0.030898, -0.031387),
            (0.04139104, -0.02879246758),
            (0.061874, -0.055686),
            (0.01541000003596937, -0.0018458490599999999),
            (0.0051577, -0.0032868999999999997),
            (0.0, -0.0),
            (0.043172456399999996, -0.0432),
            (7.929000000000002e-11, -0.004503000001225),
            (0.0, -0.0),
            (0.02275000000000016, -0.022040000000000025),
            (0.0, -0.0),
            (0.0317891, -0.0310765),
            (0.00042003, -0.0004481),
            (0.0, -0.0),
        ]


class TestTDOS(object):
    pass


class TestPDOS(object):
    pass
