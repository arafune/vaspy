#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os


from pathlib import Path
import pytest

import vaspy
from vaspy.doscar import DOSCAR


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
        assert ispin1.natom == 39
        assert len(ispin1.energies) == 301
        assert len(ispin1.tdos) == 301
        assert ispin1.tdos[0][0] == 0.0
        assert ispin1.tdos[0] == (0.0,)
        assert ispin1.tdos[-1][0] == 0.0
        assert ispin1.tdos[-13][0] == -0.2781e-16
        assert ispin1.tdos.header == "Energy\tTDOS"


class Test_ispin2_TDOS:
    def test_basic_property(self, ispin2: DOSCAR):
        assert ispin2.natom == 57
        assert len(ispin2.energies) == 31
        assert len(ispin2.tdos) == 31
        assert ispin2.tdos[0][0] == 0.0
        assert ispin2.tdos[0][0] == 0.0
        assert ispin2.tdos[-3] == (0.6508e01, -0.6492e01)
        assert ispin2.tdos.header == "Energy\tTDOS_up\tTDOS_down"


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


class TestTDOS(object):
    pass


class TestPDOS(object):
    pass
