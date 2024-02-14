#!/usr/bin/env python
"""This module provides test routines for __init__.py for vaspy module."""
import os

import vaspy


class TestInitModule:
    """Class for Test of vaspy.module.

    The functions should be tested are in __init__.py
    """

    def setup_method(self) -> None:
        self.datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"

    def test_load_in_init(self) -> None:
        """Test for __init__.py."""
        poscar = self.datadir + "POSCAR_dummy"
        procar_single = self.datadir + "PROCAR_single"
        procar_spin = self.datadir + "PROCAR_spin_dummy"
        procar_soi = self.datadir + "PROCAR_soi_dummy"
        doscar = self.datadir + "DOSCAR_dummy"
        assert isinstance(vaspy.load(poscar), vaspy.poscar.POSCAR)
        assert isinstance(vaspy.load(procar_single), vaspy.procar.PROCAR)
        assert isinstance(vaspy.load(procar_spin), vaspy.procar.PROCAR)
        assert isinstance(vaspy.load(procar_soi), vaspy.procar.PROCAR)
        assert isinstance(vaspy.load(doscar), vaspy.doscar.DOSCAR)
