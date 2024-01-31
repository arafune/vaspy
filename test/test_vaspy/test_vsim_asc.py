#! /usr/bin/env python
import os

import numpy as np

from vaspy import vsim_asc


class TestVsimASCII:
    """Class for Test of Vsim_asc."""

    def setup_method(self, method) -> None:
        """VSIM."""
        datadir = os.path.abspath(os.path.dirname(__file__)) + "/data/"
        datafile = datadir + "monolayer_hBN_phonon.ascii"
        self.hBN = vsim_asc.VSIM_ASC(datafile)

    def test_(self) -> None:
        self.hBN.system_name == "Phonopy generated file for v_sim 3.6"
        np.testing.assert_allclose([2.484999131, 0, 0], self.hBN.lattice_vectors[0])
        np.testing.assert_allclose(
            [1.242498697, 2.152072878, 0],
            self.hBN.lattice_vectors[1],
        )
        np.testing.assert_allclose([0, 0, 24.85], self.hBN.lattice_vectors[2])
        ["B", "N"] == self.hBN.atoms
        #
        #
        np.testing.assert_allclose(
            [9.354567, 17.815346, 26.162733, 32.356041, 35.832499, 39.191721],
            self.hBN.freqs,
        )
        #  mode is '3', second atom
        np.testing.assert_allclose(
            [0.000000 + -0.188937j, -0.188937 + 0.000000j, -0.000000 + 0.000000],
            self.hBN.d_vectors[3][1],
        )
        np.testing.assert_allclose(
            [-0.112619 - 0.065020j, 0.065021 - 0.112619j, 0.0 - 0.0j],
            self.hBN.d_vectors[5][1],
        )


class Test_animate_atom_phonon:
    """上手く動くようになったらテストルーチンを書いておこう)."""

    """animate_atom_phonon(position, qpt_cart, d_vector, mass=1.0,
                        n_frames=30, s_frame=0, e_frame=None,
                        magnitude=1):"""
