#!/usr/bin/env python


from vaspy import poscar


class Test_POSCAR_Head:
    def setup_method(self, method) -> None:
        self.poscar_head = poscar.PosCarHead()
        self.poscar_head.atom_types = ["Ag", "Si"]
        self.poscar_head.atomnums = [3, 5]
        self.poscar_head.system_name = "testPOSCAR"

    def test_poscar_head(self) -> None:
        assert self.poscar_head.system_name == "testPOSCAR"
        assert self.poscar_head.site_label[0] == "#0:Ag1"
