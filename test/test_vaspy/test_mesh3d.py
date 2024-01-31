#!/usr/bin/python
import os

from vaspy import mesh3d


class TestMesh3d:
    """Class for EIGENVAL class test."""

    def setup_method(self):
        os.path.abspath(os.path.dirname(__file__)) + "/data/"

    def test_check_basic_parameters(self):
        """Checking the basic parameters stored."""
        meshdata = list(range(3 * 4 * 5))
        self.simple_grid = mesh3d.Grid3D((3, 4, 5), meshdata)
