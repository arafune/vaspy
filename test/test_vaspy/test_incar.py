import os
from pathlib import Path

import numpy as np
from numpy.testing import assert_array_equal
import pytest
import vaspy
import vaspy.incar


datadir = Path(__file__).parent / "data"


@pytest.fixture
def incar0() -> vaspy.incar.INCAR:
    return vaspy.load(str(datadir / "INCAR.0"))


class TestPOSCAR(object):
    def test_incar0(self, incar0):
        assert incar0.istart == (0, True)
        assert incar0.icharg is None
        assert incar0.ismear == (1, True)
