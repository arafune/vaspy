from pathlib import Path

import numpy as np
import pytest
import vaspy
import vaspy.incar
from _pytest.fixtures import fixture
from numpy.testing import assert_array_equal

datadir = Path(__file__).parent / "data"


@pytest.fixture
def incar0() -> vaspy.incar.Incar:
    return vaspy.load(str(datadir / "INCAR.0"))


@pytest.fixture
def graphene() -> vaspy.incar.Incar:
    return vaspy.load(str(datadir / "INCAR.Graphene"))


class TestPOSCAR(object):
    def test_incar0(self, incar0: vaspy.incar.Incar) -> None:
        assert incar0["ISTART"] == (0, True)
        assert incar0["ICHARG"][1] is False
        with pytest.raises(KeyError):
            incar0["ISMER"]
        assert len(incar0) == 33

    def test_graphene(self, graphene: vaspy.incar.Incar) -> None:
        assert (
            graphene.__str__()
            == """ generic:
    SYSTEM = Graphene
    ISTART = 0
    ENCUT = 550.0
    ISMEAR = 1
    SIGMA = 0.2
    PREC = Normal
    LREAL = .FALSE.
 electron:
    EDIFF = 1.00E-06
 spin:
    ISPIN = 1
 calc:
    NCORE = 1
 output:
    LORBIT = 12
"""
        )

