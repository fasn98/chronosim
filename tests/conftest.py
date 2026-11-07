import mpmath as mp
import pytest

from chronosim.config import DPS


@pytest.fixture(autouse=True)
def full_precision():
    with mp.workdps(DPS.value):
        yield
