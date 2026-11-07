"""F7: arithmetic check of order-of-magnitude estimates. Not a test of the estimates."""

import mpmath as mp

from chronosim import estimates as e
from chronosim.provenance import Origin


TOL = mp.mpf("1e-12")  # these are a few multiplications of exact decimal inputs; only rounding remains


def close(a, b):
    return abs(a / b - 1) < TOL


def test_values_that_agree_with_the_printed_numbers():
    assert close(e.lp_over_D().value, e.PRINTED_LP_OVER_D.value)
    assert close(e.max_strain().value, e.PRINTED_MAX_STRAIN.value)
    assert close(e.planck_density().value, e.PRINTED_PLANCK_DENSITY_HAWKING.value)
    assert close(e.hawking_cutoff_dt().value, e.PRINTED_HAWKING_DT.value)


def test_maximum_density_differs_from_the_printed_number_by_ten_orders():
    assert close(e.max_density().value, mp.mpf("1e59"))
    assert close(e.max_density().value / e.PRINTED_MAX_DENSITY.value, mp.mpf("1e10"))


def test_quantum_gravity_ratio_is_D_over_l_P_near_a_mouth():
    assert close(e.quantum_gravity_ratio(1, 0, 1 / e.lp_over_D().value).value, mp.mpf("1e35"))


def test_origins_and_labels():
    for q in (e.lp_over_D(), e.max_density(), e.hawking_cutoff_dt()):
        assert q.origin is Origin.DERIVED
        assert "order of magnitude" in q.source
    assert e.PLANCK_LENGTH_CM.origin is Origin.CITED
