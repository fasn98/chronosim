"""Cross-checks of the numerical point-splitting against exact results."""

import mpmath as mp
import pytest

from chronosim.misner import closed_form_stress, hadamard_minkowski
from chronosim.pointsplit import image_part, misner_stress, rindler_part, stress_tensor
from chronosim.provenance import Origin

TOL = mp.mpf("1e-15")  # relative to the Rindler scale 1/(1440 pi^2 xi^4)


def _scale(xi):
    return 1 / (1440 * mp.pi ** 2 * mp.mpf(xi) ** 4)


def test_minkowski_vacuum_has_zero_stress():
    T = stress_tensor(lambda P, Q: hadamard_minkowski(P, Q) - hadamard_minkowski(P, Q), [0, 1, 0, 0])
    assert all(T[a][c] == 0 for a in range(4) for c in range(4))


def test_operator_on_a_known_smooth_function():
    # F = t t': d_t d_t' F = 1, so g^{cd} d_c d_d' F = -1 and
    # T_tt = 1/2 (2/3 - 1/6 * g_tt * (-1)) = 1/4,  T_xx = 1/2 (-1/6 * g_xx * (-1)) = 1/12
    T = stress_tensor(lambda P, Q: P[0] * Q[0], [0, 1, 0, 0])
    assert abs(T[0][0] - mp.mpf(1) / 4) < mp.mpf("1e-30")
    assert abs(T[1][1] - mp.mpf(1) / 12) < mp.mpf("1e-30")


def test_rindler_part_matches_the_b_to_infinity_limit():
    xi = mp.mpf("1.3")
    R = rindler_part([0, xi, 0, 0])
    expected = [-3, -1, -1, -1]
    for i in range(4):
        assert abs(R[i][i] - expected[i] * _scale(xi)) / _scale(xi) < TOL


@pytest.mark.parametrize("xi,b", [("0.7", "1"), ("1", "2pi"), ("1.5", "4"), ("1", "12")])
def test_image_sum_reproduces_li_gott_eq12(xi, b):
    xi = mp.mpf(xi)
    b = 2 * mp.pi if b == "2pi" else mp.mpf(b)
    num, cf = misner_stress(xi, b), closed_form_stress(xi, b)
    for k in ("rho", "p_xi", "p_y", "p_z"):
        assert num[k].origin is Origin.COMPUTED
        assert abs(num[k].value - cf[k].value) / _scale(xi) < TOL, k
    assert num["max_offdiag"].value / _scale(xi) < TOL
    assert abs(num["trace"].value) / _scale(xi) < TOL


def test_series_decays_like_n_to_minus_4():
    d = {}
    image_part([0, 1, 0, 0], mp.mpf(2), diagnostics=d)
    assert d["decay_spread"] < mp.mpf("1e-14")


def test_tail_is_needed_and_small():
    X, b = [0, 1, 0, 0], mp.mpf(1)
    with_tail = image_part(X, b)[0][0]
    without = image_part(X, b, tail=False)[0][0]
    rel = abs(with_tail - without) / abs(with_tail)
    assert mp.mpf("1e-6") < rel < mp.mpf("1e-4")


def test_step_must_be_smaller_than_split():
    with pytest.raises(ValueError):
        stress_tensor(lambda P, Q: 0, [0, 1, 0, 0], delta="1e-10", h="1e-10")


def test_horizon_divergence_from_numerics():
    r = misner_stress("0.01", 1)["rho"].value / misner_stress("0.02", 1)["rho"].value
    assert abs(r - 16) < mp.mpf("1e-10")
