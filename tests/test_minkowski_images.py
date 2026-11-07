"""Minkowski vacuum with images: internal consistency of the numerics.

No test here compares with Hiscock & Konkowski (1982): the article has not been
read yet. Each tolerance sits a little above the error measured at
delta = 1e-10, h = 1e-11, dps = 80, and the measured value is quoted next to it.
"""

import mpmath as mp
import pytest

from chronosim.misner import boost, hadamard_minkowski, hadamard_minkowski_image
from chronosim.pointsplit import minkowski_image_part, minkowski_images_stress
from chronosim.provenance import Origin

COMPONENTS = ("rho", "p_xi", "p_y", "p_z")


def test_image_zero_is_the_term_the_regularization_subtracts():
    X, Y = [mp.mpf("0.1"), mp.mpf("1.2"), 0, 0], [mp.mpf("-0.2"), mp.mpf("0.9"), 0, 0]
    assert hadamard_minkowski_image(0, X, Y, mp.mpf(3)) == hadamard_minkowski(X, Y)


def test_minkowski_function_is_boost_invariant():
    X = [mp.mpf("0.1"), mp.mpf("1.2"), mp.mpf("0.3"), mp.mpf("-0.2")]
    Y = [mp.mpf("-0.2"), mp.mpf("0.9"), mp.mpf("0.1"), mp.mpf("0.4")]
    ratio = hadamard_minkowski(boost(2, X), boost(2, Y)) / hadamard_minkowski(X, Y)
    assert abs(ratio - 1) < mp.mpf("1e-78")  # measured 4.3e-80


def test_outputs_are_computed():
    out = minkowski_images_stress(1, 3)
    assert all(out[k].origin is Origin.COMPUTED for k in COMPONENTS)


def test_rejects_non_positive_b():
    with pytest.raises(ValueError):
        minkowski_image_part([0, 1, 0, 0], 0)


@pytest.mark.parametrize("b", ["1", "3", "10"])
def test_tensor_is_traceless_and_has_the_conformal_form(b):
    out = minkowski_images_stress(1, mp.mpf(b))
    rho = abs(out["rho"].value)
    assert abs(out["rho"].value - 3 * out["p_xi"].value) / rho < mp.mpf("1e-19")  # measured 1.1e-20
    assert abs(out["p_xi"].value - out["p_y"].value) / rho < mp.mpf("1e-20")  # measured 8.6e-22
    assert abs(out["p_y"].value - out["p_z"].value) / rho < mp.mpf("1e-20")
    assert abs(out["trace"].value) / rho < mp.mpf("1e-21")  # measured 3.4e-23
    assert out["max_offdiag"].value / rho < mp.mpf("1e-30")  # measured 0: the n, -n pair cancels it


def test_finite_differences_converge_at_second_order():
    rho = [
        minkowski_images_stress(1, 1, delta=d, h=h)["rho"].value
        for d, h in (("1e-8", "1e-9"), ("1e-9", "1e-10"), ("1e-10", "1e-11"))
    ]
    ratio = abs(rho[0] - rho[1]) / abs(rho[1] - rho[2])
    assert abs(ratio / 100 - 1) < mp.mpf("0.01")  # measured 2.903e-16 / 2.903e-18


def test_finite_difference_error_at_the_default_step():
    a = minkowski_images_stress(1, 1)
    c = minkowski_images_stress(1, 1, delta="1e-11", h="1e-12")
    err = max(abs(a[k].value - c[k].value) for k in COMPONENTS) / abs(a["rho"].value)
    assert err < mp.mpf("1e-19")  # measured 2.9e-20


def test_scales_as_xi_to_minus_4():
    far = minkowski_images_stress(1, 1)["rho"].value
    near = minkowski_images_stress("1e-3", 1)["rho"].value * mp.mpf("1e-3") ** 4
    assert abs(near - far) / far < mp.mpf("1e-13")  # measured 2.9e-14, the (delta / xi)^2 error


@pytest.mark.parametrize("b", ["1", "3"])
def test_tensor_is_covariant_under_the_boost(b):
    xi, b, eta = mp.mpf(1), mp.mpf(b), mp.mpf("0.7")
    T0 = minkowski_image_part([0, xi, 0, 0], b)
    T = minkowski_image_part([xi * mp.sinh(eta), xi * mp.cosh(eta), 0, 0], b)
    e_eta, e_xi = (mp.cosh(eta), mp.sinh(eta)), (mp.sinh(eta), mp.cosh(eta))
    frame = lambda u, w: sum(u[a] * w[c] * T[a][c] for a in range(2) for c in range(2))
    tol = mp.mpf("1e-20") * abs(T0[0][0])  # measured at most 3.7e-21
    assert abs(frame(e_eta, e_eta) - T0[0][0]) < tol
    assert abs(frame(e_xi, e_xi) - T0[1][1]) < tol
    assert abs(frame(e_eta, e_xi)) < tol
    assert abs(T[2][2] - T0[2][2]) < tol


@pytest.mark.parametrize("b", ["1", "3", "10", "30"])
def test_decay_is_exponential_with_rate_b_not_a_power_law(b):
    d = {}
    minkowski_image_part([0, 1, 0, 0], mp.mpf(b), diagnostics=d)
    assert abs(d["decay_rate_over_b"] - 1) < mp.mpf("1e-22")  # measured at most 8.7e-24
    assert d["power_law_spread"] > mp.mpf("1e9")  # measured at least 3.8e10; n^-4 would give ~0


def test_geometric_tail_removes_the_truncation_error():
    d = {}
    X, b = [0, 1, 0, 0], mp.mpf(1)
    default = minkowski_image_part(X, b, diagnostics=d)[0][0]
    tighter = minkowski_image_part(X, b, tol="1e-40")[0][0]
    assert mp.mpf("1e-27") < d["tail_fraction"] < mp.mpf("1e-25")  # measured 4.3e-26
    assert abs(default - tighter) / tighter < mp.mpf("1e-47")  # measured 3.3e-49


def test_does_not_vanish_at_the_li_gott_value_of_b():
    # xi^4 rho is measured at 9.6e-5, against a finite-difference error near 1e-26.
    assert minkowski_images_stress(1, 2 * mp.pi)["rho"].value > mp.mpf("1e-5")
