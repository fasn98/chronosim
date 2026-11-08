"""Phase 3: consistency checks of the transcription of Echeverria, Klinkhammer & Thorne (1991).

Each test derives one printed equation from others printed in the same paper. A failure would
point at a misreading (or at an inconsistency of the paper, which is then reported, not fixed).
"""

import mpmath as mp
import pytest

from chronosim import ekt
from chronosim.provenance import Origin

PTS = [("0.2", "0.6", "0.3", "0.01", "0.001", 1), ("-0.3", "0.5", "0.8", "0.05", "0.02", -1),
       ("0.7", "0.9", "0.5", "0.2", "0.05", 1)]


@pytest.mark.parametrize("th,ph,lam,B,R,s", PTS)
def test_b13_b14_b15_are_b1_minus_b2(th, ph, lam, B, R, s):
    """B12 is (B1 - B2)/beta, so Q = Q1 - Q2, S = S1 - S2, P = P1 - P2 (p. 1098)."""
    c = ekt.appendix_b_coefficients(*(mp.mpf(x) for x in (th, ph, lam, B, R)), s).value
    for big, one, two in (("Q", "Q1", "Q2"), ("S", "S1", "S2"), ("P", "P1", "P2")):
        assert abs(c[big] - (c[one] - c[two])) < mp.mpf("1e-75")


@pytest.mark.parametrize("th,ph,lam,B,R,s", PTS)
def test_b11_follows_from_4_8_and_4_10(th, ph, lam, B, R, s):
    args = [mp.mpf(x) for x in (th, ph, lam, B, R)]
    v_410 = ekt.v1_from_lambda(args[2], args[0], args[1], args[3], args[4]).value
    v_b11 = ekt.b11_v1(*args).value
    assert abs(v_410 - v_b11) < mp.mpf("1e-75")


def test_b26_b27_b28_are_the_upper_edge_limit_of_b13_b15():
    """theta -> phi, B, R -> 0 in B13-B15 gives B26-B28."""
    ph = mp.mpf("0.7")
    c = ekt.appendix_b_coefficients(ph, ph, mp.mpf("0.4"), mp.mpf(0), mp.mpf(0), 1).value
    assert abs(c["Q"] + 2 * mp.cos(ph) ** 2) < mp.mpf("1e-75")
    assert abs(c["S"] - 2 * mp.sin(ph) ** 2) < mp.mpf("1e-75")


def test_b19_printed_forms_agree():
    th, ph = mp.mpf("0.2"), mp.mpf("0.6")
    first = 4 * mp.sin(2 * ph) * mp.cos(ph) * mp.cos(th) / (mp.sin(ph - th) * mp.sin(ph + th))
    second = 8 * mp.sin(ph) / (mp.cos(th) * (mp.tan(ph) ** 2 - mp.tan(th) ** 2))
    assert abs(first / second - 1) < mp.mpf("1e-75")


def test_b20_printed_forms_differ_by_a_factor_four_in_the_B_term():
    """Recorded, not fixed: the B coefficient of the first printed form of B20 is 1/4 of the second's.
    Which one the geometry follows is measured in R9 (outputs/report.md)."""
    th, ph, lam, B = mp.mpf("0.2"), mp.mpf("0.6"), mp.mpf("0.5"), mp.mpf("0.01")
    first = ekt.b20_first_form(th, ph, lam, B, 1).value - lam
    second = ekt.b20_second_form(th, ph, lam, B, 1).value - lam
    assert abs(second / first - 4) < mp.mpf("1e-70")


def test_4_12c_is_b20_second_form():
    th, ph, lam, B, R = (mp.mpf(x) for x in ("0.2", "0.6", "0.5", "0.01", "0.001"))
    assert abs(ekt.central_solution(th, ph, lam, B, R, 1).value[2] - ekt.b20_second_form(th, ph, lam, B, 1).value) < mp.mpf("1e-75")


def test_3_8_is_the_minimum_over_L_of_3_7():
    h, v1 = mp.mpf("0.7"), mp.mpf("2.5")
    n = lambda L: ekt.multitraversal_n(L, h, v1).value
    L_star = mp.findroot(lambda L: mp.diff(n, L), (mp.mpf("0.8"), mp.mpf(20)), solver="anderson")  # pole of eq. 3.7 at L = 0.64
    assert abs(n(L_star) / ekt.multitraversal_n_min(h, v1).value - 1) < mp.mpf("1e-40")


def test_b34_b35_are_the_two_limits_of_b32():
    """Far below and far above pi/4 (with R tan^2 << B), B32 tends to 4.14 and 4.15."""
    B, R, lam = mp.mpf("1e-2"), mp.mpf("1e-9"), mp.mpf("0.5")
    for ph, key in ((mp.mpf("0.4"), "4.14"), (mp.mpf("1.1"), "4.15")):
        b32 = ekt.b30_b32(ph, ph, lam, B, R).value["beta"]
        ref = ekt.upper_edge_414_415(ph, lam, B, R).value[key][0]
        assert abs(b32 / ref - 1) < mp.mpf("1e-3")


def test_closed_forms_are_derived():
    assert ekt.fig7b_psi(mp.mpf("0.5")).origin is Origin.DERIVED
    assert ekt.v1_range(mp.mpf("0.2"), mp.mpf("0.6"), mp.mpf("0.01"), mp.mpf("0.001")).origin is Origin.DERIVED
