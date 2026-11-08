"""Phase 3, second pass (scalpel): tighter versions of loose tolerances and the public functions and
error branches of billiards.py, ekt.py and phase3.py that the first pass left uncovered.

Every tolerance has the value measured on 2026-11-07 next to it, taken by running the code, not
copied from a comment. Fast tests only; drivers are not ported.
"""

import mpmath as mp
import pytest

from chronosim import billiards as bl
from chronosim import ekt
from chronosim import phase3 as p3

TIGHT = mp.mpf("1e-30")


# --- tolerances of the first pass, tightened ---------------------------------------------------------

def test_itinerary_window_widths_follow_B_over_two_per_traversal_within_five_percent():
    P, B = (0.5, -0.2655), 0.01
    width = {}
    for k in (2, 3, 4):
        w = bl.itinerary_windows(P, 1, k, B)
        assert len(w) == 1
        width[k] = w[0][1] - w[0][0]
    # measured: width[2] = 9.9625e-3, width[3]/width[2] = 5.0506e-3, width[4]/width[3] = 5.0506e-3
    assert abs(width[2] / B - 1) < 0.03
    assert abs(width[3] / width[2] / (B / 2) - 1) < 0.05
    assert abs(width[4] / width[3] / (B / 2) - 1) < 0.05


def test_fig7b_deviation_extrapolates_to_zero_as_B_goes_to_zero():
    """tan(psi)/v1 - 1 is linear in B with an intercept of 0 (eq. 3.1). With the deviations 0.0368192846
    (B = 0.02) and 0.0181468292 (B = 0.01) the Richardson intercept 2 d(0.01) - d(0.02) is -5.25e-4."""
    from test_ekt_reproduction import _fig7b
    d1, d2 = _fig7b("0.02"), _fig7b("0.01")
    assert abs(2 * d2 - d1) < mp.mpf("2e-3")
    assert abs(d1 / d2 - 2) < mp.mpf("0.05")  # measured 2.029


# --- ekt.py: functions without any test in the first pass ---------------------------------------------

def test_fig7c_psi_solves_eq_3_2():
    v1 = mp.mpf("0.5")
    psi = ekt.fig7c_psi(v1).value
    assert abs(mp.sin(psi) + mp.tan(psi) - 2 * v1) < mp.mpf("1e-60")
    assert 0 < psi < mp.pi / 2


def test_fig7c_event_times_are_ordered_and_antisymmetric_about_the_collision():
    t = ekt.fig7c_event_times(mp.mpf("0.4")).value
    assert all(a < b for a, b in zip(t, t[1:]))
    assert t[2] == 0
    assert abs(t[0] + t[4]) < TIGHT and abs(t[1] + t[3]) < TIGHT


def test_velocity_exchange_speed_is_one_minus_two_B():
    assert ekt.velocity_exchange_speed(mp.mpf("0.01")).value == mp.mpf("0.98")


def test_multitraversal_speed_is_v1_for_a_head_on_line():
    assert ekt.multitraversal_speed(mp.mpf(3), mp.mpf(0), mp.mpf(2)).value == 2  # h = 0: cos 0 = 1


def test_initial_line_is_h_equals_minus_B_sin_theta_and_psi_A_equals_theta_plus_phi():
    h, psi_a = ekt.initial_line(mp.mpf("0.2"), mp.mpf("0.6"), mp.mpf("0.01")).value
    assert abs(h + mp.mpf("0.01") * mp.sin(mp.mpf("0.2"))) < TIGHT
    assert abs(psi_a - mp.mpf("0.8")) < TIGHT


def test_lambda_from_speed_inverts_v1_from_lambda():
    args = [mp.mpf(x) for x in ("0.2", "0.6", "0.01", "0.001")]
    for lam in ("0", "0.3", "1"):
        v = ekt.v1_from_lambda(mp.mpf(lam), *args).value
        assert abs(ekt.lambda_from_speed(v, *args).value - mp.mpf(lam)) < mp.mpf("1e-70")


def test_v1_range_ends_are_lambda_zero_and_one():
    args = [mp.mpf(x) for x in ("0.2", "0.6", "0.01", "0.001")]
    lo, hi = ekt.v1_range(*args).value
    assert ekt.v1_from_lambda(mp.mpf(0), *args).value == lo
    assert ekt.v1_from_lambda(mp.mpf(1), *args).value == hi
    assert lo < hi


@pytest.mark.parametrize("th,ph,inside", [
    ("0.2", "0.6", True),
    ("0.7", "0.6", False),     # theta > phi: above the upper edge
    ("-0.7", "0.6", False),    # theta < B - phi: below the lower edge
    ("0.6", "0.6", False),     # the upper edge itself is excluded (strict inequality of eq. 4.7)
    ("0.2", "0.005", False),   # phi = B/2 is excluded
    ("0.2", "1.6", False),     # phi > pi/2
])
def test_dangerous_triangle_membership(th, ph, inside):
    assert bool(ekt.in_dangerous_triangle(mp.mpf(th), mp.mpf(ph), mp.mpf("0.01"))) is inside


def test_is_spurious_flags_each_listed_point_of_eq_4_6():
    pi, ph, tol = mp.pi, mp.mpf("0.6"), mp.mpf("1e-15")
    for a, b in ((0, 0), (pi, 0), (0, pi), (pi, pi), (2 * ph, 0)):
        assert ekt.is_spurious(mp.mpf(a), mp.mpf(b), ph, 1, mp.mpf("0.001"), tol), (a, b)


def test_is_spurious_flags_negative_rho_and_wrong_sign_of_beta_but_not_a_good_root():
    ph, tol, R = mp.mpf("0.6"), mp.mpf("1e-15"), mp.mpf("0.001")
    good = (mp.mpf("1.0"), mp.mpf("0.3"))
    assert not ekt.is_spurious(*good, ph, 1, R, tol)
    assert ekt.is_spurious(*good, ph, -1, R, tol)                                   # sign(beta) != s
    assert ekt.rho_of(mp.mpf("1.0"), mp.mpf("2.5")) < 0
    assert ekt.is_spurious(mp.mpf("1.0"), mp.mpf("2.5"), ph, 1, R, tol)             # rho < 0


def test_is_spurious_tolerance_is_a_distance_not_an_equality():
    ph, R = mp.mpf("0.6"), mp.mpf("0.001")
    near = (2 * ph + mp.mpf("1e-10"), mp.mpf("1e-10"))
    assert ekt.is_spurious(*near, ph, 1, R, mp.mpf("1e-9"))
    assert not ekt.is_spurious(*near, ph, 1, R, mp.mpf("1e-12"))


def test_solve_eq42_43_returns_none_at_a_singular_start_and_a_root_from_a_nearby_start():
    th, ph, B, R = (mp.mpf(x) for x in ("0.2", "0.6", "0.01", "0.001"))
    v1 = ekt.v1_from_lambda(mp.mpf("0.5"), th, ph, B, R).value
    assert ekt.solve_eq42_43(th, ph, v1, B, R, 1, (mp.mpf(0), mp.mpf(0))) is None   # measured: None
    a, b = ekt.solve_eq42_43(th, ph, v1, B, R, 1, (mp.mpf("1.2"), mp.mpf("0.001")))
    assert max(abs(x) for x in ekt.eq42_43(a, b, th, ph, v1, B, R, 1)) < TIGHT


def test_b32_has_no_real_root_when_its_discriminant_is_negative():
    """lam1 = -5 is outside [0, 1] (not a physical speed) but is a legitimate argument of the formula:
    measured p^2/4 - q = -2.0e-5, so beta is None."""
    c = ekt.b30_b32(mp.mpf("0.4"), mp.mpf("0.4"), mp.mpf(-5), mp.mpf("0.01"), mp.mpf("1e-3")).value
    assert c["beta"] is None
    c = ekt.b30_b32(mp.mpf("0.4"), mp.mpf("0.4"), mp.mpf("0.5"), mp.mpf("0.01"), mp.mpf("1e-3")).value
    assert c["beta"] is not None and c["beta"] > 0


def test_central_solution_class_II_is_the_mirror_of_class_I():
    th, ph, B, R = (mp.mpf(x) for x in ("0.2", "0.6", "0.01", "0.001"))
    bI, eI, lI = ekt.central_solution(th, ph, mp.mpf("0.3"), B, R, 1).value
    bII, eII, lII = ekt.central_solution(th, ph, mp.mpf("0.7"), B, R, -1).value  # 1 - lam1 = 0.3
    assert abs(bII + bI) < TIGHT and abs(eII + eI) < TIGHT
    assert abs((1 - lII) - lI) < TIGHT                                             # (1 - lam2') = lam2


# --- billiards.py ---------------------------------------------------------------------------------------

def test_follow_enters_mouth_2_goes_back_in_time_and_leaves_moving_away():
    legs, p, v, shift = bl.follow((1.0, -3.0), (0.0, 1.0), 0.1, 1, bl.FLOAT)
    assert [(leg[3], round(leg[2], 9)) for leg in legs] == [(1, 2.9)]   # mouth 2, distance 3 - B
    assert shift == -1
    assert abs(p[0] - 0.0) < 1e-12 and abs(p[1] + 0.1) < 1e-12           # leaves mouth 1 at the mirror point


def test_follow_returns_none_when_the_ray_misses_or_runs_out_of_mouths():
    assert bl.follow((0.5, -3.0), (0.0, 1.0), 0.1, 1, bl.FLOAT) is None   # passes between the mouths
    assert bl.follow((1.0, -3.0), (0.0, 1.0), 0.1, 2, bl.FLOAT) is None   # after one traversal it leaves for good


def test_roots_by_bracketing_rejects_a_pole_and_a_bracket_with_an_undefined_end():
    assert bl.roots_by_bracketing(lambda x: 1 / x if x != 0 else None, [-1.0, 1.0], 1e-10) == []   # jump, not a root
    assert bl.roots_by_bracketing(lambda x: None if x > 1.5 else x - 1.0, [0.5, 1.0, 2.0], 1e-12) == [1.0]


def test_roots_by_bracketing_finds_a_simple_root_to_the_requested_width():
    (r,) = bl.roots_by_bracketing(lambda x: x * x - 2, [0.0, 1.0, 2.0], 1e-12)
    assert abs(r - 2 ** 0.5) < 1e-11   # measured 1.4e-14


def test_roots_by_bracketing_finds_a_root_sitting_exactly_on_the_last_grid_node():
    """Was RED (a rising function with its zero on the last node returned []); fixed in billiards."""
    assert bl.roots_by_bracketing(lambda x: x - 2, [0.0, 1.0, 2.0], 1e-12) == [2.0]


def test_roots_by_bracketing_does_not_report_a_zero_on_an_interior_node_twice():
    """Was RED ([1.9999999999995453, 2.0] for a falling function with its zero on an interior node)."""
    assert len(bl.roots_by_bracketing(lambda x: 2 - x, [1.0, 2.0, 3.0], 1e-12)) == 1


def _flat_roots(v1_of_s, kick=1.0, present=lambda s: True):
    def fake(P, e, m, k, B, R, xm=bl.FLOAT, nu_grid=None, aim_grid=None, log=None):
        s = P[0]
        return [{"v1": v1_of_s(s), "aim": 0.0, "nu": 0.0, "kick": kick}] if present(s) else []
    return fake


def _crossings(monkeypatch, fake, target, max_split=None, s_end=1.0):
    monkeypatch.setattr(bl, "middle_roots", fake)
    p0, e = (0.0, 0.0), (1.0, 0.0)          # P = p0 + s e, so the fake sees P[0] = s
    table = [(s, fake((s, 0.0), e, 1, 1, 0.1, 0.01)) for s in (0.0, s_end)]
    return bl.crossings(table, p0, e, target, 1, 1, 0.1, 0.01, max_split=max_split)


def test_crossings_finds_the_s_where_v1_equals_the_target(monkeypatch):
    cands, unres, rej = _crossings(monkeypatch, _flat_roots(lambda s: s), 0.3)
    assert unres == [] and rej == [] and len(cands) == 1
    assert abs(cands[0][0] - 0.3) < 1e-9


def test_crossings_discards_a_crossing_that_needs_a_negative_kick(monkeypatch):
    cands, unres, rej = _crossings(monkeypatch, _flat_roots(lambda s: s, kick=-1.0), 0.3)
    assert cands == [] and unres == []
    assert len(rej) == 1 and "kick <= 0" in rej[0][1]   # B1: refused with its reason, not dropped


def test_crossings_discards_a_sign_change_across_a_pole_of_v1(monkeypatch):
    cands, _, rej = _crossings(monkeypatch, _flat_roots(lambda s: 1 / (s - 0.3)), 0.3)
    assert cands == []
    assert len(rej) == 1 and "pole" in rej[0][1]   # v1 changes sign at s = 0.3 by passing through infinity, not through the target


def test_crossings_reports_an_unlinked_interval_instead_of_dropping_it(monkeypatch):
    fake = _flat_roots(lambda s: 5.0, present=lambda s: s < 0.7)    # a root disappears at s = 0.7
    cands, unres, _ = _crossings(monkeypatch, fake, 99.0, max_split=1)
    assert cands == []
    assert unres == [("s_count_change", 0.5, 1.0, 1, 0)]


def test_crossings_with_max_split_zero_reports_the_whole_interval_without_halving(monkeypatch):
    """Was RED: 'max_split or MAX_SPLIT' turned 0 into the default 12."""
    fake = _flat_roots(lambda s: 5.0, present=lambda s: s < 0.7)
    cands, unres, _ = _crossings(monkeypatch, fake, 99.0, max_split=0)
    assert unres == [("s_count_change", 0.0, 1.0, 1, 0)]


def test_backend_threshold_is_where_k_log10_of_2_over_B_reaches_the_float_depth():
    # k log10(2/B) with B = 0.01 is 2.301 k: 4 -> 9.2 (< 11, float), 5 -> 11.5 (>= 11, mpmath)
    assert bl.backend_for(4, 0.01)[0] is bl.FLOAT
    xm, digits = bl.backend_for(5, 0.01)
    assert xm is bl.MP and digits == 32    # measured 32 = ceil(11.5) + 20


def test_visible_arc_is_none_for_a_point_inside_the_mouth_and_itinerary_windows_then_empty():
    assert bl.visible_arc((1.0, 0.0), 1, 0.1, bl.FLOAT) is None
    assert bl.itinerary_windows((1.0, 0.0), 1, 1, 0.1) == []


# --- phase3.py ---------------------------------------------------------------------------------------

def test_r2_has_one_velocity_and_one_mirror_normal_per_pair_of_velocities():
    n = 4
    r = p3.r2_collision_law(n=n)
    assert r["counts"] == {"velocity": n * n, "mirror": n * n, "other": 0}   # measured 16 / 16 / 0
    assert r["max_dev_direction"].value < mp.mpf("1e-60")                   # measured 2.1e-81


def test_r1_deviation_is_at_working_precision():
    assert p3.r1_traversal_rule(n=9)["max_dev"].value < mp.mpf("1e-70")     # measured 1.27e-80


def test_phase3_grid_is_sorted_bounded_and_contains_the_dense_window():
    g = p3.grid(0, 1, 5, [(0.5, 0.01)])
    assert g == sorted(g) and g[0] == 0.0 and g[-1] <= 1.0
    assert 0.5 in g and sum(1 for x in g if abs(x - 0.5) <= 0.01) >= p3.WINDOW_GRID.value - 1   # measured 159: 0.5 is shared with the uniform grid


def test_self_inconsistent_crossing_is_none_when_the_exit_line_is_parallel():
    # entering mouth 2 head-on along the axis: exit mouth 1 along the axis, same line as the entry
    assert p3.self_inconsistent_crossing((mp.mpf("0.9"), mp.mpf(0)), (mp.mpf(1), mp.mpf(0)), mp.mpf("0.1")) is None


def test_line_theta_phi_enters_mouth_2_on_its_circle():
    p0, e = p3.line_theta_phi(mp.mpf("0.2"), mp.mpf("0.6"), mp.mpf("0.01"))
    assert abs(mp.sqrt((p0[0] - 1) ** 2 + p0[1] ** 2) - mp.mpf("0.01")) < TIGHT
    assert abs(mp.atan2(e[1], e[0]) - mp.mpf("0.8")) < TIGHT


def test_dangerous_on_edge_only_when_asked():
    with pytest.raises(ValueError):
        p3.dangerous("0.6", "0.6", "0.5", "0.01", "0.001")   # theta == phi without on_edge
    with pytest.raises(ValueError):
        p3.dangerous("2.0", "2.0", "0.5", "0.01", "0.001", on_edge=True)   # phi > pi/2 is never allowed
