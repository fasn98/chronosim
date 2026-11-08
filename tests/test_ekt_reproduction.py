"""Phase 3: the independent geometric model (billiards.py) against Echeverria, Klinkhammer &
Thorne (1991). Small, fast versions of R1-R11; the full runs are in chronosim.run.

Tolerances on the geometry-vs-equation residuals are set from the polish precision (40 digits):
a residual above 1e-30 would mean a real disagreement, not rounding.
"""

import mpmath as mp
import pytest

from chronosim import billiards as bl
from chronosim import ekt
from chronosim import phase3 as p3
from chronosim.provenance import Origin
from chronosim.sources import ALL, EKT_1991

ROUNDING = mp.mpf("1e-30")


def test_source_is_registered_with_its_reading_status():
    assert EKT_1991 in ALL
    assert EKT_1991.doi == "10.1103/PhysRevD.44.1077"
    assert "twice" in EKT_1991.status


def test_r1_identification_reproduces_eq_2_1b():
    r = p3.r1_traversal_rule(n=9)
    assert r["max_dev"].value < mp.mpf("1e-60")
    assert r["max_dev"].origin is Origin.COMPUTED
    assert r["time_shift"] == -1  # mouth 2 -> mouth 1 goes back by T_d


def test_r2_speed_exchange_admits_only_velocity_and_mirror_exchange():
    r = p3.r2_collision_law(n=4)
    assert r["counts"]["other"] == 0
    assert r["counts"]["velocity"] > 0 and r["counts"]["mirror"] > 0
    assert r["max_dev_direction"].value < mp.mpf("1e-60")


def _fig7b(B, R="1e-7", v1="0.5"):
    with mp.workdps(40):
        grid = [-0.4 + 0.3 * i / 30 for i in range(31)]
        sols, unres = p3.solutions(p3.SYM_P0, p3.SYM_E, mp.mpf(v1), 1, 1, mp.mpf(B), mp.mpf(R), grid)
        good = [s for s in sols if s["ok"]]
        assert len(good) == 1, (sols, unres)
        psi = mp.atan(2 * abs(good[0]["sol"][0]))
        return mp.tan(psi) / mp.mpf(v1) - 1


def test_r4_fig7b_tends_to_eq_3_1_linearly_in_B():
    # measured on 2026-11-07: deviation 0.0368 (B = 0.02) and 0.0181 (B = 0.01)
    d1, d2 = _fig7b("0.02"), _fig7b("0.01")
    assert 0 < d2 < d1 < mp.mpf("0.05")
    assert abs(d1 / d2 - 2) < mp.mpf("0.1")


def test_ball_at_rest_on_the_line_of_centers_is_a_solution_for_any_position():
    for x in ("0.3", "0.6"):
        with mp.workdps(40):
            ok, fails = bl.verify((mp.mpf(x), mp.mpf(0)), p3.SYM_E, 0, 0, 0, mp.pi, 1, 1, mp.mpf("0.01"), mp.mpf("1e-6"))
        assert ok, fails


def _central_pair():
    case = {}
    with mp.workdps(40):
        th, ph, lam, B, R = (mp.mpf(x) for x in ("0.2", "0.6", "0.5", "0.01", "0.001"))
        v1 = ekt.v1_from_lambda(lam, th, ph, B, R).value
        p0, e = p3.line_theta_phi(th, ph, B)
        s0 = p3.self_inconsistent_crossing(p0, e, B)
        grid = [float(s0) - 0.01 + 0.02 * i / 60 for i in range(61)]
        sols, _ = p3.solutions(p0, e, v1, 1, 1, B, R, grid)
        case.update(th=th, ph=ph, v1=v1, B=B, R=R, sols=[s for s in sols if s["ok"]])
    return case


def test_r8_geometric_solutions_satisfy_printed_eqs_4_2_and_4_3():
    c = _central_pair()
    assert sorted("I" if s["angles"]["beta"] > 0 else "II" for s in c["sols"]) == ["I", "II"]
    with mp.workdps(40):
        for s in c["sols"]:
            a = s["angles"]
            sign = 1 if a["d"] > 0 else -1
            r42, r43 = ekt.eq42_43(a["alpha"], a["beta"], c["th"], c["ph"], c["v1"], c["B"], c["R"], sign)
            assert max(abs(r42), abs(r43)) < ROUNDING
            assert abs(ekt.d_of(a["alpha"], a["beta"], sign, c["R"]) - a["d"]) < ROUNDING      # eq. 4.4b
            assert abs(ekt.rho_of(a["alpha"], a["beta"]) * c["v1"] - a["v2"]) < ROUNDING       # eq. 4.5
            assert mp.sign(a["beta"]) == sign                                                  # eq. 4.6c
            dist = ekt.appendix_a_distances(a["alpha"], a["beta"], c["th"], c["ph"], c["B"]).value
            assert abs(dist["a"] + s["sol"][0]) < ROUNDING                                     # eq. A6: a = PQ


def test_root_of_eqs_4_2_4_3_outside_the_eq_4_6_list_has_a_negative_distance():
    """(theta, phi) = (1.0, 1.05): a root not listed in eq. 4.6, which the geometry cannot realize."""
    with mp.workdps(40):
        th, ph, lam, B, R = (mp.mpf(x) for x in ("1.0", "1.05", "0.5", "0.01", "0.001"))
        v1 = ekt.v1_from_lambda(lam, th, ph, B, R).value
        a, b = ekt.solve_eq42_43(th, ph, v1, B, R, -1, (mp.mpf("2.0468"), mp.mpf("-0.1366")))
        assert not ekt.is_spurious(a, b, ph, -1, R, mp.mpf("1e-15"))
        dist = ekt.appendix_a_distances(a, b, th, ph, B).value
        assert dist["a"] < 0 and dist["a_p"] < 0


def test_verifier_rejects_a_perturbed_solution_and_the_wrong_kick_sign():
    c = _central_pair()
    s = c["sols"][0]
    with mp.workdps(40):
        p0, e = p3.line_theta_phi(c["th"], c["ph"], c["B"])
        s_, nu, aim = s["sol"]
        ok, fails = bl.verify(p0, e, c["v1"], s_ + mp.mpf("1e-6"), nu, aim, 1, 1, c["B"], c["R"])
        assert not ok and fails
        ok, fails = bl.verify(p0, e, c["v1"], s_, nu + mp.pi, aim, 1, 1, c["B"], c["R"])
        assert not ok


def test_free_history_detects_the_dangerous_case_and_accepts_a_harmless_one():
    with mp.workdps(40):
        th, ph, lam, B, R = (mp.mpf(x) for x in ("0.2", "0.6", "0.5", "0.01", "0.001"))
        v1 = ekt.v1_from_lambda(lam, th, ph, B, R).value
        p0, e = p3.line_theta_phi(th, ph, B)
        ok, _, overlaps = bl.free_history(bl.add(p0, e, -1), e, v1, B, R)
        assert not ok and overlaps
        ok, _, overlaps = bl.free_history(p3.SYM_P0, p3.SYM_E, mp.mpf("0.5"), B, R)
        assert ok and not overlaps


def test_itinerary_windows_shrink_like_B_over_two_per_traversal():
    P, B = (0.5, -0.2655), 0.01
    w = [bl.itinerary_windows(P, 1, k, B) for k in (2, 3)]
    assert len(w[0]) == 1 and len(w[1]) == 1
    width2, width3 = w[0][0][1] - w[0][0][0], w[1][0][1] - w[1][0][0]
    assert 0.5 * B < width2 < 2 * B          # measured 0.00996
    assert 0.1 * B < width3 / width2 < 2 * B  # measured 5.05e-3, i.e. about B/2


def test_dangerous_rejects_points_outside_the_triangle():
    with pytest.raises(ValueError):
        p3.dangerous("-0.3", "0.5", "0.5", "0.3", "0.135")


def test_backend_switches_to_mpmath_for_deep_itineraries():
    assert bl.backend_for(2, 0.01)[0] is bl.FLOAT
    xm, digits = bl.backend_for(12, 0.01)
    assert xm is bl.MP and digits > 12 * 2


KINDS = {"s_count_change", "s_bisect_lost", "nu_count_change", "nu_state_lost", "nu_bisect_lost"}


def test_search_reports_every_discard_with_its_kind_or_reason():
    """Reviewer B1: nothing the bracketing cannot follow is dropped silently."""
    with mp.workdps(40):
        th, ph, lam, B, R = (mp.mpf(x) for x in ("0.2", "0.6", "0.5", "0.01", "0.001"))
        v1 = ekt.v1_from_lambda(lam, th, ph, B, R).value
        p0, e = p3.line_theta_phi(th, ph, B)
        s0 = p3.self_inconsistent_crossing(p0, e, B)
        grid = [float(s0) - 0.01 + 0.02 * i / 60 for i in range(61)]
        sols, unres = p3.solutions(p0, e, v1, 1, 1, B, R, grid)
    assert all(u[0] in KINDS for u in unres)
    assert set(p3.unresolved_by_kind(unres)) <= KINDS
    assert all(s["fails"] and all(isinstance(f, str) and f for f in s["fails"]) for s in sols if not s["ok"])


def test_nu_interval_with_a_count_change_is_logged_not_skipped():
    """A nu-interval whose aim-root count changes is split, and what is left is logged."""
    calls = []

    def fake_roots(P, nu, m, k, B, R, windows, xm, aim_grid, log=None):
        calls.append(nu)
        return [0.1] if nu < 0.3 else [0.1, 0.2]

    saved = bl._aim_roots
    bl._aim_roots = fake_roots
    try:
        log = []
        pairs = bl._nu_intervals((0.5, -0.3), 1, 1, 0.01, 0.001, [(0, 1)], bl.FLOAT, None, 4, log)
    finally:
        bl._aim_roots = saved
    assert len(log) == 1 and log[0][0] == "nu_count_change"
    lo, hi = log[0][2], log[0][3]
    assert lo < 0.3 <= hi and hi - lo <= 2 * mp.pi / 4 / 2 ** bl.MAX_SPLIT.value * 1.0001
    assert all(len(a) == len(b) for _, a, _, b in pairs)


def test_r8_rows_carry_the_eq_A6_and_sign_rule_checks():
    """Reviewer B2: the run itself checks eq. A6 and the eq. 4.6c sign rule in every R8 row."""
    c = _central_pair()
    case = {"sols": [{"ok": True, "sol": s["sol"], "angles": s["angles"]} for s in c["sols"]],
            "theta": c["th"], "phi": c["ph"], "v1": c["v1"], "B": c["B"], "R": c["R"]}
    rows = p3.r8_eq42_43(case)
    assert rows and all(r["sign_rule"] for r in rows)
    assert all(r["dev_a6"].value < ROUNDING and r["dev_a6"].origin is Origin.COMPUTED for r in rows)


def test_unresolved_fraction_on_a_hand_made_grid():
    """Grid 0, 1, 2, 3 along e = x: an s-interval [0, 0.5] with full nu, plus a nu-interval of width
    pi at s = 2 (cell [1.5, 2.5]). fine = 0.5/3 + (1 x pi)/(3 x 2 pi); conservative = (0.5 + 2)/3, since an item at a grid point
    spoils both neighbouring grid intervals."""
    sg = [0.0, 1.0, 2.0, 3.0]
    unres = [("s_count_change", 0.0, 0.5, 1, 0), ("nu_count_change", (2.0, 0.0), 0.0, mp.pi, 1, 2),
             ("nu_count_change", (2.5, 0.0), 0.0, 1.0, 1, 2)]   # off the grid: counted, no cell
    f = p3.unresolved_fraction(unres, sg, (0.0, 0.0), (1.0, 0.0))
    assert abs(f["fine"] - (0.5 / 3 + 1 / 6)) < 1e-12
    assert abs(f["conservative"] - (0.5 + 2) / 3) < 1e-12   # union of [0, 0.5], [1, 3] (both grid intervals next to s = 2) and [2, 3]
    assert f["off_grid"] == 1 and f["fine"] <= f["conservative"]


def test_refine_nu_resolves_a_fold_with_one_sign_and_keeps_one_with_both(monkeypatch):
    """Two aim roots appear at nu = 0.3. If C > 0 on both, the fold is resolved; if C changes sign
    across the two branches, the leaf stays unresolved after the allowed halvings."""
    monkeypatch.setattr(bl, "itinerary_windows", lambda *a, **k: [(0.0, 1.0)])
    monkeypatch.setattr(bl, "_aim_roots", lambda P, nu, *a, **k: [0.1] if nu < 0.3 else [0.1, 0.2])
    monkeypatch.setattr(bl, "_exchange", lambda P, nu, m, aim, *a: 1.0)
    assert bl.refine_nu((0.5, -0.3), 1, 1, 0.01, 0.001, bl.FLOAT, 0.0, 1.0, 20) == []
    monkeypatch.setattr(bl, "_exchange", lambda P, nu, m, aim, *a: 1.0 if aim < 0.15 else -1.0)
    left = bl.refine_nu((0.5, -0.3), 1, 1, 0.01, 0.001, bl.FLOAT, 0.0, 1.0, 20)
    assert len(left) == 1 and left[0][0] < 0.3 <= left[0][1] and left[0][1] - left[0][0] <= 2 ** -20
    monkeypatch.setattr(bl, "_aim_roots", lambda P, nu, *a, **k: [0.1])
    monkeypatch.setattr(bl, "_exchange", lambda P, nu, m, aim, *a: nu - 0.4)
    left = bl.refine_nu((0.5, -0.3), 1, 1, 0.01, 0.001, bl.FLOAT, 0.0, 1.0, 20)
    assert len(left) == 1 and "middle root" in left[0][2]


def test_refine_s_keeps_a_leaf_whose_new_end_logged_nu_items(monkeypatch):
    """Reviewer B1-a: a nu-item logged at a point the s-refinement evaluates must not let the leaf
    pass as resolved."""
    def fake(P, e, m, k, B, R, xm=bl.FLOAT, nu_grid=None, aim_grid=None, log=None):
        s_ = P[0]
        if log is not None and 0.4 < s_ < 0.6:
            log.append(("nu_count_change", P, 0.0, 0.1, 1, 2))
        return [{"v1": 5.0, "aim": 0.0, "nu": 0.0, "kick": 1.0}] if s_ < 0.7 else []
    monkeypatch.setattr(bl, "middle_roots", fake)
    case = {"p0": (0.0, 0.0), "e": (1.0, 0.0), "B": 0.1, "R": 0.01, "k": 1, "v1": 99.0}
    left = p3._refine_s(case, 0.0, 1.0, 5, bl.FLOAT)
    assert any("logged unresolved nu-intervals" in why for _, _, why in left)
