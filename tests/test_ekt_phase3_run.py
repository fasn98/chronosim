"""Phase 3, run/--full split and the discard log (scalpel, 2026-11-08).

Covers the code added after e7bd1d8: the logging of discards in billiards.py (reviewer B1),
unresolved_by_kind, start/collect, the quick/full task lists, report_section_quick, and the mpmath
path of solutions. Fast tests only (a few seconds each). run.main is NOT exercised: it runs phases
1 and 2 in full and writes outputs/.
"""

import mpmath as mp

from chronosim import billiards as bl
from chronosim import phase3 as p3


# --- unresolved_by_kind ---------------------------------------------------------------------------

def test_unresolved_by_kind_of_an_empty_list_is_empty():
    assert p3.unresolved_by_kind([]) == {}


def test_unresolved_by_kind_counts_each_tag_separately():
    unres = [("s_count_change", 0, 1, 1, 0), ("s_count_change", 1, 2, 2, 1), ("nu_state_lost", (0, 0), 0.1, 0.2),
             ("s_bisect_lost", 0, 1), ("nu_bisect_lost", (0, 0), 0.1, 0.2), ("nu_count_change", (0, 0), 0.1, 0.2, 1, 2)]
    assert p3.unresolved_by_kind(unres) == {"s_count_change": 2, "nu_state_lost": 1, "s_bisect_lost": 1,
                                            "nu_bisect_lost": 1, "nu_count_change": 1}


# --- B1: discards are logged ----------------------------------------------------------------------

def _flat_roots(v1_of_s, present=lambda s: True):
    def fake(P, e, m, k, B, R, xm=bl.FLOAT, nu_grid=None, aim_grid=None, log=None):
        s = P[0]
        return [{"v1": v1_of_s(s), "aim": 0.0, "nu": 0.0, "kick": 1.0}] if present(s) else []
    return fake


def test_crossings_reports_a_bisection_that_lost_its_branch_as_s_bisect_lost(monkeypatch):
    fake = _flat_roots(lambda s: s, present=lambda s: not (0.2 < s < 0.8))   # roots at both ends, none in the middle
    monkeypatch.setattr(bl, "middle_roots", fake)
    p0, e = (0.0, 0.0), (1.0, 0.0)
    table = [(s, fake((s, 0.0), e, 1, 1, 0.1, 0.01)) for s in (0.0, 1.0)]
    cands, unres, rej = bl.crossings(table, p0, e, 0.5, 1, 1, 0.1, 0.01)
    assert cands == [] and rej == []
    assert unres == [("s_bisect_lost", 0.0, 1.0)]


def _middle_roots_with_a_missing_state(monkeypatch, exchange):
    monkeypatch.setattr(bl, "_aim_roots", lambda P, nu, m, k, B, R, windows, xm, aim_grid, log=None: [0.1])
    monkeypatch.setattr(bl, "itinerary_windows", lambda *a, **kw: [(0.0, 1.0)])
    monkeypatch.setattr(bl, "_exchange", exchange)
    log = []
    out = bl.middle_roots((0.5, -0.3), (1.0, 0.0), 1, 1, 0.01, 0.001, nu_grid=4, log=log)
    return out, log


def test_middle_roots_logs_nu_state_lost_when_the_exchange_is_undefined_at_a_bracket_end(monkeypatch):
    out, log = _middle_roots_with_a_missing_state(monkeypatch, lambda *a, **kw: None)
    assert out == []
    assert log and {entry[0] for entry in log} == {"nu_state_lost"}


def test_middle_roots_without_a_log_still_returns_normally_when_the_state_is_lost(monkeypatch):
    monkeypatch.setattr(bl, "_aim_roots", lambda P, nu, m, k, B, R, windows, xm, aim_grid, log=None: [0.1])
    monkeypatch.setattr(bl, "itinerary_windows", lambda *a, **kw: [(0.0, 1.0)])
    monkeypatch.setattr(bl, "_exchange", lambda *a, **kw: None)
    assert bl.middle_roots((0.5, -0.3), (1.0, 0.0), 1, 1, 0.01, 0.001, nu_grid=4) == []


def test_log_given_to_crossings_is_the_same_list_that_comes_back_as_unresolved(monkeypatch):
    fake = _flat_roots(lambda s: 5.0, present=lambda s: s < 0.7)
    monkeypatch.setattr(bl, "middle_roots", fake)
    p0, e = (0.0, 0.0), (1.0, 0.0)
    table = [(s, fake((s, 0.0), e, 1, 1, 0.1, 0.01)) for s in (0.0, 1.0)]
    log = [("nu_state_lost", (0, 0), 0.1, 0.2)]
    _, unres, _ = bl.crossings(table, p0, e, 99.0, 1, 1, 0.1, 0.01, max_split=1, log=log)
    assert unres[0] == ("nu_state_lost", (0, 0), 0.1, 0.2)
    assert ("s_count_change", 0.5, 1.0, 1, 0) in unres


# --- N3: _circ at working precision ---------------------------------------------------------------

def test_circ_of_mpf_angles_wraps_with_mp_pi_not_the_float_pi():
    with mp.workdps(40):
        eps = mp.mpf("1e-30")
        d = bl._circ(mp.mpf(0), 2 * mp.pi - eps)
        assert abs(d - eps) < mp.mpf("1e-35")


def test_circ_of_floats_is_unchanged():
    assert abs(bl._circ(0.1, 2 * 3.141592653589793 - 0.1) - 0.2) < 1e-12


# --- solutions: the mpmath path -------------------------------------------------------------------

def test_solutions_on_the_mpmath_backend_returns_a_pair_and_restores_the_precision():
    k, B = 12, mp.mpf("0.01")
    assert bl.backend_for(k, float(B))[0] is bl.MP
    with mp.workdps(40):
        before = mp.mp.dps
        th, ph, lam, R = (mp.mpf(x) for x in ("0.2", "0.6", "0.5", "0.001"))
        p0, e = p3.line_theta_phi(th, ph, B)
        sols, unres = p3.solutions(p0, e, mp.mpf("0.5"), 1, k, B, R, [mp.mpf(-0.5), mp.mpf(-0.4)], nu_grid=3, aim_grid=3)
        assert mp.mp.dps == before
    assert isinstance(sols, list) and isinstance(unres, list)
    assert set(p3.unresolved_by_kind(unres)) <= {"s_count_change", "s_bisect_lost", "nu_count_change", "nu_state_lost", "nu_bisect_lost"}
    assert all(s["ok"] or s["fails"] for s in sols)


# --- task lists, start/collect --------------------------------------------------------------------

def test_every_quick_task_is_also_a_full_task_with_the_same_argument():
    full = p3.full_tasks()
    for t in p3.QUICK_TASKS:
        assert t in full, t


def test_full_tasks_cover_r1_to_r11():
    names = {n for n, _ in p3.full_tasks()}
    assert {"r1", "r2", "r3", "r7", "r8", "r9", "r10", "r11"} <= names


def test_quick_tasks_have_unique_names_so_collect_keys_are_unambiguous():
    names = [n for n, _ in p3.QUICK_TASKS]
    assert len(names) == len(set(names))


def test_start_then_collect_returns_results_keyed_by_task_in_order(monkeypatch):
    monkeypatch.setattr(p3, "QUICK_TASKS", (("r1", None), ("r2", None)))
    res = p3.collect(p3.start(quick=True, workers=2))
    assert list(res) == ["r1", "r2"]
    assert res["r1"][0][0] is None and res["r1"][0][1]["time_shift"] == -1
    assert res["r2"][0][1]["counts"]["other"] == 0


def test_collect_shuts_the_pool_down_even_when_a_task_fails(monkeypatch):
    monkeypatch.setattr(p3, "QUICK_TASKS", (("not_a_task", None),))
    handle = p3.start(quick=True, workers=1)
    try:
        p3.collect(handle)
        raised = False
    except ValueError:
        raised = True
    assert raised
    ex = handle[0]
    try:
        ex.submit(abs, 1)
        submitted = True
    except RuntimeError:
        submitted = False
    assert not submitted, "the executor should refuse new work after collect()"


# --- report_section_quick (formatting only, on hand-built results) --------------------------------

def _q(v, why="test"):
    return p3._c(mp.mpf(v), why)


def _quick_res(rows, physical):
    return {
        "r1": [(None, {"max_dev": _q("1e-70")})],
        "r2": [(None, {"counts": {"velocity": 3, "mirror": 2, "other": 0}})],
        "rest": [(None, [{"x": "0.3", "ok": True}])],
        "r3": [("1e-4", {"R": "1e-4", "n": 1, "dev": p3._c([mp.mpf("-1e-4")], "t"), "unresolved": 0})],
        "r8": [(("0.2", "0.6", "0.5", "0.01", "0.001"), {"rows": rows, "physical": physical, "unresolved": 2})],
        "r10": [("0.8", {"beta_I": [mp.mpf("0.001")], "b32": mp.mpf("0.001")})],
        "r11": [("0.85", {"n": 1, "beta": [mp.mpf("-0.0007")]})],
    }


def test_report_section_quick_with_no_geometric_solution_prints_dashes_not_a_crash():
    text = "\n".join(p3.report_section_quick(_quick_res([], [])))
    assert "0 soluções geométricas" in text and "sign(β) = sign(d) em 0/0" in text


def test_report_section_quick_reports_the_sign_rule_count_as_it_is():
    rows = [{"class": "I", "res": _q("1e-41"), "dev_a6": _q("1e-42"), "sign_rule": True},
            {"class": "II", "res": _q("1e-41"), "dev_a6": _q("1e-42"), "sign_rule": False}]
    text = "\n".join(p3.report_section_quick(_quick_res(rows, [(1, 1, 1)])))
    assert "sign(β) = sign(d) em 1/2" in text


def test_report_section_quick_says_it_does_not_recompute_the_full_numbers():
    text = "\n".join(p3.report_section_quick(_quick_res([], [])))
    assert "não foram recalculados" in text and "COMPUTED" in text
