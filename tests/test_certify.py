"""Interval exclusion (certify.py) and the logging of the aim-search discards (reviewer note N-a)."""

import math

import mpmath as mp
import pytest
from mpmath import iv

from chronosim import billiards as bl
from chronosim import certify as ct
from chronosim import ekt, phase3 as p3
from chronosim.provenance import CERT_EXCLUSION, Origin

EMPTY = p3.R8_CASES[2]     # (0.0, 0.4, 0.75, 0.3, 0.135)
CENTRAL = p3.R8_CASES[0]   # (0.2, 0.6, 0.5, 0.01, 0.001), three verified solutions
# One of them, from the run of phase3.dangerous(*CENTRAL) (s, nu, aim); polished again below.
CENTRAL_START = ("-0.41364648775003827368", "0.19755483169590468557", "-2.539204092560386987")


# N-a ------------------------------------------------------------------------------------------

def test_a_bracket_that_becomes_undefined_inside_the_bisection_is_logged():
    log = []
    f = lambda x: None if 0.4 < x < 0.6 else x - 0.5
    assert bl.roots_by_bracketing(f, [0.0, 1.0], 1e-12, log, "here") == []
    assert log == [("bracket_lost", "here", 0.0, 1.0)]


def test_a_sign_change_rejected_as_a_jump_is_logged():
    log = []
    step = lambda x: 1.0 if x > 0.3 else -1.0   # sign change with |f| = 1 on both sides: a jump
    assert bl.roots_by_bracketing(step, [-1.0, 1.0], 1e-10, log) == []
    assert log == [("bracket_jump", None, -1.0, 1.0)]


def test_a_root_and_an_undefined_end_are_not_logged():
    log = []
    assert bl.roots_by_bracketing(lambda x: x - 0.5, [0.0, 1.0], 1e-12, log) == [pytest.approx(0.5)]
    assert bl.roots_by_bracketing(lambda x: None if x > 1.5 else x - 1.0, [0.5, 1.0, 2.0], 1e-12, log) == [1.0]
    assert log == []


def test_aim_roots_tags_its_discards_with_the_point_and_nu(monkeypatch):
    real = bl._loop_state

    def holed(P, nu, m, a, k, B, R, xm):  # the loop is undefined in a band of aim angles
        st = real(P, nu, m, a, k, B, R, xm)
        return None if st is None or 0.4 < (a % (2 * math.pi)) / (2 * math.pi) < 0.6 else st

    monkeypatch.setattr(bl, "_loop_state", holed)
    P, log = (0.5, -1.0), []
    windows = [bl.visible_arc(P, 1, 0.1, bl.FLOAT)]
    lo, hi = windows[0][0] - windows[0][1], windows[0][0] + windows[0][1]
    bl._aim_roots(P, 0.3, 1, 1, 0.1, 0.01, [(lo, hi)], bl.FLOAT, log=log)
    assert all(item[0] in p3.AIM_KINDS and item[1] == P and item[2] == 0.3 for item in log)


def test_aim_items_count_in_the_conservative_measure_with_no_area_in_the_fine_one():
    p0, e = (0.0, 0.0), (1.0, 0.0)
    sg = [0.0, 1.0, 2.0, 3.0, 4.0]
    fr = p3.unresolved_fraction([("aim_bracket_lost", (2.0, 0.0), 0.1, 0.2, 0.3)], sg, p0, e)
    assert fr["fine"] == 0 and fr["conservative"] == pytest.approx(0.5)


def test_report_labels_the_aim_kinds():
    import inspect
    src = inspect.getsource(p3.report_section)
    assert all(k in src for k in p3.AIM_KINDS)


# The model equations of the exclusion against billiards.py -------------------------------------

def test_case_inputs_are_the_initial_line_and_v1_of_eq_4_10():
    with mp.workdps(40):
        c = ct.case_inputs(*EMPTY, xm=mp)
        th, ph, lam, B, R = (mp.mpf(x) for x in EMPTY)
        p0, e = p3.line_theta_phi(th, ph, B)
        assert c["v1"] == ekt.v1_from_lambda(lam, th, ph, B, R).value
        assert c["p0"] == p0 and c["e"] == e


@pytest.mark.parametrize("case", [EMPTY, CENTRAL])
def test_conditions_are_the_loop_quantities_of_billiards_up_to_positive_factors(case):
    with mp.workdps(40):
        c = ct.case_inputs(*case, xm=mp)
        checked = 0
        for s in ("-2.5", "-1", "-0.3", "-0.02"):
            P = bl.add(c["p0"], c["e"], mp.mpf(s))
            center, half = bl.visible_arc(P, 1, c["B"], bl.MP)
            for i in range(1, 8):
                aim = center - half + 2 * half * i / 8
                for j in range(6):
                    nu = 2 * mp.pi * j / 6 + mp.mpf("0.1")
                    st = bl._loop_state(P, nu, 1, aim, 1, c["B"], c["R"], bl.MP)
                    if st is None:
                        continue
                    q = ct.conditions(mp.mpf(s), nu, aim, c, xm=mp)
                    L0 = bl.norm(bl.add(bl.aim_point(P, 1, aim, c["B"], bl.MP), P, -1), bl.MP)
                    W = st["speed"]
                    C = bl.dot(bl.add(st["v2"], st["u"], -1), st["n"])
                    V = bl.cross(st["u"], st["n"]) - c["v1"] * bl.cross(c["e"], st["n"])
                    kick = bl.dot(bl.add(st["u"], c["e"], -c["v1"]), st["n"])
                    tol = mp.mpf(10) ** -30
                    assert abs(q["miss"] - L0 * st["miss"]) < tol
                    assert abs(q["C"] * W - L0 * C) < tol
                    assert abs(q["V"] - V) < tol and abs(q["kick"] - kick) < tol
                    assert q["dn"] < 0 and q["along"] > 0
                    checked += 1
        assert checked > 50


def test_interval_conditions_enclose_the_point_values():
    with mp.workdps(40):
        cm = ct.case_inputs(*EMPTY, xm=mp)
        qm = ct.conditions(mp.mpf("-0.7"), mp.mpf("2.1"), mp.mpf("2.9"), cm, xm=mp)
    with ct.ivprec(53):
        ci = ct.case_inputs(*EMPTY)
        qi = ct.conditions(iv.mpf("-0.7"), iv.mpf("2.1"), iv.mpf("2.9"), ci)
    for k, v in qm.items():
        assert qi[k].a <= v <= qi[k].b, k


# The exclusion is not vacuous: it does not exclude a known solution -----------------------------

def _central_solution():
    th, ph, lam, B, R = (mp.mpf(x) for x in CENTRAL)
    v1 = ekt.v1_from_lambda(lam, th, ph, B, R).value
    p0, e = p3.line_theta_phi(th, ph, B)
    sol = bl.solve_mp(p0, e, *(mp.mpf(x) for x in CENTRAL_START), v1, 1, 1, B, R)
    ok, fails = bl.verify(p0, e, v1, *sol, 1, 1, B, R)
    assert ok, fails
    return sol


def test_conditions_vanish_at_a_verified_solution():
    sol = _central_solution()
    with mp.workdps(40):
        q = ct.conditions(*sol, ct.case_inputs(*CENTRAL, xm=mp), xm=mp)
    assert max(abs(q[k]) for k in ("miss", "C", "V")) < mp.mpf(10) ** -20
    assert q["dn"] < 0 and q["along"] > 0 and q["kick"] > 0


def test_a_box_around_a_verified_solution_is_not_excluded_and_the_rest_is_reported():
    sol = [float(x) for x in _central_solution()]
    h = 1e-4
    box = [(x - h, x + h) for x in sol]
    r = ct.exclude(CENTRAL, box[0], box[1], box[2], max_boxes=3000)
    assert not r["complete"] and r["undecided"]
    assert any(all(lo <= x <= hi for (lo, hi), x in zip(b, sol)) for b in r["undecided"])
    q = ct.exclusion_quantity(r, CENTRAL, box[0], "teste")
    assert q.origin is Origin.COMPUTED and q.value == len(r["undecided"])


def test_a_box_away_from_the_solutions_is_excluded_and_labelled():
    sol = [float(x) for x in _central_solution()]
    box = ((sol[0] - 0.01, sol[0] + 0.01), (sol[1] + 1.0, sol[1] + 1.2), (sol[2] - 0.05, sol[2] + 0.05))
    r = ct.exclude(CENTRAL, *box)
    assert r["complete"] and not r["undecided"]
    q = ct.exclusion_quantity(r, CENTRAL, box[0], "teste")
    assert q.origin is Origin.CERTIFIED and q.label == CERT_EXCLUSION and "hipóteses" in q.source
    assert q.value == 0 and f"{r['boxes']} subcaixas" in q.source   # the value is the number of solutions


def test_box_budget_reports_what_is_left_as_undecided():
    r = ct.exclude(EMPTY, (-3.0, -0.00135), max_boxes=10)
    assert not r["complete"] and r["boxes"] == 10 and len(r["undecided"]) > 0


def test_outward_s_range_contains_the_decimal_ends():
    lo, hi = ct.s_range_outward("-3", "-0.00135")
    assert mp.mpf(lo) < mp.mpf("-3") and mp.mpf(hi) > mp.mpf("-0.00135")


def _chunks(left, case=EMPTY, n=2, drop=None):
    lo, hi, mid = -3.0, -0.00135, -1.0
    pieces = [(lo, mid), (mid, hi)] if n == 2 else [(lo, hi)]
    out = []
    for i, sp in enumerate(pieces):
        if i == drop:
            continue
        out.append(((i, n), {"case": case, "i": i, "n": n, "s": sp, "range": (lo, hi), "boxes": 10, "excluded_by": {"miss ≠ 0": 10},
                             "undecided": [], "n_undecided": left if i == 0 else 0, "complete": not (left and i == 0)}))
    return out


def test_report_uses_the_proof_wording_only_when_the_exclusion_is_complete():
    ref = {"cert_items": {"items": 3, "items_aim": 1, "excluded": 3, "excluded_aim": 1, "certified_roots": 0, "undecided": 0, "runs": 2, "boxes": 7}}
    done = "\n".join(p3._certification_lines(EMPTY, ref, _chunks(0)))
    assert "Provado, em aritmética intervalar, que não existe solução de uma colisão para esta condição inicial no modelo de [" in done
    assert "soluções com várias colisões não foram examinadas" in done and "Não coberto pela prova" in done
    assert "Origem: CERTIFIED" in done and "soluções na caixa: 0" in done and "Krawczyk não foi implementado" in done
    partial = "\n".join(p3._certification_lines(EMPTY, ref, _chunks(2)))
    assert "Provado" not in partial and "2 indecididas" in partial


@pytest.mark.parametrize("chunks", [_chunks(0, drop=1), _chunks(0, case=CENTRAL)])
def test_no_proof_wording_when_a_piece_is_missing_or_belongs_to_another_case(chunks):
    ref = {"cert_items": {"items": 0, "items_aim": 0, "excluded": 0, "excluded_aim": 0, "certified_roots": 0, "undecided": 0, "runs": 0, "boxes": 0}}
    text = "\n".join(p3._certification_lines(EMPTY, ref, chunks))
    assert "Provado" not in text and "Nada é afirmado" in text


def test_full_run_has_the_exclusion_chunks_of_the_empty_case():
    tasks = p3.full_tasks()
    assert [a for n, a in tasks if n == "r8c"] == [(i, p3.CERT_CHUNKS.value) for i in range(p3.CERT_CHUNKS.value)]


def test_the_section_header_does_not_deny_the_certified_line():
    import inspect
    src = inspect.getsource(p3.report_section)
    assert "Todas as grandezas desta seção são **COMPUTED**." not in src
    assert "exceto a exclusão intervalar do caso vazio" in src
