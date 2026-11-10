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
    done = "\n".join(p3._certification_lines(EMPTY, ref, _chunks(0), [_summary()]))
    assert "Provado, em aritmética intervalar, que não existe solução de uma colisão para esta condição inicial no modelo de [" in done
    assert "soluções com várias colisões não foram examinadas" in done and "Não coberto pela prova" in done
    assert "Origem: CERTIFIED" in done and "soluções na caixa: 0" in done and "Krawczyk não foi implementado" not in done
    partial = "\n".join(p3._certification_lines(EMPTY, ref, _chunks(2)))
    assert "Provado" not in partial and "2 indecididas" in partial


@pytest.mark.parametrize("chunks", [_chunks(0, drop=1), _chunks(0, case=CENTRAL)])
def test_no_proof_wording_when_a_piece_is_missing_or_belongs_to_another_case(chunks):
    ref = {"cert_items": {"items": 0, "items_aim": 0, "excluded": 0, "excluded_aim": 0, "certified_roots": 0, "undecided": 0, "runs": 0, "boxes": 0}}
    text = "\n".join(p3._certification_lines(EMPTY, ref, chunks))
    assert "Provado" not in text and "Nada é afirmado" in text


def test_certify_run_has_the_exclusion_chunks_of_the_empty_case():
    """Since 2026-11-09 (user decision (d)) the certification runs in --certify, not in --full."""
    tasks = p3.certify_tasks()
    assert [a for n, a in tasks if n == "r8c"] == [(i, p3.CERT_CHUNKS.value) for i in range(p3.CERT_CHUNKS.value)]
    assert not [a for n, a in p3.full_tasks() if n == "r8c"]


def test_the_section_header_points_to_the_certification_run():
    """--full has no CERTIFIED quantity any more; its header says where they are."""
    import inspect
    src = inspect.getsource(p3.report_section)
    assert "única fonte de grandezas CERTIFIED nesta fase, roda à parte: `python -m chronosim.run --certify`" in src
    assert "outputs/phase3_certify.md" in src



def _summary(**kw):
    x = {"spec": ("r8", 0), "boxes": 100, "excluded": 90, "undecided": 10, "unexamined": 0, "known": 3, "in_excluded": 0,
         "refined": 3, "certified": 3, "root_boxes": 3, "unmatched_roots": 0, "none": 7, "passed": True,
         "complete_near_solutions": True, "pieces": 4, "sols": [], "complete_pieces": 4, "solution_pieces": 1,
         "solution_pieces_complete": 1, "examined_fraction": 1.0, "copies": 3, "kraw_undecided": 0}
    x.update(kw)
    return x


@pytest.mark.parametrize("control, why", [
    ([], "não rodado"),
    ([_summary(in_excluded=1, passed=False)], "O controle falhou"),
    ([_summary(refined=2, complete_near_solutions=False)], "O controle está incompleto"),
])
def test_no_proof_wording_without_a_passing_positive_control(control, why):
    ref = {"cert_items": {"items": 0, "items_aim": 0, "excluded": 0, "excluded_aim": 0, "certified_roots": 0, "undecided": 0, "runs": 0, "boxes": 0}}
    text = "\n".join(p3._certification_lines(EMPTY, ref, _chunks(0), control))
    assert "Provado" not in text and why in text


def test_proof_wording_comes_with_the_control_table():
    ref = {"cert_items": {"items": 0, "items_aim": 0, "excluded": 0, "excluded_aim": 0, "certified_roots": 0, "undecided": 0, "runs": 0, "boxes": 0}}
    text = "\n".join(p3._certification_lines(EMPTY, ref, _chunks(0), [_summary()]))
    assert text.index("Provado") < text.index("Controle positivo") and "O controle passou" in text


def test_control_pieces_tile_the_range_and_halve_towards_the_mouth():
    ps = p3.control_pieces(-3.0, -1e-5)
    assert ps[0][0] == -3.0 and ps[-1][1] == -1e-5 and all(a[1] == b[0] for a, b in zip(ps, ps[1:]))
    assert all(b - a > 0 for a, b in ps)


def _piece(watch, undecided, krawczyk, s=(-1.0, -0.5), unexamined=0, unexamined_volume=0.0, holds=True, excluded=None):
    return {"spec": ("r8", 0), "s": s, "boxes": 5, "excluded_by": excluded or {"miss ≠ 0": 4}, "undecided": undecided,
            "unexamined": unexamined, "volume": 1.0, "unexamined_volume": unexamined_volume, "krawczyk": krawczyk,
            "holds_solution": holds, "watch": watch}


SOL = (-0.6, 1.0, 2.0)
SOL_BOX = ((-0.6, -0.6 + 1e-9), (1.0, 1.0 + 1e-9), (2.0, 2.0 + 1e-9))


def test_control_summary_flags_a_known_solution_in_an_excluded_box():
    w = {"sol": SOL, "copy": SOL, "excluded_in": 1, "undecided_in": [0], "unexamined_in": 0}
    piece = _piece([w], [SOL_BOX], [("root", True)])
    x = p3.control_summary(("r8", 0), [SOL], [piece])
    assert not x["passed"] and x["in_excluded"] == 1
    w["excluded_in"] = 0
    x = p3.control_summary(("r8", 0), [SOL], [piece])
    assert x["passed"] and x["complete_near_solutions"] and x["certified"] == 1 and x["unmatched_roots"] == 0 and x["excluded"] == 4


def test_a_certified_control_root_names_the_geometric_model_and_is_not_physical():
    from chronosim.provenance import CERT_PHYSICAL, CERT_ROOT_MODEL
    w = {"sol": SOL, "copy": SOL, "excluded_in": 0, "undecided_in": [0], "unexamined_in": 0}
    q = p3.control_summary(("r8", 0), [SOL], [_piece([w], [SOL_BOX], [("root", True)])])["sols"][0]["quantity"]
    assert q.origin is Origin.CERTIFIED and q.label == CERT_ROOT_MODEL != CERT_PHYSICAL


def test_a_periodic_copy_of_a_known_root_is_not_an_unmatched_root():
    box_copy = ((-0.6, -0.6 + 1e-9), (1.0 + 2 * math.pi, 1.0 + 2 * math.pi + 1e-9), (2.0, 2.0 + 1e-9))
    x = p3.control_summary(("r8", 0), [SOL], [_piece([], [box_copy], [("root", True)], holds=False)])
    assert x["root_boxes"] == 1 and x["unmatched_roots"] == 0


def test_periodic_copies_inside_the_root_box():
    assert len(ct.periodic_copies((-0.3, 6.2831852, -2.6))) == 2       # nu near 2 pi: also near 0
    assert len(ct.periodic_copies((-0.3, 3.0, -3.14159265358979))) == 2  # aim near -pi: also near pi
    assert len(ct.periodic_copies((-0.3, 3.0, 1.0))) == 1


def test_a_solution_with_two_periodic_copies_is_counted_once():
    sol = (-0.26, 2 * math.pi - 1e-7, -2.6)
    c0, c1 = [tuple(float(t) for t in q) for q in ct.periodic_copies(sol)]
    boxes = [tuple((v, v + 1e-9) for v in c) for c in (c0, c1)]
    ws = [{"sol": sol, "copy": c, "excluded_in": 0, "undecided_in": [i], "unexamined_in": 0} for i, c in enumerate((c0, c1))]
    x = p3.control_summary(("fig7b", "0.04"), [sol], [_piece(ws, boxes, [("root", True), ("root", True)])])
    assert x["known"] == 1 and x["copies"] == 2 and x["certified"] == 1 and x["root_boxes"] == 2 and x["unmatched_roots"] == 0
    ws[1]["excluded_in"] = 1                                    # a copy in an excluded box fails the control
    assert not p3.control_summary(("fig7b", "0.04"), [sol], [_piece(ws, boxes, [("root", True), ("root", True)])])["passed"]


def test_the_watch_follows_a_periodic_copy_given_inside_the_root_box():
    """A point with nu > 2 pi inside the root box is watched as given; normalized, it would fall
    outside this box and be found in no box at all."""
    p = (-0.41, 2 * math.pi + 0.05, -2.54)
    r = ct.exclude(CENTRAL, (-0.42, -0.40), (6.2, 6.4), (-2.6, -2.5), max_boxes=50, watch=[p])
    h = r["watch"][0]
    assert h["excluded_in"] or h["undecided_in"]


def test_control_lines_say_when_only_the_solution_pieces_were_completed():
    x = _summary(complete_pieces=2, pieces=19, solution_pieces=1, solution_pieces_complete=1, examined_fraction=0.4)
    text = "\n".join(p3._control_lines([x])[0])
    assert "Só os pedaços com solução foram examinados até o fim" in text and "2 de 19" in text


def test_unexamined_volume_lowers_the_examined_fraction():
    w = {"sol": SOL, "copy": SOL, "excluded_in": 0, "undecided_in": [0], "unexamined_in": 0}
    a = _piece([w], [SOL_BOX], [("root", True)])
    bpc = _piece([], [], [], s=(-0.5, -0.25), unexamined=3, unexamined_volume=0.25, holds=False)
    x = p3.control_summary(("r8", 0), [SOL], [a, bpc])
    assert abs(x["examined_fraction"] - 0.875) < 1e-12 and x["complete_pieces"] == 1 and x["solution_pieces_complete"] == 1


def test_krawczyk_certifies_a_verified_solution_and_rejects_a_box_away_from_it():
    sol = [float(x) for x in _central_solution()]
    for h in (1e-6, 1e-10):
        k = ct.krawczyk(CENTRAL, [(x - h, x + h) for x in sol])
        assert k["verdict"] == "root" and k["inequalities"]
        assert all(lo <= x <= hi for (lo, hi), x in zip(k["K"], sol))
    far = [(sol[0] - 1e-6, sol[0] + 1e-6), (sol[1] + 0.5, sol[1] + 0.5 + 2e-6), (sol[2] - 1e-6, sol[2] + 1e-6)]
    assert ct.krawczyk(CENTRAL, far)["verdict"] == "none"


def test_interval_jacobian_encloses_the_finite_difference_derivatives():
    with mp.workdps(40):
        c = ct.inputs(CENTRAL, xm=mp)
        x = [mp.mpf("-0.41"), mp.mpf("0.2"), mp.mpf("-2.54")]
        h = mp.mpf(10) ** -12
        fd = []
        for j in range(3):
            xp, xn = list(x), list(x)
            xp[j] += h
            xn[j] -= h
            fp, fn = ct.conditions(*xp, c, xm=mp), ct.conditions(*xn, c, xm=mp)
            fd.append([(fp[k] - fn[k]) / (2 * h) for k in ("miss", "C", "V")])
    w = 1e-9
    with ct.ivprec(53):
        J, _ = ct.jacobian([(float(v) - w, float(v) + w) for v in x], ct.inputs(CENTRAL))
    for j in range(3):
        for i in range(3):
            assert J[i][j].a - 1e-6 <= fd[j][i] <= J[i][j].b + 1e-6, (i, j)


def test_inputs_of_a_direct_line_are_the_decimal_strings():
    c = ct.inputs({"p0": ("0.5", "0"), "e": ("0", "1"), "v1": "0.5", "B": "0.01", "R": "1e-7"}, xm=mp)
    assert c["p0"] == (mp.mpf("0.5"), 0) and c["e"] == (0, 1) and c["v1"] == mp.mpf("0.5")


def test_watch_does_not_change_the_search():
    box = ((-0.42, -0.40), (0.1, 0.3), (-2.6, -2.5))
    a = ct.exclude(CENTRAL, *box, max_boxes=2000)
    b = ct.exclude(CENTRAL, *box, max_boxes=2000, watch=[tuple(float(v) for v in _central_solution())])
    assert a["boxes"] == b["boxes"] and a["excluded_by"] == b["excluded_by"] and a["undecided"] == b["undecided"]


# reviewer findings of 2026-11-09 (docs/handoffs/2026-11-09-reviewer-para-claude.md)

def test_a_control_case_without_known_solutions_does_not_pass():
    x = p3.control_summary(("fig7b", "0.005"), [], [_piece([], [], [], holds=False)])
    assert not x["passed"] and not x["complete_near_solutions"]
    assert "O controle falhou" in "\n".join(p3._control_lines([x])[0])


def test_a_lost_copy_is_not_hidden_by_an_unexamined_copy():
    sol = (-0.26, 2 * math.pi - 1e-7, -2.6)
    c0, c1 = [tuple(float(t) for t in q) for q in ct.periodic_copies(sol)]
    ws = [{"sol": sol, "copy": c0, "excluded_in": 0, "undecided_in": [], "unexamined_in": 1},
          {"sol": sol, "copy": c1, "excluded_in": 0, "undecided_in": [], "unexamined_in": 0}]   # in no box at all
    x = p3.control_summary(("fig7b", "0.04"), [sol], [_piece(ws, [], [], unexamined=1)])
    assert not x["passed"] and x["lost"] == 1
    assert "fora de toda caixa" in "\n".join(p3._control_lines([x])[0])


def test_no_root_in_a_box_that_holds_a_known_solution_fails_the_control():
    w = {"sol": SOL, "copy": SOL, "excluded_in": 0, "undecided_in": [0], "unexamined_in": 0}
    x = p3.control_summary(("r8", 0), [SOL], [_piece([w], [SOL_BOX], [("none", True)])])
    assert not x["passed"] and x["none_in_solution_boxes"] == 1 and x["none"] == 0


def test_krawczyk_counts_of_other_boxes_leave_out_the_solution_boxes():
    w = {"sol": SOL, "copy": SOL, "excluded_in": 0, "undecided_in": [0], "unexamined_in": 0}
    far = ((-0.9, -0.9 + 1e-9), (3.0, 3.0 + 1e-9), (0.0, 1e-9))
    x = p3.control_summary(("r8", 0), [SOL], [_piece([w], [SOL_BOX, far], [("undecided", False), ("none", True)])])
    assert x["kraw_undecided"] == 0 and x["none"] == 1 and x["passed"] and x["certified"] == 0


def test_the_certified_root_holds_the_boxes_and_states_the_hypotheses():
    w = {"sol": SOL, "copy": SOL, "excluded_in": 0, "undecided_in": [0], "unexamined_in": 0}
    q = p3.control_summary(("r8", 0), [SOL], [_piece([w], [SOL_BOX], [("root", True)])])["sols"][0]["quantity"]
    assert q.value == (SOL_BOX,)
    assert "d·n_X < 0, along > 0 e kick > 0" in q.source and "modelo geométrico de billiards.py" in q.source and "é COMPUTED" in q.source


def test_the_examined_fraction_is_never_rounded_up_to_100_percent():
    for v, want in ((0.9999692, "99,997 %"), (0.99999542, "99,9995 %"), (1.0, "100 %"), (0.98305, "98,3 %")):
        text = "\n".join(p3._control_lines([_summary(examined_fraction=v)])[0])
        assert f"| {want} |" in text, (v, text)
