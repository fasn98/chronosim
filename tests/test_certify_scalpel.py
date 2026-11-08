"""Adversarial tests of the interval exclusion (certify.py), of its use in phase3.py, of
item_span / unresolved_fraction (regression) and of the CERTIFIED rules (provenance.py).

The central invariant: a box that contains a true solution is never excluded, whatever the
budget, the side of the box the solution sits on, or the way the box is split. Every solution
below was polished and verified by billiards.solve_mp / billiards.verify in the test itself.
"""

import math
import random

import mpmath as mp
import pytest
from mpmath import iv

from chronosim import billiards as bl
from chronosim import certify as ct
from chronosim import ekt, phase3 as p3
from chronosim.provenance import (CERT_EXCLUSION, CERT_LABELS, CERT_PHYSICAL, CERT_ROOT, PHYSICAL_PROOFS, Origin, Quantity,
                                  certified, combine, physical_solution)

C0, C1, C2, C3, C4, C5 = p3.R8_CASES
TWO_PI = 2 * math.pi
NU_PAD, AIM_PAD = (getattr(x, "value", x) for x in (ct.NU_PAD, ct.AIM_PAD))

# (case, s, nu, aim): starting points found by solve_mp from scattered starts; polished below.
SOLUTIONS = [
    (C0, "-0.41364648775003841355", "0.19755483169590455834", "-2.5392040925603869066"),
    (C0, "-0.020857864927748336125", "3.5689054644025959814", "-2.6983110950412729562"),
    (C0, "-0.40816375054811242984", "3.3440407923031495829", "-2.5439834012569606911"),
    (C1, "-0.37394704976276296527", "0.036956525274129939959", "-2.435215835179755983"),
    (C3, "-0.52569286182866641394", "6.2733987961475334073", "-2.9829218359701921582"),  # nu - 2 pi = -0.0098
    (C4, "-0.056618444356847510912", "4.1427929593853546328", "-2.0926884913220521837"),
    (C4, "-0.058025720828143583021", "0.99880141039656693633", "-2.0904962574990525695"),
]


def _polished(case, s, nu, aim):
    th, ph, lam, B, R = (mp.mpf(x) for x in case)
    v1 = ekt.v1_from_lambda(lam, th, ph, B, R).value
    p0, e = p3.line_theta_phi(th, ph, B)
    sol = bl.solve_mp(p0, e, mp.mpf(s), mp.mpf(nu), mp.mpf(aim), v1, 1, 1, B, R)
    ok, fails = bl.verify(p0, e, v1, *sol, 1, 1, B, R)
    assert ok, fails
    return [float(x) for x in sol]


@pytest.fixture(scope="module")
def solutions():
    return [(case, _polished(case, *start)) for case, *start in SOLUTIONS]


def _contains(box, pt, tol=1e-12):
    return all(lo - tol <= x <= hi + tol for (lo, hi), x in zip(box, pt))


def _representatives(nu, aim):
    """All (nu + 2 pi k, aim + 2 pi j) that lie in the default search ranges."""
    out = []
    for k in range(-2, 3):
        for j in range(-2, 3):
            n_, a_ = nu + k * TWO_PI, aim + j * TWO_PI
            if -NU_PAD <= n_ <= TWO_PI + NU_PAD and -AIM_PAD <= a_ <= AIM_PAD:
                out.append((n_, a_))
    return out


# --- the search ranges cover every solution -----------------------------------------------------

def test_every_known_solution_has_a_representative_in_the_default_ranges(solutions):
    for case, (s, nu, aim) in solutions:
        assert _representatives(nu, aim), (case, nu, aim)


def test_default_ranges_cover_a_full_turn_of_nu_and_of_aim():
    assert 2 * AIM_PAD >= TWO_PI and NU_PAD > 0


def test_a_solution_whose_nu_lies_in_the_padding_is_kept_under_both_representatives(solutions):
    case, (s, nu, aim) = solutions[4]
    reps = _representatives(nu, aim)
    assert any(n_ < 0 for n_, _ in reps) and any(n_ > 6 for n_, _ in reps)


# --- a box with a true solution is never excluded -----------------------------------------------

@pytest.mark.parametrize("idx", range(len(SOLUTIONS)))
@pytest.mark.parametrize("where", ["center", "s_lo_face", "s_hi_face", "nu_lo_face", "nu_hi_face", "aim_lo_face", "aim_hi_face", "corner"])
def test_a_box_containing_a_verified_solution_is_never_excluded(solutions, idx, where):
    case, sol = solutions[idx]
    h = 2e-3
    marg = {"center": (0.5,) * 3, "s_lo_face": (1e-7, .5, .5), "s_hi_face": (1 - 1e-7, .5, .5),
            "nu_lo_face": (.5, 1e-7, .5), "nu_hi_face": (.5, 1 - 1e-7, .5),
            "aim_lo_face": (.5, .5, 1e-7), "aim_hi_face": (.5, .5, 1 - 1e-7), "corner": (1e-7, 1e-7, 1 - 1e-7)}[where]
    box = [(x - 2 * h * f, x + 2 * h * (1 - f)) for x, f in zip(sol, marg)]
    r = ct.exclude(case, *box, max_boxes=1200)
    assert not r["complete"]
    assert any(_contains(b, sol) for b in r["undecided"]), "the solution's box was dropped"


@pytest.mark.parametrize("budget", [1, 2, 7, 50, 400])
def test_budget_exhaustion_keeps_the_solution_in_the_leftover(solutions, budget):
    case, sol = solutions[0]
    box = [(x - 0.05, x + 0.05) for x in sol]
    r = ct.exclude(case, *box, max_boxes=budget)
    assert any(_contains(b, sol) for b in r["undecided"])


def test_default_ranges_with_a_narrow_s_window_keep_every_solution(solutions):
    for case, sol in solutions:
        r = ct.exclude(case, (sol[0] - 1e-3, sol[0] + 1e-3), max_boxes=1500)
        assert not r["complete"]
        reps = _representatives(sol[1], sol[2])
        assert all(any(_contains(b, (sol[0], n_, a_)) for b in r["undecided"]) for n_, a_ in reps)


def test_a_tiny_box_around_a_solution_ends_undecided_not_excluded(solutions):
    case, sol = solutions[0]
    box = [(x - 1e-6, x + 1e-6) for x in sol]
    r = ct.exclude(case, *box, max_boxes=4000)
    assert not r["complete"] and r["undecided"]
    assert any(_contains(b, sol) for b in r["undecided"])


# --- splitting never drops space ------------------------------------------------------------------

class _Width:
    def __init__(self, value):
        self.value = value


def test_splitting_tiles_the_root_box(monkeypatch):
    monkeypatch.setattr(ct, "_reason", lambda box, c: None)   # exclude nothing: every leaf is undecided
    monkeypatch.setattr(ct, "CERT_MIN_WIDTH", _Width(0.2))
    root = ((-1.0, -0.2), (0.0, 3.0), (-2.0, 1.0))
    r = ct.exclude(C2, *root, max_boxes=10**6)
    vol = lambda b: math.prod(hi - lo for lo, hi in b)
    assert math.isclose(sum(vol(b) for b in r["undecided"]), vol(root), rel_tol=1e-12)
    rng = random.Random(11)
    for _ in range(300):
        pt = [rng.uniform(*d) for d in root]
        assert any(_contains(b, pt, 0) for b in r["undecided"])


@pytest.mark.parametrize("pt", [(-1.0, 0.0, -2.0), (-0.2, 3.0, 1.0), (-0.6, 1.5, -0.5), (-1.0, 3.0, -2.0)])
def test_a_point_on_a_corner_or_on_a_split_plane_survives_a_pointwise_exact_test(monkeypatch, pt):
    """Replace the interval test by the exact test 'the box does not contain pt': pt must remain."""
    monkeypatch.setattr(ct, "_reason", lambda box, c: None if _contains(box, pt, 0) else "far")
    monkeypatch.setattr(ct, "CERT_MIN_WIDTH", _Width(1e-3))
    r = ct.exclude(C2, (-1.0, -0.2), (0.0, 3.0), (-2.0, 1.0), max_boxes=10**6)
    assert r["undecided"] and all(_contains(b, pt, 0) for b in r["undecided"])
    assert not r["complete"]


def test_complete_means_no_box_left_even_when_the_budget_ends_on_the_last_box():
    r = ct.exclude(C2, (-3.0, -0.00135), (0.0, 1.0), (-0.2, 0.2), max_boxes=1)
    assert r["complete"] and r["boxes"] == 1 and not r["undecided"]
    assert sum(r["excluded_by"].values()) == 1


def test_the_box_count_adds_up():
    r = ct.exclude(C2, (-3.0, -0.00135), max_boxes=300)
    assert r["boxes"] == 300 and sum(r["excluded_by"].values()) <= r["boxes"]
    assert len(r["undecided"]) >= 1


# --- the enclosures are sound --------------------------------------------------------------------

@pytest.mark.parametrize("case", [C0, C2, C1, C4])
def test_interval_conditions_enclose_point_values_on_random_boxes(case):
    rng = random.Random(7)
    with mp.workdps(40):
        cm = ct.case_inputs(*case, xm=mp)
    for _ in range(120):
        w = 10 ** rng.uniform(-7, 0.4)
        s0, n0, a0 = rng.uniform(-3, -0.0005), rng.uniform(-0.01, 6.3), rng.uniform(-3.15, 3.15)
        box = [(s0, s0 + w), (n0, n0 + w), (a0, a0 + w)]
        with ct.ivprec(53):
            q = ct.conditions(*(iv.mpf(b) for b in box), ct.case_inputs(*case))
        for _ in range(2):
            pt = [rng.uniform(*b) for b in box]
            with mp.workdps(40):
                qm = ct.conditions(*(mp.mpf(x) for x in pt), cm, xm=mp)
            for k, v in qm.items():
                assert q[k].a <= v <= q[k].b, (k, box, pt)


def test_interval_inputs_enclose_the_decimal_strings():
    """B, R are enclosures of the exact decimals, not of their float roundings; p0, e, v1 enclose
    the values computed in mpmath from the strings (independently of case_inputs)."""
    with mp.workdps(60):
        th, ph, lam, B, R = (mp.mpf(x) for x in C2)
        p0, e = p3.line_theta_phi(th, ph, B)
        v1 = ekt.v1_from_lambda(lam, th, ph, B, R).value
    with ct.ivprec(53):
        ci = ct.case_inputs(*C2)
    assert ci["B"].a <= B <= ci["B"].b and ci["R"].a <= R <= ci["R"].b
    assert ci["v1"].a <= v1 <= ci["v1"].b
    for k, ref in (("p0", p0), ("e", e)):
        for a, b in zip(ci[k], ref):
            assert a.a <= b <= a.b, k


def test_box_where_the_ray_origin_touches_the_entry_point_does_not_raise():
    """d = X - P contains 0: L0 = sqrt(d.d) and the divisions by L0 must degrade to wide
    intervals, not to an exception or a NaN that compares as 'excluded'."""
    a = 0.6 - math.pi    # aim of p0 on mouth 2 for C0 (s = 0 puts P at p0)
    r = ct.exclude(C0, (-1e-3, 0.0), (0.0, 6.3), (a - 0.01, a + 0.01), max_boxes=400)
    assert r["boxes"] > 0


def test_a_square_over_an_interval_straddling_zero_is_not_negative():
    with ct.ivprec(53):
        c = ct.case_inputs(*C0)
        P_box = iv.mpf([-1e-3, 0.0])
        q = ct.conditions(P_box, iv.mpf([0.0, 6.3]), iv.mpf([0.6 - math.pi - 0.01, 0.6 - math.pi + 0.01]), c)
    for k, v in q.items():
        assert not (math.isnan(float(v.a)) or math.isnan(float(v.b))), k


def test_no_undecided_box_has_a_nan_end():
    r = ct.exclude(C0, (-1e-3, 0.0), (0.0, 6.3), (0.6 - math.pi - 0.01, 0.6 - math.pi + 0.01), max_boxes=60)
    for b in r["undecided"]:
        assert all(not math.isnan(x) for side in b for x in side)


# --- precision handling ---------------------------------------------------------------------------

def test_ivprec_restores_after_an_exception():
    before = iv.prec
    with pytest.raises(RuntimeError):
        with ct.ivprec(200):
            assert iv.prec == 200
            raise RuntimeError
    assert iv.prec == before


def test_exclude_leaves_the_global_interval_precision_alone():
    saved = iv.prec
    try:
        iv.prec = 77
        ct.exclude(C2, (-3.0, -0.00135), max_boxes=5)
        assert iv.prec == 77
    finally:
        iv.prec = saved


# --- outward rounding of s and the quantity ---------------------------------------------------------

@pytest.mark.parametrize("lo,hi", [("-3", "-0.00135"), ("-0.1", "-0.1"), ("0.3", "0.30000000000000004"), ("-1e-300", "1e-300")])
def test_outward_s_range_contains_the_exact_decimals_strictly(lo, hi):
    a, b = ct.s_range_outward(lo, hi)
    with mp.workdps(60):
        assert mp.mpf(a) < mp.mpf(lo) and mp.mpf(b) > mp.mpf(hi)


def test_outward_s_range_accepts_mpf_ends():
    a, b = ct.s_range_outward(mp.mpf("-3"), mp.mpf("-0.00135"))
    assert a < -3 and b > -0.00135


def test_incomplete_result_is_never_certified_even_with_boxes_counted():
    r = {"boxes": 10, "undecided": [((0, 1),) * 3], "complete": False}
    q = ct.exclusion_quantity(r, C2, (-3.0, -0.1), "x")
    assert q.origin is Origin.COMPUTED and q.label == "" and q.value == 1


def test_complete_quantity_names_the_case_the_range_and_the_hypotheses():
    r = {"boxes": 10, "undecided": [], "complete": True}
    q = ct.exclusion_quantity(r, C2, (-3.0, -0.1), "x")
    assert q.origin is Origin.CERTIFIED and q.label == CERT_EXCLUSION and q.value == 0 and "10 subcaixas" in q.source
    for needle in (str(C2), "-3.0", "-0.1", "53", "uma colisão", "1 travessia", "arredondamento para fora"):
        assert needle in q.source


# --- phase3 glue: chunks and items ------------------------------------------------------------------

@pytest.mark.parametrize("n", [1, 3, 8])
def test_certify_chunks_tile_the_s_range_without_gaps_or_overlap(monkeypatch, n):
    seen = []

    def fake(case, s_range, *a, **k):
        seen.append(tuple(s_range))
        return {"boxes": 1, "excluded_by": {}, "undecided": [], "complete": True}

    monkeypatch.setattr(ct, "exclude", fake)
    chunks = [p3.certify_chunk(C2, i, n) for i in range(n)]
    lo, hi = chunks[0]["range"]
    assert seen[0][0] == lo and seen[-1][1] == hi
    assert all(seen[i][1] == seen[i + 1][0] for i in range(n - 1))
    assert all(a < b for a, b in seen)
    mp_lo = mp.mpf(lo)
    assert mp_lo < -mp.mpf(p3.S_FAR.value) and hi > -mp.mpf(C2[4]) * mp.mpf(p3.S_NEAR_OVER_R.value)


def test_certify_items_boxes_cover_each_items_span_and_share_equal_spans(monkeypatch):
    calls = []

    def fake(case, s_range, *a, **k):
        calls.append(tuple(s_range))
        return {"boxes": 3, "excluded_by": {}, "undecided": [], "complete": True}

    monkeypatch.setattr(ct, "exclude", fake)
    sg = [-3.0, -2.0, -1.0, -0.5, -0.1]
    p0, e = (0.0, 0.0), (1.0, 0.0)
    items = [("s_count_change", -2.5, -1.5), ("s_bisect_lost", -1.5, -2.5), ("aim_bracket_lost", (-1.0, 0.0), 0.1, 0.2, 0.3),
             ("nu_count_change", (-1.0, 0.0), 0.1, 0.2), ("bogus", 1, 2)]
    out = p3.certify_items(C2, items, sg, p0, e)
    assert out["items"] == 5 and out["excluded"] == 4 and out["undecided"] == 1
    assert out["runs"] == len(calls) == len(set(calls)) == 2   # two S items share a span; the aim and nu items at s = -1 share another
    for u, (kind, span, verdict) in zip(items, out["verdicts"]):
        if span is None:
            assert verdict == "sem intervalo em s"
            continue
        piece, _, _ = p3.item_span(u, sg, p0, e)
        assert span[0] < piece[0] and span[1] > piece[1]


def test_certify_items_reports_an_incomplete_run_as_undecided(monkeypatch):
    monkeypatch.setattr(ct, "exclude", lambda *a, **k: {"boxes": 3, "excluded_by": {}, "undecided": [1], "complete": False})
    out = p3.certify_items(C2, [("s_count_change", -2.0, -1.0)], [-3.0, -2.0, -1.0], (0.0, 0.0), (1.0, 0.0))
    assert out["excluded"] == 0 and out["undecided"] == 1


# --- item_span and unresolved_fraction: regression --------------------------------------------------

def _old_unresolved_fraction(unres, s_grid, p0, e):
    """The function of main (0ed21ec), restricted to the kinds it knew (S and NU kinds)."""
    sg = sorted(float(x) for x in s_grid)
    span = sg[-1] - sg[0]
    mids = [(a + b) / 2 for a, b in zip(sg, sg[1:])]
    cell_of = {x: ((mids[i - 1] if i else x), (mids[i] if i < len(mids) else x)) for i, x in enumerate(sg)}
    span_of = {x: (sg[max(i - 1, 0)], sg[min(i + 1, len(sg) - 1)]) for i, x in enumerate(sg)}
    two_pi = 2 * math.pi
    area, off, pieces, s_pieces = 0.0, 0, [], []
    for u in unres:
        if u[0] in p3.S_KINDS:
            a, b = sorted((float(u[1]), float(u[2])))
            s_pieces.append((a, b))
            pieces.append((a, b))
        elif u[0] in p3.NU_KINDS:
            s_ = float(bl.dot(bl.add(u[1], p0, -1), e))
            near = min(sg, key=lambda x: abs(x - s_))
            if abs(near - s_) <= 1e-12 * (1 + abs(s_)):
                lo, hi = cell_of[near]
                area += (hi - lo) * abs(float(u[3]) - float(u[2]))
                pieces.append(span_of[near])
            else:
                off += 1
                j = max(i for i, x in enumerate(sg) if x <= s_) if s_ >= sg[0] else 0
                pieces.append((sg[j], sg[min(j + 1, len(sg) - 1)]))

    def union(ps):
        covered, end = 0.0, -math.inf
        for a, b in sorted(ps):
            if b > end:
                covered += b - max(a, end)
                end = b
        return covered

    fine = union(s_pieces) / span + area / (span * two_pi)
    return {"fine": min(fine, union(pieces) / span), "off_grid": off, "conservative": union(pieces) / span}


@pytest.mark.parametrize("seed", range(8))
def test_unresolved_fraction_is_unchanged_for_the_kinds_main_knew(seed):
    rng = random.Random(seed)
    sg = sorted(rng.uniform(-3, -0.001) for _ in range(rng.choice([2, 3, 6, 15])))
    p0, e = (0.3, -0.2), (math.cos(0.4), math.sin(0.4))
    items = []
    for _ in range(rng.randint(0, 12)):
        kind = rng.choice(p3.S_KINDS + p3.NU_KINDS)
        if kind in p3.S_KINDS:
            items.append((kind, rng.uniform(-3, -0.001), rng.uniform(-3, -0.001)))
        else:
            s_ = rng.choice(sg) if rng.random() < .6 else rng.uniform(sg[0], sg[-1])
            P = (p0[0] + s_ * e[0], p0[1] + s_ * e[1])
            a = rng.uniform(0, 6)
            items.append((kind, P, a, a + rng.uniform(0, 0.3)))
    assert p3.unresolved_fraction(items, sg, p0, e) == _old_unresolved_fraction(items, sg, p0, e)


def test_unresolved_fraction_of_the_empty_list_is_zero():
    assert p3.unresolved_fraction([], [-3.0, -1.0], (0.0, 0.0), (1.0, 0.0)) == {"fine": 0.0, "off_grid": 0, "conservative": 0.0}


def test_item_span_of_an_unknown_kind_is_none_and_costs_nothing():
    assert p3.item_span(("whatever", 1), [-2.0, -1.0], (0.0, 0.0), (1.0, 0.0)) == (None, 0.0, False)
    assert p3.unresolved_fraction([("whatever", 1)], [-2.0, -1.0], (0.0, 0.0), (1.0, 0.0))["conservative"] == 0.0


def test_item_span_of_s_items_is_sorted():
    assert p3.item_span(("s_count_change", -1.0, -2.0), [-3.0, 0.0], (0, 0), (1, 0)) == ((-2.0, -1.0), 0.0, True)


@pytest.mark.parametrize("s_,piece", [(-3.0, (-3.0, -2.0)), (-1.0, (-2.0, -0.5)), (-0.1, (-0.5, -0.1))])
def test_item_span_on_grid_covers_both_neighbouring_intervals_and_the_ends(s_, piece):
    sg = [-3.0, -2.0, -1.0, -0.5, -0.1]
    got, area, on_grid = p3.item_span(("aim_bracket_lost", (s_, 0.0), 0.0, 0.1, 0.2), sg, (0.0, 0.0), (1.0, 0.0))
    assert got == (min(piece), max(piece)) or got == piece
    assert on_grid and area == 0.0


def test_item_span_off_grid_is_the_enclosing_grid_interval_and_contains_s():
    sg = [-3.0, -2.0, -1.0]
    for s_ in (-2.5, -2.0000001, -1.0000001):
        got, _, on_grid = p3.item_span(("nu_count_change", (s_, 0.0), 0.1, 0.2), sg, (0.0, 0.0), (1.0, 0.0))
        assert not on_grid and got[0] <= s_ <= got[1] and got[0] < got[1]


def test_item_span_with_a_one_point_grid_does_not_raise():
    got, _, _ = p3.item_span(("aim_bracket_lost", (-1.0, 0.0), 0.0, 0.1, 0.2), [-1.0], (0.0, 0.0), (1.0, 0.0))
    assert got == (-1.0, -1.0)


def test_nu_area_is_only_charged_to_nu_kinds():
    sg = [-3.0, -2.0, -1.0]
    _, a_aim, _ = p3.item_span(("aim_bracket_jump", (-2.0, 0.0), 0.1, 0.2, 0.3), sg, (0.0, 0.0), (1.0, 0.0))
    _, a_nu, _ = p3.item_span(("nu_count_change", (-2.0, 0.0), 0.1, 0.2), sg, (0.0, 0.0), (1.0, 0.0))
    assert a_aim == 0.0 and a_nu == pytest.approx(1.0 * 0.1)


def test_aim_item_counts_in_the_conservative_measure_exactly_like_a_nu_item_at_the_same_point():
    sg = [-3.0, -2.0, -1.0, 0.0]
    aim = p3.unresolved_fraction([("aim_bracket_lost", (-2.0, 0.0), 0.1, 0.2, 0.3)], sg, (0.0, 0.0), (1.0, 0.0))
    nu = p3.unresolved_fraction([("nu_count_change", (-2.0, 0.0), 0.1, 0.2)], sg, (0.0, 0.0), (1.0, 0.0))
    assert aim["conservative"] == nu["conservative"] == pytest.approx(2 / 3)


# --- N-a: the bracketing log ------------------------------------------------------------------------

def test_logging_does_not_change_the_roots_found():
    f = lambda x: math.sin(5 * x) + 0.1 * x
    grid = [i / 20 for i in range(-40, 41)]
    assert bl.roots_by_bracketing(f, grid, 1e-12) == bl.roots_by_bracketing(f, grid, 1e-12, [], "w")


def test_a_log_is_optional_and_a_discard_without_a_log_is_silent():
    assert bl.roots_by_bracketing(lambda x: 1.0 if x > 0.3 else -1.0, [-1.0, 1.0], 1e-10) == []


def test_a_clean_function_logs_nothing_and_finds_every_sign_change():
    log = []
    got = bl.roots_by_bracketing(lambda x: x * (x - 1) * (x + 1), [-2.0, -0.5, 0.5, 2.0], 1e-12, log)
    assert len(got) == 3 and log == []


# --- provenance: the CERTIFIED rules -----------------------------------------------------------------

def test_physical_label_cannot_be_forged_by_building_the_quantity_directly():
    """certified() refuses CERT_PHYSICAL without both interval proofs; the constructor must not be a
    way around that rule (Quantity(..., CERTIFIED, label=CERT_PHYSICAL) with no proofs)."""
    with pytest.raises(ValueError):
        Quantity(1, Origin.CERTIFIED, "no proof at all", label=CERT_PHYSICAL)


def test_proofs_given_as_one_string_do_not_unlock_the_physical_label():
    one = "; ".join(PHYSICAL_PROOFS)
    with pytest.raises(ValueError):
        certified(1, CERT_PHYSICAL, "x", proofs=one)
    q = physical_solution(1, "x", proofs=one)
    assert q.origin is Origin.COMPUTED


def test_proofs_may_be_any_iterable_and_extra_proofs_are_allowed():
    q = certified(1, CERT_PHYSICAL, "x", proofs=(p for p in ("Krawczyk",) + PHYSICAL_PROOFS))
    assert q.label == CERT_PHYSICAL
    q = physical_solution(1, "x", proofs=iter(PHYSICAL_PROOFS))
    assert q.origin is Origin.CERTIFIED and q.label == CERT_PHYSICAL


@pytest.mark.parametrize("variant", [PHYSICAL_PROOFS[0] + " ", PHYSICAL_PROOFS[1].upper(), PHYSICAL_PROOFS[0][:-1]])
def test_a_near_miss_of_a_proof_name_does_not_count(variant):
    q = physical_solution(1, "x", proofs=(variant, PHYSICAL_PROOFS[1]) if variant != PHYSICAL_PROOFS[1].upper()
                          else (PHYSICAL_PROOFS[0], variant))
    assert q.origin is Origin.COMPUTED and "não provado" in q.source


def test_physical_solution_with_no_proof_says_nothing_was_proven():
    q = physical_solution(1, "x", proofs=())
    assert q.origin is Origin.COMPUTED and "provado: nada" in q.source and all(p in q.source for p in PHYSICAL_PROOFS)


def test_physical_solution_keeps_unit_and_value():
    q = physical_solution(2.5, "x", PHYSICAL_PROOFS, unit="m")
    assert q.value == 2.5 and q.unit == "m"
    q = physical_solution(2.5, "x", (), unit="m")
    assert q.value == 2.5 and q.unit == "m"


def test_degrade_certified_to_certified_keeps_the_label_and_changes_the_source():
    q = certified(1, CERT_ROOT, "a").degrade(Origin.CERTIFIED, "b")
    assert q.origin is Origin.CERTIFIED and q.label == CERT_ROOT and q.source == "b"


@pytest.mark.parametrize("to", [Origin.COMPUTED, Origin.HEURISTIC])
@pytest.mark.parametrize("label", CERT_LABELS)
def test_degrading_any_certified_label_drops_it(to, label):
    proofs = PHYSICAL_PROOFS if label == CERT_PHYSICAL else ()
    q = certified(1, label, "a", proofs=proofs).degrade(to, "b")
    assert q.origin is to and q.label == ""


@pytest.mark.parametrize("to", [Origin.CITED, Origin.DERIVED])
def test_certified_cannot_be_promoted_to_cited_or_derived(to):
    with pytest.raises(ValueError):
        certified(1, CERT_ROOT, "a").degrade(to, "b")


@pytest.mark.parametrize("frm", [Origin.CITED, Origin.DERIVED, Origin.COMPUTED, Origin.HEURISTIC])
def test_no_origin_degrades_or_promotes_into_certified(frm):
    with pytest.raises(ValueError):
        Quantity(1, frm, "a").degrade(Origin.CERTIFIED, "b")


def test_combine_never_returns_certified_for_any_mix():
    certs = [certified(1, CERT_ROOT, "a"), certified(1, CERT_EXCLUSION, "b"), certified(1, CERT_PHYSICAL, "c", PHYSICAL_PROOFS)]
    others = [Quantity(1, o, "z") for o in (Origin.CITED, Origin.DERIVED, Origin.COMPUTED, Origin.HEURISTIC)]
    for n in range(1, 4):
        for c in certs:
            for o in others + [None]:
                got = combine(1, [c] * n + ([o] if o else []), "s")
                assert got.origin is not Origin.CERTIFIED and got.label == ""
                assert got.origin >= Origin.COMPUTED


def test_combine_with_cited_and_derived_only_stays_below_certified():
    d, c = Quantity(1, Origin.DERIVED, "a"), Quantity(1, Origin.CITED, "b")
    assert combine(1, [c, d], "s").origin is Origin.DERIVED
    assert combine(1, [c], "s").origin is Origin.CITED


def test_combine_with_no_inputs_and_certified_floor_is_refused():
    assert combine(1, [], "s").origin is Origin.CITED
    with pytest.raises(ValueError):
        combine(1, [], "s", at_least=Origin.CERTIFIED)
    with pytest.raises(ValueError):
        combine(1, [certified(1, CERT_ROOT, "a")], "s", at_least=Origin.CERTIFIED)


def test_combine_of_certified_with_at_least_computed_or_heuristic():
    c = certified(1, CERT_ROOT, "a")
    assert combine(1, [c], "s", at_least=Origin.COMPUTED).origin is Origin.COMPUTED
    assert combine(1, [c], "s", at_least=Origin.HEURISTIC).origin is Origin.HEURISTIC


def test_combine_keeps_unit_and_does_not_carry_a_label():
    q = combine(1, [certified(1, CERT_ROOT, "a")], "s", unit="m")
    assert q.unit == "m" and q.label == ""


def test_a_label_on_a_non_certified_origin_is_refused_for_every_origin():
    for o in (Origin.CITED, Origin.DERIVED, Origin.COMPUTED, Origin.HEURISTIC):
        with pytest.raises(ValueError):
            Quantity(1, o, "x", label=CERT_ROOT)


def test_an_unknown_label_is_refused_on_construction_and_on_degrade_to_certified():
    with pytest.raises(ValueError):
        Quantity(1, Origin.CERTIFIED, "x", label="provado")
    with pytest.raises(ValueError):
        certified(1, "", "x")


def test_certified_needs_a_source_and_the_physical_error_names_the_missing_proofs():
    with pytest.raises(ValueError):
        certified(1, CERT_ROOT, "   ")
    with pytest.raises(ValueError) as e:
        certified(1, CERT_PHYSICAL, "x", proofs=(PHYSICAL_PROOFS[0],))
    assert PHYSICAL_PROOFS[1] in str(e.value)


def test_every_origin_name_used_by_the_csvs_exists():
    assert {o.name for o in Origin} == {"CITED", "DERIVED", "CERTIFIED", "COMPUTED", "HEURISTIC"}
