"""Edge cases of the positive-control round of 2026-11-09 (scalpel): exclude(scale=), periodic_copies
at the borders, control_summary with an unexamined copy, and the refusal of bare root labels."""

import math

import pytest

from chronosim import certify as ct
from chronosim import phase3 as p3
from chronosim.provenance import CERT_LABELS, CERT_ROOT, CERT_ROOT_DIR, CERT_ROOT_MODEL, CERT_ROOT_PREFIX, Origin, Quantity, certified

CENTRAL = p3.R8_CASES[0]
EMPTY = p3.R8_CASES[2]
SMALL = ((-0.42, -0.40), (6.2, 6.4), (-2.6, -2.5))   # s, nu, aim: a small box of the central case


# exclude(scale=) ------------------------------------------------------------------------------

@pytest.mark.parametrize("scale", [(0, 1, 1), (1, 0.0, 1), (1, 1, -1.0), (0, 0, 0), (-1, -1, -1)])
def test_exclude_refuses_a_scale_with_a_nonpositive_width(scale):
    with pytest.raises(ValueError):
        ct.exclude(CENTRAL, SMALL[0], SMALL[1], SMALL[2], max_boxes=5, scale=scale)


def test_exclude_refuses_a_scale_that_does_not_have_three_widths():
    with pytest.raises(ValueError):
        ct.exclude(CENTRAL, SMALL[0], SMALL[1], SMALL[2], max_boxes=5, scale=(1.0, 1.0))


def test_exclude_refuses_a_scale_with_an_extra_width():
    with pytest.raises(ValueError):
        ct.exclude(CENTRAL, SMALL[0], SMALL[1], SMALL[2], max_boxes=5, scale=(1.0, 1.0, 1.0, 1.0))


def _flat(cps):
    return [t for q in cps for t in q]


def _essence(r):
    return (r["boxes"], r["excluded_by"], r["undecided"], r["complete"], r["unexamined"])


def test_a_scale_equal_to_the_box_sides_changes_nothing():
    sides = tuple(hi - lo for lo, hi in SMALL)
    a = ct.exclude(CENTRAL, *SMALL, max_boxes=300)
    b = ct.exclude(CENTRAL, *SMALL, max_boxes=300, scale=sides)
    assert _essence(a) == _essence(b)


def test_a_scale_equal_to_the_default_root_sides_changes_nothing_on_the_default_box():
    nu, aim = ct.default_ranges()
    sides = (-0.1 - -0.3, nu[1] - nu[0], aim[1] - aim[0])
    a = ct.exclude(EMPTY, (-0.3, -0.1), max_boxes=300)
    b = ct.exclude(EMPTY, (-0.3, -0.1), max_boxes=300, scale=sides)
    assert _essence(a) == _essence(b)


def test_a_larger_scale_stops_the_subdivision_earlier():
    """The stop rule is relative to the scale: a scale 1000x larger makes every box 'small enough' sooner."""
    sol = (-0.41364648775003827, 0.19755483169590468, -2.539204092560386)
    box = tuple((x - 1e-3, x + 1e-3) for x in sol)
    sides = tuple(hi - lo for lo, hi in box)
    fine = ct.exclude(CENTRAL, *box, max_boxes=2000)
    coarse = ct.exclude(CENTRAL, *box, max_boxes=2000, scale=tuple(1e10 * w for w in sides))
    assert fine["unexamined"] > 0           # at its own scale the budget runs out first
    assert coarse["unexamined"] == 0 and coarse["boxes"] < fine["boxes"]


# periodic_copies at the borders ---------------------------------------------------------------

def test_nu_exactly_zero_also_has_its_copy_at_two_pi():
    cps = ct.periodic_copies((-0.3, 0.0, 1.0))
    assert sorted(float(q[1]) for q in cps) == pytest.approx([0.0, 2 * math.pi])
    assert all(float(q[2]) == 1.0 for q in cps)


def test_aim_exactly_minus_pi_also_has_its_copy_at_plus_pi():
    cps = ct.periodic_copies((-0.3, 3.0, -math.pi))
    assert sorted(float(q[2]) for q in cps) == pytest.approx([-math.pi, math.pi])


def test_aim_exactly_plus_pi_is_the_same_solution_as_minus_pi():
    a = sorted(tuple(float(t) for t in q) for q in ct.periodic_copies((-0.3, 3.0, math.pi)))
    b = sorted(tuple(float(t) for t in q) for q in ct.periodic_copies((-0.3, 3.0, -math.pi)))
    assert _flat(a) == pytest.approx(_flat(b), abs=1e-12)


def test_nu_exactly_two_pi_is_the_same_solution_as_nu_zero():
    a = sorted(tuple(float(t) for t in q) for q in ct.periodic_copies((-0.3, 2 * math.pi, 1.0)))
    b = sorted(tuple(float(t) for t in q) for q in ct.periodic_copies((-0.3, 0.0, 1.0)))
    assert _flat(a) == pytest.approx(_flat(b), abs=1e-12)


def test_nu_zero_and_aim_minus_pi_give_four_copies_and_every_one_is_in_the_root_box():
    cps = ct.periodic_copies((-0.3, 0.0, -math.pi))
    (n0, n1), (a0, a1) = ct.default_ranges()
    assert len(cps) == 4
    assert all(n0 <= q[1] <= n1 and a0 <= q[2] <= a1 for q in cps)
    assert len({(float(q[1]), float(q[2])) for q in cps}) == 4


def test_s_is_never_changed_by_the_copies():
    assert {float(q[0]) for q in ct.periodic_copies((-0.3, 0.0, -math.pi))} == {-0.3}


def test_a_point_far_outside_the_circle_is_normalized_before_copying():
    cps = ct.periodic_copies((-0.3, 100.0, 50.0))
    (n0, n1), (a0, a1) = ct.default_ranges()
    assert cps and all(n0 <= q[1] <= n1 and a0 <= q[2] <= a1 for q in cps)


def test_a_custom_range_without_the_copy_drops_it():
    cps = ct.periodic_copies((-0.3, 0.0, 1.0), nu_range=(-0.01, 1.0))
    assert [float(q[1]) for q in cps] == [0.0]


# control_summary with an unexamined copy ------------------------------------------------------

SOL = (-0.6, 1.0, 2.0)
BOX = ((-0.6, -0.6 + 1e-9), (1.0, 1.0 + 1e-9), (2.0, 2.0 + 1e-9))
WIDE = ((-0.6, -0.6 + 1e-3), (1.0, 1.0 + 1e-3), (2.0, 2.0 + 1e-3))


def _piece(watch, undecided, krawczyk, s=(-1.0, -0.5), unexamined=0, holds=True):
    return {"spec": ("r8", 0), "s": s, "boxes": 5, "excluded_by": {"miss ≠ 0": 4}, "undecided": undecided,
            "unexamined": unexamined, "volume": 1.0, "unexamined_volume": 0.0, "krawczyk": krawczyk,
            "holds_solution": holds, "watch": watch}


def _w(sol=SOL, excluded=0, undecided=(0,), unexamined=0):
    return {"sol": sol, "copy": sol, "excluded_in": excluded, "undecided_in": list(undecided), "unexamined_in": unexamined}


def _sol(x):
    return x["sols"][0]


def test_a_copy_in_an_unexamined_box_is_not_certified_and_has_no_quantity():
    x = p3.control_summary(("r8", 0), [SOL], [_piece([_w(undecided=(), unexamined=1)], [], [], unexamined=1)])
    assert _sol(x)["certified"] is False and _sol(x)["quantity"] is None and x["certified"] == 0


def test_a_copy_in_an_unexamined_box_is_reported_as_unexamined_not_as_a_failure():
    x = p3.control_summary(("r8", 0), [SOL], [_piece([_w(undecided=(), unexamined=1)], [], [], unexamined=1)])
    assert _sol(x)["unexamined"] is True and x["in_excluded"] == 0 and x["passed"]


def test_a_copy_in_an_unexamined_box_keeps_the_proof_wording_out():
    x = p3.control_summary(("r8", 0), [SOL], [_piece([_w(undecided=(), unexamined=1)], [], [], unexamined=1)])
    assert x["complete_near_solutions"] is False and x["refined"] == 0


def test_one_certified_copy_and_one_unexamined_copy_do_not_certify_the_solution():
    ws = [_w(undecided=(0,)), _w(undecided=(), unexamined=1)]
    x = p3.control_summary(("r8", 0), [SOL], [_piece(ws, [BOX], [("root", True)], unexamined=1)])
    assert _sol(x)["copies"] == 2 and _sol(x)["certified"] is False and _sol(x)["quantity"] is None


def test_a_copy_that_is_in_no_box_at_all_fails_the_control():
    x = p3.control_summary(("r8", 0), [SOL], [_piece([_w(undecided=())], [], [])])
    assert not x["passed"]


def test_a_solution_with_no_watch_entry_fails_the_control():
    x = p3.control_summary(("r8", 0), [SOL], [_piece([], [], [], holds=False)])
    assert not x["passed"] and _sol(x)["copies"] == 0 and _sol(x)["quantity"] is None


@pytest.mark.parametrize("verdict", [("root", False), ("undecided", False), ("none", False), ("undecided", True)])
def test_a_refined_copy_without_a_proven_root_is_not_certified(verdict):
    x = p3.control_summary(("r8", 0), [SOL], [_piece([_w()], [BOX], [verdict])])
    assert _sol(x)["certified"] is False and _sol(x)["quantity"] is None


def test_a_root_box_is_certified_by_its_krawczyk_verdict_whatever_its_width():
    """Changed after the reviewer's note of 2026-11-09: the unregistered 1e-6 threshold is gone. An
    examined undecided box already stops at ekt_cert_min_width, and K(X) inside X proves a unique root
    at any width; the certified value is the box itself."""
    x = p3.control_summary(("r8", 0), [SOL], [_piece([_w()], [WIDE], [("root", True)])])
    assert _sol(x)["certified"] is True and _sol(x)["quantity"].value == (WIDE,)


def test_the_certified_quantity_carries_the_model_label_and_the_solution_key():
    x = p3.control_summary(("r8", 0), [SOL], [_piece([_w()], [BOX], [("root", True)])])
    q = _sol(x)["quantity"]
    # the value is the box where Krawczyk proved the root; the bracketing point is named as COMPUTED
    assert q.origin is Origin.CERTIFIED and q.label == CERT_ROOT_MODEL and q.value == (BOX,)
    assert f"({_sol(x)['sol']}) é COMPUTED" in q.source


def test_the_solution_is_looked_up_normalized_so_a_shifted_input_is_still_found():
    shifted = (-0.6, 1.0 + 2 * math.pi, 2.0)
    x = p3.control_summary(("r8", 0), [shifted], [_piece([_w(sol=shifted)], [BOX], [("root", True)])])
    assert _sol(x)["copies"] == 1 and _sol(x)["certified"] is True


# refusal of bare "raiz certificada" labels ----------------------------------------------------

BARE = ["raiz certificada", CERT_ROOT_PREFIX, CERT_ROOT_PREFIX + " ", " " + CERT_ROOT_PREFIX, "raiz certificada do modelo geométrico",
        "raiz certificada do modelo geométrico (billiards.py) ", "RAIZ CERTIFICADA", "raiz certificada das eqs. 4.2-4.3",
        "raiz certificada das eqs. 4.2–4.3.", "raiz  certificada das eqs. 4.2–4.3"]


@pytest.mark.parametrize("label", BARE)
def test_every_unnamed_or_misspelled_root_label_is_refused(label):
    assert label not in CERT_LABELS
    with pytest.raises(ValueError):
        certified(1, label, "Krawczyk")
    with pytest.raises(ValueError):
        Quantity(1, Origin.CERTIFIED, "Krawczyk", label=label)


@pytest.mark.parametrize("label", [CERT_ROOT, CERT_ROOT_MODEL, CERT_ROOT_DIR])
def test_a_root_label_on_a_non_certified_origin_is_refused(label):
    for origin in (Origin.COMPUTED, Origin.HEURISTIC, Origin.DERIVED, Origin.CITED):
        with pytest.raises(ValueError):
            Quantity(1, origin, "s", label=label)


@pytest.mark.parametrize("label", BARE)
def test_a_degraded_quantity_cannot_be_given_a_bare_label_back(label):
    q = certified(1, CERT_ROOT_MODEL, "Krawczyk")
    with pytest.raises(ValueError):
        Quantity(q.value, Origin.CERTIFIED, "x", label=label)


def test_exactly_the_root_labels_start_with_the_prefix_and_each_names_a_system():
    roots = [lab for lab in CERT_LABELS if lab.startswith(CERT_ROOT_PREFIX)]
    assert sorted(roots) == sorted([CERT_ROOT, CERT_ROOT_MODEL, CERT_ROOT_DIR])
    assert all(lab != CERT_ROOT_PREFIX and lab[len(CERT_ROOT_PREFIX):].startswith(" d") for lab in roots)


def test_a_degraded_certified_root_drops_its_label():
    d = certified(1, CERT_ROOT_MODEL, "Krawczyk").degrade(Origin.COMPUTED, "x")
    assert d.origin is Origin.COMPUTED and d.label == ""
