"""Phase 3, reproduction R1-R11: the solutions Echeverria, Klinkhammer & Thorne (1991) construct,
computed with the exact geometric model of billiards.py and compared with their equations (ekt.py).

Every comparison number here is COMPUTED. Inputs are exact decimal strings, so that the geometry
and the printed equations see the same numbers.
"""

from __future__ import annotations

import math

import mpmath as mp

from . import billiards as bl
from . import ekt
from .provenance import Origin, Quantity, heuristic

S_GRID = heuristic("ekt_s_grid", 240, "points of the uniform part of the outer grid of collision positions")
S_FAR = heuristic("ekt_s_far", 3, "farthest collision position searched, in D, along the incoming line (before the traversal) or the exit line (after)")
S_NEAR_OVER_R = heuristic("ekt_s_near_over_R", "0.01", "closest collision position to the line's mouth-2 entry searched before the traversal, as a fraction of R")
S_NEAR_AFTER = heuristic("ekt_s_near_after", "1e-6", "closest collision position to the mouth-1 exit searched after the traversal, in D")
REFINE_SPLIT = heuristic("ekt_refine_max_split", 30, "further halvings of an unresolved nu- or s-interval in the adaptive refinement of the empty case")
WINDOW_GRID = heuristic("ekt_window_grid", 160, "extra collision positions in each refinement window around a known special point")


def _c(value, what):
    return Quantity(value, Origin.COMPUTED, what)


def line_theta_phi(theta, phi, B):
    """Initial line of Fig. 11: enters mouth 2 at polar angle pi + phi (point C2 + B(-cos phi, -sin phi))
    with direction theta + phi from the X axis (p. 1089: psi_A = theta + phi)."""
    th, ph, B = mp.mpf(theta), mp.mpf(phi), mp.mpf(B)
    return (1 - B * mp.cos(ph), -B * mp.sin(ph)), (mp.cos(th + ph), mp.sin(th + ph))


def self_inconsistent_crossing(p0, e, B):
    """Collision-free continuation: enter mouth 2 at p0, leave mouth 1, and cross the incoming line.
    Returns the s of the crossing on the incoming line (None if the lines do not cross ahead)."""
    E, w, _ = bl.traverse(p0, e, 1, mp.mpf(B))
    den = bl.cross(e, w)
    if den == 0:
        return None
    return bl.cross(bl.add(E, p0, -1), w) / den


def grid(lo, hi, n, windows=()):
    """Uniform grid on [lo, hi] plus dense windows (center, half-width) and a geometric approach to
    both ends."""
    lo, hi = float(lo), float(hi)
    pts = {lo + (hi - lo) * i / (n - 1) for i in range(n)}
    for c, hw in windows:
        c, hw = float(c), float(hw)
        pts |= {c - hw + 2 * hw * i / (WINDOW_GRID.value - 1) for i in range(WINDOW_GRID.value)}
    span = hi - lo
    pts |= {hi - span * 10.0 ** (-j / 4) for j in range(4, 40)}
    pts |= {lo + span * 10.0 ** (-j / 4) for j in range(4, 40)}
    return sorted(p for p in pts if lo <= p <= hi)


def solutions(p0, e, v1, m, k, B, R, s_grid, nu_grid=None, aim_grid=None, origin=None):
    """All single-collision solutions on the s grid: float bracketing, mpmath polish, verification.

    Returns (out, unresolved). ``out`` holds verified solutions and rejected candidates, each with
    its reasons; ``unresolved`` is the tagged list of everything the bracketing could not follow
    (see billiards.crossings)."""
    xm, digits = bl.backend_for(k, float(B))
    with mp.workdps(max(mp.mp.dps, digits)):
        conv = (lambda x: float(x)) if xm is bl.FLOAT else (lambda x: mp.mpf(x))
        pf, ef = (conv(p0[0]), conv(p0[1])), (conv(e[0]), conv(e[1]))
        Bf, Rf, vf = conv(B), conv(R), conv(v1)
        sg = [conv(s) for s in s_grid]
        log = []
        table = bl.sweep(pf, ef, sg, m, k, Bf, Rf, xm, nu_grid, aim_grid, log)
        cands, unresolved, rejected = bl.crossings(table, pf, ef, vf, m, k, Bf, Rf, xm, nu_grid, aim_grid, log=log)
    out = [{"candidate": c, "ok": False, "fails": [why]} for c, why in rejected]
    for c in cands:
        sol = bl.solve_mp(p0, e, *c, v1, m, k, B, R)
        if sol is None:
            out.append({"candidate": c, "ok": False, "fails": ["polish did not converge"]})
            continue
        if any(abs(sol[0] - o["sol"][0]) < mp.mpf(10) ** -20 for o in out if "sol" in o):
            continue
        ok, fails = bl.verify(p0, e, v1, *sol, m, k, B, R, origin=origin)
        out.append({"sol": sol, "ok": ok, "fails": fails, "angles": bl.paper_angles(p0, e, v1, *sol, m, k, B, R)})
    return out, unresolved


def unresolved_by_kind(unres):
    """Counts of the unresolved list by kind (s_count_change, s_bisect_lost, nu_count_change,
    nu_state_lost, nu_bisect_lost)."""
    counts = {}
    for u in unres:
        counts[u[0]] = counts.get(u[0], 0) + 1
    return counts


NU_KINDS = ("nu_count_change", "nu_state_lost", "nu_bisect_lost")
S_KINDS = ("s_count_change", "s_bisect_lost")
AIM_KINDS = ("aim_bracket_lost", "aim_bracket_jump")  # (kind, P, nu, aim_a, aim_b): at one (s, nu) point


def item_span(u, sg, p0, e):
    """The s-interval an unresolved item can spoil (the piece of the conservative measure of
    unresolved_fraction), its area in the fine measure, and whether it sits on the grid.
    ``sg`` is the sorted float grid. Returns (piece or None, area, on_grid)."""
    if u[0] in S_KINDS:
        a, b = sorted((float(u[1]), float(u[2])))
        return (a, b), 0.0, True
    if u[0] not in NU_KINDS and u[0] not in AIM_KINDS:
        return None, 0.0, False
    s_ = float(bl.dot(bl.add(u[1], p0, -1), e))
    i = min(range(len(sg)), key=lambda j: abs(sg[j] - s_))
    if abs(sg[i] - s_) <= 1e-12 * (1 + abs(s_)):
        # the linkage of v1(s) at a grid point is made on both neighbouring grid intervals
        piece = (sg[max(i - 1, 0)], sg[min(i + 1, len(sg) - 1)])
        cell = ((sg[i - 1] + sg[i]) / 2 if i else sg[i]) , ((sg[i] + sg[i + 1]) / 2 if i < len(sg) - 1 else sg[i])
        # an aim item sits at one nu: no area in the fine measure
        a = (cell[1] - cell[0]) * abs(float(u[3]) - float(u[2])) if u[0] in NU_KINDS else 0.0
        return piece, a, True
    j = max(k for k, x in enumerate(sg) if x <= s_) if s_ >= sg[0] else 0
    return (sg[j], sg[min(j + 1, len(sg) - 1)]), 0.0, False


def unresolved_fraction(unres, s_grid, p0, e):
    """Fraction of the search space Omega = [s_lo, s_hi] x [0, 2 pi) covered by unresolved intervals,
    measured two ways.

    - ``fine``: s-intervals count with their full nu range (their union along s); a nu-interval logged at a position of
      the s grid counts with that position's cell (half-way to each neighbour) times its nu width.
      A nu-interval logged off the grid (inside an s-bisection) has no cell of its own; those are
      counted in ``off_grid`` and add nothing here. An aim bracket dropped at one (s, nu) point
      (AIM_KINDS) has no nu width and adds nothing here either.
    - ``conservative``: the union, along s, of the grid intervals each unresolved item can spoil,
      with the full nu range. An unresolved nu-interval at a grid point can hide a branch of middle
      roots and so spoil the linkage of v1(s) on both neighbouring grid intervals; an off-grid item
      counts with the grid interval that contains it; so does a dropped aim bracket. This is an
      upper bound on the region where a solution could have escaped only under the bracketing
      criterion and on the declared grid (a double sign change between nodes is not covered, nor
      an aim root between the last defined node and an undefined one, which roots_by_bracketing
      does not log).

    Returns {"fine", "off_grid", "conservative"}."""
    sg = sorted(float(x) for x in s_grid)
    span = sg[-1] - sg[0]
    two_pi = 2 * math.pi
    area, off, pieces, s_pieces = 0.0, 0, [], []
    for u in unres:
        piece, a, on_grid = item_span(u, sg, p0, e)
        if piece is None:
            continue
        pieces.append(piece)
        if u[0] in S_KINDS:
            s_pieces.append(piece)
        elif on_grid:
            area += a
        else:
            off += 1

    def union(ps):
        covered, end = 0.0, -math.inf
        for a, b in sorted(ps):
            if b > end:
                covered += b - max(a, end)
                end = b
        return covered

    fine = union(s_pieces) / span + area / (span * two_pi)
    return {"fine": min(fine, union(pieces) / span), "off_grid": off, "conservative": union(pieces) / span}


def _refine_s(case, sa, sb, depth, xm):
    """Adaptive refinement of one unresolved s-interval. A leaf whose middle-root count still changes
    is resolved when v1 - target has one sign at every middle root at both ends; a leaf with equal
    counts is resolved when no matched branch crosses the target. Returns the leaves left."""
    conv = float if xm is bl.FLOAT else mp.mpf
    p0 = (conv(case["p0"][0]), conv(case["p0"][1]))
    e = (conv(case["e"][0]), conv(case["e"][1]))
    B, R, k, v1 = conv(case["B"]), conv(case["R"]), case["k"], conv(case["v1"])
    left = []

    def roots(s_, new=True):
        log = []
        rs = bl._dedupe(bl.middle_roots(bl.add(p0, e, s_), e, 1, k, B, R, xm, log=log), xm)
        return rs, bool(log) and new   # nu-items at the original ends are refined on their own (refine_nu)

    def leaf(sa, ra, sb, rb, d):
        (ra, fa), (rb, fb) = ra, rb
        if fa or fb:
            left.append((sa, sb, "the middle-root search at a refined end of this leaf logged unresolved nu-intervals"))
            return
        if len(ra) == len(rb):
            if any((a["v1"] > v1) != (b["v1"] > v1) for a, b in bl._match(ra, rb)):
                left.append((sa, sb, "a branch of v1(s) crosses the target: a candidate the main search did not follow"))
            return
        if d >= depth:
            vals = [r["v1"] - v1 for r in ra + rb]
            if vals and len({x > 0 for x in vals}) != 1:
                left.append((sa, sb, "middle-root count changes and v1 - target takes both signs at the ends"))
            return
        sm = (sa + sb) / 2
        if sm == sa or sm == sb:
            left.append((sa, sb, "number resolution reached"))
            return
        rm = roots(sm)
        leaf(sa, (ra, False), sm, rm, d + 1)
        leaf(sm, rm, sb, (rb, False), d + 1)

    leaf(conv(sa), roots(conv(sa), new=False), conv(sb), roots(conv(sb), new=False), 0)
    return left


def refine_unresolved(case, depth=None):
    """Refine every unresolved interval of a dangerous() search adaptively (REFINE_SPLIT more
    halvings), and report what is left and the fraction of the search space it covers."""
    depth = REFINE_SPLIT.value if depth is None else depth
    xm, digits = bl.backend_for(case["k"], float(case["B"]))
    conv = float if xm is bl.FLOAT else mp.mpf
    p0f = (conv(case["p0"][0]), conv(case["p0"][1]))
    ef = (conv(case["e"][0]), conv(case["e"][1]))
    before = case["unresolved"]
    after = []
    with mp.workdps(max(mp.mp.dps, digits)):
        for u in before:
            if u[0] == "nu_count_change":
                for na, nb, why in bl.refine_nu(u[1], 1, case["k"], conv(case["B"]), conv(case["R"]), xm, u[2], u[3], depth):
                    kind = why.split(":")[0] if why.startswith("aim_") else "nu_count_change"
                    after.append((kind, u[1], na, nb, why) if kind == "nu_count_change" else (kind, u[1], na, None, None, why))
            elif u[0] == "s_count_change":
                for sa, sb, why in _refine_s(case, u[1], u[2], depth, xm):
                    after.append(("s_count_change", sa, sb, why))
            else:
                after.append(u)  # lost brackets (s, nu or aim) are not refined here; they stay unresolved
    f0 = unresolved_fraction(before, case["s_grid"], p0f, ef)
    f1 = unresolved_fraction(after, case["s_grid"], p0f, ef)
    reasons = {}
    for u in after:
        key = u[-1] if isinstance(u[-1], str) else u[0]
        reasons[key] = reasons.get(key, 0) + 1
    return {"before": len(before), "before_kinds": unresolved_by_kind(before), "fraction_before": f0,
            "after": len(after), "after_kinds": unresolved_by_kind(after), "fraction_after": f1,
            "reasons_after": reasons, "depth": depth, "items_after": after}


CERT_CHUNKS = heuristic("ekt_cert_chunks", 8, "pieces of the s range of the empty case run as separate interval-exclusion tasks (parallelism only)")


def certify_items(case_args, items, s_grid, p0, e):
    """Interval exclusion (certify.exclude) on the box each unresolved item can spoil: its s-interval
    of the conservative measure (rounded outward) times the whole nu and aim ranges. Items that
    share an s-interval share the run. Returns per-item verdicts and the totals."""
    from . import certify as ct
    sg = sorted(float(x) for x in s_grid)
    runs, verdicts = {}, []
    for u in items:
        piece, _, _ = item_span(u, sg, p0, e)
        if piece is None:
            verdicts.append((u[0], None, "sem intervalo em s"))
            continue
        span = ct.s_range_outward(*piece)
        if span not in runs:
            r = ct.exclude(case_args, span)
            runs[span] = {"boxes": r["boxes"], "complete": r["complete"], "undecided": len(r["undecided"])}
        verdicts.append((u[0], span, "excluído" if runs[span]["complete"] else "indecidido"))
    count = lambda pred: sum(1 for v in verdicts if pred(v))
    return {"items": len(items), "runs": len(runs), "boxes": sum(r["boxes"] for r in runs.values()),
            "excluded": count(lambda v: v[2] == "excluído"), "undecided": count(lambda v: v[2] != "excluído"),
            "excluded_aim": count(lambda v: v[2] == "excluído" and v[0] in AIM_KINDS),
            "items_aim": count(lambda v: v[0] in AIM_KINDS), "certified_roots": 0, "verdicts": verdicts}


def certify_chunk(case_args, i, n):
    """Interval exclusion on the i-th of n equal pieces of the s range of dangerous(case_args)
    (before the traversal), times the whole nu and aim ranges."""
    from . import certify as ct
    R = mp.mpf(case_args[4])
    lo, hi = ct.s_range_outward(-S_FAR.value, -R * mp.mpf(S_NEAR_OVER_R.value))
    a, b = lo + (hi - lo) * i / n, (lo + (hi - lo) * (i + 1) / n if i < n - 1 else hi)
    r = ct.exclude(case_args, (a, b))
    return {"case": tuple(case_args), "i": i, "n": n, "s": (a, b), "boxes": r["boxes"], "excluded_by": r["excluded_by"],
            "undecided": r["undecided"][:20], "n_undecided": len(r["undecided"]), "complete": r["complete"],
            "range": (lo, hi)}


# Positive control of the exclusion (user decision of 2026-11-08, docs/phase3_plan.md): the same
# machine on cases of the same sub-model (one traversal, collision before mouth 2) whose solutions
# are known; no box that contains a known solution may be excluded. Krawczyk on what is left.

FIG7B_R, FIG7B_V1, FIG7B_LO = "1e-7", "0.5", "-1.5"   # r4_r5_r6 and symmetric(): R, v1 and the low end of s
CERT_CONTROL_BUDGET = heuristic("ekt_cert_control_budget", 300_000, "boxes per piece of s in the positive control when the piece holds no known solution (the rest is reported as not examined)")
# 20M: twice the 9,963,011 boxes the solution piece of R8 #3 needed on 2026-11-09 under the former
# 10M. A safety margin against future changes of the method (user decision of 2026-11-09), not a
# tuning of the result: with any budget that completes the piece, the outcome is the same.
CERT_CONTROL_SOL_BUDGET = heuristic("ekt_cert_control_sol_budget", 20_000_000, "boxes per piece of s in the positive control when the piece holds a known solution (the rest is reported as not examined); 2x the 9,963,011 needed by R8 #3 under the former 10M, as a margin against changes of the method, not a tuning of the result")


def control_case(spec):
    """(case for certify, outward s range, label) of a positive-control case: an R8 case of the
    dangerous triangle, or the symmetric trajectory of Fig. 7(b) (one traversal)."""
    from . import certify as ct
    kind, a = spec
    if kind == "r8":
        case = R8_CASES[a]
        R = mp.mpf(case[4])
        return case, ct.s_range_outward(-S_FAR.value, -R * mp.mpf(S_NEAR_OVER_R.value)), f"R8 ({'; '.join(x.replace('.', ',') for x in case)})"
    case = {"p0": ("0.5", "0"), "e": ("0", "1"), "v1": FIG7B_V1, "B": a, "R": FIG7B_R}
    return case, ct.s_range_outward(FIG7B_LO, -mp.mpf(FIG7B_R) / 10), f"Fig. 7(b), B = {a.replace('.', ',')}, R = 1e-7, v₁ = 0,5"


def control_known(spec):
    """Verified single-collision solutions (s, nu, aim) of a control case, from the bracketing search."""
    kind, a = spec
    if kind == "r8":
        sols = [x for x in dangerous(*R8_CASES[a])["sols"] if x["ok"]]
    else:
        sols = symmetric(a, FIG7B_R, FIG7B_V1, 1)[0]
    return [tuple(x["sol"]) for x in sols]


def control_pieces(lo, hi):
    """Pieces of [lo, hi] (lo < hi < 0) halving |s| from |lo| towards |hi|: the work concentrates
    near the mouth."""
    e = [lo]
    while e[-1] / 2 < hi:
        e.append(e[-1] / 2)
    e.append(hi)
    return list(zip(e, e[1:]))


def control_piece(spec, a, b, sols):
    """Exclusion on one piece of s of a control case, at the resolution of the empty-case proof
    (sides measured against 1/CERT_CHUNKS of the s range and the whole nu and aim ranges), with every periodic copy of the known solutions
    in the box watched, then the Krawczyk test on every undecided box the machine examined (boxes
    left on the stack by the budget are not tested). A piece that holds a known solution gets
    ekt_cert_control_sol_budget boxes; the others get ekt_cert_control_budget."""
    from . import certify as ct
    case, (lo, hi), _ = control_case(spec)
    inside = [x for x in sols if a <= ct.normalize(x)[0] <= b]
    nu, aim = ct.default_ranges()
    # the resolution of the empty-case proof: there each piece is 1/CERT_CHUNKS of the s range
    scale = ((hi - lo) / CERT_CHUNKS.value, nu[1] - nu[0], aim[1] - aim[0])
    copies = [(k, q) for k, x in enumerate(inside) for q in ct.periodic_copies(x)]
    budget = (CERT_CONTROL_SOL_BUDGET if inside else CERT_CONTROL_BUDGET).value
    r = ct.exclude(case, (a, b), max_boxes=budget, watch=[q for _, q in copies], scale=scale)
    n_ex = len(r["undecided"]) - r["unexamined"]
    examined, left = r["undecided"][:n_ex], r["undecided"][n_ex:]

    def vol(bx):
        return math.prod(hi - lo for lo, hi in bx)

    kraw = [ct.krawczyk(case, box) for box in examined]
    watch = [{"sol": inside[k], "copy": tuple(float(t) for t in q), "excluded_in": len(h["excluded_in"]),
              "undecided_in": [i for i in h["undecided_in"] if i < n_ex], "unexamined_in": sum(1 for i in h["undecided_in"] if i >= n_ex)}
             for (k, q), h in zip(copies, r["watch"])]
    return {"spec": spec, "s": (a, b), "boxes": r["boxes"], "excluded_by": r["excluded_by"], "undecided": examined,
            "unexamined": r["unexamined"], "volume": vol(((a, b), nu, aim)), "unexamined_volume": sum(vol(bx) for bx in left),
            "krawczyk": [(k["verdict"], k["inequalities"]) for k in kraw], "holds_solution": bool(inside), "watch": watch}


def control_summary(spec, known, pieces):
    """Per case: pieces of s examined to the end, fraction of the volume examined, boxes excluded,
    undecided and never examined, and, per known solution (counted once; its periodic copies in the
    box are listed with it), the piece of s, the box and the Krawczyk verdict. ``passed`` is False when
    the case has no known solution, or as soon as a copy of a known solution lies in an excluded box or
    in no box at all, or Krawczyk says "no root" in a box that holds one;
    ``complete_near_solutions`` says whether, besides, every copy ended in an examined undecided box
    (those stop at ekt_cert_min_width). A certified copy gives a CERTIFIED quantity labelled CERT_ROOT_MODEL; it is not a
    physical solution (no interval proof of the absence of extra collisions)."""
    from . import certify as ct
    from .provenance import CERT_ROOT_MODEL, certified
    two_pi = 2 * math.pi
    boxes = sum(p["boxes"] for p in pieces)
    undecided = sum(len(p["undecided"]) for p in pieces)
    unexamined = sum(p["unexamined"] for p in pieces)
    excluded = sum(n for p in pieces for n in p["excluded_by"].values())
    volume = sum(p.get("volume", 0) for p in pieces)
    examined_fraction = 1 - sum(p.get("unexamined_volume", 0) for p in pieces) / volume if volume else (1.0 if not unexamined else 0.0)
    complete_pieces = sum(1 for p in pieces if p["unexamined"] == 0)
    sol_pieces = [p for p in pieces if p.get("holds_solution")]
    roots = [(p, i) for p in pieces for i, (v, ineq) in enumerate(p["krawczyk"]) if v == "root" and ineq]
    # Krawczyk verdicts on the boxes that hold no known solution (copy), apart from those that do
    holding = [{i for w in p["watch"] for i in w["undecided_in"]} for p in pieces]
    other = [(v, ineq) for p, h in zip(pieces, holding) for i, (v, ineq) in enumerate(p["krawczyk"]) if i not in h]
    none = sum(1 for v, _ in other if v == "none")
    kraw_undecided = sum(1 for v, _ in other if v == "undecided")

    def matches(box, x):  # a root box holds a known solution, or one of its copies at nu +- 2 pi, aim +- 2 pi
        s_, nu, aim = (float(t) for t in ct.normalize(x))
        return any(box[0][0] <= s_ <= box[0][1] and box[1][0] <= nu + i * two_pi <= box[1][1] and box[2][0] <= aim + j * two_pi <= box[2][1]
                   for i in (-1, 0, 1) for j in (-1, 0, 1))

    sols = []
    for x in known:
        key = tuple(float(t) for t in ct.normalize(x))
        cps = [(p, w) for p in pieces for w in p["watch"] if tuple(float(t) for t in ct.normalize(w["sol"])) == key]
        sol_boxes = [p["undecided"][i] for p, w in cps for i in w["undecided_in"]]
        widths = [max(hi - lo for lo, hi in bx[1:]) for bx in sol_boxes]
        verdicts = [p["krawczyk"][i] for p, w in cps for i in w["undecided_in"]]
        # every copy in an examined undecided box: those stop at ekt_cert_min_width (relative)
        refined = bool(cps) and all(w["undecided_in"] for _, w in cps)
        lost = sum(1 for _, w in cps if not w["undecided_in"] and not w["unexamined_in"] and not w["excluded_in"])  # in no box at all
        in_excluded = sum(w["excluded_in"] for _, w in cps)
        none_here = sum(1 for v, _ in verdicts if v == "none")  # "no root" where a known solution lies: a contradiction
        cert = refined and all(any(p["krawczyk"][i] == ("root", True) for i in w["undecided_in"]) for p, w in cps)
        boxes_txt = "; ".join("[" + ", ".join(f"[{lo!r}, {hi!r}]" for lo, hi in bx) + "]" for bx in sol_boxes)
        q = (certified(tuple(sol_boxes), CERT_ROOT_MODEL,
                       f"Krawczyk (K(X) no interior de X) em cada caixa (s, ν, mira) {boxes_txt}: uma única raiz de (miss, C, V) em cada uma, "
                       f"com d·n_X < 0, along > 0 e kick > 0 em toda a caixa; cópias periódicas da mesma solução; hipóteses: "
                       + ct.HYPOTHESES.format(s=f"[{_f(mp.mpf(cps[0][0]['s'][0]), 6)}; {_f(mp.mpf(cps[0][0]['s'][1]), 6)}]", prec=ct.CERT_PREC.value)
                       + f"; controle positivo de {spec}. O valor são as caixas; o ponto do bracketing ({key}) é COMPUTED")
             if cert else None)
        sols.append({"sol": key, "piece": cps[0][0]["s"] if cps else None, "copies": len(cps), "in_excluded": in_excluded,
                     "lost": lost, "none_here": none_here,
                     "unexamined": any(w["unexamined_in"] for _, w in cps), "refined": refined, "certified": cert, "quantity": q,
                     "width": max(widths) if widths else None, "verdicts": sorted(set(v for v, _ in verdicts))})
    root_boxes = [p["undecided"][i] for p, i in roots]
    unmatched = [bx for bx in root_boxes if not any(matches(bx, x) for x in known)]
    passed = (len(known) > 0 and all(x["copies"] for x in sols) and sum(x["in_excluded"] + x["lost"] + x["none_here"] for x in sols) == 0)
    return {"spec": spec, "boxes": boxes, "excluded": excluded, "undecided": undecided, "unexamined": unexamined,
            "pieces": len(pieces), "complete_pieces": complete_pieces, "solution_pieces": len(sol_pieces),
            "solution_pieces_complete": sum(1 for p in sol_pieces if p["unexamined"] == 0), "examined_fraction": examined_fraction,
            "known": len(known), "copies": sum(x["copies"] for x in sols), "in_excluded": sum(x["in_excluded"] for x in sols),
            "lost": sum(x["lost"] for x in sols), "none_in_solution_boxes": sum(x["none_here"] for x in sols),
            "refined": sum(x["refined"] for x in sols), "certified": sum(x["certified"] for x in sols),
            "root_boxes": len(root_boxes), "unmatched_roots": len(unmatched), "none": none, "kraw_undecided": kraw_undecided,
            "passed": passed, "complete_near_solutions": passed and all(x["refined"] for x in sols),
            "sols": sols}


# R1 --------------------------------------------------------------------------------------------

def r1_traversal_rule(n=25):
    """Exit direction of the identification model against eq. 2.1b on a grid of (phi, psi_in)."""
    worst = mp.mpf(0)
    B = mp.mpf("0.1")
    for i in range(n):
        phi = -mp.pi / 2 + mp.pi * (i + mp.mpf(1) / 2) / n
        X = (1 - B * mp.cos(phi), -B * mp.sin(phi))  # Fig. 4: angular location phi on mouth 2
        normal_in = (mp.cos(phi), mp.sin(phi))
        for j in range(n):
            psi_in = phi - mp.pi / 2 + mp.pi * (j + mp.mpf(1) / 2) / n  # entering: v . normal_in > 0
            v = (mp.cos(psi_in), mp.sin(psi_in))
            if bl.dot(v, normal_in) <= 0:
                continue
            _, w, dt = bl.traverse(X, v, 1, B)
            got = mp.atan2(w[1], w[0])
            want = ekt.traversal_psi_out(psi_in, phi).value
            dev = abs((got - want + mp.pi) % (2 * mp.pi) - mp.pi)
            worst = max(worst, dev)
    return {"max_dev": _c(worst, "max |psi_out(model) - eq. 2.1b| over the grid, rad"), "time_shift": dt}


# R2 --------------------------------------------------------------------------------------------

def r2_collision_law(n=12):
    """Elastic collisions of equal balls with the speed exchange of eqs. 2.2a-b imposed: every
    contact normal that satisfies it is one of the two cases 2.3 or 2.4."""
    counts = {"velocity": 0, "mirror": 0, "other": 0}
    worst = mp.mpf(0)
    for i in range(n):
        for j in range(n):
            a1 = 2 * mp.pi * i / n + mp.mpf("0.1")
            a2 = 2 * mp.pi * j / n + mp.mpf("0.37")
            v1 = (mp.mpf("0.8") * mp.cos(a1), mp.mpf("0.8") * mp.sin(a1))
            v2 = (mp.mpf("1.3") * mp.cos(a2), mp.mpf("1.3") * mp.sin(a2))

            def f(nu):
                n_ = (mp.cos(nu), mp.sin(nu))
                k_ = bl.dot(bl.add(v2, v1, -1), n_)
                v1p = bl.add(v1, n_, k_)
                return bl.dot(v1p, v1p) - bl.dot(v2, v2)

            nus = [mp.pi * t / 400 for t in range(401)]  # n and -n give the same collision
            for nu in bl.roots_by_bracketing(f, nus, mp.mpf(10) ** (-mp.mp.dps + 10)):
                nu = mp.findroot(f, nu)
                n_ = (mp.cos(nu), mp.sin(nu))
                k_ = bl.dot(bl.add(v2, v1, -1), n_)
                v1p = bl.add(v1, n_, k_)
                ex = bl.norm(bl.add(v1p, v2, -1), bl.MP)
                u = bl.unit(bl.add(v1, v2), bl.MP)
                refl = bl.add(bl.scale(u, 2 * bl.dot(v2, u)), v2, -1)  # v2 reflected in the line along v1 + v2
                mi = bl.norm(bl.add(v1p, refl, -1), bl.MP)
                tol = mp.mpf(10) ** (-mp.mp.dps + 15)
                if ex < tol:
                    counts["velocity"] += 1
                    dirn = bl.unit(bl.add(v2, v1, -1), bl.MP)  # eq. 2.3b direction
                    worst = max(worst, abs(abs(bl.dot(dirn, n_)) - 1))
                elif mi < tol:
                    counts["mirror"] += 1
                    worst = max(worst, abs(abs(bl.dot(u, n_)) - 1))  # eq. 2.4b direction
                else:
                    counts["other"] += 1
    return {"counts": counts, "max_dev_direction": _c(worst, "max deviation of the line of centers from eqs. 2.3b/2.4b")}


# R3, R4, R5, R6: the symmetric trajectory of Figs. 6 and 7 -------------------------------------

SYM_P0, SYM_E = (mp.mpf("0.5"), mp.mpf(0)), (mp.mpf(0), mp.mpf(1))


def symmetric(B, R, v1, k, lo="-1.5"):
    """Single-collision solutions with k traversals for the trajectory perpendicular to the line of
    centers through its midpoint (Figs. 6 and 7). s is the height of the collision."""
    windows = [(0, 6 * float(R))]
    sg = grid(mp.mpf(lo), -float(R) / 10, S_GRID.value, windows)
    sols, unresolved = solutions(SYM_P0, SYM_E, mp.mpf(v1), 1, k, mp.mpf(B), mp.mpf(R), sg)
    return [s for s in sols if s["ok"]], unresolved, [s for s in sols if not s["ok"]]


def r3_velocity_exchange(B="0.01", Rs=("1e-4", "1e-5", "1e-6"), v1="0.5"):
    """Fig. 6(b): the velocity-exchange solution; its loop speed against v1' = 1 - 2B (p. 1083)."""
    out = []
    for R in Rs:
        sols, unres, _ = symmetric(B, R, v1, 1)
        vx = [s for s in sols if abs(s["sol"][0]) < 10 * mp.mpf(R)]
        speeds = [s["angles"]["v1p"] for s in vx]
        target = ekt.velocity_exchange_speed(mp.mpf(B)).value
        out.append({"R": R, "n": len(vx), "dev": _c([sp - target for sp in speeds], "loop speed - (1 - 2B)"),
                    "unresolved": len(unres)})
    return out


def _psi_of(sol_entry):
    """Angle psi of Fig. 7 from the measured collision height: tan psi = |y| / (1/2)."""
    return mp.atan(2 * abs(sol_entry["sol"][0]))


def r4_r5_r6(Bs=("0.04", "0.02", "0.01", "0.005"), R="1e-7", v1="0.5"):
    """Figs. 7(b), 7(c), 7(d): mirror-exchange solutions with 1, 2, 3 traversals."""
    v1 = mp.mpf(v1)
    rows = []
    for B in Bs:
        row = {"B": B}
        for k in (1, 2, 3):
            sols, unres, _ = symmetric(B, R, v1, k)
            mirror = [s for s in sols if abs(s["sol"][0]) > 100 * mp.mpf(R)]
            row[k] = {"n": len(mirror), "unresolved": len(unres),
                      "y": [s["sol"][0] for s in mirror], "psi": [_psi_of(s) for s in mirror], "sols": mirror}
        if row[1]["n"] == 1:
            psi = row[1]["psi"][0]
            row["dev31"] = _c(mp.tan(psi) / v1 - 1, "tan(psi)/v1 - 1 against eq. 3.1")
        if row[2]["n"] == 1:
            psi = row[2]["psi"][0]
            row["dev32"] = _c((mp.sin(psi) + mp.tan(psi)) / (2 * v1) - 1, "(sin psi + tan psi)/(2 v1) - 1 against eq. 3.2")
            row["times"] = fig7c_times(row[2]["sols"][0], B, R, v1, psi)
        rows.append(row)
    return rows


def fig7c_times(entry, B, R, v1, psi):
    """Event times of the two-traversal solution from the verified history, against p. 1084."""
    sol = entry["sol"]
    with mp.workdps(bl.POLISH_DPS.value):
        h = bl.history(SYM_P0, SYM_E, v1, *sol, 1, 2, mp.mpf(B), mp.mpf(R))
        seg = {s[5]: s for s in h["segments"]}
        measured = [seg["return"][3], seg["loop1"][3], mp.mpf(0), seg["loop1"][4], seg["loop0"][4]]
    printed = ekt.fig7c_event_times(psi).value
    return {"measured": measured, "printed": printed,
            "max_dev": _c(max(abs(a - b) for a, b in zip(measured, printed)), "max |T(model) - T(p. 1084)|")}


def r4_rest(B="0.01", R="1e-6", xs=("0.3", "0.5", "0.7")):
    """v1 -> 0 limit of Fig. 7(b) (p. 1084): a ball at rest on the line of centers is knocked into
    mouth 2 and comes back to rest; the collision time is arbitrary. The verifier is run with the
    ball at rest at several positions; the loop speed is compared with the printed v2 = 1 (B << 1)."""
    out = []
    for x in xs:
        p0 = (mp.mpf(x), mp.mpf(0))
        with mp.workdps(bl.POLISH_DPS.value):
            ok, fails = bl.verify(p0, SYM_E, 0, 0, 0, mp.pi, 1, 1, mp.mpf(B), mp.mpf(R))
            st = bl._loop_state(p0, mp.mpf(0), 1, mp.pi, 1, mp.mpf(B), mp.mpf(R), bl.MP)
        out.append({"x": x, "ok": ok, "fails": fails, "speed": _c(st["speed"], "loop speed of the ball at rest, knocked along the line of centers")})
    return out


# R7: many traversals (Sec. III B) -------------------------------------------------------------

def r7_multitraversal(B="0.05", R="1e-4", v1="3", psi_A="0.6", h="0.5", ns=range(1, 7), lo="-50"):
    """Sec. III B: one solution for each number n of traversals; L(n) against eq. 3.7 and the
    smallest n found against eq. 3.8. The initial line has direction psi_A and passes at signed
    distance h from the center of mouth 2; L is measured from the point of closest approach."""
    B, R, v1, psi_A, h = (mp.mpf(x) for x in (B, R, v1, psi_A, h))
    e = (mp.cos(psi_A), mp.sin(psi_A))
    p0 = (1 - h * mp.sin(psi_A), h * mp.cos(psi_A))  # closest approach to C2 = (1, 0)
    rows = []
    for n in ns:
        sg = grid(mp.mpf(lo), -mp.mpf("0.5"), S_GRID.value)
        sols, unres = solutions(p0, e, v1, 1, n, B, R, sg)
        good = [s for s in sols if s["ok"]]
        Ls = [-s["sol"][0] for s in good]
        frac = unresolved_fraction(unres, sg, (float(p0[0]), float(p0[1])), (float(e[0]), float(e[1])))
        rows.append({"n": n, "L": Ls, "unresolved": len(unres), "rejected": [s["fails"] for s in sols if not s["ok"]],
                     "unresolved_kinds": unresolved_by_kind(unres), "fraction": frac,
                     "n_37": [ekt.multitraversal_n(L, h, v1).value for L in Ls]})
    return {"rows": rows, "n_min_38": ekt.multitraversal_n_min(h, v1).value}


# R8, R9: the dangerous triangle against eqs. 4.2-4.3 and 4.12-4.13 ------------------------------

def dangerous(theta, phi, lam1, B, R, k=1, extra_windows=(), on_edge=False):
    th, ph, lam, B, R = (mp.mpf(x) for x in (theta, phi, lam1, B, R))
    if not (ekt.in_dangerous_triangle(th, ph, B) or (on_edge and th == ph and B / 2 < ph < mp.pi / 2)):
        raise ValueError("(theta, phi) outside the dangerous triangle of eq. 4.7")
    v1 = ekt.v1_from_lambda(lam, th, ph, B, R).value
    p0, e = line_theta_phi(th, ph, B)
    s0 = self_inconsistent_crossing(p0, e, B)
    windows = list(extra_windows) + ([(s0, 40 * float(R))] if s0 is not None and s0 < 0 else [])
    s_range = (-S_FAR.value, -float(R) * float(S_NEAR_OVER_R.value))
    sg = grid(mp.mpf(s_range[0]), s_range[1], S_GRID.value, windows)
    sols, unres = solutions(p0, e, v1, 1, k, B, R, sg)
    return {"v1": v1, "p0": p0, "e": e, "s0": s0, "sols": sols, "unresolved": unres, "theta": th, "phi": ph, "B": B, "R": R, "lam1": lam,
            "s_range": s_range, "s_grid": sg, "k": k}


def r8_eq42_43(case):
    """Residuals of the printed eqs. 4.2-4.3 at every verified geometric solution, plus eq. 4.4b
    (d) and eq. 4.5 (v2 = rho v1) against the measured values, and the spurious tests of 4.6."""
    rows = []
    for s_ in case["sols"]:
        if not s_["ok"]:
            continue
        a = s_["angles"]
        sign = 1 if a["d"] > 0 else -1
        with mp.workdps(bl.POLISH_DPS.value):
            r42, r43 = ekt.eq42_43(a["alpha"], a["beta"], case["theta"], case["phi"], case["v1"], case["B"], case["R"], sign)
            rows.append({
                "s": s_["sol"][0], "class": "I" if a["beta"] > 0 else "II", "alpha": a["alpha"], "beta": a["beta"],
                "eps": a["alpha"] - 2 * case["phi"], "d": a["d"],
                "res": _c(max(abs(r42), abs(r43)), "max |residual| of eqs. 4.2, 4.3 at the geometric solution"),
                "dev_d": _c(abs(ekt.d_of(a["alpha"], a["beta"], sign, case["R"]) - a["d"]), "|d(eq. 4.4b) - d(model)|"),
                "dev_v2": _c(abs(ekt.rho_of(a["alpha"], a["beta"]) * case["v1"] - a["v2"]), "|rho v1 (eq. 4.5) - v2(model)|"),
                "sign_rule": (mp.sign(a["beta"]) == sign),  # eq. 4.6c: sign(beta) = sign(d)
                "dev_a6": _c(abs(ekt.appendix_a_distances(a["alpha"], a["beta"], case["theta"], case["phi"], case["B"]).value["a"] + s_["sol"][0]),
                             "|a(eq. A6) - PQ(model)|, PQ = -s"),
                "lam2": ekt.lambda_from_speed(a["v2"], case["theta"], case["phi"], case["B"], case["R"]).value,
            })
    return rows


def eq42_43_roots(case, n_alpha=48):
    """Roots of the printed eqs. 4.2-4.3, from starting points on an alpha grid times a beta set
    that is dense near beta = 0 on a log scale (class I and II deflections range from ~B R to ~1).
    Roots listed in eq. 4.6 are dropped. The others are split by the sign of the Appendix A
    distances a, a', c (eqs. A6, A7, A19): all positive, or some negative (no ball can travel a
    negative distance along its path)."""
    ph, R = case["phi"], case["R"]
    physical, negative = [], []
    with mp.workdps(bl.POLISH_DPS.value):
        betas = [mp.mpf(0)] + [sg * mp.mpf(10) ** (-k / mp.mpf(2)) for k in range(0, 17) for sg in (1, -1)]
        alphas = [mp.mpf("0.02") + (mp.pi - mp.mpf("0.04")) * i / (n_alpha - 1) for i in range(n_alpha)]
        alphas += [2 * ph + case["R"] * j for j in (-30, -10, -3, -1, 1, 3, 10, 30)]
        for sign in (1, -1):
            for a0 in alphas:
                for b0 in betas:
                    r = ekt.solve_eq42_43(case["theta"], ph, case["v1"], case["B"], R, sign, (a0, b0))
                    if r is None:
                        continue
                    a, b = r
                    if not (0 < a < mp.pi) or abs(b) >= mp.pi / 2:
                        continue
                    if max(abs(x) for x in ekt.eq42_43(a, b, case["theta"], ph, case["v1"], case["B"], R, sign)) > mp.mpf(10) ** (-mp.mp.dps + 15):
                        continue
                    if ekt.is_spurious(a, b, ph, sign, R, mp.mpf(10) ** -15):
                        continue
                    if any(abs(a - f[0]) < mp.mpf(10) ** -15 and abs(b - f[1]) < mp.mpf(10) ** -15 for f in physical + negative):
                        continue
                    dist = ekt.appendix_a_distances(a, b, case["theta"], ph, case["B"]).value
                    target = physical if min(dist["a"], dist["a_p"], dist["c"]) > 0 else negative
                    target.append((a, b, sign, dist))
    return physical, negative


def r9_central(theta="0.2", phi="0.6", lam1="0.5", Bs=("0.04", "0.02", "0.01", "0.005"), ratio="0.1"):
    """Central region: geometric (beta, eps, lam2) against eqs. 4.12-4.13, and the two forms of B20."""
    rows = []
    for B in Bs:
        R = mp.mpf(B) * mp.mpf(ratio)
        case = dangerous(theta, phi, lam1, B, R)
        for s_ in case["sols"]:
            if not s_["ok"]:
                continue
            a = s_["angles"]
            if abs(a["alpha"] - 2 * case["phi"]) > mp.mpf("0.1"):
                continue  # not a central-region (near self-inconsistent) solution
            sign = 1 if a["beta"] > 0 else -1
            beta_p, eps_p, lam2_p = ekt.central_solution(case["theta"], case["phi"], case["lam1"], case["B"], case["R"], sign).value
            lam2 = ekt.lambda_from_speed(a["v2"], case["theta"], case["phi"], case["B"], case["R"]).value
            sigma = (1 - sign) // 2
            b20a = ekt.b20_first_form(case["theta"], case["phi"], case["lam1"], case["B"], sign).value
            rows.append({"B": B, "class": "I" if sign == 1 else "II",
                         "dev_beta": _c(a["beta"] / beta_p - 1, "beta(model)/beta(4.12a or 4.13a) - 1"),
                         "dev_eps": _c((a["alpha"] - 2 * case["phi"]) / eps_p - 1, "eps(model)/eps(4.12b or 4.13b) - 1"),
                         "dev_lam2": _c((lam2 - sigma) / (lam2_p - sigma) - 1, "(lam2 - sigma)(model) / (4.12c or 4.13c) - 1"),
                         "dev_b20_first": _c((lam2 - sigma) / (b20a - sigma) - 1, "(lam2 - sigma)(model) / (B20 first form) - 1")})
    return rows


# R10, R11: the upper edge theta = phi (Fig. 14) -------------------------------------------------

def r10_upper_edge(phis, B="0.01", R="1e-5", lam1="0.5"):
    """Fig. 14: class-I deflection beta on the upper edge theta = phi, against B32 and 4.14/4.15."""
    rows = []
    for ph in phis:
        case = dangerous(ph, ph, lam1, B, R, on_edge=True)
        good = [s for s in case["sols"] if s["ok"]]
        cls1 = [s for s in good if s["angles"]["beta"] > 0]
        cls2 = [s for s in good if s["angles"]["beta"] < 0]
        b32 = ekt.b30_b32(case["theta"], case["phi"], case["lam1"], case["B"], case["R"], 1).value
        u = ekt.upper_edge_414_415(case["phi"], case["lam1"], case["B"], case["R"]).value
        rows.append({"phi": ph, "beta_I": [s["angles"]["beta"] for s in cls1], "s_I": [s["sol"][0] for s in cls1],
                     "beta_II": [s["angles"]["beta"] for s in cls2], "b32": b32["beta"], "b414": u["4.14"][0], "b415": u["4.15"][0],
                     "unresolved": len(case["unresolved"]), "rejected": [s["fails"] for s in case["sols"] if not s["ok"]]})
    return rows


def after_traversal_setup(theta, phi, lam1, B, R):
    """The collision-free exit line from mouth 1 (point, direction), v1 and the s grid used after
    the traversal."""
    th, ph, lam, Bm, Rm = (mp.mpf(x) for x in (theta, phi, lam1, B, R))
    v1 = ekt.v1_from_lambda(lam, th, ph, Bm, Rm).value
    a0, e0 = line_theta_phi(th, ph, Bm)
    E, w, _ = bl.traverse(a0, e0, 1, Bm)
    sg = grid(mp.mpf(S_NEAR_AFTER.value), mp.mpf(S_FAR.value), S_GRID.value)
    return {"E": E, "w": w, "v1": v1, "s_grid": sg, "origin": (a0, e0), "B": Bm, "R": Rm}


def after_traversal(theta, phi, lam1, B, R, k=1):
    """Single-collision solutions whose first encounter with the collision is on the collision-free
    exit line from mouth 1 (after one traversal), with k more traversals between the two visits."""
    a = after_traversal_setup(theta, phi, lam1, B, R)
    return solutions(a["E"], a["w"], a["v1"], 1, k, a["B"], a["R"], a["s_grid"], origin=a["origin"])


def r11_after_traversal(phis, B="0.01", R="1e-5", lam1="0.5"):
    """Lacuna 16 (p. 1091): on the upper edge near phi ~ pi/4, look for a class-II solution whose
    first encounter with the collision is after one traversal. The incoming line is then the
    collision-free exit from mouth 1; between the two visits there is one more traversal."""
    rows = []
    for ph in phis:
        sols, unres = after_traversal(ph, ph, lam1, B, R)
        good = [s for s in sols if s["ok"]]
        rows.append({"phi": ph, "n": len(good), "beta": [s["angles"]["beta"] for s in good],
                     "s": [s["sol"][0] for s in good], "unresolved": len(unres),
                     "rejected": [s["fails"] for s in sols if not s["ok"]]})
    return rows


# Driver for run.py: R1-R11 as independent tasks in worker processes ------------------------------

R8_CASES = (
    ("0.2", "0.6", "0.5", "0.01", "0.001"),     # central region
    ("0.2", "0.6", "0.25", "0.2", "0.05"),      # large mouth
    ("0.0", "0.4", "0.75", "0.3", "0.135"),     # large mouth and ball, R/B = 0.45
    ("0.05", "0.1", "0.5", "0.01", "0.001"),    # near the left corner
    ("1.0", "1.05", "0.5", "0.01", "0.001"),    # near the upper edge
    ("-0.78", "0.9", "0.5", "0.1", "0.02"),     # near the lower left edge
)
R9_BS = ("0.04", "0.02", "0.01", "0.005")
REFINE_FACTOR = heuristic("ekt_refine_factor", 2, "grid refinement factor of the stabilization search (empty case)")
R46_BS = ("0.04", "0.02", "0.01", "0.005")
# positive control of the exclusion (see control_case): the five R8 cases with solutions, and Fig. 7(b)
CONTROL_CASES = tuple(("r8", i) for i in (0, 1, 3, 4, 5)) + tuple(("fig7b", B) for B in R46_BS)
R10_PHIS = ("0.05", "0.1", "0.2", "0.3", "0.5", "0.7", "0.75", "0.78", "0.8", "0.82", "0.9", "1.1", "1.3", "1.5")
R11_PHIS = ("0.7", "0.75", "0.78", "0.79", "0.82", "0.85")


def _init_worker():
    mp.mp.dps = bl.POLISH_DPS.value


def _task(name, arg):
    if name == "r1":
        return r1_traversal_rule()
    if name == "r2":
        return r2_collision_law()
    if name == "rest":
        return r4_rest()
    if name == "r3":
        return r3_velocity_exchange(Rs=(arg,))[0]
    if name == "r46":
        rows = r4_r5_r6(Bs=(arg,))
        for k in (1, 2, 3):
            rows[0][k].pop("sols")
        return rows[0]
    if name == "r7":
        return r7_multitraversal(ns=(arg,))
    if name == "r7b":
        return r7_multitraversal(v1="2", h="1", ns=(arg,))
    if name == "r8x":
        kind, k, refined, lam_override = arg
        th, ph, lam, B, R = R8_CASES[2]
        lam = lam_override or lam
        global S_GRID
        saved = (bl.NU_GRID, bl.THETA_GRID, S_GRID)
        try:
            if refined:  # grids multiplied by REFINE_FACTOR for this task only
                f = REFINE_FACTOR.value
                bl.NU_GRID = Quantity(f * bl.NU_GRID.value, Origin.HEURISTIC, "ekt_nu_grid x ekt_refine_factor")
                bl.THETA_GRID = Quantity(f * bl.THETA_GRID.value, Origin.HEURISTIC, "ekt_theta_grid x ekt_refine_factor")
                S_GRID = Quantity(f * S_GRID.value, Origin.HEURISTIC, "ekt_s_grid x ekt_refine_factor")
            if kind == "before":
                case = dangerous(th, ph, lam, B, R, k=k)
                sols, unres = case["sols"], case["unresolved"]
            else:
                sols, unres = after_traversal(th, ph, lam, B, R, k=k)
        finally:
            bl.NU_GRID, bl.THETA_GRID, S_GRID = saved
        s_range = (-S_FAR.value, -float(R) * float(S_NEAR_OVER_R.value)) if kind == "before" else (float(S_NEAR_AFTER.value), S_FAR.value)
        if kind == "before":
            frac = unresolved_fraction(unres, case["s_grid"], tuple(float(x) for x in case["p0"]), tuple(float(x) for x in case["e"]))
        else:
            a = after_traversal_setup(th, ph, lam, B, R)
            frac = unresolved_fraction(unres, a["s_grid"], tuple(float(x) for x in a["E"]), tuple(float(x) for x in a["w"]))
        return {"case": R8_CASES[2], "kind": kind, "k": k, "refined": refined, "lam1": lam, "ok": len([x for x in sols if x["ok"]]),
                "rejected": [x["fails"] for x in sols if not x["ok"]], "unresolved": len(unres),
                "unresolved_kinds": unresolved_by_kind(unres), "s_range": s_range, "fraction": frac}
    if name == "r8":
        case = dangerous(*arg)
        phys, neg = eq42_43_roots(case)
        free_ok, _, overlaps = bl.free_history(bl.add(case["p0"], case["e"], -1), case["e"], case["v1"], case["B"], case["R"])
        return {"args": arg, "v1": case["v1"], "free_consistent": free_ok, "free_min_distance": min(o[2] for o in overlaps) if overlaps else None, "rows": r8_eq42_43(case), "physical": [p[:3] for p in phys],
                "negative": [(a, b, s, d) for a, b, s, d in neg], "rejected": [s["fails"] for s in case["sols"] if not s["ok"]],
                "unresolved": len(case["unresolved"]), "unresolved_kinds": unresolved_by_kind(case["unresolved"]), "s_range": case["s_range"],
                "fraction": unresolved_fraction(case["unresolved"], case["s_grid"], tuple(float(x) for x in case["p0"]), tuple(float(x) for x in case["e"]))}
    if name == "r8r":
        case = dangerous(*arg)
        ref = refine_unresolved(case)
        cert = certify_items(arg, ref["items_after"], case["s_grid"], tuple(float(x) for x in case["p0"]), tuple(float(x) for x in case["e"]))
        ref.pop("items_after")
        return {"case": arg, "ok": len([x for x in case["sols"] if x["ok"]]), "rejected": [x["fails"] for x in case["sols"] if not x["ok"]],
                **ref, "cert_items": cert}
    if name == "r8c":
        return certify_chunk(R8_CASES[2], *arg)
    if name == "r8pk":
        return control_known(arg)
    if name == "r8p":
        return control_piece(*arg)
    if name == "r9":
        return r9_central(Bs=(arg,))
    if name == "r10":
        return r10_upper_edge((arg,))[0]
    if name == "r10x":
        ph = mp.mpf(arg)
        case = {"theta": ph, "phi": ph, "B": mp.mpf("0.01"), "R": mp.mpf("1e-5"), "lam1": mp.mpf("0.5")}
        case["v1"] = ekt.v1_from_lambda(case["lam1"], ph, ph, case["B"], case["R"]).value
        phys, neg = eq42_43_roots(case)
        b = ekt.b30_b32(ph, ph, case["lam1"], case["B"], case["R"], -1).value
        disc = b["p"] ** 2 / 4 - b["q"]
        return {"physical": [p_[:3] for p_ in phys], "negative": neg,
                "b29_class2": [-b["p"] / 2 + sg * mp.sqrt(disc) for sg in (1, -1)] if disc >= 0 else []}
    if name == "r11":
        return r11_after_traversal((arg,))[0]
    raise ValueError(name)


QUICK_TASKS = (
    ("r1", None), ("r2", None), ("rest", None), ("r3", "1e-4"),
    ("r8", R8_CASES[0]),   # central region: geometry against eqs. 4.2-4.3, with A6 and the sign rule
    ("r10", "0.8"), ("r11", "0.85"),
)


def full_tasks():
    tasks = [("r1", None), ("r2", None), ("rest", None)]
    tasks += [("r3", R) for R in ("1e-4", "1e-5", "1e-6")]
    tasks += [("r46", B) for B in R46_BS]
    tasks += [("r7", n) for n in range(1, 7)]
    tasks += [("r7b", n) for n in range(1, 7)]
    tasks += [("r8", c) for c in R8_CASES]
    tasks += [("r8x", (kind, k, False, None)) for kind, k in (("before", 2), ("before", 3), ("after", 1), ("after", 2))]
    tasks += [("r8x", ("before", k, True, None)) for k in (1, 2, 3)]
    tasks += [("r8x", ("before", 1, False, lam)) for lam in ("0.1", "0.25", "0.5", "0.9")]
    tasks += [("r8r", R8_CASES[2])]
    tasks += [("r8c", (i, CERT_CHUNKS.value)) for i in range(CERT_CHUNKS.value)]
    tasks += [("r8pk", spec) for spec in CONTROL_CASES]   # each one submits its r8p pieces when done (start)
    tasks += [("r9", B) for B in R9_BS]
    tasks += [("r10", ph) for ph in R10_PHIS]
    tasks += [("r10x", ph) for ph in ("0.3", "0.9", "1.1")]
    tasks += [("r11", ph) for ph in R11_PHIS]
    return tasks


def start(quick=False, workers=None):
    """Submit the phase 3 tasks (QUICK_TASKS, or R1-R11 in full) to worker processes and return a
    handle for collect(); the caller can do other work meanwhile. Deterministic grids only; each
    task is a pure function of its arguments."""
    from concurrent.futures import ProcessPoolExecutor
    import threading
    tasks = list(QUICK_TASKS) if quick else full_tasks()
    ex = ProcessPoolExecutor(max_workers=workers, initializer=_init_worker)
    follow, lock, submitted = [], threading.Lock(), []

    def pieces_of(spec, ready):  # positive control: the exclusion pieces of a case, once its known solutions exist
        def submit(fut):
            try:
                if fut.exception() is None:  # otherwise collect() raises it
                    _, s_range, _ = control_case(spec)
                    with lock:
                        for a, b in control_pieces(*s_range):
                            arg = (spec, a, b, fut.result())
                            follow.append((("r8p", arg), ex.submit(_task, "r8p", arg)))
            finally:
                ready.set()
        return submit

    futures = []
    for n, a in tasks:
        f = ex.submit(_task, n, a)
        if n == "r8pk":
            ready = threading.Event()
            submitted.append(ready)
            f.add_done_callback(pieces_of(a, ready))
        futures.append(f)
    return ex, tasks, futures, (follow, submitted)


def collect(handle):
    """Results of start(), as a dict of (argument, result) lists keyed by task, in task order; the
    pieces of the positive control (r8p) follow, in submission order."""
    ex, tasks, futures, (follow, submitted) = handle
    try:
        done = [f.result() for f in futures]
        for ready in submitted:  # every callback has submitted its pieces
            ready.wait()
        later = [(t, f.result()) for t, f in list(follow)]
    finally:
        ex.shutdown()
    res = {}
    for (n, a), r in zip(tasks, done):
        res.setdefault(n, []).append((a, r))
    for (n, a), r in later:
        res.setdefault(n, []).append((a, r))
    return res


def reproduce(workers=None, quick=False):
    """Run R1-R11 (or the quick subset) in parallel worker processes."""
    return collect(start(quick, workers))


# Report ----------------------------------------------------------------------------------------

def _f(x, n=3):
    """Number in Portuguese notation (decimal comma)."""
    return mp.nstr(x, n).replace(".", ",")


def _pct(x):
    return _f(100 * mp.mpf(x), 3) + "%"


def _order(pairs):
    """Least-squares slope of log|dev| against log B."""
    xs = [mp.log(B) for B, _ in pairs]
    ys = [mp.log(abs(d)) for _, d in pairs]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)


def control_summaries(res):
    """control_summary of every positive-control case of a --full result."""
    known = dict(res.get("r8pk", []))
    pieces = {}
    for _, r in res.get("r8p", []):
        pieces.setdefault(r["spec"], []).append(r)
    return [control_summary(spec, known[spec], pieces.get(spec, [])) for spec in CONTROL_CASES if spec in known]


def _control_lines(summaries):
    """The positive control, as report lines, and whether it passed (every copy of every known
    solution in an examined undecided box of the minimum width, none in an excluded box)."""
    if not summaries:
        return ["    - **Controle positivo:** não rodado nesta execução; sem ele a frase de prova não é usada."], False
    from .provenance import CERT_ROOT_MODEL
    passed = all(x["passed"] for x in summaries)
    refined = all(x["complete_near_solutions"] for x in summaries)
    def pct(v):  # never rounded up to 100 %: as many digits as it takes to show what is missing
        if v >= 1:
            return "100 %"
        d = 3
        while f"{100 * v:.{d}g}" == "100" and d < 15:
            d += 1
        return f"{100 * v:.{d}g} %".replace(".", ",")
    out = ["    - **Controle positivo** (a mesma máquina, com o mesmo código, a mesma precisão e as mesmas hipóteses, só com os parâmetros do caso; o espaço de busca de cada caso é dividido em pedaços de s; critério: nenhuma caixa que contenha uma solução conhecida, ou uma cópia periódica dela, pode ser excluída):",
           "",
           "      | caso | pedaços de s examinados até o fim | com solução: examinados até o fim | volume examinado | caixas examinadas | excluídas | indecididas (largura mínima) | não examinadas (orçamento) | soluções conhecidas | em caixa excluída | Krawczyk: raiz única | Krawczyk nas caixas sem solução conhecida: sem raiz / indecidido / raiz sem solução conhecida |",
           "      |---|---|---|---|---|---|---|---|---|---|---|---|"]
    for x in summaries:
        _, _, label = control_case(x["spec"])
        out.append(f"      | {label} | {x['complete_pieces']} de {x['pieces']} | {x['solution_pieces_complete']} de {x['solution_pieces']} | {pct(x['examined_fraction'])} | {x['boxes']} | {x['excluded']} | {x['undecided']} | {x['unexamined']} | {x['known']} | {x['in_excluded']} | {x['certified']} de {x['known']} | {x['none']} / {x['kraw_undecided']} / {x['unmatched_roots']} |")
    out += ["", "      Onde estão as soluções conhecidas (cada solução contada uma vez):", "",
            "      | caso | solução (s; ν; mira) | pedaço de s | cópias periódicas na caixa | em caixa excluída | caixa indecidida (maior lado em ν, mira) | Krawczyk |",
            "      |---|---|---|---|---|---|---|"]
    for x in summaries:
        _, _, label = control_case(x["spec"])
        for y in x["sols"]:
            sol = "; ".join(_f(mp.mpf(t), 8) for t in y["sol"])
            piece = f"[{_f(mp.mpf(y['piece'][0]), 4)}; {_f(mp.mpf(y['piece'][1]), 4)}]" if y["piece"] else "nenhum"
            box = (_f(mp.mpf(y["width"]), 2) if y["refined"] else ("não examinada (orçamento)" if y["unexamined"] else "nenhuma")) if y["copies"] else "nenhuma"
            kr = f"raiz única em cada cópia, {CERT_ROOT_MODEL} (CERTIFIED)" if y["certified"] else ("; ".join(y["verdicts"]) or "não testada")
            out.append(f"      | {label} | {sol} | {piece} | {y['copies']} | {y['in_excluded']} | {box} | {kr} |")
    out.append("")
    if not passed:
        why = []
        if any(x["known"] == 0 for x in summaries):
            why.append("caso sem solução conhecida: " + "; ".join(control_case(x["spec"])[2] for x in summaries if x["known"] == 0))
        if sum(x["in_excluded"] for x in summaries):
            why.append(f"{sum(x['in_excluded'] for x in summaries)} cópia(s) de solução conhecida em caixa excluída")
        if sum(x.get("lost", 0) for x in summaries):
            why.append(f"{sum(x['lost'] for x in summaries)} cópia(s) de solução conhecida fora de toda caixa")
        if sum(x.get("none_in_solution_boxes", 0) for x in summaries):
            why.append(f"Krawczyk diz “sem raiz” em {sum(x['none_in_solution_boxes'] for x in summaries)} caixa(s) que contêm solução conhecida")
        out.append(f"      - **O controle falhou** ({'; '.join(why) or 'solução conhecida sem caixa'}). A prova do caso vazio é inválida, e a frase de prova não é usada.")
    elif not refined:
        out.append("      - Nenhuma solução conhecida caiu em caixa excluída, mas nem todas foram examinadas até a largura mínima. O controle está incompleto, e a frase de prova não é usada.")
    else:
        out.append(f"      - **O controle passou:** nenhuma das {sum(x['known'] for x in summaries)} soluções conhecidas ({sum(x['copies'] for x in summaries)} pontos, contando as cópias periódicas) caiu em caixa excluída; todas ficaram em caixas indecididas de largura mínima, e o Krawczyk provou raiz única em {sum(x['certified'] for x in summaries)} delas, em todas as cópias.")
    partial = [control_case(x["spec"])[2] for x in summaries if x["solution_pieces_complete"] == x["solution_pieces"] and x["complete_pieces"] < x["pieces"]]
    if partial:
        out.append(f"      - **Só os pedaços com solução foram examinados até o fim** em: {'; '.join(partial)}. Nos outros pedaços desses casos, o orçamento acabou antes, e o que sobrou está na coluna “não examinadas”; ali nada é afirmado, nem exclusão nem raiz.")
    from . import certify as ct
    (n0, n1), (a0, a1) = ct.default_ranges()
    out.append(f"      - **Cópias periódicas:** a caixa cobre ν ∈ [{_f(mp.mpf(n0), 4)}; {_f(mp.mpf(n1), 6)}] e mira ∈ [{_f(mp.mpf(a0), 4)}; {_f(mp.mpf(a1), 4)}], além de [0; 2π) e [−π; π). Uma solução com ν perto de 0 ou de 2π, ou com mira perto de ±π, aparece duas vezes na caixa (quatro, se as duas coisas valem); cada cópia é vigiada e testada, mas a solução é contada uma vez.")
    out.append(f"      - **Rótulo:** cada solução certificada recebe a origem CERTIFIED com o rótulo “{CERT_ROOT_MODEL}”. É raiz de (miss, C, V) com as três desigualdades válidas na caixa, e não “solução física certificada”: a ausência de colisões extras não foi provada em intervalos, e vem do verificador (COMPUTED).")
    out.append("      - **Origem de s:** no R8, s = 0 é o ponto em que a linha inicial entra na boca 2. Na Fig. 7(b), s é a altura acima da linha dos centros: s = 0 fica a meia distância das bocas. A solução da Fig. 7(b) com |s| < R é a de troca de velocidades da Fig. 6(b) do artigo (p. 1083), longe das bocas (`docs/phase3_borda.md`).")
    out.append("      - **Caixas examinadas** inclui as que foram divididas ao meio; cada uma delas termina excluída, dividida, indecidida na largura mínima ou não examinada.")
    out.append("      - **Krawczyk:** foi aplicado só às caixas indecididas que a máquina examinou. As caixas que ficaram na pilha por orçamento não foram testadas.")
    out.append("      - **Orçamento:** pedaços sem solução conhecida param em `ekt_cert_control_budget` caixas, e pedaços com solução em `ekt_cert_control_sol_budget`.")
    return out, passed and refined


def _certification_lines(arg, ref, chunks, control=None):
    """Report lines of the interval exclusion of the empty case (certify.py). The proof sentence is
    used only when the exclusion is complete and the positive control passed."""
    from . import certify as ct
    out = []
    control_out, control_ok = _control_lines(control or [])
    ci = ref.get("cert_items")
    if ci:
        old = ci["items"] - ci["items_aim"]
        out.append(f"  - **Certificação intervalar** (`certify.py`, mpmath.iv a {ct.CERT_PREC.value} bits, arredondamento para fora; enunciado e hipóteses em `docs/phase3_certification.md`).")
        out.append(f"    - **Itens não resolvidos depois do refinamento:** {ci['items']}, sendo {old} do mesmo tipo dos que já eram contados e {ci['items_aim']} colchetes da busca da mira descartados (registrados a partir desta versão). "
                   f"Resultado: {ci['excluded']} excluídos ({ci['excluded'] - ci['excluded_aim']} dos {old}; {ci['excluded_aim']} dos {ci['items_aim']}), {ci['undecided']} indecididos; o teste de Krawczyk não foi aplicado a esses itens.")
        out.append(f"    - Caixa de cada item: o intervalo em s que ele pode estragar (o da medida conservadora, arredondado para fora) × ν e mira inteiros. {ci['runs']} caixas distintas, {ci['boxes']} subcaixas examinadas.")
    if chunks:
        rs = sorted((r for _, r in chunks), key=lambda r: r["i"])
        n = rs[0]["n"]
        lo, hi = rs[0]["range"]
        # the claim needs every piece, of this case, and pieces that tile [lo, hi] without gaps
        whole = ([r["i"] for r in rs] == list(range(n)) and all(r["n"] == n and r["range"] == (lo, hi) and tuple(r["case"]) == tuple(arg) for r in rs)
                 and rs[0]["s"][0] == lo and rs[-1]["s"][1] == hi and all(a["s"][1] == b["s"][0] for a, b in zip(rs, rs[1:])))
        boxes = sum(r["boxes"] for r in rs)
        left = sum(r["n_undecided"] for r in rs)
        by = {}
        for r in rs:
            for k, m in r["excluded_by"].items():
                by[k] = by.get(k, 0) + m
        q = ct.exclusion_quantity({"complete": whole and left == 0, "boxes": boxes, "undecided": [None] * left},
                                  tuple(arg), (lo, hi), "caso vazio, Ω inteiro com uma travessia")
        hyp = ct.HYPOTHESES.format(s=f"[{_f(mp.mpf(lo), 6)}; {_f(mp.mpf(hi), 6)}]", prec=ct.CERT_PREC.value)
        reasons = "; ".join(f"{k}: {m}" for k, m in sorted(by.items(), key=lambda t: -t[1]))
        if not whole:
            out.append(f"    - **Ω inteiro, com uma travessia:** os pedaços da exclusão não cobrem o intervalo de s deste caso ({len(rs)} de {n}, ou de outra condição inicial). Nada é afirmado.")
        elif q.origin is Origin.CERTIFIED and not control_ok:
            out.append(f"    - **Ω inteiro, com uma travessia:** {boxes} subcaixas em {n} pedaços de s, todas excluídas ({reasons}). A frase de prova não é usada, porque o controle positivo abaixo não passou ou está incompleto.")
            out += control_out
        elif q.origin is Origin.CERTIFIED:
            out.append(f"    - **Ω inteiro, com uma travessia:** {boxes} subcaixas em {n} pedaços de s, todas excluídas ({reasons}).")
            out.append(f"    - **Provado, em aritmética intervalar, que não existe solução de uma colisão para esta condição inicial no modelo de [{hyp}]; soluções com várias colisões não foram examinadas.** "
                       f"Origem: {q.origin.name}, rótulo “{q.label}”; soluções na caixa: {q.value}.")
            out += control_out
            out.append("    - **Não coberto pela prova:** duas ou mais travessias entre as visitas, colisão depois de uma travessia, e posições de colisão fora desse intervalo de s. Para esses vale o que a busca por bracketing diz acima e nas variações abaixo (COMPUTED).")
            out.append("    - Nenhuma caixa ficou indecidida, então nenhuma raiz precisou do teste de Krawczyk neste caso; o Krawczyk aparece no controle positivo acima.")
        else:
            out.append(f"    - **Ω inteiro, com uma travessia:** {boxes} subcaixas em {n} pedaços de s; {left} indecididas ({q.origin.name}). A exclusão não está completa, e nada é afirmado além da contagem.")
    return out


def report_section(res):
    """Markdown lines for outputs/report.md. Every number here is COMPUTED."""
    L = []
    add = L.append
    add("")
    add("## Fase 3 · Bolas de bilhar (Echeverria, Klinkhammer & Thorne, 1991): reprodução R1–R11")
    add("")
    add("Todas as grandezas desta seção são **COMPUTED**, exceto as marcadas como CERTIFIED: a exclusão intervalar do caso vazio e as raízes do controle positivo certificadas pelo Krawczyk (R8, caso vazio), cada uma com o seu rótulo.")
    add("")
    add("- **Modelo independente** (`chronosim/billiards.py`): o modelo exato que o artigo descreve em palavras e figuras.")
    add("  - Bocas de raio B identificadas por reflexão no plano x = ½ (Fig. 2).")
    add("  - Defasagem T_d: quem entra na boca 2 sai da boca 1 em T − 1.")
    add("  - Travessia em linha reta vista da garganta.")
    add("  - Colisão elástica de esferas iguais, na forma geral.")
    add("  - Não usa nenhuma equação do artigo.")
    add("- **Equações do artigo** (`chronosim/ekt.py`): transcritas duas vezes; a segunda passada foi cega. Registro em `docs/equation_transcription.md`.")
    add("- **Unidades:** D = T_d = 1.")
    add("- **Busca:** bracketing unidimensional aninhado, em grades determinísticas, com polimento em mpmath (40 dígitos).")
    add("  - Cada solução passa por um verificador que reconstrói a história inteira e confere as leis do modelo, incluindo a ausência de qualquer outra sobreposição entre encarnações.")
    add(f"  - \"Intervalos não resolvidos\" são tudo o que o bracketing não conseguiu seguir: trechos de s ou de ν em que o número de raízes mudou e não pôde ser ligado após {bl.MAX_SPLIT.value} bisseções, bisseções (em s ou em ν) que perderam o ramo, e colchetes da busca da mira descartados (f indefinida dentro da bisseção, ou troca de sinal rejeitada como salto). Eles são contados, nunca descartados em silêncio, com uma exceção: intervalos da grade da mira com uma ponta indefinida (borda do domínio do laço) não são colchetes e não são registrados.")
    add("  - Cruzamentos rejeitados (v₁(s) passando pelo alvo num polo, ou impulso normal ≤ 0) e candidatos recusados pelo verificador entram como \"candidatos rejeitados\", com o motivo.")
    add("")
    add("**O que estas contagens são.** São soluções com **uma colisão** e o itinerário indicado, encontradas sob as grades `ekt_*` (HEURISTIC). Não são a multiplicidade do problema, e não estabelecem nenhuma afirmação universal.")
    add("")

    r1 = res["r1"][0][1]
    add(f"- **R1, eq. 2.1b.** O modelo de identificação reproduz ψ_out = ψ_in − 2φ com desvio máximo de {_f(r1['max_dev'].value)} rad numa grade de 25 × 25 em (φ, ψ_in). É uma derivação da regra, usada como checagem da transcrição.")
    r2 = res["r2"][0][1]
    c = r2["counts"]
    add(f"- **R2, eqs. 2.2–2.4.** Foram testadas 144 colisões elásticas (pares v₁, v₂), impondo a troca de módulos das eqs. 2.2a–b.")
    add(f"  - Toda normal de contato que satisfaz a condição é uma das duas da eq. 2.3 ou 2.4: troca de velocidades {c['velocity']} vezes, troca espelhada {c['mirror']} vezes, outra {c['other']} vezes.")
    add(f"  - A direção da linha dos centros desvia das eqs. 2.3b/2.4b no máximo {_f(r2['max_dev_direction'].value)}.")
    add("")

    add("### R3, R4–R6: a trajetória simétrica das Figs. 6 e 7")
    add("")
    add("A trajetória é perpendicular à linha dos centros e passa pelo ponto médio, com v₁ = 0,5.")
    add("")
    for _, r in res["r3"]:
        dev = r["dev"].value
        add(f"- **R3, Fig. 6b, troca de velocidades** (B = 0,01, R = {r['R']}): {r['n']} solução. A rapidez entre colisões menos (1 − 2B) dá {', '.join(_f(d) for d in dev)}. Houve {r['unresolved']} intervalos não resolvidos.")
    ratios = [r["dev"].value[0] / mp.mpf(R) for R, r in res["r3"] if r["dev"].value]
    add(f"  - Desvio / R = {', '.join(_f(x, 4) for x in ratios)}: o desvio é proporcional a R, e v₁′ = 1 − 2B (p. 1083) é o limite R → 0.")
    rest = res["rest"][0][1]
    add(f"- **Limite v₁ → 0 (p. 1084).** Com a bola em repouso sobre a linha dos centros em x = {'; '.join(r['x'].replace('.', ',') for r in rest)}, o verificador aceita a solução: {'; '.join('sim' if r['ok'] else 'não' for r in rest)}.")
    add(f"  - A rapidez do laço é {_f(rest[0]['speed'].value, 8)}; menos (1 − 2B − 2R), dá {_f(rest[0]['speed'].value - (1 - 2 * mp.mpf('0.01') - 2 * mp.mpf('1e-6')), 3)}. O texto dá v₂ = 1, no limite B ≪ 1.")
    add("  - Nada na história depende do instante da colisão, que é o contínuo de soluções de que fala o artigo.")
    add("")
    add("| B | soluções com 1, 2 e 3 travessias | tan ψ / v₁ − 1 (eq. 3.1) | (sin ψ + tan ψ)/(2v₁) − 1 (eq. 3.2) | instantes da p. 1084, maior desvio | intervalos não resolvidos (1/2/3) |")
    add("|---|---|---|---|---|---|")
    for B, r in res["r46"]:
        ys = "; ".join(f"y = {_f(r[k]['y'][0], 9)}" if r[k]["n"] == 1 else f"{r[k]['n']} soluções" for k in (1, 2, 3))
        add(f"| {B.replace('.', ',')} | {ys} | {_f(r['dev31'].value, 4) if 'dev31' in r else '—'} | {_f(r['dev32'].value, 4) if 'dev32' in r else '—'} | "
            f"{_f(r['times']['max_dev'].value, 3) if 'times' in r else '—'} | {r[1]['unresolved']}/{r[2]['unresolved']}/{r[3]['unresolved']} |")
    add("")
    r46 = [r for _, r in res["r46"]]
    o31 = _order([(mp.mpf(r["B"]), r["dev31"].value) for r in r46 if "dev31" in r])
    o32 = _order([(mp.mpf(r["B"]), r["dev32"].value) for r in r46 if "dev32" in r])
    add(f"ψ é medido na geometria: tan ψ = |y|/(1/2), onde y é a altura da colisão. Ajustando log|desvio| contra log B, a ordem medida é {_f(o31, 3)} para a eq. 3.1 e {_f(o32, 3)} para a eq. 3.2. As duas são, portanto, o limite B → 0 da geometria. Para a Fig. 7d (três travessias), o artigo não imprime valor (\"exercise for the reader\"), e a coluna dá o y medido.")
    add("")

    add("### R7: muitas travessias (seção III B)")
    add("")
    for key, label in (("r7", "v₁ = 3, ψ_A = 0,6, h = 0,5"), ("r7b", "v₁ = 2, ψ_A = 0,6, h = 1")):
        rows = [r for _, r in res[key]]
        nmin = rows[0]["n_min_38"]
        add(f"**{label}**, B = 0,05, R = 10⁻⁴, colisão a até 50 D. O n_min da eq. 3.8 é {_f(nmin, 4)}.")
        add("")
        add("| n | L medido | n da eq. 3.7 no L medido | intervalos não resolvidos | fração do espaço de busca, conservadora (fina) |")
        add("|---|---|---|---|---|")
        for r in rows:
            row = r["rows"][0]
            fr = row["fraction"]
            add(f"| {row['n']} | {'; '.join(_f(x, 7) for x in row['L']) or 'nenhuma solução encontrada'} | {'; '.join(_f(x, 5) for x in row['n_37']) or '—'} | {row['unresolved']} | {_pct(fr['conservative'])} ({_pct(fr['fine'])}) |")
        add("")
    first = [r["rows"][0] for _, r in res["r7"]]
    offs = [row["n_37"][0] - row["n"] for row in first if row["n_37"]]
    found_n = [row["n"] for row in first if row["L"]]
    if offs:
        add(f"- **Primeiro conjunto:** soluções encontradas para n = {', '.join(map(str, found_n))}. A eq. 3.7, avaliada no L medido, dá n mais {_f(min(offs), 3)} a {_f(max(offs), 3)} (n ≫ 1, L ≫ 1 é o regime dela). O menor n encontrado é {min(found_n)}; a eq. 3.8, também assintótica, dá {_f(res['r7'][0][1]['n_min_38'], 4)}.")
    else:
        add("- **Primeiro conjunto:** nenhuma solução encontrada sob este critério.")
    second = [r["rows"][0] for _, r in res["r7b"]]
    add(f"- **Segundo conjunto:** soluções encontradas só para n = {', '.join(str(row['n']) for row in second if row['L']) or 'nenhum'}, com {sum(row['unresolved'] for row in second)} intervalos não resolvidos no total, sem refinamento; a fração do espaço de busca que cobrem está na tabela, pelas medidas definidas em R8 (caso vazio). A busca é **inconclusiva** aqui: a ausência das outras soluções não é afirmada.")
    add("")

    add("### R8: árbitro das eqs. 4.2–4.3 contra a geometria independente")
    add("")
    add("Para cada solução geométrica verificada (colisão antes de uma travessia), α, β e d são medidos na geometria e postos nas eqs. 4.2 e 4.3 como impressas. Em separado, as eqs. 4.2–4.3 são resolvidas diretamente, por findroot a partir de uma grade de pontos de partida, e suas raízes são confrontadas com as geométricas.")
    add("")
    add("| (θ, φ, λ₁, B, R) | movimento livre consistente? | soluções geométricas | raízes físicas das eqs. 4.2–4.3 | coincidem | raízes com distância negativa | maior resíduo das eqs. 4.2–4.3 | maior desvio de d (4.4b) e de v₂ (4.5) | sign(β) = sign(d) (4.6c) | maior desvio de a (A6) contra PQ medido |")
    add("|---|---|---|---|---|---|---|---|---|---|")
    for arg, r in res["r8"]:
        rows = r["rows"]
        matched = sum(1 for (a, b, sg) in r["physical"] if any(abs(a - x["alpha"]) < mp.mpf("1e-12") and abs(b - x["beta"]) < mp.mpf("1e-12") for x in rows))
        res_max = max((x["res"].value for x in rows), default=None)
        dmax = max((max(x["dev_d"].value, x["dev_v2"].value) for x in rows), default=None)
        classes = f" ({', '.join(x['class'] for x in rows)})" if rows else ""
        a6max = max((x["dev_a6"].value for x in rows), default=None)
        signs = ("sim" if all(x["sign_rule"] for x in rows) else "não") if rows else "—"
        add(f"| ({'; '.join(a.replace('.', ',') for a in arg)}) | {'sim' if r['free_consistent'] else 'não'} | {len(rows)}{classes} | "
            f"{len(r['physical'])} | {matched}/{len(r['physical'])} | {len(r['negative'])} | {_f(res_max) if res_max is not None else '—'} | {_f(dmax) if dmax is not None else '—'} | "
            f"{signs} | {_f(a6max) if a6max is not None else '—'} |")
    add("")
    add("As raízes das eqs. 4.2–4.3 são procuradas só no domínio 0 < α < π, |β| < π/2.")
    add("")
    with_sols = [(arg, r) for arg, r in res["r8"] if r["rows"]]
    all_match = all(
        len(r["rows"]) == len(r["physical"])
        and all(any(abs(a - x["alpha"]) < mp.mpf("1e-12") and abs(b - x["beta"]) < mp.mpf("1e-12") for x in r["rows"]) for (a, b, sg) in r["physical"])
        for _, r in with_sols)
    rows_all = [x for _, r in with_sols for x in r["rows"]]
    worst = max([x["res"].value for x in rows_all] + [x["dev_d"].value for x in rows_all] + [x["dev_v2"].value for x in rows_all] + [x["dev_a6"].value for x in rows_all], default=None)
    signs_ok = all(x["sign_rule"] for x in rows_all)
    add(f"- **Árbitro.** Nos {len(with_sols)} casos com solução geométrica, as raízes físicas das eqs. 4.2–4.3 e as soluções geométricas "
        + ("são as mesmas, uma a uma." if all_match else "**não** coincidem em todos os casos (ver a tabela)."))
    if worst is not None:
        add(f"  - O maior resíduo ou desvio entre resíduo das eqs. 4.2–4.3, d (4.4b), v₂ (4.5) e a (A6) é {_f(worst)}.")
    if rows_all:
        add(f"  - A regra de sinal sign(β) = sign(d) da eq. 4.6c vale em {'todas' if signs_ok else 'nem todas'} as {len(rows_all)} soluções. A eq. A6 é comparada com a distância PQ medida em cada uma.")
    if all_match and signs_ok and worst is not None and worst < mp.mpf("1e-30"):
        add("  - Isso confere a transcrição das eqs. 4.2, 4.3, 4.4a, 4.4b, 4.5 e A6, e o modelo da travessia, nestes casos.")
    negs = [(arg, r["negative"]) for arg, r in res["r8"] if r["negative"]]
    where = "; ".join(f"(θ, φ) = ({arg[0].replace('.', ',')}; {arg[1].replace('.', ',')}): {len(n)}" for arg, n in negs) or "nenhum caso"
    add(f"- **Raízes que a eq. 4.6 não lista,** com alguma das distâncias a, a′ ou c do Apêndice A negativa: {where}.")
    add("  - A geometria não a realiza, porque nenhuma bola percorre distância negativa ao longo do caminho.")
    add("  - A lista da eq. 4.6 cobre as raízes introduzidas pelas manipulações do Apêndice A (p. 1097). Esta tem outra origem.")
    empty = [(arg, r) for arg, r in res["r8"] if not r["rows"]]
    kind_label = {"s_count_change": "s, contagem mudou", "s_bisect_lost": "s, bisseção perdeu o ramo", "nu_count_change": "ν, contagem mudou",
                  "nu_state_lost": "ν, estado perdido", "nu_bisect_lost": "ν, bisseção perdeu o ramo",
                  "aim_bracket_lost": "mira, colchete perdido (f indefinida na bisseção)", "aim_bracket_jump": "mira, colchete rejeitado como salto"}

    def kinds(kd):
        return "; ".join(f"{kind_label.get(k, k)}: {n}" for k, n in sorted(kd.items())) or "nenhum"

    pct = _pct

    def cover(fr):
        return f"cobrindo {pct(fr['conservative'])} do espaço de busca pela medida conservadora ({pct(fr['fine'])} pela medida fina)"

    for arg, r in empty:
        lo, hi = r["s_range"]
        ref = next((x for _, x in res.get("r8r", []) if tuple(x["case"]) == tuple(arg)), None)
        n_left, fr = (ref["after"], ref["fraction_after"]) if ref else (r["unresolved"], r["fraction"])
        add(f"- **Caso ({'; '.join(a.replace('.', ',') for a in arg)}):** a trajetória é perigosa, porque o movimento livre se sobrepõe a si mesmo (distância mínima {_f(r['free_min_distance'])} < 2R). "
            f"Com uma colisão antes da travessia, 1 travessia entre as visitas e as grades `ekt_*` padrão: **nenhuma solução encontrada na região resolvida; {n_left} intervalos não resolvidos, {cover(fr)}.**")
        add("  - **Espaço de busca:** Ω = [s_min; s_max] × [0; 2π) em (posição da colisão s, ângulo ν da normal de contato); a mira é resolvida dentro de cada ponto de Ω.")
        add("  - **Medida fina:** soma das células (s × largura em ν) dos intervalos não resolvidos. Intervalos em ν registrados fora da grade de s (dentro de uma bisseção em s) não têm célula e não entram nela.")
        add("  - **Medida conservadora:** fração de s coberta pela união dos intervalos da grade que cada item não resolvido pode estragar, com ν inteiro. Um intervalo em ν não resolvido num ponto da grade pode esconder um ramo de raízes e estragar a ligação de v₁(s) nos dois intervalos vizinhos; um item fora da grade conta com o intervalo que o contém.")
        add("  - A medida conservadora é cota superior da região em que uma solução poderia ter escapado **só sob o critério de bracketing e na grade declarada**. Ela não cobre uma troca de sinal dupla entre dois nós. Os colchetes da busca da mira descartados (f indefinida na bisseção, ou salto) são registrados e contados como itens no ponto (s, ν) em que ocorreram. Não são registrados os intervalos da grade da mira com uma ponta indefinida (a borda do domínio do laço, along ≤ 0): uma raiz entre o último nó definido e essa borda não aparece em contagem nenhuma.")
        if ref:
            reasons = "; ".join(f"{n}: {why}" for why, n in sorted(ref["reasons_after"].items(), key=lambda t: -t[1]))
            add(f"  - **Refinamento adaptativo** (até `ekt_refine_max_split` = {ref['depth']} halvings a mais em cada intervalo): de {ref['before']} intervalos não resolvidos ({kinds(ref['before_kinds'])}), "
                f"cobrindo {pct(ref['fraction_before']['conservative'])} pela medida conservadora, sobram {ref['after']} ({kinds(ref['after_kinds'])}); {ref['fraction_after']['off_grid']} deles estão fora da grade de s e não entram na medida fina. Motivos dos que sobram: {reasons}.")
            add("  - Critério de resolvido: num trecho em que a contagem de raízes muda, a grandeza que decide a solução (a troca normal C, em ν; v₁ menos o alvo, em s) tem um único sinal em todas as raízes das duas pontas. É o critério de bracketing da busca principal, não uma prova.")
            add(f"  - Intervalos com bisseção perdida não são refinados e continuam contados. {ref['ok']} soluções verificadas nesta busca; {len(ref['rejected'])} candidatos rejeitados.")
            for line in _certification_lines(arg, ref, res.get("r8c", []), control_summaries(res)):
                add(line)
        add(f"  - Posição da colisão s (distância com sinal ao longo da linha, a partir da entrada na boca 2; s < 0 antes dela) em [{_f(lo, 3)}; {_f(hi, 3)}] antes da travessia, e em [{_f(mp.mpf(S_NEAR_AFTER.value), 3)}; {_f(S_FAR.value, 3)}] depois dela. A grade de s tem `ekt_s_grid` pontos uniformes mais aproximações geométricas às pontas; ν e a mira usam `ekt_nu_grid` e `ekt_theta_grid`.")
        add(f"  - Mesma busca, antes do refinamento: {len(r['physical'])} raízes físicas das eqs. 4.2–4.3, "
            f"{len(r['rejected'])} candidatos geométricos rejeitados{(' (' + '; '.join(sorted(set(f[0] for f in r['rejected'] if f))) + ')') if r['rejected'] else ''}, "
            f"{r['unresolved']} intervalos não resolvidos ({kinds(r['unresolved_kinds'])}), {cover(r['fraction'])}.")
        add("  - Variações, sem refinamento:")
        for _, x in res["r8x"]:
            if tuple(x["case"]) != tuple(arg):
                continue
            kind = "Colisão antes da travessia" if x["kind"] == "before" else "Colisão depois de uma travessia"
            trav = "1 travessia" if x["k"] == 1 else f"{x['k']} travessias"
            extra = (", grades dobradas (ν, mira e s)" if x["refined"] else "") + (f", λ₁ = {x['lam1'].replace('.', ',')}" if x["lam1"] != arg[2] else "")
            add(f"    - {kind}, {trav} entre as visitas{extra}: {x['ok']} soluções, {len(x['rejected'])} candidatos rejeitados, "
                f"{x['unresolved']} intervalos não resolvidos ({kinds(x['unresolved_kinds'])}), {cover(x['fraction'])}.")
        add("  - Isto **não** é declarado paradoxo nem multiplicidade zero. Soluções com várias colisões, que o artigo menciona (p. 1080), estão fora do modelo buscado; a grade é finita; e os intervalos não resolvidos acima são lugares onde uma solução poderia ter escapado. É o primeiro alvo do mapa (pergunta Q1).")
    add("")

    add("### R9: região central do triângulo, eqs. 4.12–4.13 e as duas formas da eq. B20")
    add("")
    add("Parâmetros: (θ, φ) = (0,2; 0,6), λ₁ = ½, R = 0,1 B.")
    add("")
    add("| B | classe | β / (4.12a ou 4.13a) − 1 | ε / (4.12b ou 4.13b) − 1 | (λ₂ − σ) / (4.12c ou 4.13c) − 1 | (λ₂ − σ) / (B20, 1ª forma) − 1 |")
    add("|---|---|---|---|---|---|")
    for B, rows in res["r9"]:
        for r in rows:
            add(f"| {B.replace('.', ',')} | {r['class']} | {_f(r['dev_beta'].value, 4)} | {_f(r['dev_eps'].value, 4)} | {_f(r['dev_lam2'].value, 4)} | {_f(r['dev_b20_first'].value, 4)} |")
    add("")
    r9 = [(mp.mpf(B), r) for B, rows in res["r9"] for r in rows if r["class"] == "I"]
    ob, oe = _order([(B, r["dev_beta"].value) for B, r in r9]), _order([(B, r["dev_eps"].value) for B, r in r9])
    ol2, ol1 = _order([(B, r["dev_lam2"].value) for B, r in r9]), _order([(B, r["dev_b20_first"].value) for B, r in r9])
    add(f"- **β e ε (classe I):** os desvios das eqs. 4.12a–b caem com ordem {_f(ob, 3)} e {_f(oe, 3)} em B. As formas analíticas são o limite B → 0 da geometria.")
    add(f"- **λ₂ − σ:** com a segunda forma da eq. B20 (igual à 4.12c), o desvio cai com ordem {_f(ol2, 3)}, então o termo de ordem B confere. Com a primeira forma impressa, cai com ordem {_f(ol1, 3)}, então o coeficiente de ordem B dela não confere com a geometria.")
    add("- **As duas formas impressas da eq. B20 diferem por um fator 4 no termo em B** (ver `docs/equation_transcription.md`). A geometria segue a segunda forma.")
    add("")

    add("### R10: borda superior θ = φ (Fig. 14), B = 10⁻², R = 10⁻⁵, λ₁ = ½")
    add("")
    add("| φ | β classe I (geometria) | eq. B32 | eq. 4.14 | eq. 4.15 | classe II (geometria) |")
    add("|---|---|---|---|---|---|")
    for ph, r in res["r10"]:
        add(f"| {ph.replace('.', ',')} | {'; '.join(_f(b, 5) for b in r['beta_I']) or 'nenhuma'} | {_f(r['b32'], 5) if r['b32'] is not None else '—'} | {_f(r['b414'], 4)} | {_f(r['b415'], 4)} | {'; '.join(_f(b, 4) for b in r['beta_II']) or 'nenhuma'} |")
    add("")
    near = [abs(r["beta_I"][0] / r["b32"] - 1) for ph, r in res["r10"] if r["beta_I"] and r["b32"] and mp.mpf("0.7") <= mp.mpf(ph) <= mp.mpf("0.9")]
    add(f"- **Classe I.** Para φ entre 0,7 e 0,9, a geometria acompanha a eq. B32 com desvio relativo máximo de {_f(max(near), 2) if near else '— (sem solução nesse trecho)'}. Longe de π/4 ela se aproxima das eqs. 4.14 e 4.15, cada uma do seu lado. Tudo isso num único B = 10⁻²: nenhum limite B → 0 foi medido para B32, 4.14 ou 4.15.")
    add("  - Perto de φ → 0, β volta a crescer, e a eq. B32 deixa de acompanhar. O artigo atribui isso à quebra da expansão do Apêndice B perto do canto (p. 1092), e é onde a Fig. 14 mostra a parte numérica.")
    add("  - A forma da curva confere, a olho, com a Fig. 14: mínimo perto de φ = 0,1, máximo de cerca de 2,6 × 10⁻³ perto de φ = 0,5, e cerca de 2 × 10⁻⁵ em φ → π/2. A comparação com a figura é visual; não há números impressos.")
    add("- **Classe II.** Com a colisão antes da travessia, a geometria não encontra solução de classe II em nenhum φ da borda.")
    for ph, x in res["r10x"]:
        c2 = [n for n in x["negative"] if n[2] == -1]
        add(f"  - φ = {ph.replace('.', ',')}: as eqs. 4.2–4.3 têm {len([p_ for p_ in x['physical'] if p_[2] == -1])} raízes físicas de classe II e {len(c2)} de classe II com distância negativa "
            f"(β = {'; '.join(_f(n[1], 4) for n in c2) or '—'}, a = {'; '.join(_f(n[3]['a'], 3) for n in c2) or '—'}). A eq. B29 com s = −1 dá β = {'; '.join(_f(b, 4) for b in x['b29_class2']) or 'nenhuma raiz real'}.")
    add("  - Onde a eq. B29 prevê classe II (φ > π/4), as raízes correspondentes das eqs. 4.2–4.3 têm a = PQ < 0: a colisão ficaria além da boca 2 ao longo do caminho, e a geometria não as realiza com a colisão antes da travessia.")
    add("  - O artigo suspeita (p. 1091, lacuna 16) que a classe II que falta exista com a colisão depois da travessia. R11 procura essa solução.")
    add("")

    add("### R11: classe II com a colisão depois de uma travessia (lacuna 16, p. 1091)")
    add("")
    add("| φ (borda superior) | soluções | β | posição após a saída da boca 1 | intervalos não resolvidos |")
    add("|---|---|---|---|---|")
    for ph, r in res["r11"]:
        add(f"| {ph.replace('.', ',')} | {r['n']} | {'; '.join(_f(b, 4) for b in r['beta']) or '—'} | {'; '.join(_f(x, 4) for x in r['s']) or '—'} | {r['unresolved']} |")
    add("")
    add("- O primeiro encontro com a colisão fica sobre a linha que sai livremente da boca 1, e há mais uma travessia entre as visitas.")
    found = [ph for ph, r in res["r11"] if r["n"]]
    empty_phi = [ph for ph, r in res["r11"] if not r["n"]]
    add(f"- Soluções verificadas aparecem em φ = {', '.join(p.replace('.', ',') for p in found) or 'nenhum'}. Em φ = {', '.join(p.replace('.', ',') for p in empty_phi) or 'nenhum'}, nenhuma solução foi encontrada sob este critério.")
    add("- Só foi buscado um tipo de \"colisão depois da travessia\": uma travessia antes e uma entre as visitas. O artigo não especifica o que tem em mente.")
    add("")

    add("### O que a reprodução estabelece e o que não estabelece")
    add("")
    add("- ✅ **Estabelece:**")
    add("  - As eqs. 3.1, 3.2 (Figs. 6 e 7) e 4.12–4.13 são o limite B → 0 de um modelo geométrico exato e independente: os desvios caem com ordem ≈ 1 em B (R4–R6, R9).")
    add("  - As eqs. 4.14, 4.15 e B32 (Fig. 14) foram comparadas só em B = 10⁻²: a classe I acompanha a B32 perto de φ = π/4 (R10). Isso não é um teste de limite.")
    add("  - Na seção III B (Fig. 9), há uma solução para cada n de 1 a 6 no primeiro conjunto de R7, com a eq. 3.7 conferindo só assintoticamente. O segundo conjunto é inconclusivo.")
    add("  - As eqs. 4.2–4.3 são satisfeitas pela geometria no arredondamento de 40 dígitos, nos casos de R8 com solução (ver a tabela de R8).")
    add("- ❌ **Não estabelece:**")
    add("  - Nenhuma afirmação universal. As contagens valem para as grades, os itinerários e as janelas declarados.")
    add("  - Nada sobre soluções com várias colisões.")
    add("  - **A fase é clássica.** Ela testa a alternativa de autoconsistência de Novikov, não o mecanismo quântico de Hawking.")
    add("  - O mapa sistemático e a certificação intervalar ainda não foram feitos.")
    add("  - O teste de escala que confirmaria a observação deste projeto, de que só v₁T_d/D importa (a observação não é do artigo), ainda não foi feito.")
    add("")
    return L


def report_section_quick(res):
    """Markdown lines for the default (quick) run: a small, representative subset of R1-R11. The
    full reproduction is written by ``python -m chronosim.run --full`` to outputs/phase3_full.md."""
    L = []
    add = L.append
    add("")
    add("## Fase 3 · Bolas de bilhar (Echeverria, Klinkhammer & Thorne, 1991): execução reduzida")
    add("")
    add("Esta é a execução padrão, com um subconjunto pequeno de R1–R11. A reprodução completa está em `outputs/phase3_full.md`, escrita pela última execução de `python -m chronosim.run --full`; os números de lá não foram recalculados nesta execução. Todas as grandezas são **COMPUTED**.")
    add("")
    r1 = res["r1"][0][1]
    add(f"- **R1, eq. 2.1b:** desvio máximo {_f(r1['max_dev'].value)} rad.")
    c = res["r2"][0][1]["counts"]
    add(f"- **R2, eqs. 2.2–2.4:** troca de velocidades {c['velocity']}, troca espelhada {c['mirror']}, outra {c['other']}.")
    rest = res["rest"][0][1]
    add(f"- **Bola em repouso (p. 1084):** verificador aceita em x = {'; '.join(r['x'].replace('.', ',') for r in rest)}: {'; '.join('sim' if r['ok'] else 'não' for r in rest)}.")
    r3 = res["r3"][0][1]
    add(f"- **R3, Fig. 6b** (B = 0,01, R = {r3['R']}): {r3['n']} solução; rapidez menos (1 − 2B): {', '.join(_f(d) for d in r3['dev'].value)}; {r3['unresolved']} intervalos não resolvidos.")
    arg, r = res["r8"][0]
    rows = r["rows"]
    worst = max([x["res"].value for x in rows] + [x["dev_a6"].value for x in rows], default=None)
    add(f"- **R8, caso central ({'; '.join(a.replace('.', ',') for a in arg)}):** {len(rows)} soluções geométricas ({', '.join(x['class'] for x in rows) or '—'}), "
        f"{len(r['physical'])} raízes físicas das eqs. 4.2–4.3; maior resíduo das eqs. 4.2–4.3 ou desvio da A6: {_f(worst) if worst is not None else '—'}; "
        f"sign(β) = sign(d) em {sum(1 for x in rows if x['sign_rule'])}/{len(rows)}; {r['unresolved']} intervalos não resolvidos.")
    ph, r10 = res["r10"][0]
    add(f"- **R10, borda superior, φ = {ph.replace('.', ',')}:** β classe I = {'; '.join(_f(b, 5) for b in r10['beta_I']) or 'nenhuma'}; eq. B32 = {_f(r10['b32'], 5) if r10['b32'] is not None else '—'} (B = 10⁻²).")
    ph, r11 = res["r11"][0]
    add(f"- **R11, classe II depois da travessia, φ = {ph.replace('.', ',')}:** {r11['n']} soluções; β = {'; '.join(_f(b, 4) for b in r11['beta']) or '—'}.")
    add("")
    add("**Não estabelece** nada além destes casos; o que a fase estabelece está em `outputs/phase3_full.md`.")
    add("")
    return L
