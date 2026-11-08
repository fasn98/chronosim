"""Phase 3: billiard ball in the eternal-time-machine wormhole spacetime, exact geometric model.

This module is independent of every equation of Echeverria, Klinkhammer & Thorne (1991): it uses
only the model they state in words and figures, with finite mouth radius B and ball radius R and
no small-parameter expansion. Their equations (2.1-2.4, 3.x, 4.x, Appendix B) are the targets the
results here are compared with, in tests and in run.py.

Model (EKT 1991, Sec. I B p. 1079, Sec. II A p. 1081, Fig. 2), coplanar, in units D = T_d = 1:

- Flat space. The mouths are circles of radius B centered at C1 = (0, 0) and C2 = (1, 0).
- Points of the two mouths that are mirror images in the plane x = 1/2 are identified (Fig. 2).
- A center that enters mouth 2 at time T leaves mouth 1 at time T - 1; one that enters mouth 1
  leaves mouth 2 at T + 1 (p. 1081: "enters mouth 2 and exits from mouth 1, thereby traveling
  backward in time").
- The center moves through the throat as a point would, on a straight line as seen on the throat
  with unchanged speed (p. 1081). Our reading of "straight line on the throat" with the
  identification of Fig. 2: the velocity is reflected in the tangent line at the entry point and
  then mapped by the mirror x -> 1 - x. That this reproduces their Eq. 2.1b is a test, not an input.
- The ball is a hard sphere of radius R; a self-collision is elastic between equal masses: the
  velocity components along the line of centers are exchanged, the others are kept.

Everything produced here is COMPUTED. Floats are used for the grid searches; mpmath for polishing
and for the verifier.
"""

from __future__ import annotations

import math

import mpmath as mp

from .provenance import heuristic

# Grid-search settings (refined until counts are stable; see docs/phase3_plan.md).
THETA_GRID = heuristic("ekt_theta_grid", 64, "grid points on the visible arc of a mouth in the inner (aim) search")
NU_GRID = heuristic("ekt_nu_grid", 48, "grid points of the contact-normal angle in the middle search")
POLISH_DPS = heuristic("ekt_polish_dps", 40, "working precision when a bracketed root is polished with findroot")
VERIFY_TOL = heuristic("ekt_verify_tol", mp.mpf("1e-25"), "tolerance of the verifier (distances, speeds, residuals) at polish precision")
MAX_SPLIT = heuristic("ekt_max_split", 12, "halvings of an s-interval whose root count changes before it is reported as unresolved")
JUMP_RATIO = heuristic("ekt_jump_ratio", 1e-6, "a bracketed sign change is a root only if the refined |f| falls below this fraction of the end values (rejects jumps)")
BISECT_REL_TOL = heuristic("ekt_bisect_rel_tol", 1e-13, "bisection stops at this fraction of the bracket width (floats); with mpmath, 10^-(dps-8)")
FLOAT_DEPTH = heuristic("ekt_float_depth", 11, "floats are used while k log10(2/B) is below this; mpmath above (itinerary windows ~ (B/2)^k)")
FUTURE_TRAVERSALS = heuristic("ekt_future_traversals", 6, "traversals followed after the collision when the verifier builds the history")

C = ((0, 0), (1, 0))  # mouth centers, index 0 = mouth 1, index 1 = mouth 2


class Backend:
    def __init__(self, mod):
        self.sqrt, self.sin, self.cos, self.atan2, self.pi = mod.sqrt, mod.sin, mod.cos, mod.atan2, mod.pi


FLOAT = Backend(math)
MP = Backend(mp)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def add(a, b, k=1):
    return (a[0] + k * b[0], a[1] + k * b[1])


def scale(a, k):
    return (a[0] * k, a[1] * k)


def norm(a, xm):
    return xm.sqrt(dot(a, a))


def unit(a, xm):
    n = norm(a, xm)
    return (a[0] / n, a[1] / n)


def first_entry(p, v, B, xm, skip=None):
    """Earliest t > 0 at which p + v t enters a mouth: (t, mouth index) or None.

    ``skip`` is the mouth just exited (its exit root is t = 0)."""
    best = None
    a = dot(v, v)
    for m in (0, 1):
        if m == skip:
            continue
        pc = (p[0] - C[m][0], p[1] - C[m][1])
        b = 2 * dot(v, pc)
        c = dot(pc, pc) - B * B
        if b >= 0:
            continue
        disc = b * b - 4 * a * c
        if disc < 0:
            continue
        t = (-b - xm.sqrt(disc)) / (2 * a)
        if t > 0 and (best is None or t < best[0]):
            best = (t, m)
    return best


def traverse(X, v, m, B):
    """Exit point, exit velocity and time shift for a center entering mouth m at X with velocity v."""
    n = ((X[0] - C[m][0]) / B, (X[1] - C[m][1]) / B)
    vr = add(v, n, -2 * dot(v, n))
    return (1 - X[0], X[1]), (-vr[0], vr[1]), (-1 if m == 1 else 1)


def follow(p, v, B, k, xm, skip=None):
    """Free flight from p with velocity v through exactly k traversals.

    Returns (legs, exit_point, exit_velocity, total_time_shift) where legs is the list of
    (start, velocity, duration, mouth entered) of the k segments that end at a mouth, or None if
    the ray stops entering mouths before k traversals."""
    legs, shift = [], 0
    for _ in range(k):
        hit = first_entry(p, v, B, xm, skip)
        if hit is None:
            return None
        t, m = hit
        X = add(p, v, t)
        legs.append((p, v, t, m))
        p, v, dt = traverse(X, v, m, B)
        shift += dt
        skip = 1 - m
    return legs, p, v, shift


def aim_point(P, m, angle, B, xm):
    """Point of mouth m at polar angle ``angle``, if a straight line from P enters m there first."""
    X = (C[m][0] + B * xm.cos(angle), C[m][1] + B * xm.sin(angle))
    d = add(X, P, -1)
    hit = first_entry(P, d, B, xm)
    if hit is None or hit[1] != m or abs(hit[0] - 1) > 1e-9 * (1 + abs(hit[0])):
        return None
    return X


def visible_arc(P, m, B, xm):
    """Polar angles (center, half-width) of the arc of mouth m seen from P."""
    rel = (P[0] - C[m][0], P[1] - C[m][1])
    r = norm(rel, xm)
    if r <= B:
        return None
    return xm.atan2(rel[1], rel[0]), mp.acos(B / r) if xm is MP else math.acos(B / r)


def loop(P, Q, m, angle, k, B, xm):
    """The young ball leaves P aiming at polar angle ``angle`` of mouth m, makes k traversals and
    heads for Q. Returns dict with the transverse miss at Q, the loop length, the net time shift
    and the unit directions, or None if the itinerary does not exist."""
    X = aim_point(P, m, angle, B, xm)
    if X is None:
        return None
    first = add(X, P, -1)
    L0 = norm(first, xm)
    dir0 = scale(first, 1 / L0)
    E, w, shift = traverse(X, dir0, m, B)
    length = L0
    if k > 1:
        rest = follow(E, w, B, k - 1, xm, skip=1 - m)
        if rest is None:
            return None
        legs, E, w, s2 = rest
        length += sum(t * norm(vv, xm) for _, vv, t, _ in legs)
        shift += s2
    w = unit(w, xm)
    to_q = add(Q, E, -1)
    along = dot(to_q, w)
    if along <= 0 or shift >= 0:
        return None
    return {"miss": cross(w, to_q), "length": length + along, "shift": shift, "dir0": dir0, "dir2": w, "exit": E}


def roots_by_bracketing(f, xs, tol):
    """Sign changes of f on the grid xs, refined by bisection to width ``tol``. Values of f may be
    None (undefined); a bracket needs both ends defined. A sign change across a jump is rejected:
    the refined end values must be small compared with the bracket's original end values."""
    vals = [f(x) for x in xs]
    out = [x for x, v in zip(xs, vals) if v is not None and v == 0]  # zeros on a node, each once
    for i in range(len(xs) - 1):
        fa, fb = vals[i], vals[i + 1]
        if fa is None or fb is None or fa == 0 or fb == 0 or (fa > 0) == (fb > 0):
            continue
        a, b = xs[i], xs[i + 1]
        scale_ab = max(abs(fa), abs(fb))
        ok = True
        while abs(b - a) > tol:
            mid = (a + b) / 2
            if mid == a or mid == b:  # resolution of the number type reached
                break
            fm = f(mid)
            if fm is None:
                ok = False
                break
            if fm == 0:
                a = b = mid
                fa = fb = fm
                break
            if (fm > 0) == (fa > 0):
                a, fa = mid, fm
            else:
                b, fb = mid, fm
        if ok and min(abs(fa), abs(fb)) <= JUMP_RATIO.value * scale_ab:
            out.append((a + b) / 2)
    return sorted(out)


def backend_for(k, B):
    """Floats while the itinerary windows (about (B/2)^k wide) stay far above double resolution;
    otherwise mpmath with enough digits. Returns (backend, digits)."""
    depth = k * math.log10(2 / B)
    if depth < FLOAT_DEPTH.value:
        return FLOAT, 15
    return MP, int(math.ceil(depth)) + 20


def _stage(P, m, aim, j, B, xm):
    """State after the first j traversals of the ray from P aimed at polar angle ``aim`` of mouth m."""
    X = aim_point(P, m, aim, B, xm)
    if X is None:
        return None
    d = unit(add(X, P, -1), xm)
    E, w, _ = traverse(X, d, m, B)
    last = m  # mouth last entered; the ray now leaves through the other one, 1 - last
    for _ in range(j - 1):
        hit = first_entry(E, w, B, xm, skip=1 - last)
        if hit is None:
            return None
        t, last = hit
        E, w, _ = traverse(add(E, w, t), w, last, B)
    return E, w, last


def _gap(E, w, mm, B, xm):
    """Continuous signed clearance of the ray (E, w) from mouth mm: negative iff the ray enters it."""
    rel = (C[mm][0] - E[0], C[mm][1] - E[1])
    wu = unit(w, xm)
    if dot(rel, wu) > 0:
        return abs(cross(wu, rel)) - B
    return norm(rel, xm) - B


def _golden_min(f, a, b, tol):
    """Minimum of f on [a, b] by golden-section search (f defined on the whole interval)."""
    r = (math.sqrt(5) - 1) / 2
    c, d = b - r * (b - a), a + r * (b - a)
    fc, fd = f(c), f(d)
    while abs(b - a) > tol and a < c < d < b:
        if fc is None or fd is None:
            return None
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - r * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + r * (b - a)
            fd = f(d)
    return (a + b) / 2


def _hidden_dips(f, xs, tol):
    """Zeros of f hidden between grid points: at every grid local minimum with f > 0, refine the
    minimum by golden section; if it is negative, bisect to both zeros around it."""
    vals = [f(x) for x in xs]
    out = []
    for i in range(1, len(xs) - 1):
        v = vals[i]
        if v is None or v <= 0 or vals[i - 1] is None or vals[i + 1] is None:
            continue
        if not (v <= vals[i - 1] and v <= vals[i + 1]):
            continue
        xmin = _golden_min(f, xs[i - 1], xs[i + 1], tol)
        if xmin is None:
            continue
        fmin = f(xmin)
        if fmin is None or fmin >= 0:
            continue
        out += roots_by_bracketing(f, [xs[i - 1], xmin], tol) + roots_by_bracketing(f, [xmin, xs[i + 1]], tol)
    return out


def itinerary_windows(P, m, k, B, xm=FLOAT, grid=None):
    """Intervals of aim angle (on mouth m, seen from P) whose ray makes at least k traversals.

    A ray that leaves through one mouth cannot re-enter it (the mouths are convex), so the next
    entry, if any, is into the mouth it last entered: the itinerary repeats the same mouth and each traversal adds one condition, gap < 0, whose boundaries are found
    by bisection on a grid inside the parent interval."""
    arc = visible_arc(P, m, B, xm)
    if arc is None:
        return []
    center, half = arc
    N = grid or THETA_GRID.value
    lo, hi = center - half * (1 - 1e-12), center + half * (1 - 1e-12)

    def refine(lo, hi, j):
        if j == k:
            return [(lo, hi)]
        width = hi - lo
        xs = [lo + width * i / (N - 1) for i in range(N)]

        def g(a):
            st = _stage(P, m, a, j, B, xm)
            return None if st is None else _gap(st[0], st[1], st[2], B, xm)

        tol = _tol(width, xm)
        edges = roots_by_bracketing(g, xs, tol) + _hidden_dips(g, xs, tol)
        cuts = sorted([lo, hi] + edges)
        out = []
        for a, b in zip(cuts, cuts[1:]):
            mid = (a + b) / 2
            gm = g(mid)
            if gm is not None and gm < 0 and b > a:
                out += refine(a, b, j + 1)
        return out

    return refine(lo, hi, 1)


# ---------------------------------------------------------------------------------------------
# Single-collision solutions, by nested one-dimensional bracketing.
#
# The young ball comes in along the line p0 + s e (unit e) with speed v1. At the collision its
# center is at P = p0 + s e and the old ball's center is at Q = P - 2R n, with n = (cos nu, sin nu)
# the direction r1 - r2. After the collision the young ball leaves P aiming at polar angle ``aim``
# of mouth m, makes k traversals, and must arrive at Q as the old ball at the collision time.
#
#   inner  (aim): the loop passes through Q                         -> miss = 0
#   speed       : loop length / speed = net backward time travel     (explicit)
#   middle (nu) : normal components exchanged, (v2 - v1') . n = 0    -> C = 0
#   output      : v1, from (v1' - v1 e) parallel to n
#   outer  (s)  : v1(s) = v1*
# ---------------------------------------------------------------------------------------------


def _loop_state(P, nu, m, aim, k, B, R, xm):
    n = (xm.cos(nu), xm.sin(nu))
    Q = add(P, n, -2 * R)
    lp = loop(P, Q, m, aim, k, B, xm)
    if lp is None:
        return None
    w = lp["length"] / (-lp["shift"])
    lp.update(n=n, Q=Q, speed=w, u=scale(lp["dir0"], w), v2=scale(lp["dir2"], w))
    return lp


def _tol(width, xm):
    return width * (BISECT_REL_TOL.value if xm is FLOAT else mp.mpf(10) ** (-mp.mp.dps + 8))


def _aim_roots(P, nu, m, k, B, R, windows, xm, grid=None):
    N = grid or THETA_GRID.value
    out = []
    for lo, hi in windows:
        xs = [lo + (hi - lo) * i / (N - 1) for i in range(N)]

        def miss(a):
            st = _loop_state(P, nu, m, a, k, B, R, xm)
            return None if st is None else st["miss"]

        out += roots_by_bracketing(miss, xs, _tol(hi - lo, xm))
    return out


def _exchange(P, nu, m, aim, k, B, R, xm):
    st = _loop_state(P, nu, m, aim, k, B, R, xm)
    return None if st is None else dot(add(st["v2"], st["u"], -1), st["n"])


def _nu_intervals(P, m, k, B, R, windows, xm, aim_grid, N, log):
    """Pairs of neighbouring nu points with equal aim-root counts. A grid interval whose count
    changes is halved up to MAX_SPLIT times; what still changes is logged as unresolved."""
    def roots(nu):
        return _aim_roots(P, nu, m, k, B, R, windows, xm, aim_grid)

    out = []

    def split(na, ra, nb, rb, depth):
        if len(ra) == len(rb):
            out.append((na, ra, nb, rb))
            return
        if depth >= MAX_SPLIT.value:
            if log is not None:
                log.append(("nu_count_change", P, na, nb, len(ra), len(rb)))
            return
        nm = (na + nb) / 2
        rm = roots(nm)
        split(na, ra, nm, rm, depth + 1)
        split(nm, rm, nb, rb, depth + 1)

    nus = [2 * xm.pi * i / N for i in range(N + 1)]
    aims = [roots(nu) for nu in nus]
    for i in range(N):
        split(nus[i], aims[i], nus[i + 1], aims[i + 1], 0)
    return out


def refine_nu(P, m, k, B, R, xm, na, nb, depth, aim_grid=None):
    """Adaptive refinement of one unresolved nu-interval [na, nb] at the collision point P.

    The interval is halved up to ``depth`` more times. A leaf whose aim-root count still changes is
    resolved when the exchange condition C has one sign at every aim root at both of its ends: under
    the bracketing criterion of the main search, no C = 0 can then hide in it. A leaf with equal
    counts is resolved when no matched branch changes sign. Returns the leaves still unresolved, as
    (na, nb, reason)."""
    windows = itinerary_windows(P, m, k, B, xm, aim_grid)
    left = []

    def cvals(nu, roots):
        return [_exchange(P, nu, m, a, k, B, R, xm) for a in roots]

    def one_sign(vals):
        return vals and all(v is not None for v in vals) and len({v > 0 for v in vals}) == 1

    def leaf(na, ra, nb, rb, d):
        if len(ra) == len(rb):
            ca, cb = cvals(na, ra), cvals(nb, rb)
            if any(x is None or y is None or (x > 0) != (y > 0) for x, y in zip(ca, cb)):
                left.append((na, nb, "a matched branch changes sign: a middle root the main search did not follow"))
            return
        elif d >= depth:
            vals = cvals(na, ra) + cvals(nb, rb)
            if not (one_sign(vals) or (not ra and not rb)):
                left.append((na, nb, "aim-root count changes and C takes both signs at the ends"))
            return
        nm = (na + nb) / 2
        if nm == na or nm == nb:
            left.append((na, nb, "number resolution reached"))
            return
        rm = _aim_roots(P, nm, m, k, B, R, windows, xm, aim_grid)
        leaf(na, ra, nm, rm, d + 1)
        leaf(nm, rm, nb, rb, d + 1)

    if not windows:
        return left
    roots = lambda nu: _aim_roots(P, nu, m, k, B, R, windows, xm, aim_grid)
    leaf(na, roots(na), nb, roots(nb), 0)
    return left


def middle_roots(P, e, m, k, B, R, xm=FLOAT, nu_grid=None, aim_grid=None, log=None):
    """All (nu, aim) with miss = 0 and C = 0 at the collision point P, found by bracketing on a nu
    grid, with the aim roots linked across grid points by sorted order. Each result carries v1
    (the incoming speed this solution needs) and the kick (normal velocity transferred, > 0 for a
    physical collision). Every interval or bracket the search cannot follow is appended to ``log``
    (when given) instead of being dropped silently."""
    windows = itinerary_windows(P, m, k, B, xm, aim_grid)
    if not windows:
        return []
    N = nu_grid or NU_GRID.value
    found = []
    for nua, A, nub, Bn in _nu_intervals(P, m, k, B, R, windows, xm, aim_grid, N, log):
        for j in range(len(A)):
            ca = _exchange(P, nua, m, A[j], k, B, R, xm)
            cb = _exchange(P, nub, m, Bn[j], k, B, R, xm)
            if ca is None or cb is None:
                if log is not None:
                    log.append(("nu_state_lost", P, nua, nub))
                continue
            if (ca > 0) == (cb > 0):
                continue
            lo, hi, flo, aim = nua, nub, ca, A[j]
            ok = True
            while hi - lo > _tol(2 * xm.pi, xm):
                mid = (lo + hi) / 2
                if mid == lo or mid == hi:
                    break
                roots = _aim_roots(P, mid, m, k, B, R, windows, xm, aim_grid)
                if not roots:
                    ok = False
                    break
                aim = min(roots, key=lambda r: abs(r - aim))
                fm = _exchange(P, mid, m, aim, k, B, R, xm)
                if fm is None:
                    ok = False
                    break
                if (fm > 0) == (flo > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            nu = (lo + hi) / 2
            st = _loop_state(P, nu, m, aim, k, B, R, xm) if ok else None
            den = cross(e, st["n"]) if st is not None else 0
            if den == 0:
                if log is not None:
                    log.append(("nu_bisect_lost", P, nua, nub))
                continue
            v1 = cross(st["u"], st["n"]) / den
            kick = dot(add(st["u"], e, -v1), st["n"])
            found.append({"nu": nu, "aim": aim, "v1": v1, "kick": kick, "speed": st["speed"]})
    return found


def solve_mp(p0, e, s, nu, aim, v1, m, k, B, R, dps=None):
    """Polish a bracketed solution (s, nu, aim) at fixed incoming speed v1 with mpmath findroot on
    the three exact conditions (miss, C, v1). Returns the polished triple or None."""
    digits = max(dps or POLISH_DPS.value, backend_for(k, float(B))[1] + 20)
    with mp.workdps(digits):
        p0 = (mp.mpf(p0[0]), mp.mpf(p0[1]))
        e = (mp.mpf(e[0]), mp.mpf(e[1]))
        B, R, v1 = mp.mpf(B), mp.mpf(R), mp.mpf(v1)

        def F(s_, nu_, a_):
            P = add(p0, e, s_)
            st = _loop_state(P, nu_, m, a_, k, B, R, MP)
            if st is None:
                raise ValueError("itinerary lost")
            return [st["miss"], dot(add(st["v2"], st["u"], -1), st["n"]),
                    cross(st["u"], st["n"]) - v1 * cross(e, st["n"])]

        try:
            sol = mp.findroot(F, (mp.mpf(s), mp.mpf(nu), mp.mpf(aim)), tol=mp.mpf(10) ** (-(digits - 8)))
        except (ValueError, ZeroDivisionError):
            return None
        return tuple(sol)


def _circ(a, b):
    two_pi = 2 * mp.pi if isinstance(a, mp.mpf) or isinstance(b, mp.mpf) else 2 * math.pi
    d = (a - b) % two_pi
    return min(d, two_pi - d)


def _dedupe(roots, xm=FLOAT):
    tol = _tol(1, xm) * 1e4
    out = []
    for r in roots:
        if not any(_circ(r["nu"], o["nu"]) < tol and abs(r["aim"] - o["aim"]) < tol for o in out):
            out.append(r)
    return sorted(out, key=lambda r: (r["aim"], r["nu"] % (2 * math.pi)))


def _match(A, Bs):
    """Pair two equally long root lists by nearest (aim, nu)."""
    used, pairs = set(), []
    for a in A:
        j = min((j for j in range(len(Bs)) if j not in used),
                key=lambda j: abs(Bs[j]["aim"] - a["aim"]) + _circ(Bs[j]["nu"], a["nu"]))
        used.add(j)
        pairs.append((a, Bs[j]))
    return pairs


def sweep(p0, e, s_grid, m, k, B, R, xm=FLOAT, nu_grid=None, aim_grid=None, log=None):
    """Middle roots at every collision position of s_grid (the v1-independent part of the search)."""
    return [(s, _dedupe(middle_roots(add(p0, e, s), e, m, k, B, R, xm, nu_grid, aim_grid, log), xm)) for s in s_grid]


def crossings(table, p0, e, v1, m, k, B, R, xm=FLOAT, nu_grid=None, aim_grid=None, max_split=None, log=None):
    """Collision positions where some branch of v1(s) crosses the target v1, refined by bisection.

    Returns (candidates, unresolved, rejected). Candidates are (s, nu, aim) float triples.
    Unresolved lists, tagged by kind, everything the search could not follow: s-intervals where the
    number of middle roots changed and could not be linked after ``max_split`` halvings, s-bisections
    that lost their branch, plus whatever middle_roots logged (``log``, nu-intervals and nu-brackets).
    Rejected lists crossings found and refused, with the reason (a pole of v1(s), or kick <= 0).
    Nothing is dropped silently."""
    cands, rejected = [], []
    unresolved = log if log is not None else []
    max_split = MAX_SPLIT.value if max_split is None else max_split

    def roots_at(s):
        return _dedupe(middle_roots(add(p0, e, s), e, m, k, B, R, xm, nu_grid, aim_grid, unresolved), xm)

    def scan(sa, ra, sb, rb, depth):
        if len(ra) != len(rb):
            if depth >= max_split:
                unresolved.append(("s_count_change", sa, sb, len(ra), len(rb)))
                return
            sm = (sa + sb) / 2
            rm = roots_at(sm)
            scan(sa, ra, sm, rm, depth + 1)
            scan(sm, rm, sb, rb, depth + 1)
            return
        for a, b in _match(ra, rb):
            fa, fb = a["v1"] - v1, b["v1"] - v1
            if (fa > 0) == (fb > 0):
                continue
            lo, hi, rlo, rhi, flo = sa, sb, a, b, fa
            ok = True
            while hi - lo > _tol(1 + abs(lo), xm):
                sm = (lo + hi) / 2
                if sm == lo or sm == hi:
                    break
                rm = roots_at(sm)
                if not rm:
                    ok = False
                    break
                pick = min(rm, key=lambda r: abs(r["aim"] - rlo["aim"]) + _circ(r["nu"], rlo["nu"]))
                fm = pick["v1"] - v1
                if (fm > 0) == (flo > 0):
                    lo, rlo, flo = sm, pick, fm
                else:
                    hi, rhi = sm, pick
            if not ok:
                unresolved.append(("s_bisect_lost", sa, sb))
                continue
            tol_v = JUMP_RATIO.value * (1 + abs(v1))  # a v1 crossing across a pole of v1(s) is rejected
            if abs(rlo["v1"] - v1) >= tol_v or abs(rhi["v1"] - v1) >= tol_v:
                rejected.append(((lo, rlo["nu"], rlo["aim"]), "v1(s) crosses the target across a pole, not a root"))
            elif rlo["kick"] <= 0:
                rejected.append(((lo, rlo["nu"], rlo["aim"]), "kick <= 0: the balls would separate, not collide"))
            else:
                cands.append((lo, rlo["nu"], rlo["aim"]))

    for (sa, ra), (sb, rb) in zip(table, table[1:]):
        scan(sa, ra, sb, rb, 0)
    return cands, unresolved, rejected


# ---------------------------------------------------------------------------------------------
# Verifier: rebuild the whole history of a candidate and check every law of motion, without
# using the solver's equations.
# ---------------------------------------------------------------------------------------------


def _min_distance(p, v, t0, q, w, t1, lo, hi):
    """Minimum over t in [lo, hi] of |(p + v (t - t0)) - (q + w (t - t1))| (lo/hi may be +-inf)."""
    r0 = (p[0] - v[0] * t0 - q[0] + w[0] * t1, p[1] - v[1] * t0 - q[1] + w[1] * t1)
    vr = (v[0] - w[0], v[1] - w[1])
    a = dot(vr, vr)
    ts = [t for t in (lo, hi) if t not in (mp.inf, -mp.inf)]
    if a > 0:
        tstar = -dot(r0, vr) / a
        if lo <= tstar <= hi:
            ts.append(tstar)
    if not ts:
        return mp.sqrt(dot(r0, r0))
    return min(mp.sqrt(dot(add(r0, vr, t), add(r0, vr, t))) for t in ts)


def history(p0, e, v1, s, nu, aim, m, k, B, R, origin=None):
    """All straight segments of the ball's history for a single-collision solution, at the current
    mpmath precision. A segment is (point, velocity, t_point, t_start, t_end, label): the center is
    at point + velocity (t - t_point) for t_start <= t <= t_end. The collision is at t = 0.

    ``origin`` = (a0, e0): the incoming line p0 + s e is not the initial trajectory but the exit
    from mouth 1 of the initial line a0 + t e0, which entered mouth 2 at a0 (collision after one
    traversal)."""
    p0, e = (mp.mpf(p0[0]), mp.mpf(p0[1])), (mp.mpf(e[0]), mp.mpf(e[1]))
    B, R, v1 = mp.mpf(B), mp.mpf(R), mp.mpf(v1)
    P = add(p0, e, s)
    st = _loop_state(P, mp.mpf(nu), m, mp.mpf(aim), k, B, R, MP)
    if st is None:
        return None
    tc = mp.mpf(0)
    if origin is None:
        segs = [(P, scale(e, v1), tc, -mp.inf, tc, "incoming")]
    else:
        a0, e0 = (mp.mpf(origin[0][0]), mp.mpf(origin[0][1])), (mp.mpf(origin[1][0]), mp.mpf(origin[1][1]))
        t_exit = tc - mp.mpf(s) / v1
        segs = [(a0, scale(e0, v1), t_exit + 1, -mp.inf, t_exit + 1, "initial"),
                (P, scale(e, v1), tc, t_exit, tc, "incoming")]
    t = tc
    p, v = P, st["u"]
    skip = None
    for i in range(k):
        hit = first_entry(p, v, B, MP, skip)
        if hit is None:
            return None
        dt, mm = hit
        segs.append((p, v, t, t, t + dt, f"loop{i}"))
        X = add(p, v, dt)
        p, v, shift = traverse(X, v, mm, B)
        t = t + dt + shift
        skip = 1 - mm
    segs.append((p, v, t, t, tc, "return"))
    v2 = v
    n = st["n"]
    v2p = add(v2, n, -dot(add(v2, scale(e, v1), -1), n))
    p, v, t, skip = st["Q"], v2p, tc, None
    for i in range(FUTURE_TRAVERSALS.value):
        hit = first_entry(p, v, B, MP, skip)
        if hit is None:
            break
        dt, mm = hit
        segs.append((p, v, t, t, t + dt, f"out{i}"))
        X = add(p, v, dt)
        p, v, shift = traverse(X, v, mm, B)
        t, skip = t + dt + shift, 1 - mm
    else:
        return {"segments": segs, "state": st, "truncated": True}
    segs.append((p, v, t, t, mp.inf, "outgoing"))
    return {"segments": segs, "state": st, "truncated": False}


def verify(p0, e, v1, s, nu, aim, m, k, B, R, dps=None, origin=None):
    """Check a polished candidate against every law of the model. Returns (ok, list of failures)."""
    with mp.workdps(max(dps or POLISH_DPS.value, backend_for(k, float(B))[1] + 20)):
        tol = VERIFY_TOL.value
        h = history(p0, e, v1, s, nu, aim, m, k, B, R, origin)
        if h is None:
            return False, ["itinerary does not exist"]
        segs, st = h["segments"], h["state"]
        first = 1 if origin is None else 2
        R_, v1_ = mp.mpf(R), mp.mpf(v1)
        e_ = (mp.mpf(e[0]), mp.mpf(e[1]))
        fails = []
        P = segs[first - 1][0]
        ret = segs[first + k]
        Q = add(ret[0], ret[1], ret[4] - ret[2])
        if abs(mp.sqrt(dot(add(P, Q, -1), add(P, Q, -1))) - 2 * R_) > tol:
            fails.append("old ball does not arrive at distance 2R at the collision time")
        n = unit(add(P, Q, -1), MP)
        u, v2 = segs[first][1], ret[1]
        rel = dot(add(v2, scale(e_, v1_), -1), n)
        if rel <= 0:
            fails.append("balls are not approaching along the line of centers")
        if mp.sqrt(dot(add(u, add(scale(e_, v1_), n, rel), -1), add(u, add(scale(e_, v1_), n, rel), -1))) > tol:
            fails.append("outgoing velocity is not the elastic-collision result")
        mid = scale(add(P, Q), mp.mpf(1) / 2)
        for c in C:
            if mp.sqrt((mid[0] - c[0]) ** 2 + (mid[1] - c[1]) ** 2) <= mp.mpf(B):
                fails.append("contact point inside a mouth (excluded by the model, p. 1081 item i)")
        if origin is None:
            back = first_entry(P, scale(e_, -1), mp.mpf(B), MP)
            if back is not None:
                fails.append("incoming line, extended to the past, comes out of a mouth")
        else:
            a0, e0 = segs[0][0], segs[0][1]
            if first_entry(a0, scale(e0, -1), mp.mpf(B), MP) is not None:
                fails.append("initial line, extended to the past, comes out of a mouth")
            p0_ = segs[1][0]
            hit = first_entry(add(p0_, segs[1][1], segs[1][3] - segs[1][2]), segs[1][1], mp.mpf(B), MP, skip=0)
            if hit is not None and hit[0] < -segs[1][3]:
                fails.append("incoming leg enters a mouth before the collision")
        fwd = first_entry(ret[0], ret[1], mp.mpf(B), MP, skip=None)
        if fwd is not None and fwd[0] < (ret[4] - ret[2]):
            fails.append("return leg enters a mouth before reaching the collision")
        touching = {("incoming", "loop0"), ("return", "out0"), ("return", "outgoing")}
        for i in range(len(segs)):
            for j in range(i + 1, len(segs)):
                a, b = segs[i], segs[j]
                if (a[5], b[5]) in touching:
                    continue
                lo, hi = max(a[3], b[3]), min(a[4], b[4])
                if lo > hi:
                    continue
                dmin = _min_distance(a[0], a[1], a[2], b[0], b[1], b[2], lo, hi)
                if dmin < 2 * R_ - tol:
                    fails.append(f"overlap between {a[5]} and {b[5]} (distance {mp.nstr(dmin, 5)})")
        if h["truncated"]:
            fails.append(f"outgoing ball still traversing after {FUTURE_TRAVERSALS.value} traversals (not checked further)")
        return (not fails), fails


def paper_angles(p0, e, v1, s, nu, aim, m, k, B, R):
    """The quantities of EKT Fig. 11/17 for a solution, measured on the geometry:
    alpha = arg(v1) - arg(v2) (angle between the incoming balls), beta = arg(v1) - arg(v1')
    (deflection, positive to the right), d = signed distance the young ball travels past the
    intersection of the two incoming lines, and the speeds v1', v2."""
    with mp.workdps(max(POLISH_DPS.value, backend_for(k, float(B))[1] + 20)):
        p0, e = (mp.mpf(p0[0]), mp.mpf(p0[1])), (mp.mpf(e[0]), mp.mpf(e[1]))
        P = add(p0, e, s)
        st = _loop_state(P, nu, m, aim, k, mp.mpf(B), mp.mpf(R), MP)
        a1, a2, a1p = mp.atan2(e[1], e[0]), mp.atan2(st["v2"][1], st["v2"][0]), mp.atan2(st["u"][1], st["u"][0])
        wrap = lambda x: (x + mp.pi) % (2 * mp.pi) - mp.pi
        alpha = (a1 - a2) % (2 * mp.pi)
        beta = wrap(a1 - a1p)
        d2 = unit(st["v2"], MP)
        X_t = cross(add(st["Q"], P, -1), d2) / cross(e, d2)  # P + X_t e is on the old ball's line
        return {"alpha": alpha, "beta": beta, "d": -X_t, "v1p": st["speed"], "v2": st["speed"], "P": P, "Q": st["Q"]}


def free_history(p0, e, v1, B, R, traversals=None):
    """Collision-free history of the ball that comes in along p0 + t e (p0 on the incoming line,
    reached at t = 0) with speed v1. Returns (consistent, segments, overlaps): consistent means no
    two incarnations ever come closer than 2R, i.e. the collision-free motion is itself a solution."""
    with mp.workdps(POLISH_DPS.value):
        p, v = (mp.mpf(p0[0]), mp.mpf(p0[1])), scale((mp.mpf(e[0]), mp.mpf(e[1])), mp.mpf(v1))
        B, R = mp.mpf(B), mp.mpf(R)
        segs, t, skip, start = [], mp.mpf(0), None, -mp.inf
        for i in range(traversals or FUTURE_TRAVERSALS.value):
            hit = first_entry(p, v, B, MP, skip)
            if hit is None:
                break
            dt, m = hit
            segs.append((p, v, t, start, t + dt, f"free{i}"))
            X = add(p, v, dt)
            p, v, shift = traverse(X, v, m, B)
            t, start, skip = t + dt + shift, t + dt + shift, 1 - m
        segs.append((p, v, t, start, mp.inf, "free_out"))
        overlaps = []
        for i in range(len(segs)):
            for j in range(i + 1, len(segs)):
                a, b = segs[i], segs[j]
                lo, hi = max(a[3], b[3]), min(a[4], b[4])
                if lo > hi:
                    continue
                dmin = _min_distance(a[0], a[1], a[2], b[0], b[1], b[2], lo, hi)
                if dmin < 2 * R:
                    overlaps.append((a[5], b[5], dmin))
        return (not overlaps), segs, overlaps
