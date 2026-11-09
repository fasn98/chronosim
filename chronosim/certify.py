"""Phase 3: interval exclusion of single-collision solutions (mpmath.iv, outward rounding).

What is proven, and under which hypotheses (docs/phase3_certification.md has the full statement):

- **Model:** the exact geometric model of billiards.py (EKT 1991, Sec. I B p. 1079, Sec. II A
  p. 1081, Fig. 2), coplanar, D = T_d = 1, mouths of radius B at C1 = (0, 0) and C2 = (1, 0).
- **Family examined:** one self-collision; the young ball comes in along the initial line of
  Fig. 11, p0 + s e (phase3.line_theta_phi), and collides at P = p0 + s e with s in [s_lo, s_hi],
  before reaching mouth 2 (s < 0); after the collision it travels straight to mouth 2, makes
  exactly one traversal (out of mouth 1, one time unit earlier) and comes back as the old ball,
  which is at Q = P - 2R n at the collision, n = (cos nu, sin nu).
- **Unknowns:** (s, nu, aim), aim the polar angle of the entry point X = C2 + B (cos aim, sin aim).
- **Conditions** (the same as billiards._loop_state and billiards.solve_mp, written so they can
  be evaluated on intervals):
    miss = cross(w, Q - E) = 0           the loop passes through Q
    C    = (v2 - u) . n = 0              the normal components are exchanged
    V    = cross(u, n) - v1 cross(e, n) = 0   the incoming speed is v1 of eq. 4.10
  with the necessary inequalities d . n_X < 0 (the ray enters mouth 2 at X, transversally),
  along = (Q - E) . w > 0 (the old ball has left mouth 1 before the collision) and
  kick = (u - v1 e) . n > 0 (the impulse pushes the balls apart).
- **Test:** a box is excluded when the interval enclosure of one of miss, C, V does not contain
  0, or one inequality fails on the whole box. A box that is not excluded is halved along its
  widest side (relative to the side of the root box) until it is excluded or narrower than
  CERT_MIN_WIDTH; what is left is reported as undecided, never dropped. If nothing is left, no
  point of the root box satisfies the conditions. That is the whole claim: other families (more
  traversals, collision after a traversal, s outside the range, several collisions) are not
  examined.

The superset aim in [-AIM_PAD, AIM_PAD], nu in [-NU_PAD, 2 pi + NU_PAD] is used; aims that are not
visible from P are excluded by the inequality d . n_X < 0, not by the search.
"""

from __future__ import annotations

import math
from contextlib import contextmanager

import mpmath as mp
from mpmath import iv

from .provenance import CERT_EXCLUSION, Origin, Quantity, certified, heuristic

CERT_PREC = heuristic("ekt_cert_prec", 53, "bits of the interval arithmetic (mpmath.iv) of the exclusion test")
CERT_MIN_WIDTH = heuristic("ekt_cert_min_width", 1e-9, "relative width below which a box that is not excluded is reported as undecided")
CERT_MAX_BOXES = heuristic("ekt_cert_max_boxes", 20_000_000, "boxes examined by one exclusion run before the rest is reported as undecided")
AIM_PAD = heuristic("ekt_cert_aim_pad", 3.15, "aim range [-pad, pad] of the exclusion box (> pi: the circle with margin; any superset is valid)")
NU_PAD = heuristic("ekt_cert_nu_pad", 0.01, "nu range [-pad, 2 pi + pad] of the exclusion box (any superset is valid)")


@contextmanager
def ivprec(bits):
    """Working precision of mpmath.iv for a block (the interval context has no workprec)."""
    saved = iv.prec
    iv.prec = bits
    try:
        yield
    finally:
        iv.prec = saved


def inputs(case, xm=iv):
    """case_inputs for either form of a case: the (theta, phi, lam1, B, R) tuple of the dangerous
    triangle, or a dict {"p0": (x, y), "e": (x, y), "v1": v, "B": B, "R": R} of decimal strings for a
    line given directly (the symmetric trajectory of Figs. 6 and 7: p0 = (0.5, 0), e = (0, 1))."""
    if isinstance(case, dict):
        f = xm.mpf
        return {"p0": (f(case["p0"][0]), f(case["p0"][1])), "e": (f(case["e"][0]), f(case["e"][1])),
                "v1": f(case["v1"]), "B": f(case["B"]), "R": f(case["R"])}
    return case_inputs(*case, xm=xm)


def case_inputs(theta, phi, lam1, B, R, xm=iv):
    """The initial condition as numbers of ``xm`` (mpmath.iv or mpmath): p0 and e of the initial line
    (phase3.line_theta_phi), v1 of eqs. 4.8 and 4.10 (ekt.v1_range, ekt.v1_from_lambda), B and R.
    With iv the inputs are enclosures of the exact values of the decimal strings."""
    th, ph, lam, B, R = (xm.mpf(x) for x in (theta, phi, lam1, B, R))
    base = xm.cos(th) / xm.cos(ph) * (1 - 2 * B * xm.cos(ph))          # eq. 4.8
    lo, hi = base - 2 * R / xm.cos(ph), base + 2 * R / xm.cos(ph)
    v1 = lo + lam * (hi - lo)                                           # eq. 4.10
    p0 = (1 - B * xm.cos(ph), -B * xm.sin(ph))
    e = (xm.cos(th + ph), xm.sin(th + ph))
    return {"p0": p0, "e": e, "v1": v1, "B": B, "R": R}


def conditions(s, nu, aim, c, xm=iv):
    """miss, C and V up to positive factors (L0 = |X - P| and L0 W), and the three inequality
    quantities, at (s, nu, aim) for the inputs c of case_inputs. Every operation is one of + - * /,
    sqrt, sin, cos of ``xm``, so with mpmath.iv the results enclose the exact values.

    Returned: {"miss": L0 miss, "C": L0 (w - dir0) . n, "V": V, "dn": d . n_X, "along": L0 along,
    "kick": kick}. On the set where along > 0 the speed W = L0 + along is positive, so
    C = 0 iff (w - dir0) . n = 0, and the signs of miss and along are those of billiards."""
    p0, e, v1, B, R = c["p0"], c["e"], c["v1"], c["B"], c["R"]
    P = (p0[0] + s * e[0], p0[1] + s * e[1])
    n = (xm.cos(nu), xm.sin(nu))
    Q = (P[0] - 2 * R * n[0], P[1] - 2 * R * n[1])
    nx = (xm.cos(aim), xm.sin(aim))
    X = (1 + B * nx[0], B * nx[1])                       # on mouth 2, center C2 = (1, 0)
    d = (X[0] - P[0], X[1] - P[1])                       # L0 dir0
    dn = d[0] * nx[0] + d[1] * nx[1]
    L0 = xm.sqrt(d[0] ** 2 + d[1] ** 2)                  # ** 2, not x * x: an interval square is >= 0
    r = (d[0] - 2 * dn * nx[0], d[1] - 2 * dn * nx[1])   # reflected in the tangent at X (billiards.traverse)
    wd = (-r[0], r[1])                                   # then mirrored x -> 1 - x: L0 w
    E = (1 - X[0], X[1])                                 # exit point on mouth 1
    tq = (Q[0] - E[0], Q[1] - E[1])
    along = tq[0] * wd[0] + tq[1] * wd[1]
    miss = wd[0] * tq[1] - wd[1] * tq[0]
    Cn = (wd[0] - d[0]) * n[0] + (wd[1] - d[1]) * n[1]
    W = L0 + along / L0                                  # loop length / time shift (shift = -1)
    V = W * (d[0] * n[1] - d[1] * n[0]) / L0 - v1 * (e[0] * n[1] - e[1] * n[0])
    kick = W * (d[0] * n[0] + d[1] * n[1]) / L0 - v1 * (e[0] * n[0] + e[1] * n[1])
    return {"miss": miss, "C": Cn, "V": V, "dn": dn, "along": along, "kick": kick}


def _reason(box, c):
    s, nu, aim = (iv.mpf([lo, hi]) for lo, hi in box)
    q = conditions(s, nu, aim, c)
    if q["dn"].a >= 0:
        return "entrada na boca 2 impossível (d·n_X ≥ 0)"
    if q["along"].b <= 0:
        return "bola velha ainda não saiu da boca 1 (along ≤ 0)"
    if q["kick"].b <= 0:
        return "impulso ≤ 0"
    for name in ("miss", "C", "V"):
        f = q[name]
        if f.a > 0 or f.b < 0:
            return {"miss": "miss ≠ 0", "C": "C ≠ 0", "V": "V ≠ 0"}[name]
    return None


def normalize(point):
    """(s, nu, aim) with nu in [0, 2 pi) and aim in [-pi, pi), inside the default root box."""
    s, nu, aim = (mp.mpf(x) for x in point)
    return (s, nu % (2 * mp.pi), (aim + mp.pi) % (2 * mp.pi) - mp.pi)


def _inside(point, box):
    return all(lo <= x <= hi for x, (lo, hi) in zip(point, box))


def default_ranges():
    """(nu_range, aim_range) of the default root box (supersets of the circle)."""
    return (-NU_PAD.value, 2 * math.pi + NU_PAD.value), (-AIM_PAD.value, AIM_PAD.value)


def periodic_copies(point, nu_range=None, aim_range=None):
    """The point, normalized, and its copies at nu +- 2 pi and aim +- 2 pi that fall in the nu and
    aim ranges of the root box. They are the same solution, and each is a root of the conditions
    in the box, so a positive control watches all of them."""
    dnu, daim = default_ranges()
    nu_range, aim_range = nu_range or dnu, aim_range or daim
    s, nu, aim = normalize(point)
    out = []
    for i in (0, -1, 1):
        for j in (0, -1, 1):
            q = (s, nu + i * 2 * mp.pi, aim + j * 2 * mp.pi)
            if nu_range[0] <= q[1] <= nu_range[1] and aim_range[0] <= q[2] <= aim_range[1]:
                out.append(q)
    return out


def exclude(case, s_range, nu_range=None, aim_range=None, max_boxes=None, watch=(), scale=None):
    """Branch and bound over the box s_range x nu_range x aim_range (floats). ``case`` is a case of
    inputs(). Returns {"boxes", "excluded_by", "undecided", "complete", "watch"}: ``complete`` is True
    when every box was excluded.

    ``watch`` (positive control) is a list of points (s, nu, aim) known to be solutions; it does not
    change the search. A point whose nu and aim lie in the root box is watched as given (so that
    periodic copies can be watched); any other point is normalized first. For each point, "watch" lists the excluded boxes that contain it (each one a
    failure of the test: an excluded box cannot contain a solution) and the indices of the undecided
    boxes that contain it.

    ``scale`` (three widths, default: the sides of the root box) is what the sides are measured
    against when the widest side is chosen and when a box is small enough to stop
    (ekt_cert_min_width). A run on a piece of a larger box passes the scale of that box, so that the
    piece is resolved as finely as the whole, no more and no less."""
    dnu, daim = default_ranges()
    nu_range, aim_range = nu_range or dnu, aim_range or daim
    root = (tuple(s_range), tuple(nu_range), tuple(aim_range))
    if min(hi - lo for lo, hi in root) <= 0:
        raise ValueError(f"every side of the root box must have positive width: {root}")
    widths = list(scale) if scale else [hi - lo for lo, hi in root]
    if len(widths) != 3:
        raise ValueError(f"scale needs three widths (s, nu, aim): {widths}")
    if min(widths) <= 0:
        raise ValueError(f"every scale width must be positive: {widths}")
    max_boxes = max_boxes or CERT_MAX_BOXES.value
    points = [tuple(mp.mpf(x) for x in p) for p in watch]
    points = [p if _inside(p[1:], root[1:]) else normalize(p) for p in points]
    hits = [{"excluded_in": [], "undecided_in": []} for _ in points]
    with ivprec(CERT_PREC.value):
        c = inputs(case)
        stack, n, by, undecided, unexamined = [root], 0, {}, [], 0
        while stack:
            if n >= max_boxes:
                unexamined = len(stack)  # never examined: reported as undecided, never as excluded
                undecided += stack
                break
            box = stack.pop()
            n += 1
            why = _reason(box, c)
            if why:
                by[why] = by.get(why, 0) + 1
                for p, h in zip(points, hits):
                    if _inside(p, box):
                        h["excluded_in"].append((box, why))
                continue
            rel = [(hi - lo) / w for (lo, hi), w in zip(box, widths)]
            i = max(range(3), key=lambda j: rel[j])
            lo, hi = box[i]
            mid = (lo + hi) / 2
            if rel[i] < CERT_MIN_WIDTH.value or not lo < mid < hi:
                undecided.append(box)
                continue
            for half in ((lo, mid), (mid, hi)):
                b = list(box)
                b[i] = half
                stack.append(tuple(b))
    for p, h in zip(points, hits):
        h["undecided_in"] = [i for i, b in enumerate(undecided) if _inside(p, b)]
    return {"boxes": n, "excluded_by": by, "undecided": undecided, "complete": not undecided, "watch": hits,
            "unexamined": unexamined}


def s_range_outward(lo, hi):
    """Float s-range that contains [lo, hi] (decimal or mpf ends), rounded outward."""
    lo_f, hi_f = float(mp.mpf(lo)), float(mp.mpf(hi))
    return (math.nextafter(lo_f, -math.inf), math.nextafter(hi_f, math.inf))


HYPOTHESES = ("uma colisão; a bola jovem colide na linha inicial (Fig. 11) antes da boca 2, em s ∈ {s}; "
              "exatamente 1 travessia entre as duas visitas (entra na boca 2 e sai da boca 1 uma unidade de tempo antes); "
              "entrada transversal na boca 2; a bola velha já fora da boca 1 no instante da colisão; impulso normal > 0; "
              "modelo geométrico de billiards.py (EKT 1991, Sec. II A, p. 1081); mpmath.iv a {prec} bits, arredondamento para fora")


def exclusion_quantity(result, case, s_range, what):
    """CERT_EXCLUSION when the run is complete: the value is the number of solutions in the box, 0,
    and the source states the box, the equations, the hypotheses and the number of sub-boxes.
    Otherwise a COMPUTED count of the undecided boxes, with no claim about solutions."""
    hyp = HYPOTHESES.format(s=f"[{s_range[0]!r}; {s_range[1]!r}]", prec=CERT_PREC.value)
    if result["complete"]:
        src = (f"{what}: nenhuma solução de (miss, C, V) = 0 na caixa, para (θ, φ, λ₁, B, R) = {case}; "
               f"{result['boxes']} subcaixas, todas excluídas; hipóteses: {hyp}")
        return certified(0, CERT_EXCLUSION, src)
    return Quantity(len(result["undecided"]), Origin.COMPUTED, f"{what}: caixas indecididas da exclusão intervalar; {hyp}")


# Krawczyk test --------------------------------------------------------------------------------
#
# For F = (miss, C, V) of conditions() on a box X with midpoint x~, any real matrix Y and the
# interval Jacobian J(X):  K(X) = x~ - Y F(x~) + (I - Y J(X)) (X - x~).
# Theorem (Krawczyk 1969; Moore 1977): if K(X) lies in the interior of X, F has exactly one zero
# in X; if K(X) and X are disjoint, F has none. J(X) is enclosed by forward-mode automatic
# differentiation on intervals (the class _D below), so the derivatives are those of the very
# expressions of conditions(). The scaled F of conditions() has the same zeros as the conditions
# of the model wherever L0 > 0 and along > 0, which is checked on X with the inequalities.


class _D:
    """Interval value with its gradient in (s, nu, aim)."""

    __slots__ = ("v", "d")

    def __init__(self, v, d=(0, 0, 0)):
        self.v, self.d = v, tuple(d)

    @staticmethod
    def lift(x):
        return x if isinstance(x, _D) else _D(x)

    def __add__(self, o):
        o = _D.lift(o)
        return _D(self.v + o.v, (a + b for a, b in zip(self.d, o.d)))

    __radd__ = __add__

    def __neg__(self):
        return _D(-self.v, (-a for a in self.d))

    def __sub__(self, o):
        return self + (-_D.lift(o))

    def __rsub__(self, o):
        return _D.lift(o) + (-self)

    def __mul__(self, o):
        o = _D.lift(o)
        return _D(self.v * o.v, (a * o.v + self.v * b for a, b in zip(self.d, o.d)))

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = _D.lift(o)
        q = self.v / o.v
        return _D(q, ((a - q * b) / o.v for a, b in zip(self.d, o.d)))

    def __rtruediv__(self, o):
        return _D.lift(o) / self

    def __pow__(self, k):
        if k != 2:
            raise ValueError("only squares are used")
        return _D(self.v ** 2, (2 * self.v * a for a in self.d))


class _DX:
    """The functions conditions() needs, on _D numbers over mpmath.iv."""

    @staticmethod
    def cos(x):
        return _D(iv.cos(x.v), (-iv.sin(x.v) * a for a in x.d))

    @staticmethod
    def sin(x):
        return _D(iv.sin(x.v), (iv.cos(x.v) * a for a in x.d))

    @staticmethod
    def sqrt(x):
        r = iv.sqrt(x.v)
        return _D(r, (a / (2 * r) for a in x.d))


def jacobian(box, c):
    """Interval enclosure of the Jacobian of (miss, C, V) over the box, and the inequality values."""
    xs = [_D(iv.mpf([lo, hi]), tuple(1 if j == i else 0 for j in range(3))) for i, (lo, hi) in enumerate(box)]
    cd = {k: (tuple(_D(x) for x in v) if isinstance(v, tuple) else _D(v)) for k, v in c.items()}  # iv cannot add a _D
    q = conditions(*xs, cd, xm=_DX)
    return [list(q[k].d) for k in ("miss", "C", "V")], {k: q[k].v for k in ("dn", "along", "kick")}


def krawczyk(case, box):
    """Krawczyk test of (miss, C, V) = 0 on a box (three (lo, hi) float pairs). Returns
    {"verdict": "root" | "none" | "undecided", "K": the K(X) intervals as (lo, hi), "inequalities":
    True when d.n_X < 0, along > 0 and kick > 0 hold on the whole box}."""
    with ivprec(CERT_PREC.value):
        c = inputs(case)
        X = [iv.mpf([lo, hi]) for lo, hi in box]
        mid = [(lo + hi) / 2 for lo, hi in box]
        Fm = conditions(*(iv.mpf(m) for m in mid), c)
        F = [Fm[k] for k in ("miss", "C", "V")]
        J, ineq = jacobian(box, c)
        Jmid = mp.matrix([[float(J[i][j].mid) for j in range(3)] for i in range(3)])
        try:
            Yf = Jmid ** -1
        except ZeroDivisionError:
            return {"verdict": "undecided", "K": None, "inequalities": False}
        Y = [[iv.mpf(float(Yf[i, j])) for j in range(3)] for i in range(3)]
        K = []
        for i in range(3):
            k = iv.mpf(mid[i]) - sum(Y[i][j] * F[j] for j in range(3))
            for j in range(3):
                m = (1 if i == j else 0) - sum(Y[i][l] * J[l][j] for l in range(3))
                k = k + m * (X[j] - mid[j])
            K.append(k)
        ok_ineq = ineq["dn"].b < 0 and ineq["along"].a > 0 and ineq["kick"].a > 0
        if all(lo < k.a and k.b < hi for k, (lo, hi) in zip(K, box)):
            verdict = "root"
        elif any(k.b < lo or k.a > hi for k, (lo, hi) in zip(K, box)):
            verdict = "none"
        else:
            verdict = "undecided"
        return {"verdict": verdict, "K": [(float(k.a), float(k.b)) for k in K], "inequalities": ok_ineq}
