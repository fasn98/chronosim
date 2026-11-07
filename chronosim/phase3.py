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
    """All single-collision solutions on the s grid: float bracketing, mpmath polish, verification."""
    xm, _ = bl.backend_for(k, float(B))
    conv = (lambda x: float(x)) if xm is bl.FLOAT else (lambda x: mp.mpf(x))
    pf, ef = (conv(p0[0]), conv(p0[1])), (conv(e[0]), conv(e[1]))
    Bf, Rf, vf = conv(B), conv(R), conv(v1)
    sg = [conv(s) for s in s_grid]
    table = bl.sweep(pf, ef, sg, m, k, Bf, Rf, xm, nu_grid, aim_grid)
    cands, unresolved = bl.crossings(table, pf, ef, vf, m, k, Bf, Rf, xm, nu_grid, aim_grid)
    out = []
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
        rows.append({"n": n, "L": Ls, "unresolved": len(unres), "rejected": [s["fails"] for s in sols if not s["ok"]],
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
    sg = grid(mp.mpf(-3), -float(R) / 100, S_GRID.value, windows)
    sols, unres = solutions(p0, e, v1, 1, k, B, R, sg)
    return {"v1": v1, "p0": p0, "e": e, "s0": s0, "sols": sols, "unresolved": unres, "theta": th, "phi": ph, "B": B, "R": R, "lam1": lam}


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
                "sign_rule": (mp.sign(a["beta"]) == sign),
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


def after_traversal(theta, phi, lam1, B, R, k=1):
    """Single-collision solutions whose first encounter with the collision is on the collision-free
    exit line from mouth 1 (after one traversal), with k more traversals between the two visits."""
    th, ph, lam, Bm, Rm = (mp.mpf(x) for x in (theta, phi, lam1, B, R))
    v1 = ekt.v1_from_lambda(lam, th, ph, Bm, Rm).value
    a0, e0 = line_theta_phi(th, ph, Bm)
    E, w, _ = bl.traverse(a0, e0, 1, Bm)
    sg = grid(mp.mpf("1e-6"), mp.mpf(3), S_GRID.value)
    return solutions(E, w, v1, 1, k, Bm, Rm, sg, origin=(a0, e0))


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
        return {"kind": kind, "k": k, "refined": refined, "lam1": lam, "ok": len([x for x in sols if x["ok"]]),
                "rejected": [x["fails"] for x in sols if not x["ok"]], "unresolved": len(unres)}
    if name == "r8":
        case = dangerous(*arg)
        phys, neg = eq42_43_roots(case)
        free_ok, _, overlaps = bl.free_history(bl.add(case["p0"], case["e"], -1), case["e"], case["v1"], case["B"], case["R"])
        return {"args": arg, "v1": case["v1"], "free_consistent": free_ok, "free_min_distance": min(o[2] for o in overlaps) if overlaps else None, "rows": r8_eq42_43(case), "physical": [p[:3] for p in phys],
                "negative": [(a, b, s, d) for a, b, s, d in neg], "rejected": [s["fails"] for s in case["sols"] if not s["ok"]],
                "unresolved": len(case["unresolved"])}
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


def reproduce(workers=None):
    """Run R1-R11 in parallel worker processes (deterministic grids only; each task is a pure function
    of its arguments). Returns a dict of results keyed by task."""
    from concurrent.futures import ProcessPoolExecutor
    tasks = [("r1", None), ("r2", None), ("rest", None)]
    tasks += [("r3", R) for R in ("1e-4", "1e-5", "1e-6")]
    tasks += [("r46", B) for B in R46_BS]
    tasks += [("r7", n) for n in range(1, 7)]
    tasks += [("r7b", n) for n in range(1, 7)]
    tasks += [("r8", c) for c in R8_CASES]
    tasks += [("r8x", (kind, k, False, None)) for kind, k in (("before", 2), ("before", 3), ("after", 1), ("after", 2))]
    tasks += [("r8x", ("before", k, True, None)) for k in (1, 2, 3)]
    tasks += [("r8x", ("before", 1, False, lam)) for lam in ("0.1", "0.25", "0.5", "0.9")]
    tasks += [("r9", B) for B in R9_BS]
    tasks += [("r10", ph) for ph in R10_PHIS]
    tasks += [("r10x", ph) for ph in ("0.3", "0.9", "1.1")]
    tasks += [("r11", ph) for ph in R11_PHIS]
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker) as ex:
        futures = [ex.submit(_task, n, a) for n, a in tasks]
        done = [f.result() for f in futures]
    res = {}
    for (n, a), r in zip(tasks, done):
        res.setdefault(n, []).append((a, r))
    return res


# Report ----------------------------------------------------------------------------------------

def _f(x, n=3):
    """Number in Portuguese notation (decimal comma)."""
    return mp.nstr(x, n).replace(".", ",")


def _order(pairs):
    """Least-squares slope of log|dev| against log B."""
    xs = [mp.log(B) for B, _ in pairs]
    ys = [mp.log(abs(d)) for _, d in pairs]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)


def report_section(res):
    """Markdown lines for outputs/report.md. Every number here is COMPUTED."""
    L = []
    add = L.append
    add("")
    add("## Fase 3 · Bolas de bilhar (Echeverria, Klinkhammer & Thorne, 1991): reprodução R1–R11")
    add("")
    add("Todas as grandezas desta seção são **COMPUTED**.")
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
    add("  - \"Intervalos não resolvidos\" são trechos da grade em que o número de raízes mudou e não pôde ser ligado após 12 bisseções. Eles são contados, nunca descartados.")
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
        add("| n | L medido | n da eq. 3.7 no L medido | intervalos não resolvidos |")
        add("|---|---|---|---|")
        for r in rows:
            row = r["rows"][0]
            add(f"| {row['n']} | {'; '.join(_f(x, 7) for x in row['L']) or 'nenhuma solução encontrada'} | {'; '.join(_f(x, 5) for x in row['n_37']) or '—'} | {row['unresolved']} |")
        add("")
    first = [r["rows"][0] for _, r in res["r7"]]
    offs = [row["n_37"][0] - row["n"] for row in first if row["n_37"]]
    found_n = [row["n"] for row in first if row["L"]]
    add(f"- **Primeiro conjunto:** soluções encontradas para n = {', '.join(map(str, found_n))}. A eq. 3.7, avaliada no L medido, dá n mais {_f(min(offs), 3)} a {_f(max(offs), 3)}, então o desvio relativo cai aproximadamente como 1/n. A eq. 3.7 vale para n ≫ 1, L ≫ 1, e este conjunto está no começo desse regime. O menor n encontrado é {min(found_n)}; a eq. 3.8, também assintótica, dá {_f(res['r7'][0][1]['n_min_38'], 4)}.")
    second = [r["rows"][0] for _, r in res["r7b"]]
    add(f"- **Segundo conjunto:** soluções encontradas só para n = {', '.join(str(row['n']) for row in second if row['L']) or 'nenhum'}, com {sum(row['unresolved'] for row in second)} intervalos não resolvidos no total. A busca é **inconclusiva** aqui: a ausência das outras soluções não é afirmada.")
    add("")

    add("### R8: árbitro das eqs. 4.2–4.3 contra a geometria independente")
    add("")
    add("Para cada solução geométrica verificada (colisão antes de uma travessia), α, β e d são medidos na geometria e postos nas eqs. 4.2 e 4.3 como impressas. Em separado, as eqs. 4.2–4.3 são resolvidas diretamente, por findroot a partir de uma grade de pontos de partida, e suas raízes são confrontadas com as geométricas.")
    add("")
    add("| (θ, φ, λ₁, B, R) | movimento livre consistente? | soluções geométricas | raízes físicas das eqs. 4.2–4.3 | coincidem | raízes com distância negativa | maior resíduo das eqs. 4.2–4.3 | maior desvio de d (4.4b) e de v₂ (4.5) |")
    add("|---|---|---|---|---|---|---|---|")
    for arg, r in res["r8"]:
        rows = r["rows"]
        matched = sum(1 for (a, b, sg) in r["physical"] if any(abs(a - x["alpha"]) < mp.mpf("1e-12") and abs(b - x["beta"]) < mp.mpf("1e-12") for x in rows))
        res_max = max((x["res"].value for x in rows), default=None)
        dmax = max((max(x["dev_d"].value, x["dev_v2"].value) for x in rows), default=None)
        classes = f" ({', '.join(x['class'] for x in rows)})" if rows else ""
        add(f"| ({'; '.join(a.replace('.', ',') for a in arg)}) | {'sim' if r['free_consistent'] else 'não'} | {len(rows)}{classes} | "
            f"{len(r['physical'])} | {matched}/{len(r['physical'])} | {len(r['negative'])} | {_f(res_max) if res_max is not None else '—'} | {_f(dmax) if dmax is not None else '—'} |")
    add("")
    add("- **Árbitro.** Em todo caso perigoso, as raízes físicas das eqs. 4.2–4.3 e as soluções geométricas são as mesmas, uma a uma.")
    add("  - Os resíduos e desvios ficam no arredondamento de 40 dígitos.")
    add("  - Também se cumprem a regra de sinal sign(β) = sign(d) = s (eq. 4.6c) e a distância a = PQ da eq. A6, que coincide com a posição medida da colisão.")
    add("  - Isso confere a transcrição das eqs. 4.2, 4.3, 4.4a, 4.4b, 4.5 e A6, e o modelo da travessia.")
    negs = [(arg, r["negative"]) for arg, r in res["r8"] if r["negative"]]
    where = "; ".join(f"(θ, φ) = ({arg[0].replace('.', ',')}; {arg[1].replace('.', ',')}): {len(n)}" for arg, n in negs) or "nenhum caso"
    add(f"- **Raízes que a eq. 4.6 não lista,** com alguma das distâncias a, a′ ou c do Apêndice A negativa: {where}.")
    add("  - A geometria não a realiza, porque nenhuma bola percorre distância negativa ao longo do caminho.")
    add("  - A lista da eq. 4.6 cobre as raízes introduzidas pelas manipulações do Apêndice A (p. 1097). Esta tem outra origem.")
    empty = [(arg, r) for arg, r in res["r8"] if not r["rows"]]
    for arg, r in empty:
        add(f"- **Caso ({'; '.join(a.replace('.', ',') for a in arg)}):** a trajetória é perigosa, porque o movimento livre se sobrepõe a si mesmo (distância mínima {_f(r['free_min_distance'])} < 2R). Mesmo assim, **nenhuma solução de uma colisão foi encontrada** sob os critérios abaixo:")
        add("  - Colisão antes da travessia, com 1 travessia entre as visitas: grades `ekt_*` padrão, sem raiz das eqs. 4.2–4.3 e sem candidato geométrico.")
        for _, x in res["r8x"]:
            kind = "Colisão antes da travessia" if x["kind"] == "before" else "Colisão depois de uma travessia"
            trav = "1 travessia" if x["k"] == 1 else f"{x['k']} travessias"
            extra = (", grades dobradas (ν, mira e s)" if x["refined"] else "") + (f", λ₁ = {x['lam1'].replace('.', ',')}" if x["lam1"] != arg[2] else "")
            add(f"  - {kind}, {trav} entre as visitas{extra}: {x['ok']} soluções, {len(x['rejected'])} candidatos rejeitados, {x['unresolved']} intervalos não resolvidos.")
        add("  - Isto **não** é declarado paradoxo nem multiplicidade zero. Soluções com várias colisões, que o artigo menciona (p. 1080), estão fora do modelo buscado, e a grade é finita. É o primeiro alvo do mapa (pergunta Q1).")
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
    add(f"- **Classe I.** Para φ entre 0,7 e 0,9, a geometria acompanha a eq. B32 com desvio relativo máximo de {_f(max(near), 2)}. Longe de π/4 ela tende às eqs. 4.14 e 4.15, cada uma do seu lado.")
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
    add("  - As soluções que o artigo constrói (Figs. 6 e 7; eqs. 3.1, 3.2, 4.2–4.3, 4.12–4.15, B32) saem de um modelo geométrico exato e independente, como limites B → 0 onde o artigo usa B ≪ 1.")
    add("  - Na seção III B (Fig. 9), há uma solução para cada n de 1 a 6 no primeiro conjunto de R7, com a eq. 3.7 conferindo só assintoticamente. O segundo conjunto é inconclusivo.")
    add("  - As eqs. 4.2–4.3 são satisfeitas pela geometria no arredondamento de 40 dígitos, em todos os casos com solução.")
    add("- ❌ **Não estabelece:**")
    add("  - Nenhuma afirmação universal. As contagens valem para as grades, os itinerários e as janelas declarados.")
    add("  - Nada sobre soluções com várias colisões.")
    add("  - **A fase é clássica.** Ela testa a alternativa de autoconsistência de Novikov, não o mecanismo quântico de Hawking.")
    add("  - O mapa sistemático e a certificação intervalar ainda não foram feitos.")
    add("  - O teste de escala que confirmaria a observação deste projeto, de que só v₁T_d/D importa (a observação não é do artigo), ainda não foi feito.")
    add("")
    return L
