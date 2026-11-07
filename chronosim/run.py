"""Phase 1: validate Li & Gott (1998) eq. 12 against the eq. 8 image sum.
Phase 2 (in progress): Kim & Thorne (1991), arbiters of eq. 4a and items F1 to F7.

    python -m chronosim.run            # writes outputs/
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import mpmath as mp

from . import estimates
from . import kim_thorne as kt
from . import thermal
from . import unruh
from . import wormhole as wh
from .config import DPS
from .misner import closed_form_stress, hadamard_image, self_consistent_b
from .pointsplit import misner_stress, rindler_part
from .provenance import reported_heuristics

COMPONENTS = ("rho", "p_xi", "p_y", "p_z")
GRID_XI = ("0.5", "1", "2")
GRID_B = ("1", "3", "2pi", "10", "30")
SCALING_XI = ("1e-3", "1e-2", "1e-1", "1")
SCAN_B = ("1", "2", "3", "4", "5", "2pi", "7", "8", "10", "15", "20")
ARBITER_XI = ("0.5", "1.3", "2")
F1_BETA = ("0.1", "0.5", "0.9")
F1_N = (1, 2, 3, 5, 10)
F1_EPS = ("1e-2", "1e-4", "1e-6", "1e-8")
F3_BETA = ("0.1", "0.5", "0.9")
F3_B = ("0.01", "0.1")
F3_X = ("0.25", "0.5")
F3_DELTA = ("1e-2", "1e-4", "1e-6")
F3_RHO = ("1e-1", "1e-2", "1e-3")
F3_RHO_DELTA = ("1e-4", "1e-6")
SLOPE_DELTA = ("1e-3", "1e-4", "1e-5", "1e-6")
THERMAL_T = ("0.3", "1", "3")
MODE_POINTS = (("0.1", "0.2"), ("0.3", "0.5"), ("0.7", "1.5"))
MODE_N = 400
F4_BETA = ("0.1", "0.5", "0.9")
F4_B = ("0.01", "0.1")
F4_THETA = ("0.2", "0.7", "1.3")
F4_DELTA = ("1e-2", "1e-4", "1e-6")
F2_B = ("1e-2", "1e-4", "1e-6")
F2_N = (1, 2, 3, 5)
F5_BETA = ("0.1", "0.5", "0.9")
F5_N = (1, 2, 5, 20)
F6_UNITS = {"D": "1000", "b": "10", "X": "250"}  # Planck units, l_P = m_P = 1
F6_DT54 = ("1", "0.1", "0.01")
F6_THICKNESS = ("1e-2", "1e-4", "1e-6")
FIG_DELTA = ("1e-2", "3e-3", "1e-3", "3e-4", "1e-4", "3e-5", "1e-5")


def _num(s):
    return 2 * mp.pi if s == "2pi" else mp.mpf(s)


def validate(rows_out):
    """Compare numerical and closed-form components on the grid. Returns worst error."""
    worst = mp.mpf(0)
    for xs in GRID_XI:
        for bs in GRID_B:
            xi, b = _num(xs), _num(bs)
            num, cf = misner_stress(xi, b), closed_form_stress(xi, b)
            scale = 1 / (1440 * mp.pi ** 2 * xi ** 4)  # size of the Rindler term; nonzero even at b = 2 pi
            err = max(abs(num[k].value - cf[k].value) for k in COMPONENTS) / scale
            worst = max(worst, err)
            rows_out.append({
                "xi": xs, "b": bs,
                "rho_computed": mp.nstr(num["rho"].value, 20),
                "rho_closed_form": mp.nstr(cf["rho"].value, 20),
                "max_err_over_scale": mp.nstr(err, 3),
                "trace_over_scale": mp.nstr(abs(num["trace"].value) / scale, 3),
                "max_offdiag_over_scale": mp.nstr(num["max_offdiag"].value / scale, 3),
                "kummer_tail_fraction": mp.nstr(num["tail_fraction"].value, 3),
                "decay_spread": mp.nstr(num["decay_spread"].value, 3),
            })
    return worst


def scaling(b=mp.mpf(1)):
    """log-log slope of |rho| versus xi from COMPUTED points."""
    pts = [(mp.log(_num(x)), mp.log(abs(misner_stress(_num(x), b)["rho"].value))) for x in SCALING_XI]
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    slope = sum((p[0] - mx) * (p[1] - my) for p in pts) / sum((p[0] - mx) ** 2 for p in pts)
    return slope, pts


def arbiter():
    """Kim & Thorne eqs. 3, 4 and 20 on the Rindler vacuum, against the phase 1 operator."""
    rows = []
    F = wh.regularized_with_eq3(lambda P, Q: hadamard_image(0, P, Q, 1))
    for xs in ARBITER_XI:
        xi = mp.mpf(xs)
        scale = 1 / (1440 * mp.pi ** 2 * xi ** 4)
        T, R = wh.eq4_stress(F, [0, xi, 0, 0]), rindler_part([0, xi, 0, 0])
        rows.append({
            "xi": xs,
            "rho_eq4_over_scale": T[0][0] / scale,
            "rho_phase1_over_scale": R[0][0] / scale,
            "worst_ratio_minus_2": max((T[i][i] / R[i][i] - 2 for i in range(4)), key=abs),
        })
    return rows


def f1_geometry():
    """Covering space built from eqs. 25 and 28, against eqs. 31, 33, 34, 37 and 38."""
    D = mp.mpf(1)
    rows, exact = [], {"eq33": mp.mpf(0), "eq34": mp.mpf(0), "eq37": mp.mpf(0), "fixed_point": mp.mpf(0), "eq31": mp.mpf(0)}
    for bs in F1_BETA:
        beta = mp.mpf(bs)
        xi = kt.doppler_xi(beta)
        fixed = D / (1 - xi)
        T, X = wh.throat_map((fixed, fixed), D, beta)
        exact["fixed_point"] = max(exact["fixed_point"], abs(T - fixed), abs(X - fixed))
        event = (mp.mpf("0.123"), mp.mpf("0.456"))
        U0, V0 = kt.null_coords(kt.to_covering(event, D, xi))
        for N in F1_N:
            U, V = kt.null_coords(kt.to_covering(wh.copy_event(event, N, D, beta), D, xi))
            Un, Vn = kt.copy_n(U0, V0, N, xi)
            exact["eq33"] = max(exact["eq33"], abs(U / Un - 1), abs(V / Vn - 1))
            for es in F1_EPS:
                eps = mp.mpf(es)
                x, xp = (-eps, 2 * eps), (-2 * eps, eps + eps / 3)
                sigma = wh.exact_sigma_N(N, x, xp, D, beta).value
                Ux, Vx = kt.null_coords(kt.to_covering(wh.copy_event(x, N, D, beta), D, xi))
                Up, Vp = kt.null_coords(kt.to_covering(xp, D, xi))
                exact["eq34"] = max(exact["eq34"], abs((Vx - Vp) * (Ux - Up) / sigma - 1))
                exact["eq37"] = max(exact["eq37"], abs((Ux - Up) / kt.delta_U(N, x, xp, xi) - 1))
                rows.append({
                    "beta": bs, "N": N, "eps_over_D": es,
                    "sigma_exact": mp.nstr(sigma, 20),
                    "sigma_eq38": mp.nstr(kt.sigma_N(N, x, xp, D, beta).value, 20),
                    "eq38_rel_err": mp.nstr(abs(kt.sigma_N(N, x, xp, D, beta).value / sigma - 1), 4),
                })
        for n in (1, 2, 3):
            for tau in ("-1", "0.7"):
                Tn, Xn = kt.to_covering(wh.copy_event((mp.mpf(tau), mp.mpf(0)), n, D, beta), D, xi)
                c = xi ** n * D / (1 - xi)
                coefficient = (1 + xi ** (2 * n)) / (1 - xi ** (2 * n))
                exact["eq31"] = max(exact["eq31"], abs((Tn + c) / (Xn + c) / coefficient + 1))
    return rows, exact


def thermal_arbiter():
    """Thermal state in Minkowski space, with both operators, against Stefan-Boltzmann."""
    rows = []
    for Ts in THERMAL_T:
        T = mp.mpf(Ts)
        ref = unruh.thermal_density_closed(T)
        num = thermal.thermal_stress(T)
        finer = thermal.thermal_stress(T, split_ratio="1e-11", step_ratio="1e-12")
        for name in ("eq4a", "with_half"):
            d = num[name]
            rho = d["rho"].value
            rows.append({
                "T": Ts, "operator": name,
                "rho": rho,
                "rho_over_stefan_boltzmann": rho / ref,
                "worst_3p_over_rho_minus_1": max((3 * q.value / rho - 1 for q in d["p"]), key=abs),
                "trace_over_rho": d["trace"].value / rho,
                "max_offdiag_over_rho": d["max_offdiag"].value / abs(rho),
                "fd_err": abs(rho / finer[name]["rho"].value - 1),
                "decay_spread": d["decay_spread"].value,
                "tail_fraction": d["tail_fraction"].value,
            })
    return rows


def thermal_modes():
    """The image sum against the mode integral with Bose occupation numbers, at separated points."""
    beta = mp.mpf(1)
    rows = []
    for ts, rs in MODE_POINTS:
        t, r = mp.mpf(ts), mp.mpf(rs)
        partial = thermal.hadamard_regularized([t, r, 0, 0], [0, 0, 0, 0], beta, MODE_N)
        tail = (mp.zeta(2) - sum(mp.mpf(n) ** -2 for n in range(1, MODE_N + 1))) / (mp.pi ** 2 * beta ** 2)
        modes = thermal.hadamard_mode_integral(t, r, beta).value
        rows.append({"t": ts, "r": rs, "images": partial + tail, "modes": modes,
                     "rel_diff": abs((partial + tail) / modes - 1)})
    return rows


def unruh_arbiter():
    """Rindler vacuum from each operator, against the value from the Unruh temperature."""
    rows = []
    F = wh.regularized_with_eq3(lambda P, Q: hadamard_image(0, P, Q, 1))
    for xs in ARBITER_XI:
        xi = mp.mpf(xs)
        ref = unruh.rindler_density(xi).value
        rows.append({
            "xi": xs,
            "ref_over_scale": ref * 1440 * mp.pi ** 2 * xi ** 4,
            "with_half_over_ref": rindler_part([0, xi, 0, 0])[0][0] / ref,
            "eq4a_over_ref": wh.eq4_stress(F, [0, xi, 0, 0])[0][0] / ref,
        })
    return rows


def _f3_row(bs, bb, Xs, ds, rs):
    """One F3 point. rho = rs * X, with T a distance delta before H_1 in lambda_1 units."""
    D, b, X, delta = mp.mpf(1), mp.mpf(bb), mp.mpf(Xs), mp.mpf(ds)
    xi = kt.doppler_xi(mp.mpf(bs))
    rho = mp.mpf(rs) * X
    T = X - delta + (rho ** 2 / (2 * X) + rho ** 2 * xi / (2 * (D - X))) / (1 - xi)
    num = wh.stress_between_numeric(T, X, rho, D, b, bs, delta)
    finer = wh.stress_between_numeric(T, X, rho, D, b, bs, delta, split_ratio="1e-11", step_ratio="1e-12")
    S = kt.stress_between_kk(T, X, rho, D, b, xi, num["n_terms"]).value
    comps = ("TT", "TX", "XX")
    return {
        "beta": bs, "b_over_D": bb, "X_over_D": Xs, "delta_over_D": ds, "rho_over_X": rs,
        "eq67_kk": mp.nstr(S, 15),
        "eq4a_TT_over_eq67": mp.nstr(num["eq4a"]["TT"].value / S, 15),
        "with_half_TT_over_eq67": mp.nstr(num["with_half"]["TT"].value / S, 15),
        "eq4a_worst_dev": mp.nstr(max(abs(num["eq4a"][c].value / S - 1) for c in comps), 3),
        "with_half_worst_dev_from_half": mp.nstr(max(abs(num["with_half"][c].value / S - mp.mpf(1) / 2) for c in comps), 3),
        "rest_over_eq67": mp.nstr(num["eq4a"]["rest"].value / abs(S), 3),
        "fd_err": mp.nstr(abs(num["eq4a"]["TT"].value / finer["eq4a"]["TT"].value - 1), 3),
        "n_terms": num["n_terms"],
    }


def f3_stress():
    """T^mu^nu between the mouths from eq. 66, with both operators, against eq. 67."""
    rows = [_f3_row(bs, bb, Xs, ds, "0") for bs in F3_BETA for bb in F3_B for Xs in F3_X for ds in F3_DELTA]
    rows += [_f3_row("0.5", "0.01", "0.25", ds, rs) for rs in F3_RHO for ds in F3_RHO_DELTA]
    return rows


def f3_slope():
    """log-log slope of |T^TT| against the distance to H_N on the axis, for each operator."""
    D, b, X = mp.mpf(1), mp.mpf("0.01"), mp.mpf("0.25")
    out = {}
    for name in ("eq4a", "with_half"):
        pts = [
            (mp.log(mp.mpf(d)), mp.log(abs(wh.stress_between_numeric(X - mp.mpf(d), X, mp.mpf(0), D, b, "0.5", mp.mpf(d))[name]["TT"].value)))
            for d in SLOPE_DELTA
        ]
        mx = sum(p[0] for p in pts) / len(pts)
        my = sum(p[1] for p in pts) / len(pts)
        out[name] = sum((p[0] - mx) * (p[1] - my) for p in pts) / sum((p[0] - mx) ** 2 for p in pts)
    return out


def _f4_point(bs, bb, ths, ds, phis="0"):
    """A point a distance delta = ds * b before H_1 on the throat, in units of lambda_1."""
    D, b, theta, phi = mp.mpf(1), mp.mpf(bb), mp.mpf(ths), mp.mpf(phis)
    xi = kt.doppler_xi(mp.mpf(bs))
    delta = mp.mpf(ds) * b
    T = kt.hypersurface_time_throat(1, theta, b, xi) - delta
    return D, b, theta, phi, xi, delta, T


def _f4_dev(num, D, b, theta, phi, xi, T, variant, name):
    target = mp.mpf(1) if name == "eq4a" else mp.mpf(1) / 2
    E = kt.stress_throat(T, theta, phi, D, b, xi, num["n_terms"], variant)
    scale = max(abs(E[m][c].value) for m in range(4) for c in range(4))
    return max(abs(num[name][m][c].value - target * E[m][c].value) for m in range(4) for c in range(4)) / (target * scale)


def _f4_row(bs, bb, ths, ds, phis="0"):
    D, b, theta, phi, xi, delta, T = _f4_point(bs, bb, ths, ds, phis)
    num = wh.stress_throat_numeric(kt.throat_point(T, theta, phi, b), D, b, bs, delta)
    finer = wh.stress_throat_numeric(kt.throat_point(T, theta, phi, b), D, b, bs, delta, split_ratio="1e-11", step_ratio="1e-12")
    E = kt.stress_throat(T, theta, phi, D, b, xi, num["n_terms"], "transverse")
    return {
        "beta": bs, "b_over_D": bb, "theta": ths, "phi": phis, "delta_over_b": ds,
        "eq51_TT": mp.nstr(E[0][0].value, 15),
        "eq4a_TT_over_eq51": mp.nstr(num["eq4a"][0][0].value / E[0][0].value, 15),
        "with_half_TT_over_eq51": mp.nstr(num["with_half"][0][0].value / E[0][0].value, 15),
        "eq4a_dev_transverse": mp.nstr(_f4_dev(num, D, b, theta, phi, xi, T, "transverse", "eq4a"), 4),
        "with_half_dev_transverse": mp.nstr(_f4_dev(num, D, b, theta, phi, xi, T, "transverse", "with_half"), 4),
        "eq4a_dev_literal": mp.nstr(_f4_dev(num, D, b, theta, phi, xi, T, "literal", "eq4a"), 4),
        "with_half_dev_literal": mp.nstr(_f4_dev(num, D, b, theta, phi, xi, T, "literal", "with_half"), 4),
        "fd_err": mp.nstr(abs(num["eq4a"][0][0].value / finer["eq4a"][0][0].value - 1), 3),
        "n_terms": num["n_terms"],
    }


def f4_stress():
    """T^mu^nu at the throat from eq. 50, with both operators and both denominators of eq. 52, against eq. 51."""
    rows = [_f4_row(bs, bb, ths, ds) for bs in F4_BETA for bb in F4_B for ths in F4_THETA for ds in F4_DELTA]
    rows += [_f4_row("0.5", "0.01", "0.7", ds, "0.6") for ds in F4_DELTA]
    return rows


def f4_slope():
    """log-log slope of |T^TT| against the distance to H_1 on the throat, for each operator."""
    out = {}
    for name in ("eq4a", "with_half"):
        pts = []
        for ds in ("1e-3", "1e-4", "1e-5", "1e-6"):
            D, b, theta, phi, xi, delta, T = _f4_point("0.5", "0.01", "0.7", ds)
            num = wh.stress_throat_numeric(kt.throat_point(T, theta, phi, b), D, b, "0.5", delta)
            pts.append((mp.log(delta), mp.log(abs(num[name][0][0].value))))
        mx = sum(p[0] for p in pts) / len(pts)
        my = sum(p[1] for p in pts) / len(pts)
        out[name] = sum((p[0] - mx) * (p[1] - my) for p in pts) / sum((p[0] - mx) ** 2 for p in pts)
    return out


def f2_optics():
    """Ray-transfer amplitudes against eqs. 47 and 64 (via eq. 24), and the exact lapse against eq. 57."""
    D, X = mp.mpf(1), mp.mpf("0.25")
    xi = kt.doppler_xi(mp.mpf("0.5"))
    rows = []
    for bb in F2_B:
        b = mp.mpf(bb)
        rows.append({
            "b_over_D": bb,
            "throat": [wh.ray_amplitude_throat(N, D, b).value / kt.psi_throat(N, D, b) - 1 for N in F2_N],
            "between": [wh.ray_amplitude_between(N, X, D, b).value
                        / (kt.vleck_between(N, X, D, b, xi) / kt.zeta_between(N, X, D, xi)) - 1 for N in F2_N],
        })
    lapse = max(abs(wh.exact_lapse_on_horizon(N, X, D, bs) / kt.zeta_between(N, X, D, kt.doppler_xi(mp.mpf(bs))) - 1)
                for bs in F5_BETA for N in F2_N)
    return rows, lapse


def f5_hypersurfaces():
    """Roots of the collapsed lambda_N against eqs. 43 and 62, and N -> infinity against eqs. 44a and 63."""
    D, b, X = mp.mpf(1), mp.mpf("0.01"), mp.mpf("0.25")
    out = {"eq43": mp.mpf(0), "eq62": mp.mpf(0), "eq44a": [], "eq63": []}
    for bs in F5_BETA:
        xi = kt.doppler_xi(mp.mpf(bs))
        for th in ("0.2", "1.3"):
            theta = mp.mpf(th)
            for N in F5_N:
                out["eq43"] = max(out["eq43"], abs(wh.root_HN_throat(N, theta, b, bs).value / kt.time_HN_throat(N, theta, b, xi) - 1))
        for rs in ("1e-2", "1e-1"):
            rho = mp.mpf(rs) * X
            for N in F5_N:
                T = kt.time_HN_between(N, X, rho, D, xi)
                out["eq62"] = max(out["eq62"], abs((wh.root_HN_between(N, X, rho, D, bs).value - T) / (T - X)))
        theta, rho = mp.mpf("0.7"), mp.mpf("0.1") * X
        out["eq44a"].append((bs, abs(wh.root_HN_throat(60, theta, b, bs).value / kt.time_Hminus_throat(theta, b) - 1)))
        out["eq63"].append((bs, abs((wh.root_HN_between(60, X, rho, D, bs).value - kt.time_Hminus_between(X, rho)) / (rho ** 2 / (2 * X)))))
    return out


def f6_momentum():
    """Quadratures of eqs. 51 and 67 against eqs. 53, 54, 68 and 69, in Planck units."""
    D, b, X = (mp.mpf(F6_UNITS[k]) for k in ("D", "b", "X"))
    bs = "0.5"
    xi = kt.doppler_xi(mp.mpf(bs))
    out = {}
    out["eq69"] = max(
        abs(wh.momentum_per_area_numeric(N, X, mp.mpf(rs) * X, D, b, b_s).value
            / kt.momentum_per_area_N(N, X, D, b, kt.doppler_xi(mp.mpf(b_s))).value - 1)
        for b_s in F5_BETA for N in (1, 3) for rs in ("0", "1e-3", "1e-2")
    )
    numeric = wh.momentum_plane_numeric(mp.mpf(1), X, D, b, bs, orders=[1, 2]).value
    printed = sum((b / (2 * D)) ** N * (1 + 4 * xi ** N + xi ** (2 * N)) / (1 - xi ** N) ** 2 for N in (1, 2)) / (-6 * mp.pi)
    out["eq68"] = numeric / printed - 1
    single = -(b / (2 * D)) * (1 + 4 * xi + xi ** 2) / (1 - xi) ** 2 / (6 * mp.pi)
    out["eq54"] = [
        (dt, wh.momentum_throat_numeric(mp.mpf(dt), D, b, bs, lambda th: kt.time_Hminus_throat(th, b), orders=[1]).value
         / (single / mp.mpf(dt)) - 1)
        for dt in F6_DT54
    ]
    values = [
        wh.momentum_throat_numeric(mp.mpf(1), D, b, bs, lambda th: kt.time_HN_throat(1, th, b, xi), orders=[1], thickness=t).value
        for t in F6_THICKNESS
    ]
    out["eq53_ratio"] = [(t, v / kt.momentum_throat_N(1, mp.mpf(1), mp.mpf(t) * b, D, b, xi).value) for t, v in zip(F6_THICKNESS, values)]
    coefficient = kt.momentum_throat_N(1, mp.mpf(1), b / mp.e, D, b, xi).value  # ln(b/a) = 1
    out["eq53_slope"] = (values[2] - values[1]) / (mp.log(mp.mpf(F6_THICKNESS[1])) - mp.log(mp.mpf(F6_THICKNESS[2]))) / coefficient - 1
    return out


def f7_estimates():
    """Side-by-side arithmetic of the order-of-magnitude estimates."""
    return [
        ("l_P / D (D = 1 m)", "l_P / D", estimates.PRINTED_LP_OVER_D, estimates.lp_over_D()),
        ("(δL/L)_max, eq. 13, N − 1 = 0", "l_P / D", estimates.PRINTED_MAX_STRAIN, estimates.max_strain()),
        ("densidade de Planck m_P/l_P³", "m_P / l_P³", estimates.PRINTED_PLANCK_DENSITY_HAWKING, estimates.planck_density()),
        ("T⁰⁰_max, eq. 14", "(l_P/D)(m_P/l_P³)", estimates.PRINTED_MAX_DENSITY, estimates.max_density()),
        ("Δt do corte de Hawking (δg^VP = 1, eq. 9c)", "l_P² / D, em l_P", estimates.PRINTED_HAWKING_DT, estimates.hawking_cutoff_dt()),
    ]


def figures_phase2(out: Path):
    """Phase 2 figures: the divergence with both operators, and the two proposed cutoffs."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    blue, orange, ink, muted = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
    D, b, X = mp.mpf(1), mp.mpf("0.01"), mp.mpf("0.25")
    xi = kt.doppler_xi(mp.mpf("0.5"))

    # 3. Divergence near H_1 between the mouths, on the axis.
    fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
    deltas = [mp.mpf(d) for d in FIG_DELTA]
    nums = [wh.stress_between_numeric(X - d, X, mp.mpf(0), D, b, "0.5", d) for d in deltas]
    grid = np.logspace(-5, -2, 200)
    n_terms = nums[0]["n_terms"]
    curve = [float(abs(kt.stress_between_kk(X - mp.mpf(g), X, mp.mpf(0), D, b, xi, n_terms).value)) for g in grid]
    ax.loglog(grid, curve, color=muted, lw=1.2, ls="--", label="eq. 67 impressa (DERIVED)")
    ax.loglog([float(d) for d in deltas], [float(abs(n["eq4a"]["TT"].value)) for n in nums], "o", color=blue, ms=8, mfc="none", mew=2,
              label="eq. 4a impressa, sem ½ (COMPUTED)")
    ax.loglog([float(d) for d in deltas], [float(abs(n["with_half"]["TT"].value)) for n in nums], "s", color=orange, ms=8, mfc="none", mew=2,
              label="com ½, convenção do árbitro térmico (COMPUTED)")
    ax.set_xlabel("Δ  (distância a H₁, unidades de D)", color=ink)
    ax.set_ylabel("|T^TT|", color=ink)
    ax.set_title("Kim & Thorne: divergência Δ⁻³; os dois operadores diferem por 2", color=ink, fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.25, which="both")
    fig.savefig(out / "fig3_kt_divergence.png", dpi=130)
    plt.close(fig)

    # 4. Metric perturbation, eq. 9c with N - 1 = 0 and D = 1 m, and the two cutoffs.
    fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
    lpD = float(estimates.lp_over_D().value)
    dt = np.logspace(-40, 3, 300)
    ax.loglog(dt, lpD / dt, color=blue, lw=2, label="δg^VP ~ (l_P/D)(l_P/Δt), eq. 9c")
    ax.axvline(1, color=ink, lw=1.2, ls="--")
    ax.text(0.3, 1e-15, "corte de Kim & Thorne\nΔt ~ l_P", color=ink, fontsize=8, ha="right")
    ax.axvline(float(estimates.hawking_cutoff_dt().value), color=ink, lw=1.2, ls=":")
    ax.text(3e-35, 1e-30, "corte de Hawking\nDΔt ~ l_P²", color=ink, fontsize=8)
    ax.axhline(1, color=muted, lw=0.8)
    ax.set_xlabel("Δt / l_P  (tempo próprio até H₁)", color=ink)
    ax.set_ylabel("δg^VP  (ordem de grandeza)", color=ink)
    ax.set_title("D = 1 m: δg^VP ~ 10⁻³⁵ no corte de Kim & Thorne, ~ 1 no de Hawking", color=ink, fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    ax.grid(alpha=0.25, which="both")
    fig.savefig(out / "fig4_kt_cutoffs.png", dpi=130)
    plt.close(fig)


def figures(out: Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    # 1. Divergence at the chronology horizon
    fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
    xs = np.logspace(-3, 0.5, 200)
    for bs, color in (("1", "#1f4e79"), ("3", "#2e86c1"), ("10", "#c0392b")):
        b = _num(bs)
        C = float(abs(3 * ((2 * mp.pi / b) ** 4 - 1) / (1440 * mp.pi ** 2)))
        ax.loglog(xs, C / xs ** 4, color=color, lw=1.4, label=f"b = {bs}  (eq. 12, DERIVED)")
        pts = [_num(x) for x in SCALING_XI]
        ax.loglog([float(p) for p in pts], [float(abs(misner_stress(p, b)["rho"].value)) for p in pts],
                  "o", color=color, mfc="none", ms=7)
    ax.plot([], [], "o", color="gray", mfc="none", label="soma de imagens (COMPUTED)")
    ax.set_xlabel("ξ  (distância ao horizonte da cronologia)")
    ax.set_ylabel("|ρ|  densidade de energia do vácuo")
    ax.set_title("Espaço de Misner: |ρ| ∝ ξ⁻⁴ perto do horizonte")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.savefig(out / "fig1_horizon_divergence.png", dpi=130)
    plt.close(fig)

    # 2. The self-consistent value b = 2 pi
    fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
    bb = np.linspace(0.8, 20, 400)
    ax.plot(bb, 3 * ((2 * np.pi / bb) ** 4 - 1) / (1440 * np.pi ** 2), color="#1f4e79", lw=1.4,
            label="ξ⁴ρ, eq. 12 (DERIVED)")
    pts = [_num(s) for s in SCAN_B]
    ax.plot([float(p) for p in pts], [float(misner_stress(1, p)["rho"].value) for p in pts], "o",
            color="#c0392b", mfc="none", ms=7, label="soma de imagens, ξ = 1 (COMPUTED)")
    ax.axhline(0, color="black", lw=0.8)
    ax.axvline(float(self_consistent_b().value), color="#27ae60", ls="--", lw=1, label="b = 2π (Li & Gott)")
    ax.set_yscale("symlog", linthresh=1e-6)
    ax.set_xlabel("b  (rapidez do boost que identifica os pontos)")
    ax.set_ylabel("ξ⁴ ρ")
    ax.set_title("A energia do vácuo se anula em b = 2π")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.savefig(out / "fig2_self_consistent_b.png", dpi=130)
    plt.close(fig)


def report(out: Path, rows, worst, slope, arb, f1_rows, f1_exact, unruh_rows, f3_rows, f3_slopes, thermal_rows, mode_rows, f4_rows, f4_slopes, f2, f5, f6, f7):
    h = reported_heuristics()
    lines = [
        "# chronosim · relatório de validação",
        "",
        "## Fase 1",
        "",
        "Comparação entre a energia do vácuo **calculada** a partir da soma de imagens "
        "(Li & Gott 1998, eq. 8) com point-splitting (eq. 11) e a **fórmula fechada** (eq. 12), "
        "para um campo escalar conformemente acoplado no espaço de Misner. Unidades naturais (ħ = c = 1).",
        "",
        "### Resultado",
        "",
        f"- Maior diferença no grid ({len(rows)} pontos, 4 componentes cada), relativa à escala "
        f"1/(1440π²ξ⁴): **{mp.nstr(worst, 3)}** [COMPUTED]",
        f"- Inclinação de log|ρ| contra log ξ (b = 1, ξ de 10⁻³ a 1): **{mp.nstr(slope, 15)}** "
        "[COMPUTED]; a eq. 12 prevê exatamente −4",
        "",
        "| ξ | b | ρ calculado | ρ eq. 12 | erro / escala | traço / escala | cauda de Kummer |",
        "|---|---|---|---|---|---|---|",
        *[
            f"| {r['xi']} | {r['b']} | {r['rho_computed']} | {r['rho_closed_form']} | "
            f"{r['max_err_over_scale']} | {r['trace_over_scale']} | {r['kummer_tail_fraction']} |"
            for r in rows
        ],
        "",
        "### O que a fase 1 estabelece e o que não estabelece",
        "",
        "- **Estabelece:** a eq. 12 de Li & Gott decorre das eqs. 8 e 11, reproduzida de forma "
        "independente até a precisão numérica, incluindo a divergência ξ⁻⁴ no horizonte (o mecanismo "
        "invocado por Hawking, 1992) e o anulamento em b = 2π (o contra-exemplo de Li & Gott).",
        "- **Não estabelece:** que viagens no tempo sejam possíveis ou impossíveis. O cálculo é "
        "semiclássico (campo quântico em espaço-tempo fixo, sem reação gravitacional) e Kay, Radzikowski "
        "& Wald (1997) provaram que a função de dois pontos não pode ser Hadamard no horizonte de "
        "Cauchy; a questão depende de uma teoria de gravidade quântica.",
        "",
        "## Fase 2 (em andamento): Kim & Thorne (1991)",
        "",
        "Plano em `docs/phase2_plan.md`, transcrições em `docs/equation_transcription.md`. Até aqui: testes "
        "árbitros do operador e itens F1 a F7. Figuras: `fig3_kt_divergence.png` e `fig4_kt_cutoffs.png`. "
        "Checagem de literatura do fator 2: `docs/literature_check_factor2.md`.",
        "",
        "### Árbitro principal: estado térmico no espaço de Minkowski",
        "",
        "Não depende de nenhum artigo fechado. A função de Hadamard é a soma de imagens em tempo imaginário "
        "(período β = 1/T), derivada da condição KMS em `docs/thermal_images.md`, sem o termo n = 0, que é o "
        "vácuo de Minkowski. A referência é Stefan-Boltzmann para um escalar sem massa, π²T⁴/30, conferido pela "
        "integral de Bose. Valores COMPUTED.",
        "",
        "| T | operador | ρ / (π²T⁴/30) | maior desvio de 3p/ρ em relação a 1 | traço / ρ | erro DF |",
        "|---|---|---|---|---|---|",
        *[
            f"| {r['T']} | {'eq. 4a impressa' if r['operator'] == 'eq4a' else 'com ½'} | "
            f"{mp.nstr(r['rho_over_stefan_boltzmann'], 22)} | {mp.nstr(r['worst_3p_over_rho_minus_1'], 3)} | "
            f"{mp.nstr(r['trace_over_rho'], 3)} | {mp.nstr(r['fd_err'], 3)} |"
            for r in thermal_rows
        ],
        "",
        "Os valores coincidem entre as temperaturas porque δ e h escalam com β.",
        "",
        "**Resultado:** o operador com ½ reproduz Stefan-Boltzmann, e a eq. 4a impressa dá o dobro. O estado "
        "tem pressão ρ/3 e traço zero com os dois operadores. É uma inconsistência de convenção entre a eq. 4a "
        "impressa e o valor térmico.",
        "",
        "A soma de imagens também foi comparada com a integral de modos com ocupação de Bose, em pontos "
        f"separados (β = 1, {MODE_N} pares e o termo dominante da cauda):",
        "",
        "| t | r | imagens | modos | diferença relativa |",
        "|---|---|---|---|---|",
        *[
            f"| {r['t']} | {r['r']} | {mp.nstr(r['images'], 15)} | {mp.nstr(r['modes'], 15)} | {mp.nstr(r['rel_diff'], 3)} |"
            for r in mode_rows
        ],
        "",
        "A diferença vem do corte da soma de imagens. Birrell & Davies (1982) fica como referência pendente: a "
        "página não pôde ser conferida, e a derivação está escrita por extenso em `docs/`.",
        "",
        "### Comparação direta: eq. 4a de Kim & Thorne contra a eq. 11 de Li & Gott",
        "",
        "Operador das eqs. 4a e 4b de Kim & Thorne (p. 3933), como impressas, escrito sem usar o código da "
        "fase 1, aplicado à função de Rindler (Li & Gott, eq. 6) regularizada pela eq. 3 (p. 3933). A "
        "referência é o operador da eq. 11 de Li & Gott (arXiv gr-qc/9711074, p. 2), validado na fase 1. "
        "Escala: 1/(1440π²ξ⁴). Todos os valores são COMPUTED.",
        "",
        "| ξ | ρ eq. 4a / escala | ρ eq. 11 / escala | maior desvio da razão em relação a 2 |",
        "|---|---|---|---|",
        *[
            f"| {r['xi']} | {mp.nstr(r['rho_eq4_over_scale'], 20)} | {mp.nstr(r['rho_phase1_over_scale'], 20)} | "
            f"{mp.nstr(r['worst_ratio_minus_2'], 3)} |"
            for r in arb
        ],
        "",
        "**Resultado: dá o dobro.** Os dois artigos normalizam a função de Hadamard da mesma forma, "
        "1/(4π²σ) com σ igual à metade do quadrado da separação (Kim & Thorne, eqs. 3, 17 e 20; Li & Gott, "
        "eq. 9). A eq. 4a não tem o fator ½ da eq. 11, e nenhuma definição em volta dela o compensa. "
        "Inconsistência entre os dois artigos, nas convenções como impressas.",
        "",
        "### Árbitro secundário: efeito Unruh",
        "",
        "Referência: ρ = −π²T⁴/30 com T = 1/(2πξ), ou seja, −3/(1440π²ξ⁴). **Ressalva:** a relação \"vácuo de "
        "Rindler = −térmico\" não foi lida em Candelas & Deutsch (1977). Volume e página foram confirmados só no "
        "Crossref, e o artigo tem acesso fechado.",
        "",
        "| ξ | referência / escala | operador com ½ / referência | eq. 4a impressa / referência |",
        "|---|---|---|---|",
        *[
            f"| {r['xi']} | {mp.nstr(r['ref_over_scale'], 20)} | {mp.nstr(r['with_half_over_ref'], 20)} | "
            f"{mp.nstr(r['eq4a_over_ref'], 20)} |"
            for r in unruh_rows
        ],
        "",
        "**Resultado:** a versão com ½ concorda com o valor de Unruh, e a eq. 4a impressa dá o dobro, como no "
        "árbitro térmico. É uma inconsistência de convenção entre a eq. 4a impressa e o valor de Unruh.",
        "",
        "### F1: geometria do espaço de recobrimento",
        "",
        "Recobrimento construído por iteração do mapa de Poincaré das eqs. 25 e 28, sem usar as eqs. 29–38. "
        "Pontos a uma distância ε da garganta, com D = 1. Todos os valores são COMPUTED.",
        "",
        "| Afirmação do artigo | Maior desvio medido |",
        "|---|---|",
        f"| Origem do recobrimento em T = X = D/(1−ξ) (eq. 29) | {mp.nstr(f1_exact['fixed_point'], 3)} |",
        f"| Cópias U_n = ξ⁻ⁿU₀, V_n = ξⁿV₀ (eq. 33) | {mp.nstr(f1_exact['eq33'], 3)} |",
        f"| σ_N = ΔV·ΔU (eq. 34) | {mp.nstr(f1_exact['eq34'], 3)} |",
        f"| ΔU (eq. 37) | {mp.nstr(f1_exact['eq37'], 3)} |",
        "",
        "**Eq. 31 (p. 3940):** a eq. 31 como impressa é inconsistente com a Fig. 5 e com as eqs. 25, 28 e 29; a "
        "versão com o coeficiente de sinal oposto é consistente com as três. Nenhum resultado posterior "
        "depende desta equação.",
        "",
        "- **Leitura:** coeficiente positivo (1+ξ²ⁿ)/(1−ξ²ⁿ), nas duas transcrições e na verificação humana.",
        "- **Eqs. 25, 28 e 29:** a razão medida entre (T̃+c)/(X̃+c) e o coeficiente impresso é −1, com desvio "
        f"máximo de {mp.nstr(f1_exact['eq31'], 3)} (n = 1, 2, 3).",
        "- **Fig. 5:** as cópias n ≥ 1 são desenhadas com inclinação negativa, e a eq. 31 como impressa dá "
        "inclinação positiva para 0 < ξ < 1.",
        "- **Posição de T_0 e cruzamentos com o horizonte:** a construção reproduz X̃ = −D/(1−ξ) em n = 0, "
        "X̃ = −ξⁿD/(1−ξ) nos cruzamentos e os espaçamentos D, Dξ, Dξ² da Fig. 5. **Estes dois fatos não "
        "distinguem o sinal:** valem também para o coeficiente positivo, como um teste verifica. Só a "
        "inclinação distingue.",
        "",
        "Erro relativo da eq. 38 contra o intervalo exato. Ele cai linearmente com ε. O texto apresenta a eq. 35 "
        "como aproximação (\"very nearly the same\", p. 3941), justificada pela proximidade da garganta e do "
        "horizonte, mas **não fala em ordem em ε**. A primeira ordem é resultado desta medição e da expressão "
        "exata do erro, não uma afirmação do artigo.",
        "",
        "| β | N | ε = 1e-2 | ε = 1e-4 | ε = 1e-6 | ε = 1e-8 |",
        "|---|---|---|---|---|---|",
        *[
            f"| {bs} | {N} | " + " | ".join(
                r["eq38_rel_err"] for r in f1_rows if r["beta"] == bs and r["N"] == N
            ) + " |"
            for bs in F1_BETA for N in F1_N
        ],
        "",
        "### F3: ⟨T^μν⟩ entre as bocas",
        "",
        "Diferenças finitas aplicadas à eq. 66, com λ_N da eq. 61 e prefator lido literalmente, sem usar a eq. 67, "
        "com os dois operadores em paralelo: a eq. 4a impressa (`eq4_stress`, separação em y) e a eq. 4a com o "
        "fator ½ (`stress_tensor` da fase 1, separação em z). Δ é a distância a H_1, em unidades de λ_1. "
        "D = 1. Tudo é COMPUTED, salvo a coluna da eq. 67, que é DERIVED. A tabela completa está em "
        "`outputs/kim_thorne_f3.csv`.",
        "",
        "| β | b/D | X/D | ρ/X | Δ/D | sem ½ / eq. 67, T^TT | com ½ / eq. 67, T^TT | maior desvio sem ½ (TT, TX, XX) | demais componentes / eq. 67 |",
        "|---|---|---|---|---|---|---|---|---|",
        *[
            f"| {r['beta']} | {r['b_over_D']} | {r['X_over_D']} | {r['rho_over_X']} | {r['delta_over_D']} | "
            f"{r['eq4a_TT_over_eq67']} | {r['with_half_TT_over_eq67']} | {r['eq4a_worst_dev']} | {r['rest_over_eq67']} |"
            for r in f3_rows if r["X_over_D"] == "0.25"
        ],
        "",
        "- **Maior desvio do operador sem ½ em relação à eq. 67:** sobre o eixo, "
        f"{max(mp.mpf(r['eq4a_worst_dev']) for r in f3_rows if r['rho_over_X'] == '0')} (β = 0,9, Δ/D = 1e-2); "
        f"fora do eixo, {max(mp.mpf(r['eq4a_worst_dev']) for r in f3_rows if r['rho_over_X'] != '0')} (ρ/X = 0,1).",
        f"- **Maior erro de diferenças finitas em T^TT:** {max(mp.mpf(r['fd_err']) for r in f3_rows)}.",
        f"- **Inclinação de log|T^TT| contra log Δ, sobre o eixo (β = 0,5, b/D = 0,01, X/D = 0,25):** "
        f"−3 com desvio de {mp.nstr(f3_slopes['eq4a'] + 3, 3)} sem ½ e {mp.nstr(f3_slopes['with_half'] + 3, 3)} com ½.",
        "",
        "**Resultado:** Os resultados impressos de Kim & Thorne (eq. 67) são reproduzidos pelo operador da eq. 4a "
        "como impresso, sem o fator ½. Portanto o fator é carregado de forma consistente até o resultado final, "
        "e os valores impressos ficam 2 vezes acima dos obtidos com a convenção que reproduz o árbitro térmico. A "
        "taxa de divergência (Δ⁻³) e a troca de sinal na hipersuperfície não dependem desse fator.",
        "",
        "**Literatura** (`docs/literature_check_factor2.md`): não foi encontrada errata nem quem tenha apontado o "
        "fator. Klinkhammer (1992, tese, eqs. 3.3 e 3.8) e Visser (hep-th/9303023, eqs. 10 e 34) usam o operador "
        "com ½ e chegam à metade da eq. 67 (identidade verificada em `tests/test_literature_factor2.py`), mas "
        "declaram concordância com Kim & Thorne sem comentar o fator.",
        "",
        "O desvio residual tem a forma dos termos que o artigo despreza. Sobre o eixo, ele cai linearmente com Δ "
        "em X/D = 0,25. Em X = D/2, T^TT e T^TX ficam no nível do erro numérico e só T^XX se desvia, "
        "quadraticamente, porque a derivada primeira do prefator da eq. 66 se anula ali. Fora do eixo, o desvio "
        "em kk cresce como (ρ/X)², e as demais componentes crescem como ρ/X, a correção fracionária que o texto "
        "diz ignorar (p. 3945).",
        "",
        "### F4: ⟨T^μν⟩ na garganta",
        "",
        "Diferenças finitas aplicadas à eq. 50, com λ_N da eq. 46, sem usar a eq. 51, com os dois operadores em "
        "paralelo, em pontos da eq. 42 a uma distância Δ antes de H_1, em unidades de λ_1. A comparação usa as "
        "16 componentes do tensor, e o desvio é relativo à maior componente da eq. 51. Para a eq. 52, roda a "
        "forma impressa √(X²+Y²) e a forma √(Y²+Z²), o ρ da eq. 39. Tudo é COMPUTED, salvo as colunas da eq. 51, que são DERIVED. A tabela completa está em "
        "`outputs/kim_thorne_f4.csv`.",
        "",
        "| β | b/D | θ | φ | Δ/b | sem ½ / eq. 51, T^TT | com ½ / eq. 51, T^TT | desvio com √(Y²+Z²) | desvio com √(X²+Y²) |",
        "|---|---|---|---|---|---|---|---|---|",
        *[
            f"| {r['beta']} | {r['b_over_D']} | {r['theta']} | {r['phi']} | {r['delta_over_b']} | "
            f"{r['eq4a_TT_over_eq51']} | {r['with_half_TT_over_eq51']} | {r['eq4a_dev_transverse']} | {r['eq4a_dev_literal']} |"
            for r in f4_rows if r["b_over_D"] == "0.01" and (r["beta"] == "0.5" or r["delta_over_b"] == "1e-4")
        ],
        "",
        f"- **Maior desvio do operador sem ½, com √(Y²+Z²), em todo o grid:** "
        f"{max(mp.mpf(r['eq4a_dev_transverse']) for r in f4_rows)} (Δ/b = 1e-2); em Δ/b = 1e-6, "
        f"{max(mp.mpf(r['eq4a_dev_transverse']) for r in f4_rows if r['delta_over_b'] == '1e-6')}.",
        f"- **Maior desvio com a leitura literal √(X²+Y²):** {max(mp.mpf(r['eq4a_dev_literal']) for r in f4_rows)}; em "
        f"Δ/b = 1e-6, {max(mp.mpf(r['eq4a_dev_literal']) for r in f4_rows if r['delta_over_b'] == '1e-6')}. "
        "Não cai com Δ.",
        f"- **Maior erro de diferenças finitas em T^TT:** {max(mp.mpf(r['fd_err']) for r in f4_rows)}.",
        f"- **Inclinação de log|T^TT| contra log Δ (β = 0,5, b/D = 0,01, θ = 0,7):** −3 com desvio de "
        f"{mp.nstr(f4_slopes['eq4a'] + 3, 3)} sem ½ e {mp.nstr(f4_slopes['with_half'] + 3, 3)} com ½.",
        "",
        "**Resultado:** Os resultados impressos de Kim & Thorne (eq. 51) são reproduzidos pelo operador da eq. 4a "
        "como impresso, sem o fator ½. Portanto o fator é carregado de forma consistente até o resultado final, "
        "e os valores impressos ficam 2 vezes acima dos obtidos com a convenção que reproduz o árbitro térmico. A "
        "taxa de divergência (Δ⁻³) e a troca de sinal na hipersuperfície não dependem desse fator.",
        "",
        "**Eq. 52 (conferida):** a eq. 52 como impressa tem radicando X²+Y²; isso é inconsistente com a igualdade "
        "ρ̂ = ∂/∂ρ da própria equação, com a definição de ρ (eq. 39), com a descrição de ρ̂ como vetor unitário "
        "(p. 3943) e com a eq. 51, que só é reproduzida com √(Y²+Z²). O cálculo usa √(Y²+Z²). Com a forma "
        "impressa, a componente transversal de ρ̂ na garganta, para φ = 0, é cos(θ/2), e o desvio cresce com θ "
        "sem cair com Δ.",
        "",
        "O desvio com √(Y²+Z²) cai linearmente com Δ/b. Ele é o termo com a segunda derivada de λ_N nas direções "
        "transversais, que a eq. 51 não contém: a expressão analítica desse termo reproduz o resíduo em YY, YZ e "
        "ZZ, e as demais componentes ficam no nível do erro numérico. Por isso T^TT coincide com a eq. 51 nos 15 "
        "dígitos da tabela.",
        "",
        "### F2: desfocalização por óptica geométrica",
        "",
        "A garganta é tratada como a lente fina divergente de distância focal b/2 que o texto descreve (p. 3942), e "
        "a amplitude de uma fonte pontual é 1/B, com B o elemento da matriz de transferência. Isso usa a imagem "
        "do artigo, mas não as eqs. 47, 48 e 64. A ligação com Δ_N^(1/2) é a eq. 24. COMPUTED.",
        "",
        "| b/D | garganta, ψ/eq. 47 − 1 (N = 1, 2, 3, 5) | entre as bocas, X/D = 0,25, ψ/(eq. 64/ζ_N) − 1 (N = 1, 2, 3, 5) |",
        "|---|---|---|",
        *[
            f"| {r['b_over_D']} | {', '.join(mp.nstr(v, 4) for v in r['throat'])} | {', '.join(mp.nstr(v, 4) for v in r['between'])} |"
            for r in f2[0]
        ],
        "",
        f"As eqs. 47, 48 e 64 são o termo dominante em b/D: a correção é linear em b/D e cresce com N. A eq. 57 "
        f"coincide com o lapso exato de V obtido do espaço de recobrimento sobre o horizonte, com desvio máximo de "
        f"{mp.nstr(f2[1], 3)}.",
        "",
        "### F5: hipersuperfícies polarizadas e horizonte de Cauchy",
        "",
        f"Raiz numérica de λ_N(x, x) = 0 (eqs. 46 e 61): coincide com a eq. 43 até {mp.nstr(f5['eq43'], 3)} e com a eq. 62 "
        f"até {mp.nstr(f5['eq62'], 3)} (relativo à distância de H_N a T = X). O limite N → ∞, avaliado em N = 60, se "
        "aproxima das eqs. 44a e 63 como ξ^N:",
        "",
        "| β | desvio de N = 60 contra a eq. 44a | desvio contra a eq. 63 |",
        "|---|---|---|",
        *[f"| {bs} | {mp.nstr(v, 3)} | {mp.nstr(w2, 3)} |" for (bs, v), (_, w2) in zip(f5["eq44a"], f5["eq63"])],
        "",
        "### F6: quadrimomento",
        "",
        "Quadratura numérica das eqs. 51 e 67 (validadas em F3 e F4) nos domínios escritos no artigo, em unidades "
        f"de Planck (l_P = m_P = 1), com D = {F6_UNITS['D']}, b = {F6_UNITS['b']}, X = {F6_UNITS['X']} e β = 0,5. Por integrar "
        "as eqs. 51 e 67 impressas, estes resultados carregam a convenção sem ½; com ½, cada um cai à metade. "
        "COMPUTED.",
        "",
        f"- **Eq. 69:** desvio máximo de {mp.nstr(f6['eq69'], 3)} em β = 0,1, 0,5, 0,9, N = 1 e 3 e ρ/X = 0, 1e-3, 1e-2; independe de ρ, "
        "como o texto afirma.",
        f"- **Eq. 68** (N = 1 e 2, ΔT = l_P): desvio de {mp.nstr(f6['eq68'], 3)}.",
        "- **Eq. 54** (N = 1): a integral numérica se aproxima da eq. 54 quando ΔT/b diminui: "
        + "; ".join(f"ΔT/b = {mp.nstr(mp.mpf(dt) / mp.mpf(F6_UNITS['b']), 2)}: desvio {mp.nstr(v, 4)}" for dt, v in f6["eq54"])
        + ". A eq. 54 é o termo dominante; o termo em tan θ (cortado em a/b = 1e-6) e o fator cos θ somam a correção.",
        "- **Eq. 53** (N = 1, ΔT = l_P): o coeficiente de ln(b/a) da integral numérica difere do da eq. 53 por "
        f"{mp.nstr(f6['eq53_slope'], 3)}. A razão total fica abaixo de 1 pelo termo finito, que a eq. 53 não inclui: "
        + "; ".join(f"a/b = {t}: {mp.nstr(v, 6)}" for t, v in f6["eq53_ratio"]) + ".",
        "",
        "### F7: conferência de aritmética de estimativas de ordem de grandeza",
        "",
        "**Isto confere só a aritmética entre os valores de entrada impressos e os números impressos, não se as "
        "estimativas estão corretas.** Entradas: l_P ~ 10⁻³³ cm (p. 3930), m_P ~ 10⁻⁵ g (p. 3943), D = 1 m (p. 3935).",
        "",
        "| Grandeza | Fórmula | Valor impresso | Recalculado |",
        "|---|---|---|---|",
        *[f"| {name} | {formula} | {mp.nstr(printed.value, 3)} ({printed.source}) | {mp.nstr(value.value, 3)} |"
          for name, formula, printed, value in f7],
        "",
        "A densidade máxima impressa depois da eq. 14 é 10⁴⁹ g/cm³; com a fórmula da eq. 14 e as entradas impressas, "
        "o valor recalculado é 10⁵⁹ g/cm³. As demais quantidades conferem. A densidade de Planck de 10⁹⁴ g/cm³ que "
        "Hawking imprime (1992, p. 603) coincide com m_P/l_P³ a partir das entradas de Kim & Thorne.",
        "",
        "### O que esta parte estabelece e o que não estabelece",
        "",
        "- **Estabelece:** as eqs. 33, 34 e 37 são exatas para a geometria das eqs. 25 e 28, e a eq. 38 é "
        "uma aproximação cujo erro relativo é proporcional a ε/D, sobre o eixo de simetria.",
        "- **Estabelece também:** as eqs. 67 e 51 decorrem das eqs. 61 e 66 e das eqs. 46 e 50 pelo operador da "
        "eq. 4a impressa, com os desvios medidos acima, e a divergência é Δ⁻³. As eqs. 43, 44a, 62, 63 e 69 "
        "decorrem das expressões de partida; as eqs. 47, 48, 53, 54, 64 e 68 são os termos dominantes dos cálculos "
        "numéricos correspondentes.",
        "- **Não estabelece:** que as eqs. 46, 50, 61 e 66 estejam certas: elas são entrada do cálculo, em óptica "
        "geométrica. A ancoragem pelos árbitros foi feita no estado térmico e no vácuo de Rindler, não no "
        "buraco de minhoca. O que os árbitros mostram é qual convenção reproduz esses valores. A leitura da eq. 52 "
        "aguarda conferência humana.",
        "- **Disputa com Hawking:** Hawking (1992, p. 603 e seção VI) contesta a conclusão de Kim & Thorne. "
        "Ele argumenta que um eventual corte deveria ser em distância invariante, e não no tempo próprio, "
        "que depende do observador. A `fig4_kt_cutoffs.png` mostra os dois cortes lado a lado. Esta fase "
        "reproduz o cálculo de Kim & Thorne e não decide quem estava certo.",
        "- **F7** confere aritmética de estimativas de ordem de grandeza; não testa as estimativas.",
        "",
        "## Constantes numéricas escolhidas por nós (HEURISTIC)",
        "",
        *[f"- `{k}` = {mp.nstr(v.value, 5)}: {v.source}" for k, v in h.items()],
        "",
    ]
    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="chronosim phase 1: Misner space validation")
    p.add_argument("--out", default="outputs")
    p.add_argument("--no-figures", action="store_true")
    args = p.parse_args(argv)
    mp.mp.dps = DPS.value  # all comparisons at full working precision
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    worst = validate(rows)
    with open(out / "validation.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    slope, _ = scaling()
    if not args.no_figures:
        figures(out)
    arb = arbiter()
    f1_rows, f1_exact = f1_geometry()
    with open(out / "kim_thorne_f1.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(f1_rows[0]))
        w.writeheader()
        w.writerows(f1_rows)
    unruh_rows = unruh_arbiter()
    f3_rows = f3_stress()
    with open(out / "kim_thorne_f3.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(f3_rows[0]))
        w.writeheader()
        w.writerows(f3_rows)
    f3_slopes = f3_slope()
    thermal_rows = thermal_arbiter()
    mode_rows = thermal_modes()
    f4_rows = f4_stress()
    with open(out / "kim_thorne_f4.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(f4_rows[0]))
        w.writeheader()
        w.writerows(f4_rows)
    f4_slopes = f4_slope()
    f2 = f2_optics()
    f5 = f5_hypersurfaces()
    f6 = f6_momentum()
    f7 = f7_estimates()
    if not args.no_figures:
        figures_phase2(out)
    report(out, rows, worst, slope, arb, f1_rows, f1_exact, unruh_rows, f3_rows, f3_slopes, thermal_rows, mode_rows, f4_rows, f4_slopes, f2, f5, f6, f7)

    print(f"pior erro relativo à escala: {mp.nstr(worst, 3)} [COMPUTED]")
    print(f"inclinação perto do horizonte: {mp.nstr(slope, 15)} [COMPUTED] (eq. 12: -4)")
    print("fase 2, árbitro: eq. 4a de Kim & Thorne / eq. 11 de Li & Gott no vácuo de Rindler, maior desvio de 2: "
          f"{mp.nstr(max((abs(r['worst_ratio_minus_2']) for r in arb)), 3)} [COMPUTED]")
    print(f"fase 2, F1: eq. 33 {mp.nstr(f1_exact['eq33'], 3)}, eq. 34 {mp.nstr(f1_exact['eq34'], 3)}, "
          f"eq. 37 {mp.nstr(f1_exact['eq37'], 3)} [COMPUTED]")
    print(f"fase 2, F1: eq. 31 como impressa difere por um sinal, desvio de -1: {mp.nstr(f1_exact['eq31'], 3)} [COMPUTED]")
    print("fase 2, árbitro térmico, T = 1: com ½ / Stefan-Boltzmann = "
          f"{mp.nstr(thermal_rows[3]['rho_over_stefan_boltzmann'], 22)}; eq. 4a / Stefan-Boltzmann = "
          f"{mp.nstr(thermal_rows[2]['rho_over_stefan_boltzmann'], 22)}")
    print("fase 2, Unruh: com ½ / referência = "
          f"{mp.nstr(unruh_rows[0]['with_half_over_ref'], 20)}; eq. 4a / referência = {mp.nstr(unruh_rows[0]['eq4a_over_ref'], 20)} (xi = 0.5)")
    print(f"fase 2, F3: maior desvio sem ½ contra a eq. 67: {max(mp.mpf(r['eq4a_worst_dev']) for r in f3_rows)} [COMPUTED]")
    print(f"fase 2, F3: inclinação -3 com desvio {mp.nstr(f3_slopes['eq4a'] + 3, 3)} [COMPUTED]")
    print(f"fase 2, F4: maior desvio sem ½ contra a eq. 51, com sqrt(Y^2+Z^2): {max(mp.mpf(r['eq4a_dev_transverse']) for r in f4_rows)}; "
          f"com sqrt(X^2+Y^2): {max(mp.mpf(r['eq4a_dev_literal']) for r in f4_rows)} [COMPUTED]")
    print(f"fase 2, F4: inclinação -3 com desvio {mp.nstr(f4_slopes['eq4a'] + 3, 3)} [COMPUTED]")
    print(f"fase 2, F2: garganta, b/D = 1e-6, N = 5: {mp.nstr(f2[0][-1]['throat'][-1], 3)}; eq. 57 contra o lapso exato: {mp.nstr(f2[1], 3)}")
    print(f"fase 2, F5: eq. 43 {mp.nstr(f5['eq43'], 3)}, eq. 62 {mp.nstr(f5['eq62'], 3)} [COMPUTED]")
    print(f"fase 2, F6: eq. 69 {mp.nstr(f6['eq69'], 3)}, eq. 68 {mp.nstr(f6['eq68'], 3)}, eq. 53 (ln) {mp.nstr(f6['eq53_slope'], 3)}, "
          f"eq. 54 em ΔT/b = 1e-3: {mp.nstr(f6['eq54'][-1][1], 3)} [COMPUTED]")
    print(f"fase 2, F7: T00_max impresso 1e49 g/cm^3, recalculado {mp.nstr(f7[3][3].value, 3)} g/cm^3 (ordem de grandeza)")
    print("constantes HEURISTIC usadas:")
    for k, v in reported_heuristics().items():
        print(f"  {k} = {mp.nstr(v.value, 5)}  ({v.source})")
    print(f"relatório: {out / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
