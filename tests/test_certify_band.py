"""The direction machine near mouth 2 (certify.conditions_dir), the mouth band of the empty case and
the --certify split (user decisions (a)-(d) of 2026-11-09, docs/phase3_borda.md)."""

import math

import mpmath as mp
import pytest
from mpmath import iv

from chronosim import billiards as bl
from chronosim import certify as ct
from chronosim import phase3 as p3

EMPTY = p3.R8_CASES[2]
CONTROL = p3.BAND_CONTROLS[0]


def _aim_of(case, s, alpha):
    """Aim (polar angle on mouth 2) of the first point of the ray P + L (cos alpha, sin alpha)."""
    c = ct.inputs(case, xm=mp)
    q = ct.conditions_dir(s, 0, alpha, c, xm=mp)
    P = (c["p0"][0] + s * c["e"][0], c["p0"][1] + s * c["e"][1])
    pcs = (P[0] - 1, P[1])
    disc = q["b"] ** 2 - q["cc"]
    L0 = -q["b"] - mp.sqrt(disc)
    X = (P[0] + L0 * mp.cos(alpha), P[1] + L0 * mp.sin(alpha))
    return mp.atan2(X[1], X[0] - 1), pcs


@pytest.mark.parametrize("case", [EMPTY, CONTROL[0], p3.R8_CASES[0]])
def test_direction_conditions_are_the_loop_quantities_of_billiards(case):
    """miss, C x W and V of conditions_dir against billiards._loop_state, on a grid of the band."""
    with mp.workdps(50):
        c = ct.inputs(case, xm=mp)
        R = mp.mpf(case[4])
        ph = mp.mpf(case[1])
        checked = 0
        for s in (-R / 100, -R / 1000, -mp.mpf("1e-9") * R):
            for a in range(1, 8):
                alpha = ph + (a - 4) * mp.mpf("0.35")
                for k in range(6):
                    nu = 2 * mp.pi * k / 6 + mp.mpf("0.1")
                    q = ct.conditions_dir(s, nu, alpha, c, xm=mp)
                    if not (q["b"] < 0 and q["b"] ** 2 > q["cc"]):
                        continue
                    aim, _ = _aim_of(case, s, alpha)
                    P = (c["p0"][0] + s * c["e"][0], c["p0"][1] + s * c["e"][1])
                    st = bl._loop_state(P, nu, 1, aim, 1, c["B"], c["R"], bl.MP)
                    if st is None:
                        continue
                    assert abs(q["miss"] - st["miss"]) < mp.mpf(10) ** -30
                    W = st["speed"]
                    assert abs(W * q["C"] - bl.dot(bl.add(st["v2"], st["u"], -1), st["n"])) < mp.mpf(10) ** -30
                    assert abs(q["V"] - (bl.cross(st["u"], st["n"]) - c["v1"] * bl.cross(c["e"], st["n"]))) < mp.mpf(10) ** -30
                    checked += 1
        assert checked >= 20


def test_direction_conditions_are_finite_at_s_zero_and_L0_vanishes():
    c = ct.inputs(EMPTY, xm=mp)
    ph = mp.mpf(EMPTY[1])
    q = ct.conditions_dir(mp.mpf(0), mp.mpf(1), ph, c, xm=mp)
    assert q["cc"] == 0 and all(mp.isfinite(q[k]) for k in ("miss", "C", "V", "along", "kick"))
    q2 = ct.conditions_dir(mp.mpf("-1e-30"), mp.mpf(1), ph, c, xm=mp)
    assert abs(q2["miss"] - q["miss"]) < mp.mpf(10) ** -20


@pytest.mark.parametrize("box", [
    ((-1e-6, 0.0), (6.29, 6.2932), (1.9707, 1.9709)),               # tangent direction at s = 0
    ((-1e-4, -1e-5), (1.0, 1.1), (0.3, 0.5)),
    ((-1.35e-3, -1e-3), (4.0, 4.5), (2.0, 2.3)),
])
def test_interval_enclosures_contain_the_point_values(box):
    with ct.ivprec(53):
        q = ct.conditions_dir(*(iv.mpf(list(x)) for x in box), ct.inputs(EMPTY))
    c = ct.inputs(EMPTY, xm=mp)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                pt = [lo + (hi - lo) * t / 2 for (lo, hi), t in zip(box, (i, j, k))]
                p = ct.conditions_dir(*(mp.mpf(x) for x in pt), c, xm=mp)
                if not (p["b"] < 0 and p["b"] ** 2 >= p["cc"]):
                    continue  # outside the family: the enclosure need not contain it
                for name in ("miss", "C", "V", "dn", "along", "kick"):
                    assert q[name].a <= p[name] <= q[name].b, (name, pt)


def test_the_direction_machine_needs_a_line_that_starts_on_mouth_2():
    c = ct.inputs({"p0": ("0.5", "0"), "e": ("0", "1"), "v1": "0.5", "B": "0.01", "R": "1e-7"}, xm=mp)
    with pytest.raises(ValueError):
        ct.conditions_dir(mp.mpf(-1), 0, 0, c, xm=mp)


def _control_root():
    case, (s0, nu0, aim0) = CONTROL
    with mp.workdps(40):
        c = ct.inputs(case, xm=mp)
        sol = bl.solve_mp(c["p0"], c["e"], mp.mpf(s0), mp.mpf(nu0), mp.mpf(aim0), c["v1"], 1, 1, c["B"], c["R"])
        return case, (float(sol[0]), float(sol[1]), float(p3._alpha_of(case, sol)))


def test_krawczyk_certifies_a_band_control_root_on_the_direction_machine():
    case, x = _control_root()
    for h in (1e-8, 1e-11):
        k = ct.krawczyk(case, [(v - h, v + h) for v in x], machine="dir")
        assert k["verdict"] == "root" and k["inequalities"]
    far = [(x[0] - 1e-8, x[0] + 1e-8), (x[1] + 0.3, x[1] + 0.3 + 2e-8), (x[2] - 1e-8, x[2] + 1e-8)]
    assert ct.krawczyk(case, far, machine="dir")["verdict"] in ("none", "undecided")


def test_the_direction_machine_does_not_exclude_a_box_around_a_control_root():
    case, x = _control_root()
    r = ct.exclude(case, (x[0] - 1e-6, x[0] + 1e-6), (x[1] - 1e-4, x[1] + 1e-4), (x[2] - 1e-4, x[2] + 1e-4), machine="dir", max_boxes=3000, watch=[x])
    assert not r["watch"][0]["excluded_in"] and r["watch"][0]["undecided_in"]


def test_the_band_of_the_empty_case_near_s_zero_is_excluded():
    r = ct.exclude(EMPTY, (-1e-6, 0.0), scale=p3.band_scale(EMPTY), machine="dir", max_boxes=50_000)
    assert r["complete"]


def test_certification_runs_only_in_certify():
    full = {n for n, _ in p3.full_tasks()}
    cert = {n for n, _ in p3.certify_tasks()}
    assert not full & {"r8c", "r8pk", "r8p", "r8rc", "r8b", "r8bx", "r8bp"}
    assert {"r8rc", "r8c", "r8pk", "r8b", "r8bx", "r8bp"} <= cert
    assert [a for n, a in p3.certify_tasks() if n == "r8b"] == [(i, p3.CERT_CHUNKS.value) for i in range(p3.CERT_CHUNKS.value)]


def test_full_and_certify_are_exclusive():
    from chronosim import run
    with pytest.raises(SystemExit):
        run.main(["--full", "--certify", "--out", "/nonexistent"])


def _band_res(complete=True, agree=True, ctl_ok=True):
    lo, hi = p3.band_range(EMPTY)
    chunks = []
    n = p3.CERT_CHUNKS.value
    for i in range(n):
        a, b = lo + (hi - lo) * i / n, (lo + (hi - lo) * (i + 1) / n if i < n - 1 else hi)
        chunks.append(((i, n), {"case": tuple(EMPTY), "i": i, "n": n, "s": (a, b), "boxes": 10, "excluded_by": {"miss ≠ 0": 10},
                                "undecided": [], "n_undecided": 0 if complete else 1, "at_zero": 0, "complete": complete, "unexamined": 0, "range": (lo, hi)}))
    xs = [("aim", {"machine": "aim", "s": (-0.01, -0.00135), "boxes": 5, "complete": True, "n_undecided": 0, "unexamined": 0}),
          ("dir", {"machine": "dir", "s": (-0.01, -0.00135), "boxes": 5, "complete": agree, "n_undecided": 0 if agree else 2, "unexamined": 0})]
    box = ((-1e-3, -1e-3 + 1e-12), (4.37, 4.37 + 1e-9), (2.12, 2.12 + 1e-9))
    ctl = [(0, {"i": 0, "case": CONTROL[0], "sol": (-1e-3, 4.37, 2.12), "verify": (False, ["contact point inside a mouth (excluded by the model, p. 1081 item i)"]),
                "inequalities": {"dn": -0.33, "along": 0.76, "kick": 3.0}, "boxes": 100, "complete_search": True, "undecided": 1, "unexamined": 0,
                "in_excluded": 0 if ctl_ok else 1, "lost": 0, "copies": 1, "certified": ctl_ok, "solution_boxes": [box],
                "other": {"root": 0, "none": 0, "undecided": 0}})]
    return {"r8b": chunks, "r8bx": xs, "r8bp": ctl}


def test_band_proof_sentence_needs_exclusion_cross_check_and_control():
    lines, proved = p3._band_lines(_band_res())
    assert proved and any("Provado" in x for x in lines) and any("As duas concordam" in x for x in lines)
    assert any("não é afirmada" in x for x in lines)                         # main proof not given here
    assert any("cobre s ∈ [−3; 0)" in x for x in p3._band_lines(_band_res(), main_proved=True)[0])
    for kw in ({"complete": False}, {"agree": False}, {"ctl_ok": False}):
        lines, proved = p3._band_lines(_band_res(**kw))
        assert not proved and not any("**Provado" in x for x in lines), kw


def test_band_controls_are_labelled_as_machine_roots_not_physical_histories():
    lines, _ = p3._band_lines(_band_res())
    text = "\n".join(lines)
    assert "Não são histórias físicas" in text and "raiz certificada das equações da máquina por direção" in text
    assert "raiz certificada do modelo geométrico" not in text                 # not the model label: the verifier rejects them
    assert "testa a máquina contra raízes das próprias equações, não contra soluções físicas" in text
    assert "contact point inside a mouth" in text                              # the reasons come from the run


def test_the_verifier_sentence_follows_the_run_not_a_fixed_text():
    res = _band_res()
    res["r8bp"][0][1]["verify"] = (True, [])
    text = "\n".join(p3._band_lines(res)[0])
    assert "Não são histórias físicas" not in text and "O verificador independente aprova 1 de 1" in text


def test_krawczyk_on_the_direction_machine_requires_the_ray_to_reach_mouth_2_ahead():
    case, x = _control_root()
    k = ct.krawczyk(case, [(v - 1e-9, v + 1e-9) for v in x], machine="dir")
    assert k["verdict"] == "root" and k["inequalities"]
    with ct.ivprec(53):
        _, ineq = ct.jacobian([(v - 1e-9, v + 1e-9) for v in x], ct.inputs(case), machine="dir")
    assert ineq["b"].b < 0 and ineq["disc"].a > 0


def test_alpha_range_refuses_a_line_given_directly():
    with pytest.raises(ValueError):
        ct.alpha_range({"p0": ("0.5", "0"), "e": ("0", "1"), "v1": "0.5", "B": "0.01", "R": "1e-7"})
