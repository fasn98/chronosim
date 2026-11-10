"""The numbers of docs/phase3_borda.md, section 4 (item (i) of EKT p. 1081 and the mouth band):
what the band controls violate, the arc of contact normals item (i) forbids in the empty case, and
the grid showing that item (i) plus the machine inequalities leave points of the band."""

import mpmath as mp
import pytest

from chronosim import billiards as bl
from chronosim import certify as ct
from chronosim import phase3 as p3

EMPTY = p3.R8_CASES[2]


def _dist(X, cx):
    return mp.sqrt((X[0] - cx) ** 2 + X[1] ** 2)


# s, |P - C2| - B, |M - C2| - B, |Q - C2| - B, n . r_hat (table 4.2, 4 significant digits)
TABLE_42 = [("-1e-3", "8.780e-4", "-0.1213", "-0.2010", "0.9456"),
            ("-1e-5", "8.776e-6", "-0.1212", "-0.1981", "0.9414"),
            ("-1e-7", "8.776e-8", "-0.1212", "-0.1981", "0.9414")]


@pytest.mark.parametrize("i", range(3))
def test_the_band_controls_violate_item_i_and_have_the_old_ball_inside_mouth_2(i):
    case, (s0, nu0, aim0) = p3.BAND_CONTROLS[i]
    with mp.workdps(60):
        c = ct.inputs(case, xm=mp)
        s, nu, aim = bl.solve_mp(c["p0"], c["e"], mp.mpf(s0), mp.mpf(nu0), mp.mpf(aim0), c["v1"], 1, 1, c["B"], c["R"])
        ok, fails = bl.verify(c["p0"], c["e"], c["v1"], s, nu, aim, 1, 1, case[3], case[4])
        assert not ok and sorted(fails) == sorted(["contact point inside a mouth (excluded by the model, p. 1081 item i)",
                                                   "return leg enters a mouth before reaching the collision"])
        B, R = c["B"], c["R"]
        P = (c["p0"][0] + s * c["e"][0], c["p0"][1] + s * c["e"][1])
        n = (mp.cos(nu), mp.sin(nu))
        M, Q = (P[0] - R * n[0], P[1] - R * n[1]), (P[0] - 2 * R * n[0], P[1] - 2 * R * n[1])
        rh = ((P[0] - 1) / _dist(P, 1), P[1] / _dist(P, 1))
        got = (s, _dist(P, 1) - B, _dist(M, 1) - B, _dist(Q, 1) - B, n[0] * rh[0] + n[1] * rh[1])
        for g, want in zip(got, TABLE_42[i]):
            assert mp.nstr(g, 4) == mp.nstr(mp.mpf(want), 4)


def test_item_i_forbids_an_arc_of_about_153_to_154_degrees_in_the_empty_band():
    """Empty case (theta = 0: the line is radial, rho - B = |s|): the contact point is outside mouth 2
    iff n . r_hat < (rho^2 + R^2 - B^2) / (2 R rho); the old ball's center iff n . r_hat <
    (rho^2 + 4R^2 - B^2) / (4 R rho), which item (i) implies. Mouth 1 is farther than 2R + R."""
    c = ct.inputs(EMPTY, xm=mp)
    B, R = c["B"], c["R"]
    for s, cm_want, cq_want, arc_want in ((p3.band_range(EMPTY)[0], "0.234", "0.453", "152.9"), (0.0, "0.225", "0.45", "154.0")):
        s = mp.mpf(s)
        P = (c["p0"][0] + s * c["e"][0], c["p0"][1] + s * c["e"][1])
        rho = _dist(P, 1)
        assert abs(rho - B - abs(s)) < mp.mpf(10) ** -12
        cm = (rho ** 2 + R ** 2 - B ** 2) / (2 * R * rho)
        cq = (rho ** 2 + 4 * R ** 2 - B ** 2) / (4 * R * rho)
        assert mp.nstr(cm, 3) == cm_want and mp.nstr(cq, 3) == cq_want and cm < cq
        assert mp.nstr(2 * mp.degrees(mp.acos(cm)), 4) == arc_want
        assert _dist(P, 0) - B > 3 * R


def test_item_i_and_the_machine_inequalities_leave_points_of_the_band():
    """Grid of section 4.3: 3 values of s x 120 of nu x 60 of alpha = 21600 points; 21120 are in the
    family (the ray reaches mouth 2: b < 0 and b^2 > c); 6413 satisfy the machine inequalities
    (dn < 0, along > 0, kick > 0); 1179 also item (i), the contact point outside both mouths."""
    with mp.workdps(30):
        c = ct.inputs(EMPTY, xm=mp)
        B, R = c["B"], c["R"]
        a0, a1 = ct.alpha_range(EMPTY)
        family = ineq = item_i = 0
        first = None
        for s in (mp.mpf("-0.00135"), mp.mpf("-1e-5"), mp.mpf("-1e-9")):
            P = (c["p0"][0] + s * c["e"][0], c["p0"][1] + s * c["e"][1])
            for i in range(120):
                nu = 2 * mp.pi * i / 120
                M = (P[0] - R * mp.cos(nu), P[1] - R * mp.sin(nu))
                for j in range(60):
                    alpha = a0 + (a1 - a0) * (j + mp.mpf(1) / 2) / 60
                    q = ct.conditions_dir(s, nu, alpha, c, xm=mp)
                    if not (q["b"] < 0 and q["b"] ** 2 > q["cc"]):
                        continue
                    family += 1
                    if q["dn"] < 0 and q["along"] > 0 and q["kick"] > 0:
                        ineq += 1
                        if _dist(M, 1) > B and _dist(M, 0) > B:
                            item_i += 1
                            first = first or (s, nu, alpha, q)
        assert (family, ineq, item_i) == (21120, 6413, 1179)
        s, nu, alpha, q = first
        assert s == mp.mpf("-0.00135") and mp.nstr(nu, 4) == "1.466" and mp.nstr(alpha, 3) == "0.637"
        assert [mp.nstr(q[k], 2) for k in ("miss", "C", "V")] == ["-0.2", "-0.73", "-0.22"]
