"""Adversarial tests of the direction machine near mouth 2, the --certify split and its plumbing
(complement of tests/test_certify_band.py; docs/phase3_borda.md, decisions (a)-(d) of 2026-11-09)."""

import itertools
import math
from concurrent.futures import Future

import mpmath as mp
import pytest
from mpmath import iv

from chronosim import certify as ct
from chronosim import phase3 as p3
from chronosim import run

EMPTY = p3.R8_CASES[2]            # (theta, phi, lam1, B, R) = (0, 0.4, 0.75, 0.3, 0.135)
CONTROL_CASE, CONTROL_START = p3.BAND_CONTROLS[0]
PHI = float(EMPTY[1])


def _grid(lo, hi, n):
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


def _in_family(p):
    return p["b"] < 0 and p["b"] ** 2 >= p["cc"]


# conditions_dir: interval enclosure of the point values -----------------------------------------

# Boxes around the tangent corner (b = 0 at alpha = phi +- pi/2) at s = 0 and next to it, plus boxes
# that straddle the corner, touch s = 0 from below, or sit at the ends of the band.
TANGENT_UP, TANGENT_DOWN = PHI + math.pi / 2, PHI - math.pi / 2
CORNER_BOXES = [
    ((-1e-9, 0.0), (0.0, 6.3), (TANGENT_UP - 1e-3, TANGENT_UP + 1e-3)),
    ((-1e-12, 0.0), (1.0, 1.001), (TANGENT_UP - 1e-9, TANGENT_UP + 1e-9)),
    ((-1e-6, 0.0), (6.29, 6.2932), (TANGENT_UP - 2e-4, TANGENT_UP + 2e-4)),
    ((-1e-6, 0.0), (2.0, 2.5), (TANGENT_UP - 1e-2, TANGENT_UP)),
    ((-1e-6, 0.0), (2.0, 2.5), (TANGENT_UP - 0.05, TANGENT_UP - 0.04)),
    ((-1.35e-3, 0.0), (0.0, 1.0), (TANGENT_UP - 0.1, TANGENT_UP + 0.01)),
    ((-1e-6, 0.0), (3.0, 3.5), (TANGENT_DOWN - 1e-3, TANGENT_DOWN + 1e-3)),
    ((-1e-6, 0.0), (3.0, 3.5), (TANGENT_DOWN, TANGENT_DOWN + 0.05)),
    ((-1.35e-3, -1e-3), (5.0, 5.5), (TANGENT_UP - 1e-4, TANGENT_UP + 1e-4)),
    ((-1e-3, 0.0), (0.0, 6.3), (PHI - 1e-3, PHI + 1e-3)),                # head-on direction
    ((-1e-5, -1e-6), (1.0, 1.01), (PHI - 1.0, PHI - 0.9)),
    ((0.0, 0.0), (2.0, 2.0), (PHI + 0.3, PHI + 0.3)),                      # a degenerate (point) box at s = 0
]


@pytest.mark.parametrize("box", CORNER_BOXES)
def test_interval_conditions_enclose_every_family_point_near_the_tangent_corner(box):
    with ct.ivprec(53):
        q = ct.conditions_dir(*(iv.mpf(list(x)) for x in box), ct.inputs(EMPTY))
    c = ct.inputs(EMPTY, xm=mp)
    with mp.workdps(30):
        for pt in itertools.product(*(_grid(lo, hi, 5) for lo, hi in box)):
            p = ct.conditions_dir(*(mp.mpf(x) for x in pt), c, xm=mp)
            if not _in_family(p):
                continue
            for name in ("miss", "C", "V", "dn", "along", "kick"):
                assert q[name].a <= p[name] <= q[name].b, (name, pt, float(p[name]), q[name])


def test_the_corner_box_at_s_zero_has_a_finite_bounded_enclosure():
    """The quotient L0 = c/(-b + sqrt) is 0/0 there; the intersection must keep it finite and >= 0."""
    box = ((-1e-9, 0.0), (1.0, 2.0), (TANGENT_UP - 1e-3, TANGENT_UP + 1e-3))
    with ct.ivprec(53):
        q = ct.conditions_dir(*(iv.mpf(list(x)) for x in box), ct.inputs(EMPTY))
    for name in ("miss", "C", "V", "along", "kick"):
        assert mp.isfinite(q[name].a) and mp.isfinite(q[name].b), name
        assert q[name].b - q[name].a < 10, name   # not a blown-up hull


def test_interval_enclosure_is_not_vacuous_on_a_small_box():
    """A box of width 1e-9 in every coordinate must give enclosures of comparable width (no division blow-up)."""
    x = (-3e-4, 2.0, 1.0)
    box = [(v - 5e-10, v + 5e-10) for v in x]
    with ct.ivprec(53):
        q = ct.conditions_dir(*(iv.mpf(list(b)) for b in box), ct.inputs(EMPTY))
    for name in ("miss", "C", "V"):
        assert q[name].b - q[name].a < 1e-6, name


def test_conditions_dir_at_a_point_is_a_thin_interval_around_the_mp_value():
    x = (mp.mpf("-4e-4"), mp.mpf("2.5"), mp.mpf("0.9"))
    with mp.workdps(40):
        ref = ct.conditions_dir(*x, ct.inputs(EMPTY, xm=mp), xm=mp)
    with ct.ivprec(53):
        q = ct.conditions_dir(*(iv.mpf([float(v), float(v)]) for v in x), ct.inputs(EMPTY))
    for name in ("miss", "C", "V"):
        assert abs(q[name].mid - ref[name]) < 1e-9, name


# _meet ---------------------------------------------------------------------------------------------

def test_meet_of_overlapping_intervals_is_their_intersection():
    r = ct._meet(iv.mpf([0, 5]), iv.mpf([2, 8]), iv.mpf([1, 4]))
    assert (r.a, r.b) == (2, 4)


def test_meet_of_nested_intervals_is_the_inner_one():
    r = ct._meet(iv.mpf([-10, 10]), iv.mpf([1, 2]))
    assert (r.a, r.b) == (1, 2)


def test_meet_of_touching_intervals_is_the_common_point():
    r = ct._meet(iv.mpf([0, 1]), iv.mpf([1, 3]))
    assert (r.a, r.b) == (1, 1)


def test_meet_of_a_single_interval_is_itself():
    r = ct._meet(iv.mpf([-1, 2]))
    assert (r.a, r.b) == (-1, 2)


def test_meet_of_disjoint_intervals_falls_back_to_a_valid_hull_that_contains_both():
    """Enclosures of the same value cannot be disjoint; if they are, the result must not be an
    inverted interval, and it must not silently drop either enclosure."""
    r = ct._meet(iv.mpf([0, 1]), iv.mpf([2, 3]))
    assert r.a <= r.b
    assert r.a <= 0 and r.b >= 3


def test_meet_with_unbounded_half_line_keeps_the_nonnegative_part():
    r = ct._meet(iv.mpf([-2, 3]), iv.mpf([0, mp.inf]))
    assert (r.a, r.b) == (0, 3)


# _nonneg -------------------------------------------------------------------------------------------

def test_nonneg_clips_the_negative_part_of_an_interval():
    assert (lambda r: (r.a, r.b))(ct._nonneg(iv.mpf([-3, 2]))) == (0, 2)
    assert (lambda r: (r.a, r.b))(ct._nonneg(iv.mpf([1, 2]))) == (1, 2)


def test_nonneg_of_a_wholly_negative_interval_is_zero_not_an_error():
    r = ct._nonneg(iv.mpf([-3, -1]))
    assert (r.a, r.b) == (0, 0)


# _reason_dir: one test per reason ------------------------------------------------------------------

def _reason(box):
    with ct.ivprec(53):
        return ct._reason_dir(box, ct.inputs(EMPTY))


def test_reason_point_inside_mouth_2():
    assert "dentro da boca" in _reason(((0.1, 0.2), (1.0, 2.0), (0.0, 1.0)))


def test_reason_ray_not_heading_to_mouth_2():
    assert "b ≥ 0" in _reason(((-0.01, -0.009), (1.0, 2.0), (-2.5, -2.4)))


def test_reason_ray_does_not_reach_mouth_2():
    assert "não atinge" in _reason(((-0.01, -0.009), (1.0, 2.0), (1.8, 1.9)))


def test_reason_kick_not_positive():
    assert "impulso" in _reason(((-0.01, -0.009), (0.0, 0.5), (-1.05, -0.92)))


def test_reason_old_ball_has_not_left_mouth_1():
    assert "along" in _reason(((-0.01, -0.009), (3.0, 3.5), (-1.05, -0.92)))


def test_reason_miss_nonzero():
    assert _reason(((-0.05, -0.0499), (1.6, 1.61), (0.48, 0.482))) == "miss ≠ 0"


def test_reason_C_nonzero():
    assert _reason(((-0.05, -0.0499), (3.2, 3.21), (0.68, 0.682))) == "C ≠ 0"


def test_reason_none_on_a_box_that_holds_a_root():
    x = _control_root()
    box = tuple((v - 1e-9, v + 1e-9) for v in x)
    with ct.ivprec(53):
        assert ct._reason_dir(box, ct.inputs(CONTROL_CASE)) is None


def _q(**kw):
    base = {"dn": (-1, -0.5), "along": (0.5, 1), "kick": (0.5, 1), "miss": (-1, 1), "C": (-1, 1), "V": (-1, 1)}
    base.update(kw)
    return {k: iv.mpf(list(v)) for k, v in base.items()}


def test_reason_from_covers_dn_and_V_and_priorities():
    """dn >= 0 and V != 0 are not reachable on a box of the empty case (scan of ~8000 boxes), so
    they are exercised on crafted enclosures."""
    assert "d·n_X" in ct._reason_from(_q(dn=(0, 1)))
    assert "d·n_X" in ct._reason_from(_q(dn=(0, 1), along=(-1, 0), kick=(-1, 0)))       # dn is checked first
    assert ct._reason_from(_q(V=(0.1, 1))) == "V ≠ 0"
    assert ct._reason_from(_q(V=(-1, -0.1))) == "V ≠ 0"
    assert ct._reason_from(_q()) is None
    assert ct._reason_from(_q(miss=(0, 1))) is None                   # an enclosure touching 0 is not excluded
    assert "along" in ct._reason_from(_q(along=(-1, 0)))               # along <= 0 exactly: excluded
    assert ct._reason_from(_q(along=(-1, 1e-12))) is None
    assert "impulso" in ct._reason_from(_q(kick=(-1, 0)))


def _scan_boxes():
    for s in [(-0.05, -0.0499), (-1e-2, -9.9e-3), (-1.3e-3, -1.2e-3), (-1e-9, -0.9e-9)]:
        for i in range(0, 40, 6):
            nu = (i * 0.16, i * 0.16 + 0.01)
            for k in range(0, 60, 3):
                a = -1.17 + k * 0.05
                yield (s, nu, (a, a + 0.002))


def test_every_reason_claimed_on_a_scan_is_true_at_the_family_points_of_the_box():
    """Soundness of _reason_dir against point values: whenever a reason is given for a box, every
    sampled family point of that box violates the corresponding condition (and is not a root)."""
    cm = ct.inputs(EMPTY, xm=mp)
    ci = None
    n_reasoned = 0
    with ct.ivprec(53):
        ci = ct.inputs(EMPTY)
        reasons = [(box, ct._reason_dir(box, ci)) for box in _scan_boxes()]
    with mp.workdps(30):
        for box, why in reasons:
            if why is None:
                continue
            n_reasoned += 1
            for pt in itertools.product(*(_grid(lo, hi, 3) for lo, hi in box)):
                p = ct.conditions_dir(*(mp.mpf(x) for x in pt), cm, xm=mp)
                if not _in_family(p):
                    continue
                if "along" in why:
                    assert p["along"] <= 0, (why, pt)
                elif "impulso" in why:
                    assert p["kick"] <= 0, (why, pt)
                elif "d·n_X" in why:
                    assert p["dn"] >= 0, (why, pt)
                elif why == "miss ≠ 0":
                    assert p["miss"] != 0, (why, pt)
                elif why == "C ≠ 0":
                    assert p["C"] != 0, (why, pt)
                else:  # "raio não atinge" / "b ≥ 0" / "dentro da boca": no family point at all
                    pytest.fail(f"{why}: a family point exists at {pt}")
    assert n_reasoned > 100


# krawczyk(machine="dir") ---------------------------------------------------------------------------

def _control_root():
    from chronosim import billiards as bl
    s0, nu0, a0 = CONTROL_START
    with mp.workdps(40):
        c = ct.inputs(CONTROL_CASE, xm=mp)
        sol = bl.solve_mp(c["p0"], c["e"], mp.mpf(s0), mp.mpf(nu0), mp.mpf(a0), c["v1"], 1, 1, c["B"], c["R"])
        return [float(sol[0]), float(sol[1]), float(p3._alpha_of(CONTROL_CASE, sol))]


@pytest.mark.parametrize("off,h", [([0, 0.3, 0], 1e-3), ([0, 0.02, 0], 1e-4), ([1e-3, 0, 0], 1e-6), ([1e-4, 0, 0], 1e-6), ([0, 0.3, 0], 1e-6)])
def test_krawczyk_dir_proves_a_box_without_a_root_has_none(off, h):
    x = _control_root()
    k = ct.krawczyk(CONTROL_CASE, [(v + o - h, v + o + h) for v, o in zip(x, off)], machine="dir")
    assert k["verdict"] == "none"
    assert any(kk[1] < lo or kk[0] > hi for kk, (lo, hi) in zip(k["K"], [(v + o - h, v + o + h) for v, o in zip(x, off)]))


@pytest.mark.parametrize("h", [1e-3, 1e-5, 1e-7, 1e-9])
def test_krawczyk_dir_never_says_none_on_a_box_containing_the_root(h):
    x = _control_root()
    for shift in (0.0, 0.4, -0.4, 0.9):
        box = [(v + shift * h - h, v + shift * h + h) for v in x]
        assert ct.krawczyk(CONTROL_CASE, box, machine="dir")["verdict"] != "none", (h, shift)


def test_krawczyk_dir_machine_and_aim_machine_differ():
    """The machine argument must reach the Jacobian: the same box has different K on the two machines."""
    x = _control_root()
    box = [(v - 1e-6, v + 1e-6) for v in x]
    k_dir = ct.krawczyk(CONTROL_CASE, box, machine="dir")
    k_aim = ct.krawczyk(CONTROL_CASE, box, machine="aim")
    assert k_dir["K"] != k_aim["K"]


def test_krawczyk_rejects_an_unknown_machine():
    with pytest.raises(KeyError):
        ct.krawczyk(CONTROL_CASE, [(0, 1)] * 3, machine="nope")


# exclude(machine="dir") with watch -----------------------------------------------------------------

def test_exclude_dir_watch_inside_an_excluded_box_is_reported_as_excluded():
    x = _control_root()
    pt = (x[0], x[1] + 0.3, x[2])                                 # not a root: a nearby non-root point
    r = ct.exclude(CONTROL_CASE, (pt[0] - 1e-6, pt[0] + 1e-6), (pt[1] - 1e-4, pt[1] + 1e-4), (pt[2] - 1e-4, pt[2] + 1e-4),
                   machine="dir", watch=[pt])
    assert r["complete"]
    assert r["watch"][0]["excluded_in"] and not r["watch"][0]["undecided_in"]


def test_exclude_dir_watch_is_not_normalized_when_outside_the_root_box():
    """With machine='dir' the watched point is used as given: a point outside the box is in no box
    (lost), instead of being wrapped by nu + 2 pi or aim + 2 pi as the aim machine does."""
    x = _control_root()
    pt = (x[0], x[1] + 2 * math.pi, x[2])                         # nu above the default range (2 pi + pad)
    r = ct.exclude(CONTROL_CASE, (x[0] - 1e-6, x[0] + 1e-6), machine="dir", watch=[pt], max_boxes=2000)
    assert not r["watch"][0]["excluded_in"] and not r["watch"][0]["undecided_in"]


def test_exclude_dir_watch_with_nu_above_2pi_in_a_narrow_box_is_left_alone():
    """The aim machine would wrap this point by 2 pi; the dir machine must not."""
    with mp.workdps(30):
        pt = (mp.mpf("-1e-4"), mp.mpf(2) * mp.pi + mp.mpf("1.5"), mp.mpf("0.5"))
    r_dir = ct.exclude(EMPTY, (-1.1e-4, -0.9e-4), (1.4, 1.6), (0.4, 0.6), machine="dir", watch=[pt], max_boxes=3000)
    assert not r_dir["watch"][0]["excluded_in"] and not r_dir["watch"][0]["undecided_in"]


def test_exclude_dir_uses_the_alpha_range_by_default():
    """The default third range of the dir machine is alpha_range(case), not the aim range [-pad, pad]."""
    lo, hi = ct.alpha_range(EMPTY)
    assert lo == pytest.approx(PHI - math.pi / 2 - ct.ALPHA_PAD.value) and hi == pytest.approx(PHI + math.pi / 2 + ct.ALPHA_PAD.value)
    # a watched point with alpha outside alpha_range is in no box even though it is inside [-pad, pad] of the aim range
    pt = (-1e-4, 2.0, hi + 0.05)
    r = ct.exclude(EMPTY, (-1.1e-4, -0.9e-4), (1.9, 2.1), watch=[pt], machine="dir", max_boxes=5000)
    assert not r["watch"][0]["excluded_in"] and not r["watch"][0]["undecided_in"]


def test_exclude_rejects_an_unknown_machine():
    with pytest.raises(KeyError):
        ct.exclude(EMPTY, (-1e-4, -1e-5), machine="nope")


def test_exclude_dir_of_a_case_that_is_not_on_mouth_2_raises():
    case = {"p0": ("0.5", "0"), "e": ("0", "1"), "v1": "0.5", "B": "0.01", "R": "1e-7"}
    with pytest.raises((ValueError, KeyError)):
        ct.exclude(case, (-0.5, -0.1), machine="dir", max_boxes=10)


def test_exclude_dir_stops_at_the_box_budget_and_says_not_complete():
    r = ct.exclude(EMPTY, (-1.35e-3, 0.0), scale=p3.band_scale(EMPTY), machine="dir", max_boxes=50)
    assert not r["complete"] and r["unexamined"] > 0


# start(mode="certify") ---------------------------------------------------------------------------

class _SyncExecutor:
    """Runs submitted work at once, in this process (the real one forks workers)."""
    def __init__(self, *a, **k):
        pass

    def submit(self, fn, *args):
        f = Future()
        try:
            f.set_result(fn(*args))
        except BaseException as e:  # noqa: BLE001
            f.set_exception(e)
        return f

    def shutdown(self, *a, **k):
        pass


@pytest.fixture
def sync_pool(monkeypatch):
    import concurrent.futures as cf
    monkeypatch.setattr(cf, "ProcessPoolExecutor", _SyncExecutor)
    ran = []
    monkeypatch.setattr(p3, "_task", lambda name, arg: ran.append((name, arg)) or {"task": name})
    return ran


def test_start_certify_submits_exactly_the_certify_task_list(sync_pool):
    handle = p3.start(mode="certify")
    assert handle[1] == p3.certify_tasks()
    assert [t for t in sync_pool if t[0] != "r8p"] == p3.certify_tasks()   # r8p pieces are follow-ups of r8pk


def test_start_certify_submits_the_positive_control_pieces_after_r8pk(sync_pool):
    res = p3.collect(p3.start(mode="certify"))
    assert len(res["r8pk"]) == len(p3.CONTROL_CASES)
    assert len(res.get("r8p", [])) > 0
    assert set(res) >= {"r8rc", "r8c", "r8pk", "r8p", "r8b", "r8bx", "r8bp"}
    assert len(res["r8b"]) == p3.CERT_CHUNKS.value and len(res["r8bp"]) == len(p3.BAND_CONTROLS)
    assert [a for a, _ in res["r8bx"]] == ["aim", "dir"]


def test_start_certify_runs_none_of_the_phase_1_to_11_tasks(sync_pool):
    res = p3.collect(p3.start(mode="certify"))
    assert not set(res) & {"r1", "r2", "rest", "r3", "r46", "r8x", "r8r", "r9", "r10", "r10x"}


def test_start_full_and_quick_no_longer_run_certification(sync_pool):
    for kw in ({"quick": True}, {}):
        res = p3.collect(p3.start(**kw))
        assert not set(res) & {"r8c", "r8pk", "r8p", "r8rc", "r8b", "r8bx", "r8bp"}, kw


def test_start_quick_flag_and_default_mode_are_the_quick_subset(sync_pool):
    assert p3.start(quick=True)[1] == list(p3.QUICK_TASKS)
    assert p3.start()[1] == p3.full_tasks()
    assert p3.start(quick=True, mode="certify")[1] == p3.certify_tasks()   # an explicit mode wins


def test_start_rejects_an_unknown_mode(sync_pool):
    with pytest.raises(KeyError):
        p3.start(mode="nope")


def test_the_r8rc_task_runs_the_interval_test_and_r8r_does_not(monkeypatch):
    calls = []
    monkeypatch.setattr(p3, "dangerous", lambda *a: {"sols": [], "s_grid": [], "p0": (mp.mpf(0), mp.mpf(0)), "e": (mp.mpf(1), mp.mpf(0))})
    monkeypatch.setattr(p3, "refine_unresolved", lambda case: {"items_after": ["x"]})
    monkeypatch.setattr(p3, "certify_items", lambda *a: calls.append(a) or "CERT")
    assert p3._task("r8r", EMPTY)["cert_items"] is None and not calls
    assert p3._task("r8rc", EMPTY)["cert_items"] == "CERT" and len(calls) == 1


# run.main(["--certify"]) --------------------------------------------------------------------------

def _fake_res():
    import sys
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
    import test_certify_band as tb
    return tb._band_res()


def test_certify_writes_only_phase3_certify_md_and_leaves_report_alone(tmp_path, monkeypatch, capsys):
    """The real run is ~100 min; the pool is replaced by a canned result of the right shape."""
    seen = {}
    monkeypatch.setattr(p3, "start", lambda **kw: seen.update(kw) or "handle")
    monkeypatch.setattr(p3, "collect", lambda h: _fake_res())
    (tmp_path / "report.md").write_text("untouched")
    assert run.main(["--certify", "--out", str(tmp_path)]) == 0
    assert seen == {"mode": "certify"}
    assert sorted(x.name for x in tmp_path.iterdir()) == ["phase3_certify.md", "report.md"]
    assert (tmp_path / "report.md").read_text() == "untouched"
    out = capsys.readouterr().out
    assert "Provado" in out and "constantes HEURISTIC usadas" in out and "ekt_cert_alpha_pad" in out


def test_certify_output_has_no_proof_sentence_when_the_band_is_incomplete(tmp_path, monkeypatch, capsys):
    import sys
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
    import test_certify_band as tb
    monkeypatch.setattr(p3, "start", lambda **kw: "handle")
    monkeypatch.setattr(p3, "collect", lambda h: tb._band_res(complete=False))
    assert run.main(["--certify", "--out", str(tmp_path)]) == 0
    assert "**Provado" not in (tmp_path / "phase3_certify.md").read_text()
    assert "Provado" not in capsys.readouterr().out


def test_certify_does_not_start_the_other_phases(tmp_path, monkeypatch):
    monkeypatch.setattr(p3, "start", lambda **kw: "handle")
    monkeypatch.setattr(p3, "collect", lambda h: _fake_res())
    monkeypatch.setattr(run, "report", lambda *a, **k: pytest.fail("--certify must not write report.md"))
    run.main(["--certify", "--out", str(tmp_path)])


def test_certify_section_with_an_empty_result_does_not_crash_and_claims_nothing():
    lines = p3.certify_section({})
    text = "\n".join(lines)
    assert "Não rodada nesta execução" in text and "**Provado" not in text


@pytest.mark.parametrize("flags", [["--certify", "--full"], ["--full", "--certify"]])
def test_certify_and_full_are_mutually_exclusive_in_either_order(flags, tmp_path):
    with pytest.raises(SystemExit):
        run.main(flags + ["--out", str(tmp_path)])
    assert not list(tmp_path.iterdir())
