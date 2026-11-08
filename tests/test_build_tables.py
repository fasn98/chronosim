"""paper/build_tables.py: the paper's numbers come from here, so its guards must really stop the
build. Git-dependent functions run in a throw-away repository; ROOT is redirected, never the real one."""

import importlib.util
import re
import subprocess
from pathlib import Path

import mpmath as mp
import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "paper" / "build_tables.py"
REAL_ROOT = SCRIPT.parent.parent


@pytest.fixture(scope="module")
def bt():
    """Import the script; it must not change the global mpmath precision (checked below)."""
    saved = mp.mp.dps
    spec = importlib.util.spec_from_file_location("build_tables_under_test", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
        mod.GLOBAL_DPS_AFTER_IMPORT = mp.mp.dps
        mod.DPS_SET_BY_SCRIPT = mod.DIGITS  # the precision the script works at, read back
    finally:
        mp.mp.dps = saved
    mod.GLOBAL_DPS_BEFORE_IMPORT = saved
    return mod


@pytest.fixture
def at40(bt):
    with mp.workdps(bt.DPS_SET_BY_SCRIPT):  # the precision the script itself runs at
        yield


# ---------------------------------------------------------------- tex()

def test_tex_keeps_twenty_digit_input_without_loss(bt, at40):
    assert bt.tex("1.0000000000000000001", sig=20) == "1.0000000000000000001"


def test_tex_twenty_digit_small_number_in_scientific_form(bt, at40):
    assert bt.tex("0.00012345678901234567891", sig=20) == "1.2345678901234567891\\times10^{-4}"


def test_tex_rounds_nineteen_nines_up_to_two(bt, at40):
    assert bt.tex("1.9999999999999999999") == "2"


def test_tex_three_significant_digits_of_a_twenty_digit_value(bt, at40):
    assert bt.tex("3.9600000000000000001e-20") == "3.96\\times10^{-20}"


def test_tex_exact_zero(bt, at40):
    assert bt.tex("0") == "0"
    assert bt.tex("0.0000000000000000000") == "0"


def test_tex_power_of_ten_has_no_unit_mantissa(bt, at40):
    assert bt.tex("1e-20") == "10^{-20}"
    assert bt.tex("-1e-20") == "-10^{-20}"


def test_tex_rounding_carry_into_the_exponent(bt, at40):
    assert bt.tex("9.9996e-5") == "10^{-4}"


def test_tex_negative_keeps_sign_in_both_forms(bt, at40):
    assert bt.tex("-0.5") == "-0.5"
    assert bt.tex("-3.96e-20") == "-3.96\\times10^{-20}"


def test_tex_integer_valued_float_has_no_trailing_point_zero(bt, at40):
    assert bt.tex("2.0") == "2"


def test_script_precision_holds_twenty_digit_csv_values(bt):
    """The script's own dps must be enough for the 20 significant digits of the CSVs."""
    assert bt.DPS_SET_BY_SCRIPT >= 25


def test_tex_keeps_twenty_digits_even_when_the_caller_runs_at_15(bt):
    """The script's precision is local: a caller at the mpmath default does not lose digits."""
    with mp.workdps(15):
        assert bt.tex("1.0000000000000000001", sig=20) == "1.0000000000000000001"


def test_importing_the_script_leaves_the_global_precision_alone(bt):
    """Scalpel note 2: the import used to set mp.mp.dps = 40 for the whole process."""
    assert bt.GLOBAL_DPS_AFTER_IMPORT == bt.GLOBAL_DPS_BEFORE_IMPORT


# ---------------------------------------------------------------- Macros

def test_macro_name_must_be_letters_only(bt):
    m = bt.Macros()
    for bad in ("Phase1", "a_b", "", "x y"):
        with pytest.raises(ValueError):
            m[bad] = "1"


def test_macro_defined_twice_is_an_error(bt):
    m = bt.Macros()
    m["Abc"] = "1"
    with pytest.raises(ValueError):
        m["Abc"] = "2"


# ---------------------------------------------------------------- docs counts

LIT = REAL_ROOT / "docs" / "literature_check_factor2.md"
EQ = REAL_ROOT / "docs" / "equation_transcription.md"


@pytest.fixture
def docs_root(tmp_path, bt, monkeypatch):
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs" / "literature_check_factor2.md").write_text(LIT.read_text(encoding="utf-8"), encoding="utf-8")
    (tmp_path / "docs" / "equation_transcription.md").write_text(EQ.read_text(encoding="utf-8"), encoding="utf-8")
    monkeypatch.setattr(bt, "ROOT", tmp_path)
    return tmp_path


def _edit(path, old, new, count=1):
    t = path.read_text(encoding="utf-8")
    assert old in t, f"fixture premise broken: {old!r} not in {path.name}"
    path.write_text(t.replace(old, new, count), encoding="utf-8")


def test_literature_counts_on_the_real_doc(bt, docs_root):
    assert bt.literature_counts() == (175, 131, 99, 44)


def test_literature_counts_stop_when_the_inspire_line_changes(bt, docs_root):
    _edit(docs_root / "docs" / "literature_check_factor2.md", "trabalhos que citam Kim & Thorne;", "works citing Kim & Thorne;")
    with pytest.raises(SystemExit, match="not found"):
        bt.literature_counts()


def test_literature_counts_stop_when_the_mention_cell_changes(bt, docs_root):
    _edit(docs_root / "docs" / "literature_check_factor2.md", "mencionam Kim & Thorne", "citam Kim & Thorne no texto")
    with pytest.raises(SystemExit, match="not found"):
        bt.literature_counts()


def test_literature_counts_stop_when_the_missing_sentence_changes(bt, docs_root):
    _edit(docs_root / "docs" / "literature_check_factor2.md", "não têm versão no arXiv", "não estão no arXiv")
    with pytest.raises(SystemExit, match="not found"):
        bt.literature_counts()


def test_literature_counts_stop_when_arithmetic_is_inconsistent(bt, docs_root):
    _edit(docs_root / "docs" / "literature_check_factor2.md", "44 dos 175", "43 dos 175")
    with pytest.raises(SystemExit, match="inconsistent"):
        bt.literature_counts()


def test_literature_counts_stop_when_total_differs_from_citing(bt, docs_root):
    _edit(docs_root / "docs" / "literature_check_factor2.md", "44 dos 175", "44 dos 174")
    with pytest.raises(SystemExit, match="inconsistent"):
        bt.literature_counts()


def test_transcribed_equations_on_the_real_doc(bt, docs_root):
    assert bt.transcribed_equations() == 52


def test_transcribed_equations_stop_when_stated_count_changes(bt, docs_root):
    _edit(docs_root / "docs" / "equation_transcription.md", "as 52 equações conferem", "as 51 equações conferem")
    with pytest.raises(SystemExit, match="disagree"):
        bt.transcribed_equations()


def test_transcribed_equations_stop_when_the_statement_disappears(bt, docs_root):
    _edit(docs_root / "docs" / "equation_transcription.md", "as 52 equações conferem", "todas conferem")
    with pytest.raises(SystemExit, match="disagree"):
        bt.transcribed_equations()


def test_transcribed_equations_stop_when_a_row_is_not_agreeing(bt, docs_root):
    p = docs_root / "docs" / "equation_transcription.md"
    t = p.read_text(encoding="utf-8")
    line = next(l for l in t.split("# Echeverria")[0].splitlines() if re.match(r"^\| \d+[a-z]? \| \d{4} \|", l))
    assert line.rstrip().endswith("| sim |")
    p.write_text(t.replace(line, line.rstrip()[: -len("| sim |")] + "| não |", 1), encoding="utf-8")
    with pytest.raises(SystemExit, match="disagree"):
        bt.transcribed_equations()


def test_transcribed_equations_stop_when_a_row_is_deleted(bt, docs_root):
    p = docs_root / "docs" / "equation_transcription.md"
    t = p.read_text(encoding="utf-8")
    line = next(l for l in t.split("# Echeverria")[0].splitlines() if re.match(r"^\| \d+[a-z]? \| \d{4} \|", l))
    p.write_text(t.replace(line + "\n", "", 1), encoding="utf-8")
    with pytest.raises(SystemExit, match="disagree"):
        bt.transcribed_equations()


def test_transcribed_equations_ignore_rows_after_the_echeverria_section(bt, docs_root):
    p = docs_root / "docs" / "equation_transcription.md"
    with p.open("a", encoding="utf-8") as fh:
        fh.write("\n| 99 | 1234 | `x` | `x` | sim |\n")
    assert bt.transcribed_equations() == 52


# ---------------------------------------------------------------- git: release() and trailers()

def _git(repo, *args, env_extra=None):
    import os
    env = {**os.environ, "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t", "GIT_COMMITTER_NAME": "t",
           "GIT_COMMITTER_EMAIL": "t@t", "GIT_CONFIG_GLOBAL": "/dev/null", "GIT_CONFIG_SYSTEM": "/dev/null"}
    r = subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
    return r.stdout.strip()


def _commit(repo, name, body=""):
    (repo / name).write_text(name, encoding="utf-8")
    _git(repo, "add", name)
    _git(repo, "commit", "-q", "-m", f"add {name}" + (f"\n\n{body}" if body else ""))
    return _git(repo, "rev-parse", "HEAD")


@pytest.fixture
def repo(tmp_path, bt, monkeypatch):
    r = tmp_path / "repo"
    r.mkdir()
    _git(r, "init", "-q", "-b", "main")
    monkeypatch.setattr(bt, "ROOT", r)
    return r


def test_release_without_tag_is_pending_when_not_required(bt, repo):
    _commit(repo, "a")
    assert bt.release(False).startswith("PENDING")


def test_release_without_tag_stops_when_required(bt, repo):
    _commit(repo, "a")
    with pytest.raises(SystemExit, match="does not exist"):
        bt.release(True)


def test_release_with_head_on_the_tag_and_clean_tree(bt, repo):
    h = _commit(repo, "a")
    _git(repo, "tag", bt.RELEASE_TAG)
    out = bt.release(False)
    assert out == f"{bt.RELEASE_TAG} (commit \\texttt{{{h[:12]}}})"


def test_release_stops_when_head_moved_past_the_tag(bt, repo):
    _commit(repo, "a")
    _git(repo, "tag", bt.RELEASE_TAG)
    _commit(repo, "b")
    with pytest.raises(SystemExit, match="not the commit of tag"):
        bt.release(False)


def test_release_head_mismatch_stops_even_without_require(bt, repo):
    _commit(repo, "a")
    _git(repo, "tag", bt.RELEASE_TAG)
    _commit(repo, "b")
    with pytest.raises(SystemExit):
        bt.release(True)


def test_release_stops_on_modified_tracked_file(bt, repo):
    _commit(repo, "a")
    _git(repo, "tag", bt.RELEASE_TAG)
    (repo / "a").write_text("changed", encoding="utf-8")
    with pytest.raises(SystemExit, match="uncommitted"):
        bt.release(False)


def test_release_stops_on_staged_change(bt, repo):
    _commit(repo, "a")
    _git(repo, "tag", bt.RELEASE_TAG)
    (repo / "a").write_text("changed", encoding="utf-8")
    _git(repo, "add", "a")
    with pytest.raises(SystemExit, match="uncommitted"):
        bt.release(False)


def test_release_ignores_untracked_files(bt, repo):
    """Documented behaviour (--untracked-files=no): generated files such as paper/generated may exist."""
    _commit(repo, "a")
    _git(repo, "tag", bt.RELEASE_TAG)
    (repo / "scratch.txt").write_text("x", encoding="utf-8")
    assert bt.release(False).startswith(bt.RELEASE_TAG)


def test_release_with_annotated_tag_resolves_to_the_commit(bt, repo):
    h = _commit(repo, "a")
    _git(repo, "tag", "-a", bt.RELEASE_TAG, "-m", "release")
    assert h[:12] in bt.release(False)


def test_release_in_a_directory_that_is_not_a_repository(bt, tmp_path, monkeypatch):
    plain = tmp_path / "plain"
    plain.mkdir()
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path))
    monkeypatch.setattr(bt, "ROOT", plain)
    assert bt.release(False).startswith("PENDING")
    with pytest.raises(SystemExit):
        bt.release(True)


CLAUDE = "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"


def test_trailers_counts_each_kind(bt, repo):
    _commit(repo, "a")
    _commit(repo, "b", CLAUDE)
    _commit(repo, "c", CLAUDE)
    assert bt.trailers() == {"plain_with": 2, "plain_without": 1, "merge_with": 0, "merge_without": 0}


def test_trailers_ignores_a_human_coauthor(bt, repo):
    _commit(repo, "a", "Co-Authored-By: Maria <maria@example.org>")
    assert bt.trailers()["plain_without"] == 1 and bt.trailers()["plain_with"] == 0


def test_trailers_ignores_claude_mentioned_only_in_the_body(bt, repo):
    _commit(repo, "a", "Claude helped with this commit.")
    assert bt.trailers()["plain_with"] == 0


def test_trailers_sees_merge_commits_apart_from_plain_ones(bt, repo):
    _commit(repo, "a")
    _git(repo, "checkout", "-q", "-b", "side")
    _commit(repo, "b", CLAUDE)
    _git(repo, "checkout", "-q", "main")
    _commit(repo, "c")
    _git(repo, "merge", "-q", "--no-ff", "-m", "merge side", "side")
    assert bt.trailers() == {"plain_with": 1, "plain_without": 2, "merge_with": 0, "merge_without": 1}


def test_trailers_merge_with_trailer_is_counted_as_merge_with(bt, repo):
    _commit(repo, "a")
    _git(repo, "checkout", "-q", "-b", "side")
    _commit(repo, "b")
    _git(repo, "checkout", "-q", "main")
    _commit(repo, "c")
    _git(repo, "merge", "-q", "--no-ff", "-m", "merge side\n\n" + CLAUDE, "side")
    assert bt.trailers()["merge_with"] == 1


def test_trailers_with_pipe_in_subject_does_not_break_the_parse(bt, repo):
    (repo / "a").write_text("a", encoding="utf-8")
    _git(repo, "add", "a")
    _git(repo, "commit", "-q", "-m", "fix a | b | c\n\n" + CLAUDE)
    assert bt.trailers()["plain_with"] == 1


def test_trailers_on_empty_repository_is_all_zero(bt, repo):
    """git log fails on a repository with no commit; the function must not invent counts."""
    assert bt.trailers() == {"plain_with": 0, "plain_without": 0, "merge_with": 0, "merge_without": 0}


def test_trailers_on_the_real_repository_never_attribute_a_merge(bt):
    """The build aborts if a merge carries the trailer; make sure the real history satisfies it."""
    assert bt.trailers()["merge_with"] == 0


# ------------------------------------- added after paper/build_tables.py changed during this review
# (rows() origin check, computed(), and the paper/generated exclusion in release()); the working-tree
# version of the script is what is tested here, not the committed one.

def test_release_ignores_a_modified_tracked_file_under_paper_generated(bt, repo):
    (repo / "paper" / "generated").mkdir(parents=True)
    (repo / "paper" / "generated" / "numbers.tex").write_text("v1", encoding="utf-8")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "gen")
    _git(repo, "tag", bt.RELEASE_TAG)
    (repo / "paper" / "generated" / "numbers.tex").write_text("v2", encoding="utf-8")
    assert bt.release(False).startswith(bt.RELEASE_TAG)


def test_release_still_stops_on_a_modified_file_outside_paper_generated(bt, repo):
    (repo / "paper" / "generated").mkdir(parents=True)
    (repo / "paper" / "generated" / "numbers.tex").write_text("v1", encoding="utf-8")
    (repo / "paper" / "other.tex").write_text("v1", encoding="utf-8")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "gen")
    _git(repo, "tag", bt.RELEASE_TAG)
    (repo / "paper" / "other.tex").write_text("v2", encoding="utf-8")
    with pytest.raises(SystemExit, match="uncommitted"):
        bt.release(False)


def _csv_root(tmp_path, bt, monkeypatch, text):
    (tmp_path / "outputs").mkdir()
    (tmp_path / "outputs" / "t.csv").write_text(text, encoding="utf-8")
    monkeypatch.setattr(bt, "OUT", tmp_path / "outputs")


def test_rows_rejects_an_unknown_origin(bt, tmp_path, monkeypatch):
    _csv_root(tmp_path, bt, monkeypatch, "a,origin\n1,GUESSED\n")
    with pytest.raises(SystemExit, match="unknown origin"):
        bt.rows("t.csv")


def test_rows_accepts_a_known_origin_with_a_suffix_after_semicolon(bt, tmp_path, monkeypatch):
    _csv_root(tmp_path, bt, monkeypatch, "a,origin\n1,COMPUTED; extra\n2,HEURISTIC\n")
    assert len(bt.rows("t.csv")) == 2


def test_rows_checks_prefixed_origin_columns_too(bt, tmp_path, monkeypatch):
    _csv_root(tmp_path, bt, monkeypatch, "a,printed_origin,recomputed_origin\n1,CITED,WHATEVER\n")
    with pytest.raises(SystemExit, match="recomputed_origin"):
        bt.rows("t.csv")


def test_computed_rejects_a_row_that_is_not_computed(bt):
    with pytest.raises(SystemExit, match="not COMPUTED"):
        bt.computed([{"origin": "COMPUTED"}, {"origin": "DERIVED"}], "t.csv")


def test_phase_one_two_pi_relative_value_is_rho_over_the_scale(bt):
    """Reviewer B1: PhaseOneTwoPiRel is |rho| / (1 / (1440 pi^2 xi^4)) at the b = 2 pi row with the
    largest |rho|, recomputed here from outputs/validation.csv, not max_err_over_scale."""
    import csv
    rows = [r for r in csv.DictReader(open(bt.OUT / "validation.csv", encoding="utf-8")) if r["b"] == "2pi"]
    worst = max(rows, key=lambda r: abs(mp.mpf(r["rho_computed"])))
    with mp.workdps(40):
        expected = abs(mp.mpf(worst["rho_computed"])) * 1440 * mp.pi ** 2 * mp.mpf(worst["xi"]) ** 4
    macros = dict(re.findall(r"\\newcommand\{\\([A-Za-z]+)\}\{(.*)\}", (bt.GEN / "numbers.tex").read_text(encoding="utf-8")))
    assert macros["PhaseOneTwoPiRel"] == bt.tex(expected)
    assert macros["PhaseOneTwoPiRel"] != bt.tex(worst["max_err_over_scale"])
