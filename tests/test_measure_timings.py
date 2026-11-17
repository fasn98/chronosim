"""paper/measure_timings.py: LaTeX escaping of system strings and the per-step timeout."""

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture()
def mt():
    spec = importlib.util.spec_from_file_location("measure_timings", ROOT / "paper" / "measure_timings.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_latex_escape_all_special_characters(mt):
    assert mt.latex_escape("\\{}$&#^_%~") == (
        "\\textbackslash{}\\{\\}\\$\\&\\#\\textasciicircum{}\\_\\%\\textasciitilde{}"
    )


def test_latex_escape_leaves_plain_text(mt):
    s = "13th Gen Intel(R) Core(TM) i5-13420H @ 2.10GHz"
    assert mt.latex_escape(s) == s


def test_write_tex_escapes_system_strings(mt, tmp_path, monkeypatch):
    monkeypatch.setattr(mt, "GEN", tmp_path)
    data = {"seconds": {"pytest": 1.0, "run": 2.0, "tables": 0.1, "pdf": 3.0},
            "cpu": "CPU_X & Co #1", "logical_cpus": 4, "ram_gib": 8.0, "os": "OS 100%",
            "kernel": "6.1~rc{1}^$", "python": "3.13.13", "measured": "2026-01-01"}
    mt.write_tex(data)
    tex = (tmp_path / "timings.tex").read_text(encoding="utf-8")
    assert "\\newcommand{\\TimeCPU}{CPU\\_X \\& Co \\#1}" in tex
    assert "\\newcommand{\\TimeOS}{OS 100\\%}" in tex
    assert "\\newcommand{\\TimeKernel}{6.1\\textasciitilde{}rc\\{1\\}\\textasciicircum{}\\$}" in tex


def test_timeout_stops_and_writes_nothing(mt, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(mt, "GEN", tmp_path)
    monkeypatch.setattr(mt, "STEP_TIMEOUT_S", 0.5)
    monkeypatch.setattr(mt, "STEPS", [("slow", "sleeper", [sys.executable, "-c", "print('partial', flush=True); import time; time.sleep(10)"])])
    monkeypatch.setattr(sys, "argv", ["measure_timings.py"])
    with pytest.raises(SystemExit, match="exceeded the 0.5 s limit per step; no timings written"):
        mt.main()
    assert not (tmp_path / "timings.json").exists()
    assert not (tmp_path / "timings.tex").exists()
    assert "partial" in capsys.readouterr().out
