"""Anti-fabrication guard.

No module of the package may reach a random-number generator: a randomly
generated or noise-perturbed value cannot enter a scientific result unnoticed.
The test walks the syntax tree of every module and fails the build otherwise.
"""

import ast
from pathlib import Path

import chronosim

PKG = Path(chronosim.__file__).parent
FORBIDDEN_MODULES = {"random", "secrets", "numpy.random", "uuid"}
FORBIDDEN_ATTRS = {"random", "rand", "randn", "randint", "default_rng", "randomstate", "seed", "urandom"}


def _violations(path: Path):
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name in FORBIDDEN_MODULES or a.name.startswith("numpy.random"):
                    found.append(f"import {a.name}")
        elif isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            if mod in FORBIDDEN_MODULES or mod.startswith("numpy.random"):
                found.append(f"from {mod} import ...")
            for a in node.names:
                if a.name.lower() in FORBIDDEN_ATTRS:
                    found.append(f"from {mod} import {a.name}")
        elif isinstance(node, ast.Attribute) and node.attr.lower() in FORBIDDEN_ATTRS:
            found.append(f".{node.attr} (line {node.lineno})")
    return found


def test_no_module_can_reach_a_random_number_generator():
    files = sorted(PKG.rglob("*.py"))
    assert files, "package not found"
    bad = {str(f.relative_to(PKG)): v for f in files if (v := _violations(f))}
    assert not bad, f"random-number generators reachable: {bad}"


def test_guard_detects_a_planted_violation(tmp_path):
    p = tmp_path / "bad.py"
    p.write_text("import numpy as np\nx = np.random.normal()\nimport mpmath\ny = mpmath.rand()\n")
    v = _violations(p)
    assert any("random" in s for s in v) and any("rand" in s for s in v)
