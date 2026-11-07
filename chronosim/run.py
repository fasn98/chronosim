"""Phase 1 run: validate Li & Gott (1998) eq. 12 against the eq. 8 image sum.

    python -m chronosim.run            # writes outputs/
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import mpmath as mp

from .config import DPS
from .misner import closed_form_stress, self_consistent_b
from .pointsplit import misner_stress
from .provenance import reported_heuristics

COMPONENTS = ("rho", "p_xi", "p_y", "p_z")
GRID_XI = ("0.5", "1", "2")
GRID_B = ("1", "3", "2pi", "10", "30")
SCALING_XI = ("1e-3", "1e-2", "1e-1", "1")
SCAN_B = ("1", "2", "3", "4", "5", "2pi", "7", "8", "10", "15", "20")


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


def report(out: Path, rows, worst, slope):
    h = reported_heuristics()
    lines = [
        "# chronosim · fase 1 · relatório de validação",
        "",
        "Comparação entre a energia do vácuo **calculada** a partir da soma de imagens "
        "(Li & Gott 1998, eq. 8) com point-splitting (eq. 11) e a **fórmula fechada** (eq. 12), "
        "para um campo escalar conformemente acoplado no espaço de Misner. Unidades naturais (ħ = c = 1).",
        "",
        "## Resultado",
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
        "## Constantes numéricas escolhidas por nós (HEURISTIC)",
        "",
        *[f"- `{k}` = {mp.nstr(v.value, 5)}: {v.source}" for k, v in h.items()],
        "",
        "## O que isto estabelece e o que não estabelece",
        "",
        "- **Estabelece:** a eq. 12 de Li & Gott decorre das eqs. 8 e 11, reproduzida de forma "
        "independente até a precisão numérica, incluindo a divergência ξ⁻⁴ no horizonte (o mecanismo "
        "invocado por Hawking, 1992) e o anulamento em b = 2π (o contra-exemplo de Li & Gott).",
        "- **Não estabelece:** que viagens no tempo sejam possíveis ou impossíveis. O cálculo é "
        "semiclássico (campo quântico em espaço-tempo fixo, sem reação gravitacional) e Kay, Radzikowski "
        "& Wald (1997) provaram que a função de dois pontos não pode ser Hadamard no horizonte de "
        "Cauchy; a questão depende de uma teoria de gravidade quântica.",
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
    report(out, rows, worst, slope)

    print(f"pior erro relativo à escala: {mp.nstr(worst, 3)} [COMPUTED]")
    print(f"inclinação perto do horizonte: {mp.nstr(slope, 15)} [COMPUTED] (eq. 12: -4)")
    print("constantes HEURISTIC usadas:")
    for k, v in reported_heuristics().items():
        print(f"  {k} = {mp.nstr(v.value, 5)}  ({v.source})")
    print(f"relatório: {out / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
