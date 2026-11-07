"""Phase 1: validate Li & Gott (1998) eq. 12 against the eq. 8 image sum.
Phase 2 (partial): the Minkowski vacuum with images, numerics only. The comparison
with Hiscock & Konkowski (1982) is pending: the article has not been read yet.

    python -m chronosim.run            # writes outputs/
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import mpmath as mp

from .config import DPS
from .misner import closed_form_stress, self_consistent_b
from .pointsplit import minkowski_images_stress, misner_stress
from .provenance import reported_heuristics

COMPONENTS = ("rho", "p_xi", "p_y", "p_z")
GRID_XI = ("0.5", "1", "2")
GRID_B = ("1", "3", "2pi", "10", "30")
SCALING_XI = ("1e-3", "1e-2", "1e-1", "1")
SCAN_B = ("1", "2", "3", "4", "5", "2pi", "7", "8", "10", "15", "20")
SCAN_VACUA_B = ("0.8", "1", "1.5", "2", "3", "4", "5", "2pi", "8", "10", "30")


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


def minkowski_grid(rows_out):
    """Minkowski vacuum with images on the phase 1 grid. Returns the worst finite-difference error."""
    worst = mp.mpf(0)
    for xs in GRID_XI:
        for bs in GRID_B:
            xi, b = _num(xs), _num(bs)
            num = minkowski_images_stress(xi, b)
            finer = minkowski_images_stress(xi, b, delta="1e-11", h="1e-12")
            rho = num["rho"].value
            err = max(abs(num[k].value - finer[k].value) for k in COMPONENTS) / abs(rho)
            worst = max(worst, err)
            rows_out.append({
                "xi": xs, "b": bs,
                "rho_computed": mp.nstr(rho, 20),
                "xi4_rho": mp.nstr(xi ** 4 * rho, 17),
                "fd_err_over_rho": mp.nstr(err, 3),
                "rho_minus_3p_over_rho": mp.nstr(abs(rho - 3 * num["p_xi"].value) / abs(rho), 3),
                "trace_over_rho": mp.nstr(abs(num["trace"].value) / abs(rho), 3),
                "max_offdiag_over_rho": mp.nstr(num["max_offdiag"].value / abs(rho), 3),
                "decay_rate_over_b_minus_1": mp.nstr(num["decay_rate_over_b"].value - 1, 3),
                "power_law_spread": mp.nstr(num["power_law_spread"].value, 3),
                "tail_fraction": mp.nstr(num["tail_fraction"].value, 3),
                "n_used": int(num["n_used"].value),
            })
    return worst


def vacua_scan():
    """xi^4 rho of both vacua at xi = 1, and their difference D, with its finite-difference error."""
    rows = []
    for bs in SCAN_VACUA_B:
        b = _num(bs)
        mink = minkowski_images_stress(1, b)["rho"].value
        finer = minkowski_images_stress(1, b, delta="1e-11", h="1e-12")["rho"].value
        lg = closed_form_stress(1, b)["rho"].value
        D = mink - lg
        rows.append({
            "b": bs,
            "minkowski": mp.nstr(mink, 15),
            "li_gott": mp.nstr(lg, 15),
            "difference": mp.nstr(D, 8),
            "fd_err": mp.nstr(abs(mink - finer), 3),
            "difference_b4_exp": mp.nstr(D * b ** 4 * mp.exp(4 * mp.pi ** 2 / b), 8),
        })
    return rows


def scaling(b=mp.mpf(1), stress=misner_stress):
    """log-log slope of |rho| versus xi from COMPUTED points."""
    pts = [(mp.log(_num(x)), mp.log(abs(stress(_num(x), b)["rho"].value))) for x in SCALING_XI]
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


def report(out: Path, rows, worst, slope, mink_rows, mink_worst, mink_slope, scan):
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
        "## Fase 2 (parcial): vácuo de Minkowski com imagens",
        "",
        "**Comparação com Hiscock & Konkowski (1982): PENDENTE.** O artigo ainda não foi lido, e nenhuma "
        "fórmula dele foi usada. Tudo abaixo é só a soma numérica G = Σ_{n≠0} G_M(X, Λⁿ X′) com o "
        "point-splitting da eq. 11. As referências a Li & Gott vêm do arXiv, com a auditoria contra o PRL "
        "também pendente.",
        "",
        "### Resultado",
        "",
        f"- Maior erro de diferenças finitas no grid ({len(mink_rows)} pontos, comparando δ = 1e-10, "
        f"h = 1e-11 com δ = 1e-11, h = 1e-12), relativo a |ρ|: **{mp.nstr(mink_worst, 3)}** [COMPUTED]",
        f"- Inclinação de log|ρ| contra log ξ (b = 1, ξ de 10⁻³ a 1): −4 com desvio de "
        f"**{mp.nstr(mink_slope + 4, 3)}** [COMPUTED]",
        "- Lei de decaimento medida: exponencial, não n⁻⁴. No último termo somado, −ln(r)/b difere de 1 por "
        f"no máximo **{mp.nstr(max(abs(mp.mpf(r['decay_rate_over_b_minus_1'])) for r in mink_rows), 3)}** [COMPUTED]",
        "",
        "| ξ | b | ρ calculado | ξ⁴ρ | erro DF / ρ | (ρ−3p) / ρ | traço / ρ | fração da cauda | N |",
        "|---|---|---|---|---|---|---|---|---|",
        *[
            f"| {r['xi']} | {r['b']} | {r['rho_computed']} | {r['xi4_rho']} | {r['fd_err_over_rho']} | "
            f"{r['rho_minus_3p_over_rho']} | {r['trace_over_rho']} | {r['tail_fraction']} | {r['n_used']} |"
            for r in mink_rows
        ],
        "",
        "### Os dois vácuos em ξ = 1",
        "",
        "D = ρ(Minkowski com imagens) − ρ(eq. 12). A coluna seguinte é uma **observação empírica**, sem "
        "derivação nem fonte: D·b⁴·exp(4π²/b) fica quase constante em b ≈ 1–3 e sobe para b maior. Em b = 0,8, "
        "D é só cerca de 3 vezes o erro de diferenças finitas, e a coluna não é confiável ali.",
        "",
        "| b | ξ⁴ρ Minkowski [COMPUTED] | ξ⁴ρ eq. 12 [DERIVED] | D [COMPUTED] | erro DF | D·b⁴·exp(4π²/b) |",
        "|---|---|---|---|---|---|",
        *[
            f"| {r['b']} | {r['minkowski']} | {r['li_gott']} | {r['difference']} | {r['fd_err']} | "
            f"{r['difference_b4_exp']} |"
            for r in scan
        ],
        "",
        "### O que esta parte estabelece e o que não estabelece",
        "",
        "- **Estabelece:** o valor numérico de ⟨T_ab⟩ para a soma de imagens de Minkowski com o operador da "
        "eq. 11, com convergência em segunda ordem, a estrutura ρ = 3p sem traço, escala ξ⁻⁴ e decaimento "
        "exponencial dos termos.",
        "- **Não estabelece:** que este resultado reproduz Hiscock & Konkowski. Essa é a comparação pendente. "
        "Também não estabelece a forma de D(b), que é só um ajuste observado, nem nada sobre viagens no tempo.",
        "",
        "## Constantes numéricas escolhidas por nós (HEURISTIC)",
        "",
        *[f"- `{k}` = {mp.nstr(v.value, 5)}: {v.source}" for k, v in h.items()],
        "",
    ]
    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="chronosim: Misner space validation")
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
    mink_rows = []
    mink_worst = minkowski_grid(mink_rows)
    with open(out / "minkowski_images.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(mink_rows[0]))
        w.writeheader()
        w.writerows(mink_rows)
    mink_slope, _ = scaling(stress=minkowski_images_stress)
    scan = vacua_scan()
    if not args.no_figures:
        figures(out)
    report(out, rows, worst, slope, mink_rows, mink_worst, mink_slope, scan)

    print(f"pior erro relativo à escala: {mp.nstr(worst, 3)} [COMPUTED]")
    print(f"inclinação perto do horizonte: {mp.nstr(slope, 15)} [COMPUTED] (eq. 12: -4)")
    print(f"fase 2 (parcial), maior erro de diferenças finitas / |rho|: {mp.nstr(mink_worst, 3)} [COMPUTED]")
    print(f"fase 2 (parcial), inclinação perto do horizonte: -4 com desvio {mp.nstr(mink_slope + 4, 3)} [COMPUTED]")
    print("fase 2: comparação com Hiscock & Konkowski PENDENTE (artigo não lido)")
    print("constantes HEURISTIC usadas:")
    for k, v in reported_heuristics().items():
        print(f"  {k} = {mp.nstr(v.value, 5)}  ({v.source})")
    print(f"relatório: {out / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
