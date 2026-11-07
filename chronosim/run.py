"""Phase 1: validate Li & Gott (1998) eq. 12 against the eq. 8 image sum.
Phase 2 (in progress): Kim & Thorne (1991), arbiter of eq. 4a and item F1.

    python -m chronosim.run            # writes outputs/
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import mpmath as mp

from . import kim_thorne as kt
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


def report(out: Path, rows, worst, slope, arb, f1_rows, f1_exact):
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
        "Plano em `docs/phase2_plan.md`, transcrições em `docs/equation_transcription.md`. Até aqui: teste "
        "árbitro do operador e item F1. F3, F4, F2, F5, F6 e F7 ainda não foram implementados.",
        "",
        "### Teste árbitro: eq. 4a de Kim & Thorne contra a eq. 11 de Li & Gott",
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
        f"**Discrepância na eq. 31 (p. 3940).** A cópia T_n da garganta, obtida das eqs. 25, 28 e 29, satisfaz "
        "a relação da eq. 31 com o sinal oposto do coeficiente impresso: a razão medida entre "
        "(T̃+c)/(X̃+c) e (1+ξ²ⁿ)/(1−ξ²ⁿ) é −1, com desvio máximo de "
        f"{mp.nstr(f1_exact['eq31'], 3)} (n = 1, 2, 3). As duas transcrições leem o sinal +. Nenhuma "
        "fórmula usada adiante depende da eq. 31.",
        "",
        "Erro relativo da eq. 38 contra o intervalo exato. Ele cai linearmente com ε, como esperado da "
        "aproximação da eq. 35:",
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
        "### O que esta parte estabelece e o que não estabelece",
        "",
        "- **Estabelece:** as eqs. 33, 34 e 37 são exatas para a geometria das eqs. 25 e 28, e a eq. 38 é "
        "uma aproximação cujo erro relativo é proporcional a ε/D, sobre o eixo de simetria.",
        "- **Não estabelece:** nada ainda sobre ⟨T_μν⟩ de Kim & Thorne. A eq. 38 é só o ponto de partida dos "
        "itens seguintes, e o caso fora do eixo (eqs. 46 e 61) não foi testado. O árbitro mostra o desacordo "
        "entre as convenções impressas, mas não diz qual dos artigos o causa.",
        "- **Disputa com Hawking:** Hawking (1992, p. 603 e seção VI) contesta a conclusão de Kim & Thorne. "
        "Ele argumenta que um eventual corte deveria ser em distância invariante, e não no tempo próprio, "
        "que depende do observador. Esta fase reproduz o cálculo de Kim & Thorne e não decide quem estava certo.",
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
    report(out, rows, worst, slope, arb, f1_rows, f1_exact)

    print(f"pior erro relativo à escala: {mp.nstr(worst, 3)} [COMPUTED]")
    print(f"inclinação perto do horizonte: {mp.nstr(slope, 15)} [COMPUTED] (eq. 12: -4)")
    print("fase 2, árbitro: eq. 4a de Kim & Thorne / eq. 11 de Li & Gott no vácuo de Rindler, maior desvio de 2: "
          f"{mp.nstr(max((abs(r['worst_ratio_minus_2']) for r in arb)), 3)} [COMPUTED]")
    print(f"fase 2, F1: eq. 33 {mp.nstr(f1_exact['eq33'], 3)}, eq. 34 {mp.nstr(f1_exact['eq34'], 3)}, "
          f"eq. 37 {mp.nstr(f1_exact['eq37'], 3)} [COMPUTED]")
    print(f"fase 2, F1: eq. 31 como impressa difere por um sinal, desvio de -1: {mp.nstr(f1_exact['eq31'], 3)} [COMPUTED]")
    print("constantes HEURISTIC usadas:")
    for k, v in reported_heuristics().items():
        print(f"  {k} = {mp.nstr(v.value, 5)}  ({v.source})")
    print(f"relatório: {out / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
