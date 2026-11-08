"""The paper uses the English figures (suffix _en); the report uses the Portuguese ones. Both must
plot exactly the same values: only the text may differ."""

import pytest

from chronosim import run


@pytest.fixture(scope="module")
def plotted(tmp_path_factory):
    out = tmp_path_factory.mktemp("figs")
    data = {**run.figures(out), **run.figures_phase2(out)}
    return out, data


@pytest.mark.parametrize("fig", ["fig1", "fig2", "fig3", "fig4"])
def test_english_and_portuguese_figures_plot_the_same_values(plotted, fig):
    _, data = plotted
    pt, en = data[fig]["pt"], data[fig]["en"]
    assert pt and any(xs for xs, _ in pt)
    assert en == pt


@pytest.mark.parametrize("name", ["fig1_horizon_divergence", "fig2_self_consistent_b", "fig3_kt_divergence", "fig4_kt_cutoffs"])
def test_both_language_files_are_written(plotted, name):
    out, _ = plotted
    assert (out / f"{name}.png").stat().st_size > 0
    assert (out / f"{name}_en.png").stat().st_size > 0
