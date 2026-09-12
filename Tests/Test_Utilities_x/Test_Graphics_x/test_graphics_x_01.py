"""Tests for TidalPy.Utilities_x.graphics_x (radial-function and interior plots).

Figures are drawn on the non-interactive Agg backend and closed after each test.
"""
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest
from matplotlib import pyplot as plt

from TidalPy.RadialSolver_x import homogeneous_love_numbers
from TidalPy.rheology_x import Maxwell
from TidalPy.Utilities_x.graphics_x import (
    BENCHMARK_YS,
    TOBIE2005_X_LIMITS,
    load_benchmark_ys,
    plot_interior,
    plot_ys,
)

N = 40
PLANET_RADIUS = 1600.0e3
RADIUS = np.linspace(0.0, PLANET_RADIUS, N)


def _fake_ys(scale=1.0):
    """A smooth complex (6, N) array standing in for a radial solution."""
    x = RADIUS / PLANET_RADIUS
    rows = [scale * (i + 1) * x**(i + 1) * (1.0 + 0.1j) for i in range(6)]
    return np.asarray(rows, dtype=np.complex128)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# =====================================================================================================================
# plot_ys
# =====================================================================================================================

def test_plot_ys_single_solution():
    figure, axes = plot_ys(_fake_ys(), RADIUS)
    assert axes.shape == (2, 3)
    assert all(len(axis.get_lines()) == 1 for axis in axes.ravel())
    assert len(figure.axes) == 6                       # no imaginary twin axes by default
    assert figure.legends == []                        # single curve: no legend
    assert axes[0, 0].get_ylabel() == "Radius [km]"
    assert axes[0, 0].get_lines()[0].get_ydata()[-1] == pytest.approx(PLANET_RADIUS / 1000.0)


def test_plot_ys_transposed_input_and_shared_radius():
    """(N, 6) arrays are accepted and one radius array can serve every solution."""
    figure, axes = plot_ys([_fake_ys().T, _fake_ys(2.0)], RADIUS, labels=("a", "b"))
    assert all(len(axis.get_lines()) == 2 for axis in axes.ravel())
    assert len(figure.legends) == 1
    assert [text.get_text() for text in figure.legends[0].get_texts()] == ["a", "b"]


def test_plot_ys_options():
    figure, axes = plot_ys(
        [_fake_ys(), _fake_ys(0.5)], [RADIUS, RADIUS], colors="k", line_styles=["-", "--"],
        depth_plot=True, planet_radius=PLANET_RADIUS, plot_imaginary=True, use_tobie_limits=True,
        y_limits=(0.0, 1600.0), figure_size=(6.0, 6.0))
    assert len(figure.axes) == 12                      # six panels plus six imaginary twins
    assert axes[0, 0].get_ylabel() == "Depth [km]"
    assert axes[0, 0].get_xlim() == TOBIE2005_X_LIMITS[0]
    assert axes[0, 1].get_xlim() == TOBIE2005_X_LIMITS[1]
    assert axes[0, 0].get_ylim() == (0.0, 1600.0)
    assert axes[0, 0].get_lines()[0].get_ydata()[0] == pytest.approx(PLANET_RADIUS / 1000.0)  # r = 0 is full depth
    assert axes[0, 1].get_lines()[1].get_linestyle() == "--"
    assert figure.get_size_inches()[0] == pytest.approx(6.0)


def test_plot_ys_explicit_x_limits():
    limits = ((0.0, 1.0), None, (-1.0, 1.0), None, None, (0.0, 2.0))
    _, axes = plot_ys(_fake_ys(), RADIUS, x_limits=limits, use_tobie_limits=True)
    assert axes[0, 0].get_xlim() == (0.0, 1.0)
    assert axes[1, 2].get_xlim() == (0.0, 2.0)
    assert axes[0, 1].get_xlim() != TOBIE2005_X_LIMITS[1]   # explicit limits win over the Tobie preset


@pytest.mark.parametrize("name", sorted(BENCHMARK_YS))
def test_load_benchmark_ys(name):
    curves = load_benchmark_ys(name)
    assert set(curves) == {"y1", "y2", "y3", "y4"}
    for series in curves.values():
        assert set(series) == set(BENCHMARK_YS[name]["markers"])
        for values, radius in series.values():
            assert values.shape == radius.shape and values.size > 5
            assert np.all(np.isfinite(values)) and np.all(np.isfinite(radius))
            assert 0.0 <= radius.min() and radius.max() <= 1.7e6      # Enceladus radii [m]


def test_load_benchmark_ys_aliases_and_errors():
    assert load_benchmark_ys("T05").keys() == load_benchmark_ys("tobie2005").keys()
    assert load_benchmark_ys("rn08").keys() == load_benchmark_ys("roberts_nimmo2008").keys()
    with pytest.raises(ValueError, match="Unknown benchmark"):
        load_benchmark_ys("nope")


def test_plot_ys_with_benchmarks():
    figure, axes = plot_ys(_fake_ys(), RADIUS, benchmarks=("tobie2005", "rn08"), labels=["model"])
    # y1..y4 carry the benchmark scatter series (3 Tobie + 2 Roberts-Nimmo); y5, y6 do not.
    assert all(len(axes.ravel()[i].collections) == 5 for i in range(4))
    assert all(len(axes.ravel()[i].collections) == 0 for i in (4, 5))
    legend_texts = [text.get_text() for text in figure.legends[0].get_texts()]
    assert legend_texts[0] == "model" and "T05 HG" in legend_texts and "RN08 LC" in legend_texts


@pytest.mark.parametrize("kwargs, match", (
    (dict(radial_solutions=[_fake_ys(), _fake_ys()], radius=[RADIUS]), "one radius array per"),
    (dict(radial_solutions=_fake_ys()[:5], radius=RADIUS), r"\(6, N\)"),
    (dict(radial_solutions=_fake_ys(), radius=RADIUS[:-1]), "radial points"),
    (dict(radial_solutions=_fake_ys(), radius=RADIUS, depth_plot=True), "planet_radius"),
    (dict(radial_solutions=_fake_ys(), radius=RADIUS, benchmarks="nope"), "Unknown benchmark"),
    (dict(radial_solutions=_fake_ys(), radius=RADIUS, x_limits=((0, 1),)), "x_limits"),
    (dict(radial_solutions=_fake_ys(), radius=RADIUS, y_limits=(0, 1, 2)), "y_limits"),
    (dict(radial_solutions=[_fake_ys(), _fake_ys()], radius=RADIUS, labels=("only one",)), "labels"),
))
def test_plot_ys_invalid_inputs(kwargs, match):
    with pytest.raises(ValueError, match=match):
        plot_ys(**kwargs)


# =====================================================================================================================
# plot_interior
# =====================================================================================================================

def _interior_arrays():
    x = RADIUS / PLANET_RADIUS
    density = 3500.0 - 500.0 * x
    gravity = 1.2 * x
    pressure = 5.0e9 * (1.0 - x**2)
    return dict(radius=RADIUS, gravity=gravity, pressure=pressure, density=density)


def test_plot_interior_two_panels():
    figure, axes = plot_interior(**_interior_arrays(), bulk_density=3300.0)
    assert len(axes) == 2
    assert len(figure.axes) == 3                       # gravity, pressure, density twin
    texts = [text.get_text() for axis in figure.axes for text in axis.texts]
    assert any("g_{s}" in text for text in texts)
    assert any("P_{0}" in text for text in texts)
    assert any("3300.0" in text for text in texts)
    assert figure.axes[0].get_xlabel().startswith("Gravity")


def test_plot_interior_full_options():
    arrays = _interior_arrays()
    shear = (5.0e10 + 1.0e9j) * np.ones(N)
    bulk = 1.0e11 * np.ones(N)
    figure, axes = plot_interior(
        **arrays, temperature=1500.0 * np.ones(N), shear_modulus=shear, bulk_modulus=bulk,
        planet_radius=PLANET_RADIUS, bulk_density=3400.0, planet_name="test world", depth_plot=True,
        use_scatter=True, annotate=False)
    assert len(axes) == 3
    assert len(figure.axes) == 6                       # 3 panels + density, temperature, imaginary twins
    assert axes[0].get_ylabel() == "Depth [km]"
    assert figure.axes[-1].get_xscale() == "log"       # positive imaginary parts use a log axis
    assert axes[2].get_legend() is not None
    assert figure._suptitle.get_text() == "Test World"
    assert all(len(axis.texts) == 0 for axis in figure.axes)
    assert len(axes[0].collections) == 1 and len(axes[0].get_lines()) == 0


def test_plot_interior_real_moduli_linear_imaginary():
    figure, axes = plot_interior(**_interior_arrays(), shear_modulus=5.0e10 * np.ones(N))
    assert len(axes) == 3
    assert len(figure.axes) == 4                       # real modulus only: no imaginary twin
    assert axes[2].get_xlabel() == "Modulus [GPa]"


@pytest.mark.parametrize("override, match", (
    (dict(density=np.ones(N - 1)), "density"),
    (dict(depth_plot=True), "planet_radius"),
    (dict(shear_modulus=np.ones(3)), "shear_modulus"),
))
def test_plot_interior_invalid_inputs(override, match):
    kwargs = _interior_arrays()
    kwargs.update(override)
    with pytest.raises(ValueError, match=match):
        plot_interior(**kwargs)


# =====================================================================================================================
# RadialSolverSolution methods
# =====================================================================================================================

@pytest.fixture(scope="module")
def solution():
    mu = Maxwell().calc_complex_modulus(60.0e9, 1.0e15, 4.1e-5)
    result = homogeneous_love_numbers(1.8216e6, 3529.0, mu, 4.1e-5, num_slices=40, solve_for=("tidal", "loading"))
    assert result.success
    return result


def test_solution_plot_ys(solution):
    figure, axes = solution.plot_ys(show_plot=False)
    assert axes.shape == (2, 3)
    assert all(len(axis.get_lines()) == 2 for axis in axes.ravel())
    assert [text.get_text() for text in figure.legends[0].get_texts()] == ["Tidal", "Loading"]
    figure, _ = solution.plot_ys(show_plot=False, plot_imaginary=True, depth_plot=True, planet_radius=1.8216e6)
    assert len(figure.axes) == 12


def test_solution_plot_interior(solution):
    figure, axes = solution.plot_interior(show_plot=False, planet_name="Io")
    assert len(axes) == 3
    assert figure._suptitle.get_text() == "Io"
    _, axes = solution.plot_interior(show_plot=False, depth_plot=True)
    assert axes[0].get_ylabel() == "Depth [km]"
