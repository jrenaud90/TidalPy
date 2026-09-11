"""Plot the six radial functions (y1..y6) of radial-solver solutions.

The radial functions describe a planet's viscoelastic-gravitational response to a unit degree-l potential:
radial displacement (y1), radial stress (y2), tangential displacement (y3), tangential stress (y4), the
gravitational potential perturbation (y5), and the potential stress (y6). `plot_ys` draws them in a
two-by-three panel figure for one or more solutions, and can overlay the published Enceladus curves of
Tobie et al. (2005) and Roberts & Nimmo (2008) for benchmarking.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

# =====================================================================================================================
# Constants
# =====================================================================================================================

DATA_DIRECTORY = Path(__file__).resolve().parent / "data"

# Published benchmark curves: name -> (data file, plotted color, per-series marker, series label).
BENCHMARK_YS: Dict[str, dict] = {
    "tobie2005": {
        "file": "T05-Data.csv",
        "color": "r",
        "label": "T05",
        "markers": {"HG": ".", "LC1": "+", "LC2": "1"},
        "reference": "Tobie, Mocquet & Sotin (2005), Icarus 177, 534-549; Enceladus models.",
    },
    "roberts_nimmo2008": {
        "file": "RN08-Data.csv",
        "color": "b",
        "label": "RN08",
        "markers": {"HG": ".", "LC": "+"},
        "reference": "Roberts & Nimmo (2008), Icarus 194, 675-689; Enceladus models.",
    },
}
BENCHMARK_ALIASES = {"t05": "tobie2005", "tobie": "tobie2005", "rn08": "roberts_nimmo2008", "roberts": "roberts_nimmo2008"}

# Axis limits used by Tobie et al. (2005) for y1..y4 (y5, y6 left automatic).
TOBIE2005_X_LIMITS: Tuple[Optional[Tuple[float, float]], ...] = (
    (0.0, 0.15), (-2100.0, 4500.0), (-0.04, 0.04), (0.0, 2000.0), None, None)

Y_TITLES = ("Radial Disp.", "Radial Stress", "Tang. Disp.", "Tang. Stress", "Grav. Potential Perturb.", "Potential Stress")
Y_UNITS = ("s$^{2}$ / m", "kg / m$^{3}$", "s$^{2}$ / m", "kg / m$^{3}$", "unitless", "1 / m")

ArrayLike = Union[np.ndarray, Sequence[float]]


# =====================================================================================================================
# Benchmark data
# =====================================================================================================================

def load_benchmark_ys(name: str) -> Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]]:
    """Load a published set of radial-function curves.

    Parameters
    ----------
    name : str
        Benchmark name: ``"tobie2005"`` (aliases ``"t05"``, ``"tobie"``) or ``"roberts_nimmo2008"``
        (aliases ``"rn08"``, ``"roberts"``).

    Returns
    -------
    dict
        ``{"y1": {"HG": (values, radius_m), ...}, "y2": ..., "y3": ..., "y4": ...}``. Each series is a pair of
        arrays: the radial-function values and the radii [m] they were digitized at (NaN padding removed).
        Only y1..y4 were published.

    Assumptions
    -----------
    - The curves were digitized from the published figures; expect digitization-level scatter.
    """
    key = BENCHMARK_ALIASES.get(name.lower(), name.lower())
    if key not in BENCHMARK_YS:
        raise ValueError(f"Unknown benchmark {name!r}; choose from {sorted(BENCHMARK_YS)}.")
    table = np.genfromtxt(DATA_DIRECTORY / BENCHMARK_YS[key]["file"], delimiter=",", names=True)
    curves: Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]] = {}
    for column in table.dtype.names:
        y_name, series, axis = column.split("_")
        if axis != "x":
            continue
        values = np.asarray(table[column], dtype=np.float64)
        radius = np.asarray(table[f"{y_name}_{series}_y"], dtype=np.float64)
        keep = np.isfinite(values) & np.isfinite(radius)
        curves.setdefault(y_name, {})[series] = (values[keep], radius[keep])
    return curves


# =====================================================================================================================
# Input handling
# =====================================================================================================================

def _as_solution_list(radial_solutions, radius) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """Normalize the solution/radius inputs to equal-length lists of ``(6, N)`` and ``(N,)`` arrays."""
    single = isinstance(radial_solutions, np.ndarray) or (
        len(radial_solutions) > 0 and np.ndim(radial_solutions[0]) == 1)
    if single:
        radial_solutions = [radial_solutions]
        radius = [radius]
    elif isinstance(radius, np.ndarray) and radius.ndim == 1:
        radius = [radius] * len(radial_solutions)
    if len(radial_solutions) != len(radius):
        raise ValueError(
            f"Provide one radius array per radial solution (got {len(radial_solutions)} solutions and "
            f"{len(radius)} radius arrays).")
    solutions_out, radius_out = [], []
    for index, (solution, radius_array) in enumerate(zip(radial_solutions, radius)):
        solution = np.asarray(solution)
        radius_array = np.asarray(radius_array, dtype=np.float64).ravel()
        if solution.ndim != 2 or 6 not in solution.shape:
            raise ValueError(
                f"Radial solution {index} must be a (6, N) or (N, 6) array; found shape {solution.shape}.")
        if solution.shape[0] != 6:
            solution = solution.T
        if solution.shape[1] != radius_array.size:
            raise ValueError(
                f"Radial solution {index} has {solution.shape[1]} radial points but its radius array has "
                f"{radius_array.size}.")
        solutions_out.append(solution)
        radius_out.append(radius_array)
    return solutions_out, radius_out


def _per_solution(values, count: int, default: Sequence, name: str) -> List:
    """Broadcast a per-solution style option (list or single value) to ``count`` entries."""
    if values is None:
        return [default[i % len(default)] for i in range(count)]
    if isinstance(values, str):
        return [values] * count
    values = list(values)
    if len(values) != count:
        raise ValueError(f"`{name}` must have one entry per radial solution ({count}); found {len(values)}.")
    return values


# =====================================================================================================================
# Plotting
# =====================================================================================================================

def plot_ys(
        radial_solutions: Union[ArrayLike, Sequence[ArrayLike]],
        radius: Union[ArrayLike, Sequence[ArrayLike]],
        labels: Optional[Sequence[str]] = None,
        colors: Optional[Union[str, Sequence[str]]] = None,
        line_styles: Optional[Union[str, Sequence[str]]] = None,
        depth_plot: bool = False,
        planet_radius: Optional[float] = None,
        plot_imaginary: bool = False,
        benchmarks: Union[str, Sequence[str]] = (),
        use_tobie_limits: bool = False,
        x_limits: Optional[Sequence[Optional[Tuple[float, float]]]] = None,
        y_limits: Optional[Tuple[float, float]] = None,
        figure_size: Tuple[float, float] = (9.5, 9.5),
        show_plot: bool = False,
        ) -> Tuple[Figure, np.ndarray]:
    """Plot the six radial functions y1..y6 of one or more radial-solver solutions.

    Parameters
    ----------
    radial_solutions : array or sequence of arrays
        Complex radial functions as a ``(6, N)`` array (``(N, 6)`` is transposed automatically), or a
        sequence of such arrays, one per solution.
    radius : array or sequence of arrays
        Radius [m] of each radial point; one array shared by every solution, or one array per solution.
    labels : sequence of str, optional
        Legend label per solution (default ``"Solution 0"``, ``"Solution 1"``, ...).
    colors, line_styles : str or sequence of str, optional
        Matplotlib color / line style per solution (a single value applies to all). Defaults follow the
        active color cycle and solid lines; imaginary parts use dotted lines.
    depth_plot : bool, default False
        Plot against depth (``planet_radius - radius``) instead of radius.
    planet_radius : float, optional
        Planet radius [m]; required when `depth_plot` is True.
    plot_imaginary : bool, default False
        Also draw the imaginary parts on a twin (top) x-axis of each panel.
    benchmarks : str or sequence of str, default ()
        Published curves to overlay (see `load_benchmark_ys`): ``"tobie2005"`` and/or ``"roberts_nimmo2008"``.
    use_tobie_limits : bool, default False
        Use the y1..y4 axis limits of Tobie et al. (2005); ignored when `x_limits` is given.
    x_limits : sequence of 6 (min, max) pairs or None, optional
        Explicit x-axis limits per panel (``None`` entries stay automatic).
    y_limits : (min, max), optional
        Radius/depth limits [km] applied to every panel.
    figure_size : (width, height), default (9.5, 9.5)
        Figure size in inches.
    show_plot : bool, default False
        Call ``matplotlib.pyplot.show()`` before returning.

    Returns
    -------
    figure : matplotlib.figure.Figure
    axes : numpy.ndarray of matplotlib.axes.Axes, shape (2, 3)
        Panels for y1, y2, y3 (top row) and y4, y5, y6 (bottom row).

    Assumptions
    -----------
    - The radial functions follow TidalPy's y1..y6 convention (unit forcing potential, SI units).
    - Radii are in meters; the axes show kilometers.
    """
    solutions, radii = _as_solution_list(radial_solutions, radius)
    count = len(solutions)
    if depth_plot and planet_radius is None:
        raise ValueError("`planet_radius` is required for a depth plot.")
    if isinstance(benchmarks, str):
        benchmarks = (benchmarks,)
    benchmark_keys = [BENCHMARK_ALIASES.get(name.lower(), name.lower()) for name in benchmarks]
    for key in benchmark_keys:
        if key not in BENCHMARK_YS:
            raise ValueError(f"Unknown benchmark {key!r}; choose from {sorted(BENCHMARK_YS)}.")
    if x_limits is not None and len(x_limits) != 6:
        raise ValueError("`x_limits` needs one (min, max) pair (or None) per radial function y1..y6.")
    if y_limits is not None and len(y_limits) != 2:
        raise ValueError("`y_limits` must be a (min, max) pair.")

    labels = _per_solution(labels, count, [f"Solution {i}" for i in range(count)], "labels")
    colors = _per_solution(colors, count, plt.rcParams["axes.prop_cycle"].by_key()["color"], "colors")
    line_styles = _per_solution(line_styles, count, ["-"], "line_styles")

    figure, axes = plt.subplots(nrows=2, ncols=3, figsize=figure_size)
    figure.subplots_adjust(wspace=0.25, hspace=0.45 if plot_imaginary else 0.3)
    panels: List[Axes] = list(axes.ravel())
    imaginary_panels: List[Axes] = [panel.twiny() for panel in panels] if plot_imaginary else []

    vertical_label = "Depth [km]" if depth_plot else "Radius [km]"
    for y_index, panel in enumerate(panels):
        panel.set(xlabel=f"$y_{{{y_index + 1}}}$ [{Y_UNITS[y_index]}]", title=Y_TITLES[y_index])
        if y_index % 3 == 0:
            panel.set_ylabel(vertical_label)
        else:
            panel.yaxis.set_ticklabels([])
    for panel in imaginary_panels:
        panel.set_xlabel("Imaginary part (dotted)", fontsize="small")

    # Solutions.
    for index, (solution, radius_array) in enumerate(zip(solutions, radii)):
        vertical = (planet_radius - radius_array) if depth_plot else radius_array
        vertical_km = vertical / 1000.0
        for y_index, panel in enumerate(panels):
            panel.plot(np.real(solution[y_index]), vertical_km, label=labels[index], c=colors[index],
                       ls=line_styles[index])
            if plot_imaginary:
                imaginary_panels[y_index].plot(np.imag(solution[y_index]), vertical_km, c=colors[index], ls=":")

    # Published benchmark curves (y1..y4 only).
    for key in benchmark_keys:
        spec = BENCHMARK_YS[key]
        for y_name, series in load_benchmark_ys(key).items():
            panel = panels[int(y_name[1:]) - 1]
            for series_name, (values, radius_m) in series.items():
                vertical = (planet_radius - radius_m) if depth_plot else radius_m
                panel.scatter(values, vertical / 1000.0, label=f"{spec['label']} {series_name}", c=spec["color"],
                              marker=spec["markers"].get(series_name, "x"), s=50)

    # Limits.
    if x_limits is not None:
        limits = x_limits
    elif use_tobie_limits:
        limits = TOBIE2005_X_LIMITS
    else:
        limits = (None,) * 6
    for panel, limit in zip(panels, limits):
        if limit is not None:
            panel.set_xlim(limit)
    if y_limits is not None:
        for panel in panels:
            panel.set_ylim(y_limits)

    # One legend for the whole figure when there is more than one curve.
    handles, legend_labels = panels[0].get_legend_handles_labels()
    if len(legend_labels) > 1:
        columns = min(len(legend_labels), 4)
        rows = (len(legend_labels) - 1) // columns + 1
        figure.legend(handles, legend_labels, loc="lower center", ncol=columns, fancybox=True,
                      bbox_to_anchor=(0.5, 0.005))
        figure.subplots_adjust(bottom=0.08 + 0.028 * rows)

    if show_plot:
        plt.show()
    return figure, axes
