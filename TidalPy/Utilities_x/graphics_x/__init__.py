"""TidalPy.Utilities_x.graphics_x - plotting helpers for the new backend.

`plot_ys` draws the six radial functions (y1..y6) of one or more radial-solver solutions, optionally
against published benchmark curves; `plot_interior` draws a planet's interior profiles (gravity, density,
pressure, optional temperature and moduli). Both return the matplotlib figure and axes so callers can
adjust them further. `RadialSolverSolution.plot_ys` and `.plot_interior` call these with the solution's
own arrays.
"""

from TidalPy.Utilities_x.graphics_x.radial_functions import (
    BENCHMARK_YS,
    TOBIE2005_X_LIMITS,
    load_benchmark_ys,
    plot_ys,
)
from TidalPy.Utilities_x.graphics_x.interior import (
    INTERIOR_PLOT_STYLE,
    plot_interior,
)

__all__ = [
    "BENCHMARK_YS",
    "INTERIOR_PLOT_STYLE",
    "TOBIE2005_X_LIMITS",
    "load_benchmark_ys",
    "plot_interior",
    "plot_ys",
]
