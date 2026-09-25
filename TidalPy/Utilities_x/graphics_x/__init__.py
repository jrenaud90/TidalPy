"""Plotting helpers for the new backend.

Every entry point returns the matplotlib figure and axes so callers can adjust them further.
"""

from TidalPy.Utilities_x.graphics_x.radial_functions import (
    BENCHMARK_YS,
    TOBIE2005_X_LIMITS,
    load_benchmark_ys,
    plot_ys,
)
from TidalPy.Utilities_x.graphics_x.interior import (
    INTERIOR_PLOT_STYLE,
    load_interior_plot_style,
    plot_interior,
)
from TidalPy.Utilities_x.graphics_x.maps import (
    MAP_PLOT_STYLE,
    MAP_PROJECTIONS,
    make_map_axes,
    plot_map,
)

__all__ = [
    "BENCHMARK_YS",
    "INTERIOR_PLOT_STYLE",
    "load_interior_plot_style",
    "MAP_PLOT_STYLE",
    "MAP_PROJECTIONS",
    "TOBIE2005_X_LIMITS",
    "load_benchmark_ys",
    "make_map_axes",
    "plot_interior",
    "plot_map",
    "plot_ys",
]
