# Utilities_x: Graphics (`graphics_x`)

_Updated: 2026-09-10_

`TidalPy.Utilities_x.graphics_x` contains plotting helpers. They are mostly matplotlib,
return the figure and axes they draw, and are what the radial-solver solution's
`plot_ys` and `plot_interior` methods call.

| Function | Draws |
|---|---|
| `plot_ys` | The six radial functions y1..y6 of one or more radial-solver solutions, with optional published benchmark curves. |
| `plot_interior` | A planet's interior profiles: gravity and density, pressure (and temperature), shear and bulk moduli. |
| `load_benchmark_ys` | The digitized Enceladus curves of Tobie et al. (2005) and Roberts & Nimmo (2008) as arrays. |

---

## Radial Functions: `plot_ys`

```python
import numpy as np
from TidalPy.RadialSolver_x import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.rheology_x import Elastic, Maxwell
from TidalPy.Utilities_x.graphics_x import plot_ys

build_data = build_rs_input_homogeneous_layers(
    1600.0e3, 2.0 * np.pi / (86400.0 * 1.37), (3500.0,), (1.2e11,), (6.7e10,), (1e30,), (1e20,),
    ("solid",), (False,), (False,), Maxwell(), Elastic(), radius_fraction_tuple=(1.0,), slice_per_layer=80)
solution = radial_solver(*build_data, degree_l=2, solve_for=("tidal", "loading"))

# Straight from the solution: one curve per solved boundary-condition type ("Tidal", "Loading").
figure, axes = solution.plot_ys(show_plot=False)

# Or from arrays: a (6, N) or (N, 6) complex array (or a list of them) plus the radius grid.
figure, axes = plot_ys(solution.result, solution.radius_array, labels=["Enceladus"],
                       benchmarks="tobie2005", use_tobie_limits=True)
```

Parameters:

| Argument | Meaning |
|---|---|
| `radial_solutions`, `radius` | One `(6, N)` array (an `(N, 6)` array is transposed) or a list of them, plus one radius array [m] shared by all or one per solution. |
| `labels`, `colors`, `line_styles` | Per-solution legend labels, colors, and line styles (a single color or style applies to all). |
| `depth_plot`, `planet_radius` | Plot against depth instead of radius (`planet_radius` [m] required). |
| `plot_imaginary` | Also draw the imaginary parts (dotted) on a twin x-axis of each panel. |
| `benchmarks` | `"tobie2005"` and/or `"roberts_nimmo2008"` (aliases `"t05"`, `"rn08"`): overlay the published y1..y4 curves. |
| `use_tobie_limits`, `x_limits`, `y_limits` | The y1..y4 axis limits used by Tobie et al. (2005); or explicit per-panel x limits; radius/depth limits [km]. |
| `figure_size`, `show_plot` | Figure size in inches; call `plt.show()` before returning (default off). |

Returns `(figure, axes)` with `axes` a `(2, 3)` array: y1, y2, y3 on the top row and y4, y5, y6 below. A
legend is added when more than one curve is drawn.

The plot is the quickest instability check for a solve: large spikes, constant sinusoids, or curves that do
not vary smoothly with radius mean the integration did not converge (see the radial solver pages).

### Benchmark Data

`load_benchmark_ys(name)` returns `{"y1": {"HG": (values, radius), ...}, ..., "y4": ...}` for the
homogeneous (`HG`) and liquid-core (`LC`, `LC1`, `LC2`) Enceladus models digitized from Tobie, Mocquet &
Sotin (2005, Icarus 177) and Roberts & Nimmo (2008, Icarus 194). The data files ship with the package
(`Utilities_x/graphics_x/data/`). `Benchmarks_x/RadialSolver/Enceladus_Tobie_Roberts.ipynb` reproduces both
figures.

---

## Interior Profiles: `plot_interior`

```python
from TidalPy.Utilities_x.graphics_x import plot_interior

# From a radial-solver solution (the EOS solve must have succeeded):
figure, axes = solution.plot_interior(show_plot=False, planet_name="Enceladus")

# Or from arrays (SI units in, km / GPa on the axes):
figure, axes = plot_interior(
    solution.radius_array, solution.gravity_array, solution.pressure_array, solution.density_array,
    shear_modulus=solution.shear_modulus_array, bulk_modulus=solution.bulk_modulus_array,
    planet_radius=solution.radius, bulk_density=solution.density_bulk, depth_plot=True)
```

Panels: gravity with density on a twin axis; pressure with an optional `temperature` twin axis; and, when
`shear_modulus` and/or `bulk_modulus` are given, the moduli in GPa (real parts solid, imaginary parts of
complex moduli dotted on a twin axis, log-scaled when they are all positive). `use_scatter` draws points
instead of lines, `annotate` (default on) labels the surface gravity, central pressure, and bulk density,
and `planet_name` becomes the title. Styling (colors, line styles, marker size, fonts, panel size) lives in
the `INTERIOR_PLOT_STYLE` dictionary; edit it in place to restyle every plot.

---

## Migrating from `TidalPy.utilities.graphics`

| Classic | New |
|---|---|
| `from TidalPy.utilities.graphics.multilayer import yplot` | `from TidalPy.Utilities_x.graphics_x import plot_ys` |
| `yplot(ys, radius, plot_tobie=True, plot_roberts=True)` | `plot_ys(ys, radius, benchmarks=("tobie2005", "roberts_nimmo2008"))` |
| `plot_imags=`, `other_xlimits=`, `other_ylimits=` | `plot_imaginary=`, `x_limits=`, `y_limits=` |
| `show_plot=True` default | `show_plot=False` default for the functions; the solution methods keep `show_plot=True` |
| `from TidalPy.utilities.graphics.planet_plot import planet_plot` | `from TidalPy.Utilities_x.graphics_x import plot_interior` |
| `planet_plot(radii, gravity, pressure, density, temperature, shear, bulk, planet_radius, bulk_density, auto_show=...)` | `plot_interior(radius, gravity, pressure, density, temperature=, shear_modulus=, bulk_modulus=, planet_radius=, bulk_density=, show_plot=)` |
| Styling from `TidalPy.config["graphics"]` | `INTERIOR_PLOT_STYLE` |
| `rs_solution.plot_ys()` (no arguments) | `solution.plot_ys(show_plot=True, **plot_kwargs)`; likewise `plot_interior` |

Behavior fixes relative to the classic helpers: invalid input raises `ValueError` (the classic code built
the exceptions without raising them), a legend appears whenever more than one curve is drawn, imaginary twin
axes exist only when requested, density is plotted in kg m⁻³ (the classic plot divided by 1000 but kept the
label), `(N, 6)` radial-function arrays from `get_radial_solution_array` are accepted, and the solution
methods raise an informative error instead of returning `None` when the solve or the EOS solve failed.
