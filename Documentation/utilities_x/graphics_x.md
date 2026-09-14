# Graphics (`Utilities_x.graphics_x`)

_Updated: 2026-09-13_

Two plotting helpers are used by the radial-solver solution's `plot_ys` and `plot_interior` methods call, and they can also be used directly on arrays.

A radial-function plot is an easy, visual instability check for a Love-number solve. Spikes, sustained oscillations, or curves that do not vary smoothly with radius mean the integration did not converge, and that shows up in a glance at the figure.

| Function | Draws |
|---|---|
| `plot_ys` | The six radial functions of one or more solutions, with optional published benchmark curves. |
| `plot_interior` | A body's interior profiles: gravity and density, pressure and optionally temperature, and the moduli. |
| `load_benchmark_ys` | The digitized Enceladus curves of Tobie et al. (2005) and Roberts and Nimmo (2008), as arrays. |

## Radial Functions

```python
import numpy as np
from TidalPy.RadialSolver_x import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.rheology_x import Elastic, Maxwell
from TidalPy.Utilities_x.graphics_x import plot_ys

build_data = build_rs_input_homogeneous_layers(
    1600.0e3, 2.0 * np.pi / (86400.0 * 1.37), (3500.0,), (1.2e11,), (6.7e10,), (1e30,), (1e20,),
    ("solid",), (False,), (False,), Maxwell(), Elastic(),
    radius_fraction_tuple=(1.0,), slice_per_layer=80)
solution = radial_solver(*build_data, degree_l=2, solve_for=("tidal", "loading"))

# From the solution: one set of curves per solved boundary-condition type.
figure, axes = solution.plot_ys(show_plot=False)

# From arrays: an (N, 6) or (6, N) complex array, or a list of them, plus the radius grid.
radial_functions = solution.get_radial_solution_array(solution.radius_array, 0)
figure, axes = plot_ys(radial_functions, solution.radius_array,
                       labels=["Enceladus"], benchmarks="tobie2005", use_tobie_limits=True)
```

Note the difference between the two forms. The solution method knows how many boundary-condition types were solved and splits them itself. The bare function takes one solution's six functions at a time, so pass `get_radial_solution_array(radius, ytype_index)` rather than the raw result array when more than one type was requested.

| Argument | Meaning |
|---|---|
| `radial_solutions`, `radius` | One `(6, N)` array, an `(N, 6)` array which is transposed for you, or a list of either, plus one radius array \[m\] shared by all or one per solution. |
| `labels`, `colors`, `line_styles` | Per-solution legend labels, colors, and line styles. A single color or style applies to all. |
| `depth_plot`, `planet_radius` | Plot against depth instead of radius; the planet radius \[m\] is then required. |
| `plot_imaginary` | Also draw the imaginary parts, dotted, on a twin axis in each panel. |
| `benchmarks` | `"tobie2005"` and `"roberts_nimmo2008"`, or the aliases `"t05"` and `"rn08"`, to overlay the published curves. |
| `use_tobie_limits`, `x_limits`, `y_limits` | The axis limits used by Tobie et al. (2005), explicit per-panel limits, or radius and depth limits in km. |
| `figure_size`, `show_plot` | Figure size in inches, and whether to call `plt.show()` before returning. The default is not to. |

The returned `axes` is a two-by-three array. The first three radial functions across the top row and the last three below. A legend appears whenever more than one curve is drawn.

### Benchmark data

`load_benchmark_ys(name)` returns a nested dict keyed first by radial function, `y1` through `y4`, then by model. Tobie et al. (2005) supplies a homogeneous model (`HG`) and two liquid-core models (`LC1`, `LC2`); Roberts and Nimmo (2008) supplies a homogeneous model (`HG`) and one liquid-core model (`LC`). Each entry is a `(values, radius)` pair. The digitized data files ship with the package, and `Benchmarks_x/RadialSolver/Enceladus_Tobie_Roberts.ipynb` reproduces both published figures.

## Interior profiles

```python
from TidalPy.Utilities_x.graphics_x import plot_interior

# From a solution; the equation-of-state solve must have succeeded.
figure, axes = solution.plot_interior(show_plot=False, planet_name="Enceladus")

# From arrays: MKS in, km and GPa on the axes.
figure, axes = plot_interior(
    solution.radius_array, solution.gravity_array, solution.pressure_array, solution.density_array,
    shear_modulus=solution.shear_modulus_array, bulk_modulus=solution.bulk_modulus_array,
    planet_radius=solution.radius, bulk_density=solution.density_bulk, depth_plot=True)
```

The panels are gravity with density on a twin axis; pressure, with an optional temperature twin axis; and, when either modulus is supplied, the moduli in GPa. Real parts are solid lines, the imaginary parts of complex moduli are dotted on a twin axis, and the modulus panel is log-scaled when every value is positive.

`use_scatter` draws points instead of lines. `annotate`, on by default, labels the surface gravity, central pressure, and bulk density. `planet_name` becomes the title. Styling, meaning colors, line styles, marker size, fonts, and panel size, lives in the `INTERIOR_PLOT_STYLE` dictionary; edit it in place to restyle every plot the module draws.
