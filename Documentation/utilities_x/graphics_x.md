# Graphics (`Utilities_x.graphics_x`)

_Updated: 2026-09-16_

Two plotting helpers back the radial-solver solution's `plot_ys` and `plot_interior` methods and can also be used directly on arrays. Two more draw surface maps of the 3D tidal fields.

A radial-function plot is a visual instability check for a Love-number solve: spikes, sustained oscillations, or curves that do not vary smoothly with radius mean the integration did not converge.

| Function | Draws |
|---|---|
| `plot_ys` | The six radial functions of one or more solutions, with optional published benchmark curves. |
| `plot_interior` | A body's interior profiles: gravity and density, pressure and optionally temperature, and the moduli. |
| `load_benchmark_ys` | The digitized Enceladus curves of Tobie et al. (2005) and Roberts and Nimmo (2008), as arrays. |
| `plot_map` | One colatitude-by-longitude slice of a 3D field, such as tidal heating or stress, as a global map. |
| `make_map_axes` | A figure of one or more map panels for `plot_map` to draw into. |

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

The returned `axes` is a two-by-three array: the first three radial functions across the top row and the last three below. A legend appears whenever more than one curve is drawn.

### Benchmark Data

`load_benchmark_ys(name)` returns a nested dict keyed first by radial function, `y1` through `y4`, then by model. Tobie et al. (2005) supplies a homogeneous model (`HG`) and two liquid-core models (`LC1`, `LC2`); Roberts and Nimmo (2008) supplies a homogeneous model (`HG`) and one liquid-core model (`LC`). Each entry is a `(values, radius)` pair. The digitized data files ship with the package, and `Benchmarks_x/RadialSolver/Enceladus_Tobie_Roberts.ipynb` reproduces both published figures.

## Interior Profiles

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

## Surface Maps

`plot_map` draws one colatitude-by-longitude slice of a 3D field as a global map, and `make_map_axes` builds a figure of map panels for it to draw into. The field is shaped `(len(colatitudes), len(longitudes))`, the axis order of the world 3D grids: one radius of `calc_3d_tides(...)["heating"]` is `heating[radius_index]`, and one stress component at one radius and time is `calc_3d_stress_strain(...)["stress"][radius_index, :, :, time_index, component_index]`. Cartopy supplies the projections when it is installed, as part of the optional `graphics` extra. Without it the matplotlib Mollweide and rectangular axes are used. No coastlines or other Natural Earth data are drawn, so nothing is downloaded.

```python
import numpy as np
from TidalPy.Utilities_x.graphics_x import make_map_axes, plot_map

longitudes = np.linspace(0.0, 2.0 * np.pi, 73)     # [rad]; the repeated 0 and 2 pi column is drawn once
colatitudes = np.linspace(0.02, np.pi - 0.02, 36)   # [rad]
colatitude_grid, longitude_grid = np.meshgrid(colatitudes, longitudes, indexing="ij")
pattern = np.sin(colatitude_grid)**2 * np.cos(2.0 * longitude_grid)   # Shaped (colatitudes, longitudes)

# One map of a signed field, with the color scale centered on zero
figure, axis = plot_map(
    longitudes,
    colatitudes,
    pattern,
    title="Degree-2 sectoral pattern",
    colorbar_label="Amplitude",
    symmetric=True)

# Two panels drawn on one color scale
figure, axes = make_map_axes(
    nrows=1,
    ncols=2)
for panel, phase in zip(axes.flat, (0.0, 0.5 * np.pi)):
    plot_map(
        longitudes,
        colatitudes,
        np.sin(colatitude_grid)**2 * np.cos(2.0 * longitude_grid - phase),
        axis=panel,
        value_limits=(-1.0, 1.0),
        colormap="RdBu_r",
        title=f"Phase {phase:.2f} rad")
```

Longitudes are wrapped onto the 360 degrees centered on the map center and sorted, so the seam of the data falls on the map edge, and a longitude that repeats after wrapping keeps its first column. Each sample is drawn as the cell around it, with the outer cells clipped to the globe, and NaN cells are left blank.

| Argument | Meaning |
|---|---|
| `longitudes`, `colatitudes`, `values` | The grid axes \[rad\], colatitude 0 at the north pole, and the field shaped `(len(colatitudes), len(longitudes))`. |
| `projection`, `central_longitude` | `"mollweide"` (the default), `"plate_carree"`, or `"robinson"`, and the longitude at the map center \[deg\]. Robinson and a nonzero center need cartopy. |
| `symmetric`, `value_limits`, `log_scale` | A color scale centered on zero with a diverging colormap, for stress and strain; explicit limits; or a log scale over the positive values, for heating. |
| `title`, `colorbar_label`, `colormap`, `colorbar`, `grid_lines` | Labels and appearance. The default colormaps come from `MAP_PLOT_STYLE`. |
| `axis`, `use_cartopy`, `show_plot` | A panel from `make_map_axes` to draw into; `True` to require cartopy or `False` to skip it, where the default uses it when installed; and whether to call `plt.show()` before returning. The default is not to. |

`make_map_axes(nrows, ncols, projection, central_longitude, use_cartopy, figure_size)` returns the figure and a `(nrows, ncols)` array of panels, which are cartopy map axes when cartopy is used. Figure size, colormaps, grid lines, fonts, and colorbar spacing live in the `MAP_PLOT_STYLE` dictionary; edit it in place to restyle every map.
