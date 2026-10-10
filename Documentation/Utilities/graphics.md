# Graphics (`Utilities.graphics`)

_Updated: 2026-10-02_

`TidalPy.Utilities.graphics` plots radial functions (`plot_ys`, with published benchmark curves from `load_benchmark_ys`), interior profiles (`plot_interior`), and surface maps of the 3D tidal fields (`plot_map`, `make_map_axes`). A radial-solver solution has its own `plot_ys` and `plot_interior` methods.

A radial-function plot is a visual check of a Love-number solve: spikes, sustained oscillations, or curves that do not vary smoothly with radius generally mean the integration did not converge.

## Radial Functions

```python
import numpy as np
from TidalPy.RadialSolver import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.Rheology import Elastic, Maxwell
from TidalPy.Utilities.graphics import plot_ys

build_data = build_rs_input_homogeneous_layers(
    1600.0e3,
    2.0 * np.pi / (86400.0 * 1.37),
    (3500.0,),
    (1.2e11,),
    (6.7e10,),
    (1e30,),
    (1e20,),
    ("solid",),
    (False,),
    (False,),
    Maxwell(),
    Elastic(),
    radius_fraction_tuple=(1.0,),
    slice_per_layer=80
)
solution = radial_solver(*build_data, degree_l=2, solve_for=("tidal", "loading"))

# From the solution: one set of curves per solved boundary-condition type.
figure, axes = solution.plot_ys(show_plot=False)

# From arrays: an (N, 6) or (6, N) complex array, or a list of them, plus the radius grid.
radius = solution.sample_radii()
radial_functions = solution.get_radial_solution_array(radius, 0)
figure, axes = plot_ys(
    radial_functions,
    radius,
    labels=["Enceladus"],
    benchmarks="tobie2005",
    use_tobie_limits=True
)
```

The solution method splits the solved boundary-condition types itself. The bare function takes one solution's six functions at a time, so pass `get_radial_solution_array(radius, ytype_index)`, not the raw result array, when more than one type was solved.

| Argument | Meaning |
|---|---|
| `radial_solutions`, `radius` | One `(6, N)` or `(N, 6)` array, or a list of them, plus one radius array \[m\] shared by all or one per solution. |
| `labels`, `colors`, `line_styles` | Per-solution legend labels, colors, and line styles. A single color or style applies to all. |
| `depth_plot`, `planet_radius` | Plot against depth instead of radius. The planet radius \[m\] is then required. |
| `plot_imaginary` | Also draw the imaginary parts, dotted, on a twin axis in each panel. |
| `benchmarks` | `"tobie2005"` and `"roberts_nimmo2008"`, or the aliases `"t05"` and `"rn08"`, to overlay the published curves. |
| `use_tobie_limits`, `x_limits`, `y_limits` | The axis limits used by Tobie et al. (2005), explicit per-panel limits, or radius and depth limits in km. |
| `figure_size`, `show_plot` | Figure size in inches, and whether to call `plt.show()` before returning (default no). |

The returned `axes` is a two-by-three array, the first three radial functions on top and the last three below, with a legend when more than one curve is drawn.

### Benchmark Data

`load_benchmark_ys(name)` returns a nested dict keyed by radial function (`y1` to `y4`), then by model, each entry a `(values, radius)` pair. Tobie et al. (2005) supplies a homogeneous model (`HG`) and two liquid-core models (`LC1`, `LC2`); Roberts and Nimmo (2008) a homogeneous (`HG`) and a liquid-core (`LC`) model. `Benchmarks/RadialSolver/Enceladus_Tobie_Roberts.ipynb` reproduces both published figures.

## Interior Profiles

```python
from TidalPy.Utilities.graphics import plot_interior

# From a solution. The equation-of-state solve must have succeeded.
figure, axes = solution.plot_interior(show_plot=False, planet_name="Enceladus")

# From arrays: MKS in, km and GPa on the axes. The solution can be evaluated at any radius, so pick the grid you want.
radius = solution.sample_radii()
figure, axes = plot_interior(
    radius,
    solution.get_gravity(radius),
    solution.get_pressure(radius),
    solution.get_density(radius),
    shear_modulus=solution.get_shear_modulus(radius),
    bulk_modulus=solution.get_bulk_modulus(radius),
    planet_radius=solution.radius,
    bulk_density=solution.density_bulk,
    depth_plot=True
)
```

The panels are gravity with density on a twin axis; pressure, with an optional temperature twin axis; and, when either modulus is given, the moduli in GPa (log-scaled when every value is positive). Real parts are solid lines and the imaginary parts of complex moduli dotted on a twin axis. The density axis leaves `DENSITY_AXIS_MARGIN` (0.1) of its range clear on each side, so a uniform-density layer does not sit on a spine.

`use_scatter` draws points instead of lines, `annotate` (on by default) labels the surface gravity, central pressure, and bulk density, and `planet_name` becomes the title. Styling (colors, line styles, marker size, fonts, panel size) lives in the `INTERIOR_PLOT_STYLE` dictionary, filled from the `[graphics.interior]` table of `TidalPy_Configs.toml` at import. `load_interior_plot_style()` reads the table again and `reload_interior_plot_style()` refills the dictionary (after a `TidalPy.reinit`, say). Edit the dictionary in place to restyle every plot.

## Surface Maps

`plot_map` draws one colatitude-by-longitude slice of a 3D field as a global map, and `make_map_axes` builds a figure of map panels for it. The field is shaped `(len(colatitudes), len(longitudes))`, the axis order of the world 3D grids: one radius of `calc_3d_tides(...)["heating"]` is `heating[radius_index]`, and one stress component at one radius and time is `calc_3d_stress_strain(...)["stress"][radius_index, :, :, time_index, component_index]`. Cartopy (in the optional `graphics` extra) supplies the projections when installed; otherwise the matplotlib Mollweide and rectangular axes are used. No coastlines or other Natural Earth data are drawn, so nothing is downloaded.

```python
import numpy as np
from TidalPy.Utilities.graphics import make_map_axes, plot_map

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
    symmetric=True
)

# Two panels drawn on one color scale
figure, axes = make_map_axes(
    nrows=1,
    ncols=2
)
for panel, phase in zip(axes.flat, (0.0, 0.5 * np.pi)):
    plot_map(
        longitudes,
        colatitudes,
        np.sin(colatitude_grid)**2 * np.cos(2.0 * longitude_grid - phase),
        axis=panel,
        value_limits=(-1.0, 1.0),
        colormap="RdBu_r",
        title=f"Phase {phase:.2f} rad"
    )
```

Longitudes are wrapped onto the 360 degrees around the map center and sorted, so the data seam falls on the map edge; a longitude that repeats after wrapping keeps its first column. Each sample is drawn as the cell around it, clipped to the globe, and NaN cells are left blank.

| Argument | Meaning |
|---|---|
| `longitudes`, `colatitudes`, `values` | The grid axes \[rad\], colatitude 0 at the north pole, and the field shaped `(len(colatitudes), len(longitudes))`. |
| `projection`, `central_longitude` | `"mollweide"` (the default), `"plate_carree"`, or `"robinson"`, and the longitude at the map center \[deg\]. Robinson and a nonzero center need cartopy. |
| `symmetric`, `value_limits`, `log_scale` | A color scale centered on zero with a diverging colormap, for stress and strain; explicit limits; or a log scale over the positive values, for heating. |
| `title`, `colorbar_label`, `colormap`, `colorbar`, `grid_lines` | Labels and appearance. The default colormaps come from `MAP_PLOT_STYLE`. |
| `axis`, `use_cartopy`, `show_plot` | A panel from `make_map_axes` to draw into; `True` to require cartopy or `False` to skip it (the default uses it when installed); and whether to call `plt.show()` before returning (default no). |

`make_map_axes(nrows, ncols, projection, central_longitude, use_cartopy, figure_size)` returns the figure and a `(nrows, ncols)` array of panels (cartopy map axes when cartopy is used). Figure size, colormaps, grid lines, fonts, and colorbar spacing live in the `MAP_PLOT_STYLE` dictionary; edit it in place to restyle every map.
