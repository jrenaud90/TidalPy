# TidalPy Configurations

_Updated: 2026-09-23_

TidalPy's settings and parameters are read when the package is first imported. They live in a configuration file in the user's documents directory, whose location varies by operating system.

**Windows**

"C:\\Users\\\<username\>\\Documents\\TidalPy\\\<TidalPy Version\>\\Config\\TidalPy_Configs.toml"

**MacOS**

"/Users/\<username\>/Documents/TidalPy/\<TidalPy Version\>/Config/TidalPy_Configs.toml"

**Linux**

"/home/\<username\>/Documents/TidalPy/\<TidalPy Version\>/Config/TidalPy_Configs.toml"

The `toml` file holds all of TidalPy's settings, with comments giving context for each one. Changes reach TidalPy the next time it is imported, so you will likely need to restart your kernel to see them.

## Original Configuration File

The first time TidalPy is loaded it copies over the contents of ["defaultc.py"](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/defaultc.py). Developers who want to change or add a configuration should edit the string in that file.

## Overriding Defaults with a New Config File

To change configurations for one project while leaving the default config file alone, provide a new configuration file:

```python
import TidalPy

print(TidalPy.configs["logging"]["write_log_to_disk"])  # by default this is false.

# Say we made a new config file (e.g., "TidalPy_Config2.toml"). It can contain any number of specific configuration changes.
# Only those present will overwrite the default configs (e.g., this logging.write_log_to_disk config).
TidalPy.reinit("TidalPy_Config2.toml")

print(TidalPy.configs["logging"]["write_log_to_disk"])  # Should now show true.

# You can also do this by providing a dictionary instead of a file.
TidalPy.reinit({"logging": {"write_log_to_disk": True}})
```

Overrides need to be performed again after a kernel restart. For permanent changes, edit the config in your documents directory (make a copy first).

## New-Backend Configuration (`TidalPy_Configs_x.toml`)

The new C++ backend (the `_x` modules) reads its settings from a second file in the same directory, `TidalPy_Configs_x.toml`, loaded to `TidalPy.config_x`. Until the two files merge in TidalPy 0.9.0, the package-wide settings on this page (logging, pathing, debugging) stay in `TidalPy_Configs.toml`, while the new backend's numerical settings, tidal defaults, world defaults, and per-material layer defaults live in `TidalPy_Configs_x.toml`. The packaged defaults are in ["defaultc_x.py"](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/defaultc_x.py).

TidalPy loads the packaged defaults first and merges your file over them, so your file only needs the values you change, and a default added in a later release reaches an existing file without regenerating it. Tables merge key by key and any other value (a list included) replaces the default whole. A physics-model table that names a different `model` than the default replaces the default table instead of merging with it, so no parameter of the default model carries over to a model that does not take it.

### Solver Defaults

The `[eos_solver]` and `[radial_solver]` sections set the defaults for every whole-planet EOS solve and every shooting-method Love-number solve: `LayeredWorld.solve_eos` and `solve_love_numbers`, the standalone `RadialSolver_x.radial_solver`, and the solves that `calc_tides` and the 3D tidal maps run internally. A call overrides only the arguments it passes, so a configuration file plus a world or system TOML fixes every numerical setting of a result.

| Key | `[eos_solver]` | `[radial_solver]` |
|---|---|---|
| `integration_method` | `"DOP853"` | `"DOP853"` |
| `rtol`, `atol` | `1.0e-10`, `1.0e-14` | `1.0e-6`, `1.0e-10` |
| `pressure_tol` | `1.0e-8` (relative to the central-pressure scale) | |
| `max_iters` | `100` | |
| `slices_per_layer` | `100` | |
| `nondimensionalize` | `true` | `true` |
| `use_kamata` | | `false` |
| `start_radius_tolerance` | | `1.0e-5` |
| `scale_rtols` | | `false` |
| `max_num_steps`, `expected_size`, `max_ram_mb` | | `500000`, `128`, `500` |

Both solves run in non-dimensional units (the planet radius, its bulk density, and $1/\sqrt{\pi G \rho}$ as the length, density, and time units), so one tolerance pair means the same thing for every planet. The packaged values come from a convergence study over the bundled worlds and synthetic homogeneous, rocky, icy-ocean, and liquid-core models at degrees 2 and 3 and periods from a day to a hundred days. DOP853 gave the most accuracy per millisecond at every tolerance on both solves: RK45 needs a hundred times tighter `rtol` for the same Love-number error, RK23 far more, and the implicit methods are slower without being more accurate here.

Tightening the EOS tolerance costs almost nothing, so it is set where the mass, moment of inertia, and surface gravity are converged to about 1e-8. The Love tolerance is the loosest pair at which the degree-2 and degree-3 Love numbers of every well-conditioned case stay within about 3e-8 (real part) and 1e-6 (imaginary part) of a reference solved a million times tighter. Each step tighter in `rtol` gains roughly a factor of ten for about a quarter more time, and a Love solve takes a fraction of a millisecond on a cached world. A dynamic liquid layer at a long forcing period is ill-conditioned at any tolerance (use a static liquid there), and an interpolated PREM-style profile is limited by its own tabulation, whose Love numbers move in the fifth digit with the slice count, rather than by the integrator.

### Numerical Settings

`[numerical]` holds the floors and tolerances the C++ code reads through its shared configuration singleton: the frequency extremes (`minimum_frequency`, `maximum_frequency`, `min_spin_orbit_diff`), the material floors (`minimum_viscosity`, `minimum_modulus`), `minimum_solid_rigidity` (the rigidity $\mu / (\bar{\rho} g R)$ below which a melt-weakened solid is solved by the radial solver as a static liquid, applied when the EOS is solved), the geometry floor `minimum_layer_thickness`, the guarded-denominator `numerical_floor`, `layer_continuity_rtol` (also how far past a layer end a profile read snaps onto the end), `max_start_radius_fraction`, `minimum_surface_rcond` (the reciprocal condition number below which a radial solve's surface boundary-condition system counts as singular and the solve fails), `frequency_match_rtol` (how close two tidal-mode frequencies must be to share one radial solve, and how small a frequency counts as zero), `minimum_nusselt` (the floor of the convection cooling model), `eos_invert_rtol` and `eos_invert_max_iters` (the density-from-pressure inversion of the Birch-Murnaghan and Vinet material models, used by any model built without its own `invert_rtol` or `invert_max_iters`), and the quadrature resolutions of the 3D tidal heating integrals (`tides_3d_latitude_nodes`, `tides_3d_longitude_nodes`, `tides_3d_radial_slices`, the defaults of the matching `calc_3d_tides` arguments). `TidalPy.constants.update_constants_x()` pushes an edited value into the C++ side without a restart.

### Graphics

`[graphics.interior]` restyles `plot_interior` (see [graphics](../utilities_x/graphics_x.md)): the matplotlib colors of each profile, the line styles and markers of real and imaginary parts, the marker size, the panel size in inches, and the label and title font sizes. The table is read when `TidalPy.Utilities_x.graphics_x` is imported.

### Layer Material Defaults

Every section of `[layers]` is keyed by a material `type`. A layer that names no `type` takes the `[layers.default]` block, a copy of `[layers.mantle_rock]`, and the same block is the fallback for factories that do not know a layer type. A layer written by `get_config_dict` or `save_to_toml` carries `type = "none"`, which applies no material defaults: the saved layer already lists every model it holds.

The `[layers.default]` model tables are also what a physics-model factory (`make_rheology`, `make_viscosity`, `make_partial_melt`, `make_cooling`, `make_radiogenics`, `make_material_eos`, and `make_tide` from `[tides]`) takes when it is called with no config at all, so a model built by hand and the same model attached by the world builder read one file. Pass an empty dict to get the model's own defaults instead. The table names one model, and only that model takes its parameters: with a `[layers.default.material.shear_viscosity]` table that names the `reference` law, `make_viscosity("reference")` reads it and `make_viscosity("constant")` is built from its own defaults, so a parameter of one model never carries over to another that reads the same key with a different meaning. The names are matched through the family's own aliases, so `make_radiogenics("isotopes")` finds a table that names `isotope`. Inside a world the layer's table is merged over its material block key by key whatever models the two name; each model is then built from its own keys, so a `fixed` radiogenics table over the `isotope` defaults takes no dataset.

### Warnings

`[warnings]` switches the Python warnings the configuration and world-building code can give, each on by default and given at most once per cause per session: `stale_worldpack_copy`, for a data-directory copy of a bundled world or data file that differs from the packaged one (see the [world pack page](../structures_x/config/worldpack.md)); `schema_version`, for a world or system file whose `schema_version` is missing or differs from this build's in its minor version (a major difference is refused, not warned about); `truncation_promotion`, for a `[tides]` truncation level that is not tabulated and is promoted to the next tabulated one; `short_degree_list`, for a `[tides]` per-degree list (`fixed_k`, `fixed_q`, `fixed_dt_s`) the tide model reads that stops short of `max_degree_l`, whose missing degrees are zero and so dissipate nothing; and `unknown_config_key`, for a key of `TidalPy_Configs_x.toml` (or of a configuration passed to `TidalPy.reinit`) that nothing reads, which is how a misspelled or outdated key shows itself. A `[layers.<type>]` block of your own is allowed; its keys are checked against the layer schema, and what its model tables hold is checked when a layer is built from it.

### Reproducing a Run

A configuration file together with a world or system TOML reproduces a result on another computer running the same TidalPy version. A world file may also pin the `[eos_solver]` and `[radial_solver]` keys its results depend on (see the [TOML schema](../structures_x/config/toml_schema.md)), in which case the file alone fixes the solver settings and the configuration supplies the rest. Every saved world, system, and configuration file starts with a comment header naming the TidalPy, SciPy, and CyRK versions that wrote it. Save the configuration in effect with `TidalPy.save_config_x`, and load it with `TidalPy.reinit`:

```python
import TidalPy

# Override settings for this session (a file path works too); only the given values change.
TidalPy.reinit(provided_config_x={"numerical": {"minimum_viscosity": 1.0e3}})

# Save the full effective configuration next to the world file.
TidalPy.save_config_x("my_run_config.toml")

# On another computer (or later): start from the saved settings.
TidalPy.reinit(provided_config_x="my_run_config.toml")

# Return to the packaged defaults merged with your TidalPy_Configs_x.toml.
TidalPy.reinit(provided_config_x="default")
```

A saved configuration begins with a comment header recording the TidalPy, SciPy, and CyRK versions that produced it. The header is a note for the reader: TidalPy does not check it when loading the file. Physical constants such as Newton's constant come from the installed SciPy, so a different SciPy release can shift results slightly even with the same configuration.

## Cleaning the Configurations Directory
If you often install different versions of TidalPy, old config files will build up in the configurations directory.

To clear them, delete the "TidalPy" directory (and all subdirectories) mentioned at the top of this page. TidalPy builds a new directory with the latest config file the next time it is imported.

> [!WARNING]
> If you made changes to the config file, back it up before deleting.

## Important Settings
Many of these settings are explained elsewhere in the documentation. A few are worth changing early on.

### Logging
`write_log_to_disk = false`

Controls whether TidalPy's log is written to a file. Off by default. Set it to `true`, or call `import TidalPy; TidalPy.log_to_file()` to write the log to disk for the current session.

`file_level = "DEBUG"`

The level of messages saved to a file, if file saving is enabled by the previous setting. The default of "DEBUG" includes every message, which is useful when learning TidalPy or debugging an issue.

`console_level = "INFO"`

The level of logging printed to the console. The default of "INFO" covers messages you likely want to see, including warnings and errors.

`print_log_notebook = false`

Controls whether the log is printed to the console when TidalPy is used in a Jupyter notebook. Off by default, because notebook cell output makes it noisy.

### Config
`save_configs_locally = false`

If true, TidalPy saves a copy of its current configurations to a file in your working directory, which records the exact state TidalPy was in for a session.

`use_cwd_for_config = false`

If true, TidalPy uses a file in the working directory as its configuration. The file must be in the working directory and named "TidalPy_Configs.toml".
