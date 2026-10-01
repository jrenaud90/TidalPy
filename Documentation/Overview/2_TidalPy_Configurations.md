# TidalPy Configurations

_Updated: 2026-09-29_

TidalPy's settings and parameters are read when the package is first imported. They live in one configuration file, `TidalPy_Configs.toml`, in the TidalPy data directory inside the user's documents directory, whose location varies by operating system.

**Windows**

"C:\\Users\\\<username\>\\Documents\\TidalPy\\\<major.minor\>.X\\Config\\TidalPy_Configs.toml"

**MacOS**

"/Users/\<username\>/Documents/TidalPy/\<major.minor\>.X/Config/TidalPy_Configs.toml"

**Linux**

"/home/\<username\>/Documents/TidalPy/\<major.minor\>.X/Config/TidalPy_Configs.toml"

The version folder holds only the major and minor version (`0.8.X` for every 0.8 release), so patch releases share one configuration. Beside `Config` it holds `Logs` (log files, when they are written there) and `Worlds` (the editable copies of the bundled worlds, see the [world pack page](../Structures/config/worldpack.md)). `TidalPy.paths.get_config_dir()`, `get_log_dir()`, and `get_worlds_dir()` return the three paths, and `get_data_dir()` the version folder.

To keep the data directory somewhere else, set the `TIDALPY_DATA_DIR` environment variable before importing TidalPy. It replaces the "TidalPy" folder in the paths above, so the version folder goes inside it: `TIDALPY_DATA_DIR=/scratch/me/tidalpy` gives "/scratch/me/tidalpy/\<major.minor\>.X/Config/TidalPy_Configs.toml".

> [!NOTE]
> Where the data directory cannot be created or written (a read-only home directory on a cluster node, in a container, or on a CI runner), TidalPy still imports. It warns once and runs without the directory. Instead, the configuration is the packaged defaults, no log file is written to `Logs`, and the bundled worlds are read from the package. Set `TIDALPY_DATA_DIR` to a writable directory to keep a configuration file there.

The file holds all of TidalPy's settings, with comments giving context for each one. The loaded settings are the dictionary `TidalPy.config`. Edits to the file reach TidalPy the next time it is imported, or immediately with `TidalPy.reinit(provided_config="default")`.

## Packaged Defaults

The packaged defaults are the string in ["defaultc.py"](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/defaultc.py) (`TidalPy.configurations.get_packaged_config()` returns them as a dict). The first time TidalPy is imported it writes them, headed by a comment naming the TidalPy, SciPy, and CyRK versions, to `TidalPy_Configs.toml`. Developers who want to change or add a configuration should edit the string in `defaultc.py`.

TidalPy loads the packaged defaults first and merges your file over them, so your file only needs the values you change, and a default added in a later release reaches an existing file without regenerating it. Tables merge key by key and any other value (a list included) replaces the default whole. A physics-model table merges the same way when it names a different `model` than the default: a key the new model does not read is ignored by it. The exception is a key whose meaning depends on the model reading it (`configurations.MODEL_SPECIFIC_KEYS`, the radiogenics `ref_time_s`), which is dropped from the default table first.

The file has these sections:

| Section | Controls |
|---|---|
| `[pathing]` | The run output directory created in the working directory when something is saved there. |
| `[logging]` | The console and file log levels and where the log file goes. |
| `[configs]` | Whether a configuration in the working directory is used, and whether the loaded configuration is saved beside the run. |
| `[numerical]` | Floors, tolerances, and quadrature resolutions read by the C++ code. |
| `[eos_solver]`, `[radial_solver]` | Defaults of every equation-of-state and Love-number solve. |
| `[tides]` | Default tidal model per world family, harmonic degrees, truncation levels, and per-degree fixed parameters. |
| `[warnings]` | Switches for the Python warnings of the configuration and world-building code. |
| `[worlds]` | Default world properties (albedo, emissivity, obliquity, spin, and the star-only values). |
| `[radiogenics.known_isotope_data]` | User-defined radiogenic isotope datasets. |
| `[graphics]` | Styling of the plotting helpers. |
| `[layers.<type>]` | Per-material layer defaults used by the world builder and the model factories. |

## Overriding Settings for a Session

To change settings for one project while leaving your configuration file alone, pass a new configuration file or a dictionary to `TidalPy.reinit`:

```python
import TidalPy

print(TidalPy.config["radial_solver"]["rtol"])  # 1e-06 unless your file changes it

# A configuration file (e.g., "TidalPy_Config2.toml") holding any number of settings. Only those present
# override the loaded configuration.
TidalPy.reinit(provided_config="TidalPy_Config2.toml")

# The same with a dictionary.
TidalPy.reinit(provided_config={"radial_solver": {"rtol": 1.0e-8}})
print(TidalPy.config["radial_solver"]["rtol"])  # 1e-08

# Return to the packaged defaults merged with your TidalPy_Configs.toml.
TidalPy.reinit(provided_config="default")
```

`TidalPy.reinit` merges the override over `TidalPy.config`, reconfigures the logger, and pushes the numerical and solver settings into the C++ code. Overrides need to be performed again after a kernel restart. For permanent changes, edit the configuration file in your documents directory (make a copy first).

## Package Setup

These three sections are read when TidalPy is imported or reinitialized.

### Pathing

`save_directory = "TidalPy-Run"` and `append_datetime = true`

The run output directory, created under the current working directory only when something is saved to it (a log file with `[logging] use_cwd`, or a copy of the configuration with `[configs] save_configs_locally`). With `append_datetime` the date and time are appended to its name.

### Logging

TidalPy logs through one C++ logger (see [Logging](../Utilities/logging.md)). The levels are `"trace"`, `"debug"`, `"info"`, `"warning"`, `"error"`, `"critical"`, and `"off"`.

`write_log_to_disk = false`

Controls whether TidalPy's log is written to a file, `TidalPy_<YYYYMMDD-HHMMSS>.log`. Off by default. Set it to `true`, or call `TidalPy.log_to_file()` to write the log to disk for the current session. Test mode (the `TIDALPY_TEST_MODE` environment variable) never writes a log file.

`use_cwd = true`

Writes the log file to the `Logs` folder of the run output directory (`[pathing]`) rather than the `Logs` folder of the TidalPy data directory.

`file_level = "debug"`

The level of messages saved to a file, if file saving is enabled. The default of `"debug"` includes nearly every message, which is useful when learning TidalPy or debugging an issue.

`console_level = "info"`

The level of logging printed to the console. The default of `"info"` covers messages you likely want to see, including warnings and errors.

`print_log_notebook = false` and `write_log_notebook = false`

Control whether the log is printed below the cell, and whether a log file is written, when TidalPy is used in a Jupyter notebook. Both are off by default, because notebook cell output makes the console log noisy.

### Configs

`save_configs_locally = false`

If true, TidalPy saves a copy of the loaded configuration to the run output directory at startup, which records the exact settings of a session.

`use_cwd_for_config = false`

If true, TidalPy merges a file named "TidalPy_Configs.toml" found in the current working directory over the loaded configuration at startup. A working directory without that file leaves the loaded configuration as it is.

## Solver Defaults

The `[eos_solver]` and `[radial_solver]` sections set the defaults for every whole-planet EOS solve and every shooting-method Love-number solve: `BaseWorld.solve_eos` and `solve_love_numbers`, the standalone `RadialSolver.radial_solver`, and the solves that `calc_tides` and the 3D tidal maps run internally. A call overrides only the arguments it passes, so a configuration file plus a world or system TOML fixes every numerical setting of a result.

| Key | `[eos_solver]` | `[radial_solver]` |
|---|---|---|
| `integration_method` | `"DOP853"` | `"DOP853"` |
| `rtol`, `atol` | `1.0e-10`, `1.0e-14` | `1.0e-6`, `1.0e-10` |
| `pressure_tol` | `1.0e-8` (relative to the central-pressure scale) | |
| `max_iters` | `100` | |
| `solve_temperature` | `true` | |
| `slices_per_layer` | `100` | |
| `nondimensionalize` | `true` | `true` |
| `use_kamata` | | `false` |
| `start_radius_tolerance` | | `1.0e-5` |
| `scale_rtols` | | `false` |
| `max_num_steps`, `expected_size`, `max_ram_mb` | | `500000`, `128`, `500` |

`solve_temperature` carries temperature and heat flow through the structure solve, so each layer's profile follows its cooling model; its viscosity and melt models see the local temperature (see [Worlds](../Structures/worlds/worlds.md)). `slices_per_layer` only sets the number of radial samples in the profile a solve reports. The Love solves and the profile getters read the solve's dense output at the exact radius.

Both solves run in non-dimensional units (the planet radius, its bulk density, and $1/\sqrt{\pi G \rho}$ as the length, density, and time units), so one tolerance pair means the same thing for every planet. The packaged values come from a convergence study over the bundled worlds and synthetic homogeneous, rocky, icy-ocean, and liquid-core models at degrees 2 and 3 and periods from a day to a hundred days. DOP853 gave the most accuracy per millisecond at every tolerance on both solves: RK45 needs a hundred times tighter `rtol` for the same Love-number error, RK23 far more, and the implicit methods are slower without being more accurate here.

Tightening the EOS tolerance costs almost nothing, so it is set where the mass, moment of inertia, and surface gravity are converged to about 1e-8. The Love tolerance is the loosest pair at which the degree-2 and degree-3 Love numbers of every well-conditioned case stay within about 3e-8 (real part) and 1e-6 (imaginary part) of a reference solved a million times tighter. Each step tighter in `rtol` gains roughly a factor of ten for about a quarter more time, and a Love solve takes a fraction of a millisecond on a cached world. A dynamic liquid layer at a long forcing period is ill-conditioned at any tolerance (use a static liquid there), and an interpolated PREM-style profile is limited by its own tabulation, whose Love numbers move in the fifth digit with the slice count, rather than by the integrator.

## Numerical Settings

`[numerical]` holds the floors and tolerances the C++ code reads through its shared configuration singleton: the frequency extremes (`minimum_frequency`, below which a tidal mode counts as static, and `maximum_frequency`), the material floor `minimum_modulus`, `minimum_solid_rigidity` (the rigidity $\mu / (\bar{\rho} g R)$ below which a melt-weakened solid is solved by the radial solver as a static liquid, applied when the EOS is solved), the geometry floor `minimum_layer_thickness`, the guarded-denominator `numerical_floor`, `layer_continuity_rtol` (also how far past a layer end a profile read snaps onto the end), `max_start_radius_fraction`, `minimum_surface_rcond` (the reciprocal condition number below which a radial solve's surface boundary-condition system counts as singular and the solve fails), `frequency_match_rtol` (how close two tidal-mode frequencies must be to share one radial solve, and how small a frequency counts as zero), `minimum_nusselt` (the floor of the convection cooling model), `maximum_eos_mass_ratio` (the largest factor by which a world's solved mass may differ from its stated mass before its EOS solve fails), `eos_invert_rtol` and `eos_invert_max_iters` (the density-from-pressure inversion of the Birch-Murnaghan and Vinet material models, used by any model built without its own `invert_rtol` or `invert_max_iters`), and the quadrature resolutions of the 3D tidal heating integrals (`tides_3d_latitude_nodes`, `tides_3d_longitude_nodes`, `tides_3d_radial_slices`, the defaults of the matching `calc_3d_tides` arguments), the fewest radii each thread takes in the analytic colatitude integral of `calc_3d_tides` and of the per-layer heating (`tides_3d_min_radii_per_thread`), and the threads of the Love-number solves and the per-layer heating in `calc_tides` (`love_solve_threads`, where 0 uses `max(num_logical_processors - 4, 1)` see [Parallel Love Solves](../RadialSolver/parallel.md)). `test_constant` is read only by the test suite. `TidalPy.constants.update_constants()` pushes an edited value into the C++ side without a restart (see [Constants](../Utilities/constants.md)).

## Tides

`[tides]` supplies the global (1D) tidal defaults a world takes when its own `[tides]` table omits a value (see the [TOML schema](../Structures/config/toml_schema.md)):

- `min_degree_l = 2` and `max_degree_l = 2`: the harmonic degrees of the mode sum (2 to 10 are supported).
- `eccentricity_trunc_lvl = 10`: the eccentricity truncation level, one of 2, 4, 6, 8, 10, 20, 50, or `"exact"`; `eccentricity_exact_tolerance = 1.0e-4` sets the mode range of `"exact"` (see [Eccentricity Functions](../Tides/Eccentricity.md)).
- `obliquity_trunc_lvl = "off"`: the obliquity truncation level, one of `"off"` (0), 2, 4, or `"gen"` (see [Obliquity Functions](../Tides/Obliquity.md)).
- `layer_tidal_heating = true`: whether `calc_tides` also resolves each layer's heating when the Love numbers come from the radial solver.
- `fixed_k`, `fixed_q`, `fixed_dt_s`: per-degree Love numbers, quality factors, and time lags \[s\] of the analytic tide models, lists indexed from $l = 2$.

`[tides.default_model]` names the default tide model per world family: `fixed_q` for stars, `fixed_dt` for gas giants, and `rheology` for terrestrial and layered worlds. A `[tides.<world type>]` table replaces the per-degree lists for that family; the packaged `[tides.star]` holds the fluid Love numbers of an $n = 3$ polytrope with a modified quality factor $Q' = 10^6$.

## Worlds

`[worlds]` holds the default world properties: `albedo = 0.3`, `emissivity = 1.0`, `obliquity_rad = 0.0`, and `spin_frequency_rad_s = 0.0`. `[worlds.star]` adds the star-only `effective_temperature_k = 5772.0` and `luminosity_w = 0.0` (zero derives the luminosity from the effective temperature). A world's own table wins over these, and they win over the class default.

## Radiogenic Isotope Datasets

`[radiogenics.known_isotope_data]` is empty by default. Each table inside it is a user-defined isotope dataset that a radiogenics model selects by name, the same way it selects the built-in `modern_day_chondritic`, `llri_and_slri`, and `bulk_silicate_earth` sets. Half lives and the reference time are in Myr and the heat production rate `hpr` in \[W kg$^{-1}$\]:

```toml
[radiogenics.known_isotope_data.my_dataset]
    ref_time = 4600.0
    [radiogenics.known_isotope_data.my_dataset.U238]
        iso_mass_fraction = 0.9928
        hpr = 9.48e-5
        half_life = 4470.0
        element_concentration = 0.012e-6
```

A layer then uses it with `[layers.<name>.radiogenics]` `model = "isotope"` and `isotopes = "my_dataset"`. A built-in dataset takes precedence over a user dataset of the same name. See [Radiogenic Models](../Radiogenics/radiogenics_models.md).

## Graphics

`[graphics.interior]` restyles `plot_interior` (see [graphics](../Utilities/graphics.md)): the matplotlib colors of each profile, the line styles and markers of real and imaginary parts, the marker size, the panel size in inches, and the label and title font sizes. The table is read when `TidalPy.Utilities.graphics` is imported.

## Layer Material Defaults

Every section of `[layers]` is keyed by a material `type`: the packaged types are `iron`, `mantle_rock`, `ice`, `hp_ice`, and `gas`. A layer that names no `type` takes the `[layers.default]` block, a copy of `[layers.mantle_rock]`, and the same block is the fallback for factories that do not know a layer type. A layer written by `get_config_dict` or `save_to_toml` carries `type = "none"`, which applies no material defaults: the saved layer already lists every model it holds.

The `[layers.default]` model tables are also what a physics-model factory (`make_rheology`, `make_viscosity`, `make_partial_melt`, `make_cooling`, `make_radiogenics`, `make_material_eos`, and `make_tide` from `[tides]`) takes when it is called with no config at all, so a model built by hand and the same model attached by the world builder read one file. Pass an empty dict to get the model's own defaults instead. A factory follows the builder's model-switch rule: a model the table does not name still takes its parameters, except the few keys another model of the family would read with a different meaning (a radiogenic dataset's `ref_time_s`, which a fixed rate would take as its own reference time). With a `[layers.default.material.shear_viscosity]` table that names the `reference` law, `make_viscosity("constant")` takes its reference viscosity too, and a key a model does not read is ignored by it. The names are matched through the family's own aliases, so `make_radiogenics("isotopes")` counts as the table's `isotope` model. Inside a world the layer's table is merged over its material block key by key whatever models the two name; each model is then built from its own keys, so a `fixed` radiogenics table over the `isotope` defaults takes no dataset.

## Warnings

`[warnings]` switches the Python warnings the configuration and world-building code can give, each on by default and given at most once per cause per session:

- `stale_worldpack_copy`: a data-directory copy of a bundled world or data file differs from the packaged one (see the [world pack page](../Structures/config/worldpack.md)).
- `schema_version`: a world or system file's `schema_version` is missing or differs from this build's in its minor version (a major difference is refused, not warned about).
- `truncation_promotion`: a `[tides]` truncation level is not tabulated and is promoted to the next tabulated one.
- `short_degree_list`: a `[tides]` per-degree list (`fixed_k`, `fixed_q`, `fixed_dt_s`) the tide model reads stops short of `max_degree_l`, so its missing degrees are zero and dissipate nothing.
- `unknown_config_key`: a key of `TidalPy_Configs.toml` (or of a configuration passed to `TidalPy.reinit`) that nothing reads, which is how a misspelled or outdated key shows itself.

A `[layers.<type>]` block of your own is allowed; its keys are checked against the layer schema, and what its model tables hold is checked when a layer is built from it.

## Reproducing a Run

A configuration file together with a world or system TOML reproduces a result on another computer running the same TidalPy version. A world file may also pin the `[eos_solver]` and `[radial_solver]` keys its results depend on (see the [TOML schema](../Structures/config/toml_schema.md)), in which case the file alone fixes the solver settings and the configuration supplies the rest. Save the configuration in effect with `TidalPy.save_config`, and load it with `TidalPy.reinit`:

```python
import TidalPy

# Override settings for this session (a file path works too); only the given values change.
TidalPy.reinit(provided_config={"radial_solver": {"rtol": 1.0e-8}})

# Save the full effective configuration next to the world file.
TidalPy.save_config("my_run_config.toml")

# On another computer (or later): start from the saved settings.
TidalPy.reinit(provided_config="my_run_config.toml")

# Return to the packaged defaults merged with your TidalPy_Configs.toml.
TidalPy.reinit(provided_config="default")
```

Every saved world, system, and configuration file begins with a comment header recording the TidalPy, SciPy, and CyRK versions that produced it. The header is a note for the reader: TidalPy does not check it when loading the file. Physical constants such as Newton's constant come from the installed SciPy, so a different SciPy release can shift results slightly even with the same configuration.

## Cleaning the Configurations Directory

Each minor version of TidalPy has its own data directory, so old ones build up if you often install different versions.

To clear them, delete the "TidalPy" directory (and all subdirectories) mentioned at the top of this page. `TidalPy.clear_data()` deletes the `Config`, `Logs`, and `Worlds` folders of the installed version after asking for confirmation. TidalPy builds a new directory with the latest config file the next time it is imported.

> [!WARNING]
> If you made changes to the config file or to the world files in `Worlds`, back them up before deleting.
