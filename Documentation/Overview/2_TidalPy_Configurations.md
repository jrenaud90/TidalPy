# TidalPy Configurations

_Updated: 2026-10-06_

TidalPy's settings and parameters are read when the package is first imported. They live in one configuration file, `TidalPy_Configs.toml`, in the TidalPy data directory inside the user's documents directory, whose location varies by operating system.

**Windows**

"C:\\Users\\\<username\>\\Documents\\TidalPy\\\<major.minor\>.X\\Config\\TidalPy_Configs.toml"

**MacOS**

"/Users/\<username\>/Documents/TidalPy/\<major.minor\>.X/Config/TidalPy_Configs.toml"

**Linux**

"/home/\<username\>/Documents/TidalPy/\<major.minor\>.X/Config/TidalPy_Configs.toml"

The version folder holds only the major and minor version (`0.8.X` for every 0.8 release), so patch releases share one configuration. Beside `Config` it holds `Logs` (log files, when they are written there), `Worlds` (the editable copies of the bundled worlds, see the [world pack page](../Structures/config/worldpack.md)), and `Materials` (the editable copies of the bundled materials, see the [MatPack page](../Material/matpack.md)). `TidalPy.paths.get_config_dir()`, `get_log_dir()`, `get_worlds_dir()`, and `get_materials_dir()` return the four paths, and `get_data_dir()` the version folder.

To keep the data directory somewhere else, set the `TIDALPY_DATA_DIR` environment variable before importing TidalPy. It replaces the "TidalPy" folder in the paths above, so the version folder goes inside it: `TIDALPY_DATA_DIR=/scratch/me/tidalpy` gives "/scratch/me/tidalpy/\<major.minor\>.X/Config/TidalPy_Configs.toml".

> [!NOTE]
> Where the data directory cannot be created or written (a read-only home directory on a cluster node, in a container, or on a CI runner), TidalPy still imports. It warns once and runs without the directory. Instead, the configuration is the packaged defaults, no log file is written to `Logs`, and the bundled worlds and materials are read from the package. Set `TIDALPY_DATA_DIR` to a writable directory to keep a configuration file there.

The file holds all of TidalPy's settings, with comments giving context for each one. The loaded settings are the dictionary `TidalPy.config`. Edits to the file reach TidalPy the next time it is imported, or immediately with `TidalPy.reinit(provided_config="default")`.

## Packaged Defaults

The packaged defaults are the string in ["defaultc.py"](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/defaultc.py) (`TidalPy.configurations.get_packaged_config()` returns them as a dict). The first time TidalPy is imported it writes them, headed by a comment naming the TidalPy, SciPy, and CyRK versions, to `TidalPy_Configs.toml`. Developers who want to change or add a configuration should edit the string in `defaultc.py`.

TidalPy loads the packaged defaults first and merges your file over them, so your file only needs the values you change, and a default added in a later release reaches an existing file without regenerating it. Tables merge key by key and any other value (a list included) replaces the default whole. A key that nothing reads is warned about (see [Warnings](#warnings)), and the per-material `[layers.<type>]` tables of older files are dropped with one warning (see [Layer Defaults](#layer-defaults)).

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
| `[radiogenics]` | The isotope dataset an isotope radiogenics model takes by default, and user-defined isotope datasets. |
| `[graphics]` | Styling of the plotting helpers. |
| `[layers]` | The material of a layer that names none. |

## Overriding Settings for a Session

To change settings for one project while leaving your configuration file alone, pass a new configuration file or a dictionary to `TidalPy.reinit`:

```python
import TidalPy

print(TidalPy.config["radial_solver"]["rtol"])  # 1e-06 unless your file changes it

# A configuration file holding any number of settings. Only those present override the loaded configuration.
# "TidalPy_Config2.toml" stands for your own file; this one is written here so that the example runs.
with open("TidalPy_Config2.toml", "w") as config_file:
    config_file.write("[radial_solver]\nrtol = 1.0e-7\n")
TidalPy.reinit(provided_config="TidalPy_Config2.toml")  # Override from the file
print(TidalPy.config["radial_solver"]["rtol"])  # 1e-07

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

`print_log_notebook = false`, `notebook_console_level = "warning"`, and `write_log_notebook = false`

Control what is printed below a cell, and whether a log file is written, when TidalPy is used in a Jupyter notebook. By default a notebook prints only the messages at `notebook_console_level` and above (the accuracy and convergence warnings), and never below `console_level`; `print_log_notebook = true` prints every message at `console_level`. No log file is written from a notebook unless `write_log_notebook` is true, because notebook cell output makes the full log noisy. `TidalPy.capture_log()` collects a block's messages as a list (see [Logging](../Utilities/logging.md)).

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
| `rtol`, `atol` | `1.0e-10`, `1.0e-14` | `3.0e-8`, `3.0e-8` |
| `pressure_tol` | `1.0e-8` (relative to the central-pressure scale) | |
| `max_iters` | `100` | |
| `solve_temperature` | `true` | |
| `max_thermal_passes`, `thermal_tol` | `12`, `1.0e-8` | |
| `slices_per_layer` | `100` | |
| `nondimensionalize` | `true` | `true` |
| `use_kamata` | | `false` |
| `start_radius_tolerance` | | `1.0e-5` |
| `scale_rtols` | | `false` |
| `max_num_steps`, `expected_size`, `max_ram_mb` | | `500000`, `128`, `500` |

`solve_temperature` carries temperature and heat flow through the structure solve, so each layer's profile follows its cooling model, and its material is evaluated at the local temperature (see [Worlds](../Structures/worlds/worlds.md#temperature-and-heat-flow)). Such a solve relaxes its thermal network against the structure in passes, one structure integration each, until the interface temperatures and heat flows change by less than `thermal_tol` (relative) between passes, or `max_thermal_passes` is reached. `slices_per_layer` only sets the number of radial samples in the profile a solve reports. The Love solves and the profile getters read the solve's dense output at the exact radius, and the integration finds the edges of the solid and liquid zones itself.

Both solves run in non-dimensional units (the planet radius, its bulk density, and $1/\sqrt{\pi G \rho}$ as the length, density, and time units), so one tolerance pair means the same thing for every planet. The packaged values come from a convergence study over the bundled worlds and synthetic homogeneous, rocky, icy-ocean, and liquid-core models at degrees 2 and 3 and periods from a day to a hundred days. DOP853 gave the most accuracy per millisecond at every tolerance on both solves: RK45 needs a hundred times tighter `rtol` for the same Love-number error, RK23 far more, and the implicit methods are slower without being more accurate here.

Tightening the EOS tolerance costs almost nothing, so it is set where the mass, moment of inertia, and surface gravity are converged to about 1e-8. The Love tolerance comes from a sweep of 54 `rtol` and `atol` pairs over the bundled worlds (degrees 2 to 4, tidal and loading, periods within a factor of 30 of each world's own) and the synthetic and profile cases, checked on cases the sweep had not seen. At `rtol = atol = 3e-8` the Love numbers are within about 4e-6 of a reference solved at `rtol = 1e-11`, and 90 percent of them within 4e-7, for the same number of integration steps as the former `1e-6` and `1e-10` with errors two to eight times smaller. A small `atol` costs steps on solution components near zero without improving the Love numbers, so `atol` equal to `rtol` is the cheaper way to tighten. Below about 1e-7, the error also depends on step placement, so neighboring tolerances can differ by a factor of a few. A Love solve takes a fraction of a millisecond on a cached world. A dynamic liquid layer at a long forcing period is ill-conditioned at any tolerance (use a static liquid there), and an interpolated PREM-style profile is limited by its own tabulation, whose Love numbers move in the fifth digit with the slice count, rather than by the integrator.

## Numerical Settings

`[numerical]` holds the floors and tolerances the C++ code reads through its shared configuration singleton: the frequency extremes (`minimum_frequency`, below which a tidal mode counts as static, and `maximum_frequency`), the material floor `minimum_modulus`, `minimum_solid_rigidity` (the post-melt rigidity $\mu / (\bar{\rho} g R)$, with the world's stated bulk density, surface gravity, and radius, at or below which a layer that can change state is liquid: the EOS solve splits the layer into solid and liquid zones there, the radial solver takes each liquid zone as a liquid, and a convecting interior that is liquid there takes the magma-ocean scaling; see [Worlds](../Structures/worlds/worlds.md#pieces-and-zones)), `minimum_complex_rigidity` (the floor of a solid zone's complex rigidity $|\mu(\omega)| / (\bar{\rho} g R)$ at the forcing frequency in a world's Love solve, reached by raising the real (elastic) part and keeping the imaginary (dissipative) part: a viscously relaxed solid forced near `minimum_frequency` has $\mu(\omega) \approx i \omega \eta$, and the solid equations, which divide by it, fail or slow down without the floor, while the zone's dissipation still vanishes with the frequency; 0 turns it off; the standalone `radial_solver` and `solve_love_numbers_supplied` take their moduli as given), `minimum_zone_fraction` (the thinnest solid or liquid zone, as a fraction of the world radius, that the radial solver takes as a layer of its own; a thinner zone takes the state of its thicker neighbor), the geometry floor `minimum_layer_thickness`, the guarded-denominator `numerical_floor`, `layer_continuity_rtol` (also how far past a layer end a profile read snaps onto the end), `max_start_radius_fraction`, `minimum_surface_rcond` (the reciprocal condition number below which a radial solve's surface boundary-condition system counts as singular and the solve fails), `frequency_match_rtol` (how close two tidal-mode frequencies must be to share one radial solve, and how small a frequency counts as zero), `minimum_nusselt` (the floor of the convection cooling model), `maximum_eos_mass_ratio` (the largest factor by which a world's solved mass may differ from its stated mass before its EOS solve fails), `eos_invert_rtol` and `eos_invert_max_iters` (the density-from-pressure inversion of the Birch-Murnaghan and Vinet equation-of-state laws, used by any law built without its own `invert_rtol` or `invert_max_iters`), and the quadrature resolutions of the 3D tidal heating integrals (`tides_3d_latitude_nodes`, `tides_3d_longitude_nodes`, `tides_3d_radial_slices`, the defaults of the matching `calc_3d_tides` arguments), the fewest radii each thread takes in the analytic colatitude integral of `calc_3d_tides` and of the per-layer heating (`tides_3d_min_radii_per_thread`), the threads of the Love-number solves and the per-layer heating in `calc_tides` (`love_solve_threads`, where 0 uses `max(num_logical_processors - 4, 1)`; see [Parallel Love Solves](../RadialSolver/parallel.md)), and the fewest unique Love solves `calc_tides` spreads over threads (`love_solve_min_parallel`; fewer run on the calling thread). `test_constant` is read only by the test suite. `TidalPy.constants.update_constants()` pushes an edited value into the C++ side without a restart (see [Constants](../Utilities/constants.md)).

## Tides

`[tides]` supplies the global (1D) tidal defaults a world takes when its own `[tides]` table omits a value (see the [TOML schema](../Structures/config/toml_schema.md)):

- `min_degree_l = 2` and `max_degree_l = 2`: the harmonic degrees of the mode sum (2 to 10 are supported).
- `eccentricity_trunc_lvl = 10`: the eccentricity truncation level, one of 2, 4, 6, 8, 10, 20, 50, or `"exact"`; `eccentricity_exact_tolerance = 1.0e-4` sets the mode range of `"exact"` (see [Eccentricity Functions](../Tides/Eccentricity.md)).
- `obliquity_trunc_lvl = "off"`: the obliquity truncation level, one of `"off"` (0), 2, 4, or `"gen"` (see [Obliquity Functions](../Tides/Obliquity.md)).
- `love_method = "radial_solver"`: where the rheology tide model takes its Love numbers, the world's interior (`"radial_solver"`) or a quasi-homogeneous method (`"homogeneous"`, `"cpl"`, `"ctl"`).
- `layer_tidal_heating = true`: whether `calc_tides` also resolves each layer's heating when the Love numbers come from the radial solver.
- `fixed_k`, `fixed_q`, `fixed_dt_s`: per-degree Love numbers, quality factors, and time lags \[s\] of the analytic tide models, lists indexed from $l = 2$.

Three keys of a world's `[tides]` table have no packaged default but are read here when set: `global_tidal_model` (the tide model of every world type, in place of `[tides.default_model]`), and `love_fixed_q` and `love_fixed_dt_s` (the quality factor and time lag \[s\] of the `cpl` and `ctl` Love methods).

`[tides.default_model]` names the default tide model per world family: `fixed_q` for stars, `fixed_dt` for gas giants, and `rheology` for terrestrial and layered worlds. A `[tides.<world type>]` table replaces the per-degree lists for that family; the packaged `[tides.star]` holds the fluid Love numbers of an $n = 3$ polytrope with a modified quality factor $Q' = 10^6$.

## Worlds

`[worlds]` holds the default world properties: `albedo = 0.3`, `emissivity = 1.0`, `obliquity_rad = 0.0`, `spin_frequency_rad_s = 0.0`, and `moment_of_inertia_factor = 0.4` (the spin model's $C/(M R^2)$, the world's moment of inertia until its EOS is solved). `[worlds.star]` adds the star-only `effective_temperature_k = 5772.0` and `luminosity_w = 0.0` (zero derives the luminosity from the effective temperature) and sets `albedo = 0.0` and `moment_of_inertia_factor = 0.0754` (an $n = 3$ polytrope). A world's own table wins over these, and they win over the class default. A world built by hand (`TerrestrialWorld(name, radius, mass)` and the other classes) takes the same defaults for every property it is not given.

## Radiogenic Isotope Datasets

`[radiogenics]` `isotopes = "modern_day_chondritic"` is the dataset an isotope radiogenics model takes when its table (or `make_radiogenics("isotope")`) names neither a dataset nor its own isotope arrays.

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

## Layer Defaults

`[layers]` holds one key, `material = "simple_rock"`: the material of a layer whose world-file table names none, as a MatPack name (`TidalPy.Material.available_materials()`) or a material table. Everything else a layer leaves out takes the `Layer` default, the simplest case: no thermal expansion or melting, no heating, the material's own rheology, and no cooling or radiogenics model. A material's values live in its MatPack file in the `Materials` folder of the data directory (`TidalPy.Material.material_info(name)["path"]`), or in the layer's own material table (see the [TOML schema](../Structures/config/toml_schema.md#the-material-table)).

A table under `[layers]` other than `material` in your file (the per-material blocks of older files, such as `[layers.mantle_rock]`) is not read: it is dropped at load with one warning naming it, under the `unknown_config_key` switch. Delete it from the file to silence the warning.

A physics-model factory called with no config takes the model's own defaults (`make_rheology("andrade")`, `make_viscosity("arrhenius")`, `make_cooling("convection")`), with two exceptions: `make_radiogenics("isotope")` takes the `[radiogenics] isotopes` dataset, and the tide factories (`make_tide`) take the per-degree lists of `[tides]`.

## Warnings

`[warnings]` switches the Python warnings the configuration and world-building code can give, each on by default and given at most once per cause per session:

- `stale_worldpack_copy`: a data-directory copy of a bundled world or data file differs from the packaged one (see the [world pack page](../Structures/config/worldpack.md)).
- `stale_matpack_copy`: the same for a bundled material (see the [MatPack page](../Material/matpack.md)).
- `schema_version`: a world, system, or material file's `schema_version` is missing or differs from this build's in its minor version (a major difference is refused, not warned about).
- `truncation_promotion`: a `[tides]` truncation level is not tabulated and is promoted to the next tabulated one.
- `short_degree_list`: a `[tides]` per-degree list (`fixed_k`, `fixed_q`, `fixed_dt_s`) the tide model reads stops short of `max_degree_l`, so its missing degrees are zero and dissipate nothing.
- `unknown_config_key`: a key of `TidalPy_Configs.toml` (or of a configuration passed to `TidalPy.reinit`) that nothing reads, which is how a misspelled or outdated key shows itself.

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

To clear them, delete the "TidalPy" directory (and all subdirectories) mentioned at the top of this page. `TidalPy.clear_data()` deletes the `Config`, `Logs`, `Worlds`, and `Materials` folders of the installed version after asking for confirmation. TidalPy builds a new directory with the latest config file the next time it is imported.

> [!WARNING]
> If you made changes to the config file, to the world files in `Worlds`, or to the material files in `Materials`, back them up before deleting.
