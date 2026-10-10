# TidalPy Configurations

_Updated: 2026-10-09_

TidalPy reads its settings from one file, `TidalPy_Configs.toml`, when the package is first imported. The file sits in the TidalPy data directory inside your documents directory:

| System | Path |
|---|---|
| Windows | "C:\\Users\\\<username\>\\Documents\\TidalPy\\\<major.minor\>.X\\Config\\TidalPy_Configs.toml" |
| MacOS | "/Users/\<username\>/Documents/TidalPy/\<major.minor\>.X/Config/TidalPy_Configs.toml" |
| Linux | "/home/\<username\>/Documents/TidalPy/\<major.minor\>.X/Config/TidalPy_Configs.toml" |

The version folder holds only the major and minor version (`0.8.X` for every 0.8 release), so patch releases share one configuration. Beside `Config` it holds `Logs`, `Worlds` (editable copies of the bundled worlds, see the [world pack page](../Structures/config/worldpack.md)), and `Materials` (editable copies of the bundled materials, see the [MatPack page](../Material/matpack.md)). `TidalPy.paths.get_config_dir()`, `get_log_dir()`, `get_worlds_dir()`, and `get_materials_dir()` return the four paths, and `get_data_dir()` the version folder.

To keep the data directory elsewhere, set the `TIDALPY_DATA_DIR` environment variable before importing TidalPy. It replaces the "TidalPy" folder: `TIDALPY_DATA_DIR=/scratch/me/tidalpy` gives "/scratch/me/tidalpy/\<major.minor\>.X/Config/TidalPy_Configs.toml".

> [!NOTE]
> Where the data directory cannot be created or written (a read-only home directory on a cluster node, a container, a CI runner), TidalPy still imports. It warns once and uses the packaged defaults, writes no log file to `Logs`, and reads the bundled worlds and materials from the package. Set `TIDALPY_DATA_DIR` to a writable directory to keep a configuration file.

The file comments each setting. The loaded settings are the dictionary `TidalPy.config`. Edits to the file take effect at the next import, or at once with `TidalPy.reinit(provided_config="default")`.

## Overriding Settings for a Session

To change settings for one project without touching your file, pass another configuration file or a dictionary to `TidalPy.reinit`. Only the values given change:

```python
import TidalPy

print(TidalPy.config["radial_solver"]["rtol"])  # 3e-08 unless your file changes it

# "TidalPy_Config2.toml" stands for your own file; this one is written here so that the example runs.
with open("TidalPy_Config2.toml", "w") as config_file:
    config_file.write("[radial_solver]\nrtol = 1.0e-7\n")
TidalPy.reinit(provided_config="TidalPy_Config2.toml")
print(TidalPy.config["radial_solver"]["rtol"])  # 1e-07

TidalPy.reinit(provided_config={"radial_solver": {"rtol": 1.0e-8}})  # The same with a dictionary
print(TidalPy.config["radial_solver"]["rtol"])  # 1e-08

TidalPy.reinit(provided_config="default")  # Packaged defaults merged with your TidalPy_Configs.toml
```

`TidalPy.reinit` merges the override over `TidalPy.config`, reconfigures the logger, and pushes the numerical and solver settings into the C++ code. Overrides are lost on a kernel restart. For permanent changes, edit your configuration file (make a copy first).

## Reproducing a Run

A configuration file plus a world or system TOML reproduces a result on another computer running the same TidalPy version. A world file may also pin the `[eos_solver]` and `[radial_solver]` keys its results depend on (see the [TOML schema](../Structures/config/toml_schema.md)); it then fixes the solver settings alone and the configuration supplies the rest.

```python
import TidalPy

TidalPy.save_config("my_run_config.toml")             # The full effective configuration
TidalPy.reinit(provided_config="my_run_config.toml")  # Later, or on another computer
```

Every saved world, system, and configuration file begins with a comment header naming the TidalPy, SciPy, and CyRK versions that wrote it. TidalPy does not check it when loading. Physical constants such as Newton's constant come from the installed SciPy, so another SciPy release can shift results slightly.

## Cleaning the Configurations Directory

Each minor version has its own data directory, so old ones build up. Delete the "TidalPy" directory above to clear them all. `TidalPy.clear_data()` deletes the `Config`, `Logs`, `Worlds`, and `Materials` folders of the installed version after asking for confirmation. The next import builds a fresh directory.

> [!WARNING]
> Back up any changes to the config file, the world files in `Worlds`, or the material files in `Materials` before deleting.

## Packaged Defaults

The packaged defaults are the string in ["defaultc.py"](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/defaultc.py) (`TidalPy.configurations.get_packaged_config()` returns them as a dict); developers who want to change or add a configuration edit that string. On first import TidalPy writes them to `TidalPy_Configs.toml` under the version comment header.

TidalPy loads the packaged defaults and merges your file over them, so your file needs only the values you change, and a default added in a later release reaches an existing file. Tables merge key by key; any other value (a list included) replaces the default whole. A key that nothing reads is warned about (see [Warnings](#warnings)), and the per-material `[layers.<type>]` tables of older files are dropped with one warning (see [Layer Defaults](#layer-defaults)).

## Package Setup

`[pathing]`, `[logging]`, and `[configs]` are read when TidalPy is imported or reinitialized.

### Pathing

- `save_directory = "TidalPy-Run"`: the run output directory, created under the working directory only when something is saved to it (a log file with `[logging] use_cwd`, or a configuration copy with `[configs] save_configs_locally`).
- `append_datetime = true`: append the date and time to its name.

### Logging

TidalPy logs through one C++ logger (see [Logging](../Utilities/logging.md)). The levels are `"trace"`, `"debug"`, `"info"`, `"warning"`, `"error"`, `"critical"`, and `"off"`.

- `write_log_to_disk = false`: write the log to `TidalPy_<YYYYMMDD-HHMMSS>.log`. `TidalPy.log_to_file()` does this for the current session. Test mode (`TIDALPY_TEST_MODE`) never writes one.
- `use_cwd = true`: put the log file in the run output directory's `Logs` folder rather than the data directory's.
- `file_level = "debug"`: the file level; `"debug"` keeps nearly every message, useful when learning or debugging.
- `console_level = "info"`: the console level, which includes warnings and errors.
- `print_log_notebook = false`: in a Jupyter notebook, print every message at `console_level` below the cell.
- `notebook_console_level = "warning"`: otherwise a notebook prints only this level and above (the accuracy and convergence warnings), never below `console_level`.
- `write_log_notebook = false`: write a log file from a notebook, where cell output makes the full log noisy.

`TidalPy.capture_log()` collects a block's messages as a list (see [Logging](../Utilities/logging.md)).

### Configs

- `save_configs_locally = false`: save a copy of the loaded configuration to the run output directory at startup.
- `use_cwd_for_config = false`: merge a "TidalPy_Configs.toml" in the working directory over the loaded configuration at startup, if one is there.

## Solver Defaults

`[eos_solver]` and `[radial_solver]` set the defaults of every whole-planet EOS solve and shooting-method Love solve: `solve_eos`, `solve_love_numbers`, the standalone `radial_solver`, and the solves inside `calc_tides` and the 3D tidal maps. A call overrides only the arguments it passes.

| Key | `[eos_solver]` | `[radial_solver]` |
|---|---|---|
| `integration_method` | `"DOP853"` | `"DOP853"` |
| `rtol`, `atol` | `1.0e-10`, `1.0e-14` | `3.0e-8`, `3.0e-8` |
| `nondimensionalize` | `true` | `true` |

`integration_method` takes `"DOP853"`, `"RK45"`, `"RK23"`, or the implicit `"BDF"`, `"LSODA"`, `"Radau"`. `nondimensionalize` integrates in non-dimensional units. The other `[eos_solver]` keys:

- `pressure_tol = 1.0e-8`: the surface-pressure mismatch, relative to the central-pressure scale, that ends the central-pressure iteration; `max_iters = 100` caps the iterations.
- `solve_temperature = true`: carry temperature and heat flow through the solve, so each layer's profile follows its cooling model and its material is evaluated at the local temperature (see [Worlds](../Structures/worlds/worlds.md#temperature-and-heat-flow)). The solve relaxes its thermal network in passes, one structure integration each, until interface temperatures and heat flows change by less than `thermal_tol = 1.0e-8` (relative) between passes, or `max_thermal_passes = 12` is reached.
- `slices_per_layer = 100`: samples per layer in the reported profile only. Love solves and profile getters read the dense output at the exact radius, and the integration finds the edges of the solid and liquid zones itself.

The other `[radial_solver]` keys:

- `starting_method = "takeuchi"`: `"takeuchi"`, `"kamata"`, `"power_series"`, or `"unity"` (see [Starting Conditions](../RadialSolver/starting_conditions.md)).
- `degree1_frame = "CE"`: the frame of degree-1 load Love numbers, `"CE"`, `"CM"`, `"CF"`, `"CL"`, or `"CH"` (Blewitt 2003). Only a degree-1 loading solve reads it; see [Degree-1 Load Love Numbers](../RadialSolver/calculating_love_numbers.md#degree-1-load-love-numbers).
- `start_radius_tolerance = 1.0e-5`: sets the automatic starting radius $R\, \tau^{1/l}$, where $\tau$ is this tolerance, capped by `[numerical] max_start_radius_fraction`.
- `scale_rtols = false`: tighten the `rtol` of the stress-like radial functions by layer type (experimental).
- `max_num_steps = 500000`, `expected_size = 128`, `max_ram_mb = 500`: the step cap, the initial storage in steps (it grows as needed), and the memory cap \[MB\].

### Accuracy of the Defaults

Both solves run in non-dimensional units (the planet radius, its bulk density, and $1/\sqrt{\pi G \rho}$ as the length, density, and time units), so one tolerance pair means the same for every planet. The defaults come from a convergence study over the bundled worlds and synthetic models:

- DOP853 gave the most accuracy per millisecond on both solves. RK45 needs a hundred times tighter `rtol` for the same Love-number error; the implicit methods are slower and no more accurate.
- The EOS tolerances converge mass, moment of inertia, and surface gravity to about 1e-8, at almost no cost.
- At `rtol = atol = 3e-8`, Love numbers are within about 4e-6 of an `rtol = 1e-11` reference (90 percent within 4e-7). A Love solve takes a fraction of a millisecond on a cached world.
- To tighten, lower `rtol` and `atol` together; a small `atol` alone costs steps without improving the Love numbers. Below about 1e-7, neighboring tolerances can differ in error by a factor of a few.

Limits: a dynamic liquid layer at a long forcing period is ill-conditioned at any tolerance (use a static liquid), and an interpolated PREM-style profile is limited by its tabulation (its Love numbers move in the fifth digit with the slice count), not by the integrator.

## Numerical Settings

`[numerical]` holds the floors and tolerances the C++ code reads. `TidalPy.constants.update_constants()` pushes an edited value to C++ without a restart (see [Constants](../Utilities/constants.md)).

- `minimum_frequency = 1.0e-16` \[rad s$^{-1}$\] is the lowest continuation frequency. A tidal mode of frequency $\omega$ takes the Love number at $\omega' = \sqrt{\omega^2 + \omega_c^2}$, with $\omega_c$ its world's continuation frequency, and its $-\mathrm{Im}(k)$ is scaled by $|\omega| / \omega'$. Well below $\omega_c$ its dissipation falls linearly to zero, well above it the mode is solved at its own frequency, and the tidal torque is smooth in the spin through a spin-orbit lock and across $\omega_c$. $\omega_c$ is `minimum_frequency`, or higher where most of the world's solid would be near-fluid, since a Love solve there is slow and returns solver noise. That higher value is the frequency at which $\omega \eta$ falls to the liquid threshold (`minimum_solid_rigidity` times $\bar{\rho} g R$), with $\eta$ the viscosity below which the softer half of the solid volume lies. `BaseWorld.calc_continuation_frequency()` returns $\omega_c$ from the last EOS solve, and a world with an analytic tide model (`fixed_q`, `fixed_dt`, ...) always uses `minimum_frequency`. The continuation is a numerical regularization, exact for a Maxwell body below its slowest relaxation and approximate otherwise (Andrade's transient creep keeps rising as $\omega^{-\alpha}$). Only a mode of exactly zero frequency is dropped, and a direct Love solve below `minimum_frequency` raises `ValueError`.
- `maximum_frequency = 1.0e8` is the largest forcing frequency a Love solve accepts, to catch a frequency not given in rad s$^{-1}$.
- `minimum_modulus = 1.0e-3`: a modulus below this is treated as zero.
- `minimum_solid_rigidity = 1.0e-6`: the post-melt rigidity $\mu / (\bar{\rho} g R)$, with the world's stated bulk density, surface gravity, and radius, at or below which a layer that can change state is liquid. The EOS solve splits the layer into solid and liquid zones there, the radial solver takes each liquid zone as a liquid, and a convecting interior that is liquid there takes the magma-ocean scaling (see [Worlds](../Structures/worlds/worlds.md#pieces-and-zones)).
- `minimum_complex_rigidity = 1.0e-9`: the floor of a solid zone's complex rigidity $|\mu(\omega)| / (\bar{\rho} g R)$ at the forcing frequency in a world's Love solve, reached by raising the real (elastic) part and keeping the imaginary (dissipative) part. A viscously relaxed solid forced near zero frequency has $\mu(\omega) \approx i \omega \eta$, and the solid equations, which divide by it, fail or slow down without the floor, while the zone's dissipation still vanishes with the frequency. 0 turns it off. The standalone `radial_solver` and `solve_love_numbers_supplied` take their moduli as given.
- `minimum_zone_fraction = 1.0e-7`: the thinnest solid or liquid zone, as a fraction of the world radius, kept as its own layer; a thinner one takes its thicker neighbor's state.
- `minimum_layer_thickness = 0.1`: the geometry floor; a thinner layer is ignored.
- `numerical_floor = 1.0e-100`: the smallest magnitude of a guarded denominator; in practice it only replaces a true zero.
- `layer_continuity_rtol = 1.0e-6`: how closely a layer's inner radius must match the previous outer radius, and how far past a layer end a profile read snaps onto the end.
- `max_start_radius_fraction = 0.90`: the largest radial-solver starting radius over the planet radius; the automatic choice is capped here and a larger supplied one refused.
- `minimum_surface_rcond = 1.0e-12`: the reciprocal condition number below which the surface boundary-condition system counts as singular and the radial solve fails.
- `frequency_match_rtol = 1.0e-9`: how close two mode frequencies must be to share a radial solve, and how small a frequency counts as zero.
- `minimum_nusselt = 1.0`: the floor of the convection cooling model's Nusselt number.
- `maximum_eos_mass_ratio = 10.0`: the largest factor between a world's solved and stated mass before its EOS solve fails.
- `eos_invert_rtol = 1.0e-13`, `eos_invert_max_iters = 60`: the density-from-pressure inversion of the Birch-Murnaghan and Vinet laws, for a law built without its own `invert_rtol` or `invert_max_iters`.
- `tides_3d_latitude_nodes = 16`, `tides_3d_longitude_nodes = 64`, `tides_3d_radial_slices = 16`: the quadrature resolutions of the 3D tidal heating integrals, the defaults of the matching `calc_3d_tides` arguments.
- `tides_3d_min_radii_per_thread = 8`: the fewest radii per thread in the colatitude integral of `calc_3d_tides` and of the per-layer heating.
- `love_solve_threads = 0`: the threads of the Love solves and per-layer heating in `calc_tides`; 0 uses `max(num_logical_processors - 4, 1)` (see [Parallel Love Solves](../RadialSolver/parallel.md)).
- `love_solve_min_parallel = 3`: the fewest unique Love solves `calc_tides` spreads over threads; fewer run on the calling thread.
- `test_constant = 42.0`: read only by the test suite.

## Tides

`[tides]` supplies the global (1D) tidal defaults a world takes when its own `[tides]` table omits a value (see the [TOML schema](../Structures/config/toml_schema.md)):

- `min_degree_l = 2` and `max_degree_l = 2`: the harmonic degrees of the mode sum (2 to 10 are supported).
- `eccentricity_trunc_lvl = 10`: one of 2, 4, 6, 8, 10, 20, 50, or `"exact"`; `eccentricity_exact_tolerance = 1.0e-4` sets the mode range of `"exact"` (see [Eccentricity Functions](../Tides/Eccentricity.md)).
- `obliquity_trunc_lvl = "off"`: one of `"off"` (0), 2, 4, or `"gen"` (see [Obliquity Functions](../Tides/Obliquity.md)).
- `love_method = "radial_solver"`: where the rheology tide model takes its Love numbers, the world's interior (`"radial_solver"`) or a quasi-homogeneous method (`"homogeneous"`, `"cpl"`, `"ctl"`).
- `layer_tidal_heating = true`: whether `calc_tides` also resolves each layer's heating when the Love numbers come from the radial solver.
- `fixed_k`, `fixed_q`, `fixed_dt_s`: per-degree Love numbers, quality factors, and time lags \[s\] of the analytic tide models, lists indexed from $l = 2$.

Three keys have no packaged default but are read when a world's `[tides]` table sets them: `global_tidal_model` (the tide model of every world type, in place of `[tides.default_model]`), and `love_fixed_q` and `love_fixed_dt_s` (the quality factor and time lag \[s\] of the `cpl` and `ctl` Love methods).

`[tides.default_model]` names the default tide model per world family: `fixed_q` for stars, `fixed_dt` for gas giants, and `rheology` for terrestrial and layered worlds. A `[tides.<world type>]` table replaces the per-degree lists for that family; the packaged `[tides.star]` holds the fluid Love numbers of an $n = 3$ polytrope with a modified quality factor $Q' = 10^6$.

## Worlds

`[worlds]` holds the default world properties: `albedo = 0.3`, `emissivity = 1.0`, `obliquity_rad = 0.0`, `spin_frequency_rad_s = 0.0`, and `moment_of_inertia_factor = 0.4` (the spin model's $C/(M R^2)$, used until the EOS is solved). `[worlds.star]` adds `effective_temperature_k = 5772.0` and `luminosity_w = 0.0` (zero derives the luminosity from the effective temperature) and sets `albedo = 0.0` and `moment_of_inertia_factor = 0.0754` (an $n = 3$ polytrope). A world's own table wins over these, and they win over the class default. A world built by hand (`TerrestrialWorld(name, radius, mass)` and the other classes) takes the same defaults.

## Evolution

`[evolution]` holds the defaults of `System.evolve` (see [Evolving a World About Its Host](../Structures/system/system.md#evolving-a-world-about-its-host)). An argument to `evolve` wins, then the system file's own `[evolution]` table, then this section.

- `evolve_thermal = true`: evolve the layer temperatures of each world with layers. A thermal run takes its surface temperature from the system's star.
- `method = "Radau"`: CyRK's implicit integrator, `"Radau"`, `"BDF"`, or `"LSODA"` (any case). Radau holds a spin-orbit lock with the longest steps.
- `semi_major_axis_rtol = 1.0e-5`: relative tolerance on the change $a/a_0 - 1$ (its absolute tolerance is $10^{-12}$).
- `eccentricity_rtol = 1.0e-4`, `eccentricity_atol = 1.0e-8`: tolerances on the change in $e$ since its reference value.
- `spin_rtol = 1.0e-3`: relative tolerance on each spin's offset from its commensurability. The absolute tolerance is a tenth of the offset at which the slowest resonant mode reaches the world's continuation frequency.
- `thermal_rtol = 1.0e-5`: relative tolerance on the layer temperatures.
- `radial_rtol = 1.0e-10`, `radial_atol = 1.0e-10`: the Love-solve tolerances during a run. The rates must be smooth at the scale of the integrator's difference steps, which the `[radial_solver]` defaults are not near a lock. A tighter `radial_atol` makes a solve of a near-fluid solid at low frequency far slower.
- `max_wall_time = inf`: a wall-clock cap \[s\]; a run that reaches it returns with `success` False.

TidalPy checks that the method is implicit, the tolerances are finite and positive, and the wall time is at least zero. A bad value in your configuration file falls back to its default with a warning, and `TidalPy.reinit` refuses one. A system file's `[evolution]` table is checked when the file is read, and an unknown key there is refused.

## Radiogenic Isotope Datasets

`[radiogenics]` `isotopes = "modern_day_chondritic"` is the dataset an isotope radiogenics model takes when its table (or `make_radiogenics("isotope")`) names neither a dataset nor its own isotope arrays.

`[radiogenics.known_isotope_data]` is empty by default. Each table in it is a user-defined dataset that a model selects by name, like the built-in `modern_day_chondritic`, `llri_and_slri`, and `bulk_silicate_earth`. Half lives and the reference time are in Myr and the heat production rate `hpr` in \[W kg$^{-1}$\]:

```toml
[radiogenics.known_isotope_data.my_dataset]
    ref_time = 4600.0
    [radiogenics.known_isotope_data.my_dataset.U238]
        iso_mass_fraction = 0.9928
        hpr = 9.48e-5
        half_life = 4470.0
        element_concentration = 0.012e-6
```

A layer uses it with `[layers.<name>.radiogenics]` `model = "isotope"` and `isotopes = "my_dataset"`. A built-in dataset wins over a user dataset of the same name. See [Radiogenic Models](../Radiogenics/radiogenics_models.md).

## Graphics

`[graphics.interior]` styles `plot_interior` (see [graphics](../Utilities/graphics.md)): the matplotlib color of each profile, the line styles and markers of real and imaginary parts, the marker size, the panel size in inches, and the label and title font sizes. It is read when `TidalPy.Utilities.graphics` is imported.

## Layer Defaults

`[layers]` holds one key, `material = "simple_rock"`: the material of a layer whose world-file table names none, as a MatPack name (`TidalPy.Material.available_materials()`) or a material table. Everything else a layer leaves out takes the `Layer` default: no thermal expansion or melting, no heating, the material's own rheology, and no cooling or radiogenics model. A material's values live in its MatPack file in the data directory's `Materials` folder (`TidalPy.Material.material_info(name)["path"]`), or in the layer's own material table (see the [TOML schema](../Structures/config/toml_schema.md#the-material-table)).

Any other table under `[layers]` (such as the per-material `[layers.mantle_rock]` of older files) is not read: it is dropped at load with one warning, under the `unknown_config_key` switch. Delete it to silence the warning.

A physics-model factory called with no config takes the model's own defaults (`make_rheology("andrade")`, `make_viscosity("arrhenius")`, `make_cooling("convection")`), except `make_radiogenics("isotope")`, which takes the `[radiogenics] isotopes` dataset, and `make_tide`, which takes the per-degree lists of `[tides]`.

## Warnings

`[warnings]` switches the Python warnings of the configuration and world-building code. Each is `true` by default and given at most once per cause per session:

- `stale_worldpack_copy`: a data-directory copy of a bundled world or data file differs from the packaged one (see the [world pack page](../Structures/config/worldpack.md)).
- `stale_matpack_copy`: the same for a bundled material (see the [MatPack page](../Material/matpack.md)).
- `schema_version`: a world, system, or material file's `schema_version` is missing or differs from this build's in its minor version (a major difference is refused).
- `truncation_promotion`: a `[tides]` truncation level is not tabulated and is promoted to the next tabulated one.
- `short_degree_list`: a per-degree list (`fixed_k`, `fixed_q`, `fixed_dt_s`) the tide model reads stops short of `max_degree_l`, so the missing degrees dissipate nothing.
- `unknown_config_key`: a key of `TidalPy_Configs.toml` (or of a configuration passed to `TidalPy.reinit`) that nothing reads, such as a misspelled or outdated key.
