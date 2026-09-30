# Migrating from TidalPy 0.7.X

_Updated: 2026-09-29_

TidalPy 0.8.0 replaced the Python, Cython, and numba code of 0.7.X and earlier with a C++ backend wrapped by Cython. The modules, classes, functions, configuration file, and logging all changed, so 0.7.X scripts need updating. This page maps the 0.7.X API onto 0.8.0 and shows how to port common workflows. The <a href="code_map.html">interactive code map</a> shows the main classes and functions of 0.8.0, the calls between them, and the purpose, inputs, and outputs of each.

The 0.7.X API will receive no new features but will get bug fixes until approx. the end of 2026. We recommend switching to the 0.8.0 API. To keep the 0.7.X API, pin the version:

```bash
pip install "TidalPy<0.8"
```

or, with conda, 

```bash
conda install -c conda-forge "tidalpy<0.8"
```

## What Changed

- The physics runs in C++, wrapped by thin Cython layers. Without numba, nothing compiles on a function's first call.
- Classes store configuration and return results from explicit `solve_*`, `get_*`, and `calc_*` calls. Changing an attribute no longer updates the world, its layers, and its orbit.
- Orbital state no longer lives on a world. A `System` holds the orbits and passes them to each tidal calculation.
- Worlds and systems are described by TOML files that carry a `schema_version`, or by the equivalent Python dict.
- Every physics family (rheology, viscosity, partial melt, cooling, radiogenics, material equations of state, tides, stellar luminosity) follows one pattern: model classes, a `make_<family>(name, config)` factory, direct functions, vectorized `calc_*` methods, and binary save and load.
- One configuration file, `TidalPy_Configs.toml`, in a data directory scoped to the minor version (`<Documents>/TidalPy/0.8.X/`).
- One logger, written in C++ with spdlog.
- The new code raises the built-in `ValueError`, `RuntimeError`, `TypeError`, and `NotImplementedError` instead of TidalPy's own exception classes.
- The required dependencies are NumPy, SciPy, matplotlib, platformdirs, toml, and CyRK. numba, dill, pathos, astropy, and astroquery are no longer used. psutil moved to the `dev` extra. The `burnman` and `julia` extras were removed.
- Importing TidalPy no longer warns about the backend change. `TidalPy.exceptions.TidalPyDeprecationWarning` still exists, so code that filters it keeps working.

## Module Map

Imports are case-sensitive on every operating system, so `import TidalPy.rheology` raises `ModuleNotFoundError` even though `TidalPy.Rheology` exists.

| 0.7.X | 0.8.0 | Documentation |
|---|---|---|
| `TidalPy.structures` (worlds, layers, `Orbit`) | `TidalPy.Structures` (worlds, layers, `System`) | [Structures](Structures/index.md), [System](Structures/system/system.md), [TOML schema](Structures/config/toml_schema.md) |
| `TidalPy.RadialSolver` | `TidalPy.RadialSolver` | [RadialSolver](RadialSolver/index.md) |
| `TidalPy.Material.eos` | `TidalPy.Material.eos` (material EOS models and the whole-planet solver) | [Material and EOS](Material/index.md) |
| `TidalPy.rheology` | `TidalPy.Rheology` | [Rheology](Rheology/index.md) |
| `TidalPy.rheology.viscosity` | `TidalPy.Viscosity` | [Viscosity](Viscosity/index.md) |
| `TidalPy.rheology.partial_melt` | `TidalPy.PartialMelt` | [Partial Melting](PartialMelt/index.md) |
| `TidalPy.cooling` | `TidalPy.Cooling` | [Cooling](Cooling/index.md) |
| `TidalPy.radiogenics` | `TidalPy.Radiogenics` | [Radiogenics](Radiogenics/index.md) |
| `TidalPy.tides` | `TidalPy.Tides` (mostly used through world methods) | [Tides](Tides/index.md) |
| `TidalPy.dynamics`, `TidalPy.orbit` | `TidalPy.Dynamics` and the `System` class | [Dynamics](Dynamics/index.md), [System](Structures/system/system.md) |
| `TidalPy.stellar` | `TidalPy.Stellar` (luminosity models attached to a star) and `System` (insolation) | [Stellar](Stellar/index.md) |
| `TidalPy.utilities` | `TidalPy.Utilities` | [Utilities](Utilities/index.md) |
| `TidalPy.logger` | `TidalPy.Utilities.logging` | [Logging](Utilities/logging.md) |
| `TidalPy.toolbox` | no replacement yet (see [Quick Tidal Dissipation](#quick-tidal-dissipation)) | |
| `TidalPy.Extending` (BurnMan), `TidalPy.output`, `TidalPy.numba_scipy` | removed | |
| `TidalPy.WorldPack` (bundled worlds) | `TidalPy.WorldPack` (bundled worlds in the new TOML schema) | [WorldPack](Structures/config/worldpack.md) |

## Exceptions

`TidalPy.exceptions` keeps `TidalPyException`, `InitializationError`, `ArgumentException`, `ConfigurationException`, `ModelException`, `UnknownModelError`, `TidalPyIntegrationException`, `SolutionFailedError`, and `TidalPyDeprecationWarning`. The other 0.7.X exception classes are removed. An `except TidalPyException` block no longer catches errors from the physics modules:

| 0.7.X | 0.8.0 |
|---|---|
| `ArgumentException` (bad arguments, array sizes, radius ordering) | `ValueError` |
| `UnknownModelError` (an unknown model, layer type, or integration method name) | `ValueError`; an unknown config key also names the closest accepted key |
| `SolutionFailedError` from `radial_solver(raise_on_fail=True)` | `SolutionFailedError` (unchanged), also raised when the solver's EOS solve fails |
| `AttributeNotSetError` (a world state that was never set) | `RuntimeError` (for example, a tidal solve on a world whose EOS is not solved) |
| a wrong argument type | `TypeError` |

## Configuration and Data Directory

TidalPy 0.8.0 keeps its settings in `TidalPy_Configs.toml` in `<Documents>/TidalPy/0.8.X/Config` (see [TidalPy Configurations](Overview/2_TidalPy_Configurations.md)). Each minor version has its own data directory, so 0.8.0 does not read the 0.7.X file. Copy over any setting you changed, using the new keys:

| 0.7.X | 0.8.0 |
|---|---|
| `[pathing]`, `[logging]`, `[configs]` | Unchanged except as listed below. The log levels `"trace"` and `"off"` are new. |
| `[logging] console_error_level` | Removed. |
| `[configs] use_cwd_for_world_dir`, `overwrite_configs` | Removed. The editable worlds are always in the data directory's `Worlds` folder. |
| `[debug]` (`extensive_logging`, `extensive_checks`) and `[numba]` | Removed, with `TidalPy.extensive_logging` and `TidalPy.extensive_checks`. |
| `[tides.modes] minimum_frequency`, `maximum_frequency` | `[numerical] minimum_frequency`, `maximum_frequency`. |
| `[tides.modes] min_spin_orbital_diff` | Removed, with `TidalPy.constants.MIN_SPIN_ORBITAL_DIFF`; `minimum_frequency` is the only zero-frequency floor. |
| `[tides.models.*] eccentricity_truncation_lvl`, `max_tidal_order_l`, `obliquity_tides_on` | `[tides] eccentricity_trunc_lvl`, `max_degree_l`, `obliquity_trunc_lvl`. |
| `[tides.models.global_approx] fixed_q`, `static_k2`, `fixed_dt`, `use_ctl` | `[tides] fixed_q`, `fixed_k`, `fixed_dt_s` (lists indexed from $l = 2$), and the model named in `[tides.default_model]`. |
| `[layers.ice]`, `[layers.rock]`, `[layers.iron]` | `[layers.ice]`, `[layers.mantle_rock]`, `[layers.iron]`, plus `hp_ice`, `gas`, and `default`, each with one table per physics model. |
| `[worlds.types.*]` | Removed. `[worlds]` now holds default world properties (albedo, emissivity, obliquity, spin). |
| `[physics.radiogenics.known_isotope_data]` | `[radiogenics.known_isotope_data]`. |
| `[graphics.planet_plots]` | `[graphics.interior]`. |
| none | `[numerical]`, `[eos_solver]`, `[radial_solver]`, and `[warnings]` are new. |

The configuration functions also changed:

- `TidalPy.reinit(provided_config_file)` is now `TidalPy.reinit(provided_config)`, which takes a file path, a dict, or `"default"` and merges it over the loaded configuration.
- `TidalPy.save_config(path)` is new. It saves the configuration in effect, which together with a world or system file reproduces a result.
- `TidalPy.world_config_dir` is replaced by `TidalPy.paths.get_worlds_dir()`.
- TidalPy warns about any key it does not read, which flags keys carried over from 0.7.X.

```python
import TidalPy

TidalPy.reinit(provided_config={"tides": {"eccentricity_trunc_lvl": 20}})   # Override for this session
TidalPy.save_config("my_run_config.toml")                                  # Record the settings of a run
TidalPy.reinit(provided_config="default")                                  # Back to your saved file
```

## Logging

0.7.X logged through Python's `logging` module (`TidalPy.logger.get_logger`). 0.8.0 logs through one C++ logger (spdlog). The C++ code writes to it directly, often with the interpreter lock released. Python code writes to it through `TidalPy.Utilities.logging`:

```python
from TidalPy.Utilities.logging import log_info, set_console_level

set_console_level("debug")                            # Show debug messages in the console for this session
log_info("Starting the Io run")                       # Written to the same console and file as TidalPy's messages
```

Handlers attached to Python's `logging` (including pytest's `caplog`) do not see these messages. The log file keeps its name, `TidalPy_<YYYYMMDD-HHMMSS>.log`, and is written only when `[logging] write_log_to_disk` is set.

## Worlds

The bundled worlds are new TOML files, and the set changed:

- `earth` is replaced by `earth_prem` (built from the PREM profile) beside `earth_simple`.
- `jupiter_simple` is new.
- `io_simple`, `triton_simple`, `55cnc`, `55cnce`, `55cnce_simple`, and `nereid_dev` are removed.

`TidalPy.Structures.available_worlds()` lists the bundled worlds. The bundled files are copied to `<Documents>/TidalPy/0.8.X/Worlds`, where they can be edited. `TidalPy.Structures.install_worldpack(force=True)` restores the packaged copies and discards those edits.

A world computes its interior in `solve_eos` and its Love numbers in `solve_love_numbers`. Its tides come from a `System` or from a direct `calc_tides` call that takes the orbital state as arguments.

| 0.7.X | 0.8.0 |
|---|---|
| `build_world(world_name, world_config)` | `build_world(source)`, where `source` is a bundled name, a TOML path, or a dict |
| `build_from_world(old_world, new_config)`, `scale_from_world(...)` | Edit `world.get_config_dict()` and pass it to `build_world` |
| `world.paint()` | `TidalPy.Utilities.graphics.plot_interior` with the world's profile getters |
| `world.save_world()` | `world.save_to_toml(path)`, or `world.save_binary(path)` for the binary format |
| `world.orbit`, `world.set_state(...)`, and the orbital setters (`semi_major_axis`, `eccentricity`, `orbital_frequency`) | The `System` that holds the world (see [Systems](#systems)) |
| `world.set_spin_frequency`, `world.set_obliquity` | Unchanged |
| `world.set_fixed_q`, `world.set_fixed_dt` | `world.set_tide_model(make_tide("fixed_q" or "fixed_dt", config))` |
| `world.tidal_heating_global` | `world.get_tidal_heating()` after `calc_tides`, or `System.calc_world_evolution` |

```python
import numpy as np

from TidalPy.Structures import build_world
from TidalPy.Utilities.graphics import plot_interior

world = build_world("earth_simple")                   # A bundled name, a TOML path, or a dict
world.solve_eos()                                     # Interior structure: density, gravity, pressure, mass
world.solve_love_numbers(
    frequency=2.0 * np.pi / (12.42 * 3600.0))         # Semidiurnal forcing [rad s-1]
print(world.love_number_k, world.love_number_h)

# Interior plot, in place of world.paint()
radius = np.linspace(0.0, world.radius, 200)          # [m]
figure, axes = plot_interior(
    radius,
    world.get_gravity(radius),
    world.get_pressure(radius),
    world.get_density(radius),
    planet_radius=world.radius,
    planet_name="Earth",
    show_plot=False)

# Save and rebuild, in place of world.save_world()
world.save_to_toml("my_earth.toml")
reloaded = build_world("my_earth.toml")

# A modified copy, in place of build_from_world()
config = world.get_config_dict()
config["name"] = "Earth-Cold-Mantle"
config["layers"]["mantle"]["temperature_k"] = 1200.0
modified = build_world(config)
```

## Systems

`System` replaces `Orbit` (`PhysicsOrbit`). The system returns the rates for the current state and holds no integrator: a time evolution integrates these rates with an integrator of your choice (demo 12 uses CyRK).

| 0.7.X | 0.8.0 |
|---|---|
| `Orbit(star, tidal_host, tidal_bodies)` | `System(name)` then one `add_world` per world |
| `add_tidal_world`, `add_tidal_host`, `add_star` | `add_world(world, tidal_host=..., is_star=..., semi_major_axis=..., eccentricity=...)` and `set_tidal_host` |
| `orbit.set_state(world, eccentricity=..., semi_major_axis=...)` | `set_eccentricity(world, e)`, `set_semi_major_axis(world, a)` |
| `calculate_orbital_derivatives(world)`, `TidalPy.dynamics` single-body functions | `calc_world_evolution(world)`: tidal heating, `dU_dM`, `dU_dw`, `dU_dO`, `da_dt`, `de_dt`, `dn_dt`, `dspin_dt`, and energy-balance terms |
| the `TidalPy.dynamics` dual-body functions | `calc_pair_evolution(world)`: both bodies raise tides on their shared orbit |
| `calculate_insolation(world)` | `calc_insolation_flux(world)` and `calc_equilibrium_temperature(world)` |

```python
from TidalPy.Structures import build_system, build_world
from TidalPy.Structures.system import System

io = build_world("io")
io.solve_eos()   # Needed by the rheology tide model
jupiter = build_world("jupiter_simple")

system = System("jovian")
system.add_world(jupiter)   # Jupiter raises Io's tides and is treated as a point mass
system.add_world(
    io,
    tidal_host=jupiter,
    semi_major_axis=4.217e8,
    eccentricity=0.0041)
io.set_spin_frequency(system.calc_orbital_frequency(io))   # Synchronous rotation

rates = system.calc_world_evolution(io)
print(rates["tidal_heating"])   # [W]
print(rates["da_dt"], rates["de_dt"], rates["dspin_dt"])

# Systems can also be built from TOML; this one is bundled
sol_system = build_system("sol_system")
print(sol_system.calc_insolation_flux("earth"))   # [W m-2], about 1361
```

## Quick Tidal Dissipation

`TidalPy.toolbox` (`quick_tidal_dissipation` and `quick_dual_body_tidal_dissipation`) has no replacement yet. `calc_world_evolution` returns the same quantities for any world in a `System`: the tidal heating, the potential derivatives, and the orbit and spin rates. The world's `love_method` sets how its Love numbers are found:

- `"radial_solver"` (the shooting method, and the default) integrates the radial equations through every layer.
- `"homogeneous"` treats each tidal layer as a homogeneous sphere of its averaged material, with no radial solve, as `quick_tidal_dissipation` did.

The bundled Io has a core, a mantle, and a thin asthenosphere that does almost all of the dissipating. The example calculates its heating and rates with both methods:

```python
from TidalPy.Structures import build_world
from TidalPy.Structures.system import System

io = build_world("io")                                # Core, mantle, and a dissipating asthenosphere
io.solve_eos()                                        # Interior structure, used by both Love methods
jupiter = build_world("jupiter_simple")               # The host acts as a point mass

system = System("jovian")
system.add_world(jupiter)
system.add_world(
    io,
    tidal_host=jupiter,
    semi_major_axis=4.217e8,
    eccentricity=0.0041)
io.set_spin_frequency(system.calc_orbital_frequency(io))   # Synchronous rotation

for love_method in ("radial_solver", "homogeneous"):
    io.set_tide_config(love_method=love_method)
    rates = system.calc_world_evolution(io)
    print(love_method, rates["tidal_heating"])       # [W]: 9.33e13 (radial solver), 1.37e14 (homogeneous)
    print(rates["dU_dM"], rates["dU_dw"], rates["dU_dO"])
    print(rates["da_dt"], rates["de_dt"], rates["dspin_dt"])   # [m s-1], [s-1], [rad s-2]
```

The radial solver resolves the layering and reproduces the 9.33e13 W that the bundled Io is calibrated to (Lainey et al. 2009), in about 2 ms per call. The homogeneous method weights the layers by volume, which only approximates how a thin weak layer deforms inside a stiffer planet. Here it overestimates the heating and every rate by 47 percent but takes under 0.1 ms. Use it for fast sweeps and first estimates, and the radial solver when the interior structure matters. See [Love numbers](Tides/love/love_numbers.md) for the other methods.

A one-layer world built from bare numbers reproduces `quick_tidal_dissipation`. The world below is its default homogeneous Maxwell body at its default degree 2 and eccentricity truncation 2. The result matches the homogeneous-sphere closed form.

```python
import numpy as np

from TidalPy.Structures import build_world
from TidalPy.Structures.system import System
from TidalPy.Tides.classes import make_tide

radius = 1.8215e6                                     # [m]
mass = 8.93e22                                        # [kg]
density = mass / (4.0 / 3.0 * np.pi * radius**3)      # [kg m-3]

target = build_world({
    "schema_version": "0.2.0",
    "name": "target",
    "type": "terrestrial",
    "radius_m": radius,
    "mass_kg": mass,
    "tides": {
        "global_tidal_model": "rheology",   # Dissipation from the layer's complex modulus
        "love_method": "homogeneous",       # Homogeneous-sphere Love numbers, no radial solve
        "max_degree_l": 2,
        "eccentricity_trunc_lvl": 2,
        "obliquity_trunc_lvl": "off"},
    "layers": {
        "interior": {
            "class": "solidliquid",
            "radius_fraction": 1.0,
            "temperature_k": 1600.0,
            "material": {
                "model": "constant",
                "reference_density_kg_m3": density,
                "shear_modulus_static_pa": 6.0e10,
                "shear_viscosity": {
                    "model": "constant",
                    "reference_viscosity_pas": 1.0e19},
                "partial_melt": {
                    "model": "off"}},
            "shear_rheology": {
                "model": "maxwell"}}}})
target.solve_eos()                     # Surface gravity and moment of inertia

host = build_world("jupiter_simple")   # Any world can be the host; it acts as a point mass
system = System("pair")
system.add_world(host)
system.add_world(
    target,
    tidal_host=host,
    semi_major_axis=4.217e8,
    eccentricity=0.0041)
target.set_spin_frequency(system.calc_orbital_frequency(target))   # Synchronous rotation

result = system.calc_world_evolution(target)
print(result["tidal_heating"])                        # [W], about 2.6e10
print(result["dU_dM"], result["dU_dw"], result["dU_dO"])
print(result["da_dt"], result["de_dt"], result["dspin_dt"])

# A constant phase lag (fixed-Q) body, like rheology="fixed_q" in quick_tidal_dissipation
target.set_tide_model(make_tide(
    "fixed_q",
    {"fixed_k": [0.3],          # k_2
     "fixed_q": [100.0]}))      # Q_2
result = system.calc_world_evolution(target)
print(result["tidal_heating"])  # [W], equals (21/2)(k_2/Q_2) G M^2 R^5 n e^2 / a^6
```

For dual-body dissipation, as in `quick_dual_body_tidal_dissipation`, make each body the other's tidal host with `system.set_tidal_host(host, target)` and call `system.calc_pair_evolution(target)`. Each body dissipates through its own tide model. A gas giant built from a file carries a `fixed_dt` model by default. The per-degree Love numbers that `quick_tidal_dissipation` returned come from `target.get_tidal_love_k(l, m, p, q)` after a `calc_tides` call, or from the closed-form functions in `TidalPy.Tides.love`.

## Radial Solver

`radial_solver` takes the same positional arguments and nearly the same keywords as in 0.7.X. The differences:

- `use_prop_matrix=True` is now `love_method="propagation_matrix"`.
- The solver settings (`integration_method`, `integration_rtol`, `integration_atol`, `expected_size`, the `eos_*` arguments, and the rest) default to `None`, which reads the `[radial_solver]` and `[eos_solver]` sections of the configuration. The packaged tolerances are tighter than the 0.7.X defaults (for example, `integration_rtol` 1e-6 and `integration_atol` 1e-10, against 1e-5 and 1e-8).
- `solve_for` is a tuple of case-insensitive strings (_e.g._, `("tidal", "loading")`).
- Invalid inputs raise `ValueError` instead of `ArgumentException` or `UnknownModelError`.
- When several boundary conditions are solved for, `k`, `h`, and `l` are complex128 arrays; after a failed solve they are complex128 NaN arrays (float64 in 0.7.X).
- `moi_factor` is now the conventional $C/(M R^2)$, 0.4 for a uniform sphere. The 0.7.X value, $C/(0.4 M R^2)$, is now the new `moi_sphere_ratio`. Code that reads `moi_factor` gets a value 2.5 times smaller.
- The solution is evaluated at any radius through dense output (`get_radial_solution(radius)`), and `plot_ys` and `plot_interior` take `show_plot` and plotting keywords.
- The input builders `build_rs_input_homogeneous_layers` and `build_rs_input_from_data` keep their argument names. Their rheology arguments take `TidalPy.Rheology` models or model names, and one model can stand in for every layer. `perform_checks` is accepted and ignored: inputs are always validated.

```python
import numpy as np

from TidalPy.RadialSolver import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.Rheology import Elastic, Maxwell

# Use a helper function to build the radial solver inputs (see "build_inputs.md" for details on these helpers)
build_data = build_rs_input_homogeneous_layers(
    6000.0e3,                                         # Planet radius [m]
    2.0 * np.pi / (86400.0 * 7.5),                    # Forcing frequency [rad s-1]
    density_tuple=(8000.0, 3300.0),
    static_bulk_modulus_tuple=(2.0e11, 1.0e11),
    static_shear_modulus_tuple=(1.0e11, 5.0e10),
    bulk_viscosity_tuple=(1.0e18, 1.0e18),
    shear_viscosity_tuple=(1.0e20, 1.0e19),
    layer_type_tuple=("solid", "solid"),
    layer_is_static_tuple=(False, False),
    layer_is_incompressible_tuple=(False, False),
    shear_rheology_model_tuple=Maxwell(),             # One model for every layer
    bulk_rheology_model_tuple=Elastic(),
    radius_fraction_tuple=(0.5, 1.0),
    slices_tuple=(20, 30))

solution = radial_solver(
    *build_data,
    degree_l=2,
    solve_for=("tidal", "loading"))                     # Both boundary conditions in one solve
print(solution.k)                                       # One complex k per boundary condition
print(solution.moi_factor, solution.moi_sphere_ratio)   # C/(M R^2) and the 0.7.X moi_factor
```

See [Calculating Love Numbers](RadialSolver/calculating_love_numbers.md) and [Helper Functions](RadialSolver/build_inputs.md).

## Rheology and Other Physics Models

The rheology, viscosity, partial-melt, cooling, radiogenics, and luminosity functions of 0.7.X are model classes in 0.8.0, built by name with a factory. Config keys carry their units (`reference_viscosity_pas`, `solidus_k`).

| 0.7.X | 0.8.0 |
|---|---|
| `find_rheology(name)` | `make_rheology(name, config)`, which returns a model instance |
| `rheology(frequency, modulus, viscosity)` | `rheology.calc_complex_modulus(modulus, viscosity, frequency)` (note the order) |
| `vectorize_frequency`, `vectorize_modulus_viscosity` | `calc_complex_modulus_vectorize_frequency`, `_vectorize_modulus`, `_vectorize_all` |
| `Voigt(args=(5.0, 0.02))`, `Andrade(args=(0.3, 1.0))`, `change_args` | Keyword parameters: `Voigt(voigt_modulus_frac=5.0, voigt_viscosity_frac=0.02)`, `Andrade(alpha=0.3, zeta=1.0)` |
| `Newton`, `SundbergCooper`; the names `voigtkelvin` and `sundbergcooper` | `Viscous`, `Sundberg`; the factory takes `newton`, `voigt-kelvin`, and `sundberg-cooper` as aliases |
| the complex compliance functions (`rheology.complex_compliance`) | removed: the models return the complex modulus, whose reciprocal is the compliance |
| `rheology.viscosity` functions (`arrhenius`, `reference`, `constant`) | `make_viscosity(name, config)` and `calc_viscosity(temperature, pressure)` |
| `rheology.partial_melt` (`spohn`, `henning`, `calculate_melt_fraction`) | `make_partial_melt(name, config)`, `calc_melt_fraction`, `calc_partial_melt` |
| `cooling` functions (`convection`, `conduction`, `off`) | `make_cooling(name, config)` and `calc_cooling`, or the direct functions `convective`, `conductive`, `cooling_off` |
| `radiogenics` functions (`isotope`, `fixed`, `off`) with times in Myr | `make_radiogenics(name, config)` and `calc_heating(time, mass)` with times in seconds; isotope sets are named datasets |
| `stellar.luminosity_from_mass` | `TidalPy.Stellar.mass_to_luminosity(mass)` or `make_luminosity("mass_to_luminosity")` |

```python
import numpy as np

from TidalPy.Radiogenics import make_radiogenics
from TidalPy.Rheology import make_rheology
from TidalPy.Viscosity import make_viscosity

rheology = make_rheology(
    "andrade",
    {"alpha": 0.25,
     "zeta": 1.0})
# Arguments are modulus [Pa], viscosity [Pa s], frequency [rad s-1]
complex_modulus = rheology.calc_complex_modulus(6.0e10, 1.0e19, 2.0e-5)
complex_moduli = rheology.calc_complex_modulus_vectorize_frequency(
    6.0e10,
    1.0e19,
    np.logspace(-8, -3, 50)   # A frequency sweep [rad s-1]
)

viscosity = make_viscosity(
    "reference",
    {"reference_viscosity_pas": 1.0e21,
     "reference_temperature_k": 1600.0,
     "molar_activation_energy_j_mol": 3.0e5})
print(viscosity.calc_viscosity(1500.0, 1.0e9))  # [Pa s] at 1500 K and 1 GPa

radiogenics = make_radiogenics(
    "isotope",
    {"isotopes": "modern_day_chondritic"})
print(radiogenics.calc_heating(radiogenics.ref_time, 1.0))   # [W] for 1 kg at the reference time
```

A world attaches these models to its layers from its TOML file or dict, so most scripts do not build them by hand. See [Rheology](Rheology/index.md), [Viscosity](Viscosity/index.md), [Partial Melting](PartialMelt/index.md), [Cooling](Cooling/index.md), [Radiogenics](Radiogenics/index.md), and [Stellar](Stellar/index.md).

## Tides

- The eccentricity functions are unsquared: `TidalPy.Tides.eccentricity_func(eccentricity, degree_l, truncation)`, with `eccentricity_squared_func` for the squares. 0.7.X tabulated the squared functions per degree and level (`eccentricity_funcs_l2_trunc10` and so on).
- As in 0.7.X, truncation level $N$ keeps every product of two eccentricity functions, and so every heating term, through $e^N$. Levels 2 to 10 and 20 keep the same terms in both versions.
- 0.8.0 tabulates levels 2, 4, 6, 8, 10, 20, and 50, plus `"exact"`. A configured level of 12 to 18 or 22 is promoted to the next tabulated level with a warning, and the direct functions raise `NotImplementedError` for it. The default level is now 10 (it was 6). `recommend_eccentricity_truncation` picks a level for a given eccentricity.
- Obliquity was on or off in 0.7.X. 0.8.0 offers `"off"`, levels 2 and 4, and the general functions `"gen"`.
- Degrees 2 to 10 are supported (2 to 7 in 0.7.X).
- Tidal modes are keyed by $(l, m, p, q)$ instead of names such as `'2o-n'`: `world.get_tidal_love_k(l, m, p, q)`.
- The grid potential functions (`tidal_potential_nsr`, `tidal_potential_obliquity_nsr`, and the others, with their `use_static` switch) and the multilayer mode collapse are replaced by `BaseWorld.calc_3d_tides` and `calc_3d_displacements` (see [3D Tidal Stress, Strain, and Heating](Tides/multilayer_3d_heating.md)). `calculate_displacements` is replaced by `TidalPy.Tides.displacement_point` and `calc_3d_displacements`.
- The `tides.love1d` helpers are in `TidalPy.Tides.love`: `calc_effective_rigidity(shear_modulus, density, gravity, radius, degree_l=2)` (the argument order changed from `effective_rigidity(shear_modulus, gravity, radius, density)`), `calc_homogeneous_love_numbers(complex_shear_modulus, density, gravity, radius, degree_l=2)` in place of `complex_love`, and `apply_fixed_q` and `apply_fixed_dt`.

```python
from TidalPy.Tides.love import calc_effective_rigidity, calc_homogeneous_love_numbers

effective_rigidity = calc_effective_rigidity(
    6.0e10,                                           # Shear modulus [Pa]
    3500.0,                                           # Density [kg m-3]
    1.8,                                              # Surface gravity [m s-2]
    1.82e6)                                           # Radius [m]
love = calc_homogeneous_love_numbers(
    6.0e10 + 1.0e8j,                                  # Complex shear modulus [Pa]
    3500.0,
    1.8,
    1.82e6)
print(effective_rigidity, love.k)
```

## Utilities and Constants

| 0.7.X | 0.8.0 |
|---|---|
| `TidalPy.constants.PI_DBL`, `NAN_DBL`, `DBL_MANT_DIG` | `TidalPy.constants.pi`, `nan`, `dbl_mant_digits` (C++: `d_PI`, `d_NAN`, `d_DBL_MANT_DIGITS`) |
| `TidalPy.constants.MIN_VISCOSITY`, `MIN_SPIN_ORBITAL_DIFF` | Removed |
| `sec2myr`, `myr2sec` with 3.154e13 s per Myr | The same functions with the exact Julian value, `TidalPy.constants.seconds_per_myr` = 3.15576e13 s (a 0.06% difference) |
| `orbital_motion2semi_a`, `semi_a2orbital_motion` | Unchanged names; invalid inputs raise `ValueError` |
| `build_nondimensional_scales(frequency, mean_radius, bulk_density)` | `build_nondimensional_scales(mean_radius, bulk_density)` |
| `utilities.graphics.multilayer.yplot` | `TidalPy.Utilities.graphics.plot_ys`: `plot_tobie` and `plot_roberts` become `benchmarks=("tobie2005", "roberts_nimmo2008")`, `plot_imags` becomes `plot_imaginary`, `other_xlimits` and `other_ylimits` become `x_limits` and `y_limits` |
| `utilities.graphics.planet_plot` | `TidalPy.Utilities.graphics.plot_interior`, with `show_plot` in place of `auto_show`; styling from `INTERIOR_PLOT_STYLE` (`[graphics.interior]`) |
| `projection_map`, `GridPlot`, `success_grid_plot` | `plot_map` and `make_map_axes` for surface maps; the grid plots are removed |
| `stellar.insolation` functions (`equilibrium_insolation_mendez`, `_williams`, `_no_eccentricity`) and `calc_equilibrium_temperature` | `System.calc_insolation_flux` (the orbit-averaged Méndez and Rivera-Valentín form) and `System.calc_equilibrium_temperature` |

## Removed With No Replacement

- `TidalPy.toolbox`, with `quick_tidal_dissipation` and `quick_dual_body_tidal_dissipation` (see [Quick Tidal Dissipation](#quick-tidal-dissipation) for the equivalent).
- The BurnMan interior builds (`TidalPy.Extending`). BurnMan remains an optional comparison package for `Benchmarks/EOS/EOS_vs_BurnMan.ipynb` only.
- Orbit averaging (`TidalPy.orbit.orbit_average` and its 3D and 4D forms).
- The multiprocessing driver (`TidalPy.utilities.multiprocessing`) and the exoplanet archive download (`TidalPy.utilities.exoplanets.get_exoplanet_data`). The Love solves release the interpreter lock, so standard thread and process pools run them in parallel (see [Parallel Love Solves](RadialSolver/parallel.md)).
- numba support (`TidalPy.numba_scipy` and the `[numba]` configuration).
- `TidalPy.output`.
- The state graph that updated a world, its layers, and its orbit when one attribute changed.
- Mode outputs keyed by mode name, and eccentricity truncation levels 12 to 18 and 22.
- The selectable insolation models; one orbit-averaged form remains.
- `calculate_temperature_frommelt` and its array form, and `calculate_mass_gravity_arrays` (the EOS solve fills mass and gravity).
- `calc_tidal_susceptibility`.

## Performance

The timings below compare TidalPy 0.7.6 and 0.8.0. Ratios change with the machine and the problem size, so treat them as rough magnitudes and measure your own workload before relying on them.

The largest gains are where 0.7.X called out to BurnMan or compiled numba kernels. 0.8.0 is slower on the standalone radial solver at tight tolerances, on two vectorized sweeps, and on scalar rheology calls.

### Where It Is Faster

| Task | 0.7.X | 0.8.0 | Change |
|---|---|---|---|
| Build a planet with its interior (Io, 3 layers) | 178 ms | 0.73 ms | 245x faster |
| Radiogenic heating, one evaluation | 0.29 us | 0.046 us | 6.4x faster |
| Orbit-averaged 3D heating map (50 x 16 x 32) | 13.0 ms | 2.12 ms | 6.1x faster |
| Global tidal heating, degrees 2 to 4, e^10 | 0.050 ms | 0.013 ms | 3.9x faster |
| Rheology, 10k complex moduli | 0.175 ms | 0.055 ms | 3.2x faster |
| Build a world from config (2 layers, no interior solve) | 1.38 ms | 0.44 ms | 3.1x faster |
| Global tidal heating, e^10 truncation | 0.026 ms | 0.0084 ms | 3.1x faster |
| Global tidal heating, e^4 truncation | 0.021 ms | 0.0071 ms | 3.0x faster |
| Global tidal heating, e^2 truncation | 0.020 ms | 0.0069 ms | 3.0x faster |
| Instantaneous 3D heating map (50 x 16 x 32 x 8 times) | 12.9 ms | 5.23 ms | 2.5x faster |
| Homogeneous Love numbers (closed form) | 0.15 us | 0.063 us | 2.4x faster |
| Convective cooling, one evaluation | 0.16 us | 0.14 us | 1.2x faster |

For planet building, 0.7.X handed the interior to BurnMan, which does mineral-physics lookups and its own root finding. 0.8.0 integrates the equation of state in C++.

The 3D maps are timed at eccentricity truncation level 2 in both versions, where both keep the same heating terms.

The global tidal heating rows use the homogeneous Love method, which solves the same problem as the 0.7.X `quick_tidal_dissipation`. They use eccentricity truncations that both versions tabulate and at which both sum the same heating terms (see [Tides](#tides)). In 0.8.0 the cost follows the number of distinct forcing frequencies, not the number of modes. Love numbers are solved once per frequency and degree. The homogeneous methods form the layer-averaged shear modulus once per frequency and share it across degrees, so adding degrees adds little. With the `radial_solver` Love method, each frequency and degree is a full radial solve, and that solve dominates the cost.

### Where It Is Slower

| Task | 0.7.X | 0.8.0 | Change |
|---|---|---|---|
| `radial_solver`, 1 layer, 10 slices | 0.30 ms | 0.45 ms | 0.68x, 1.5x slower |
| Rheology, one complex modulus | 0.060 us | 0.081 us | 0.74x, 1.4x slower |
| Radiogenic heating, 10k times | 0.085 ms | 0.109 ms | 0.78x, 1.3x slower |
| Convective cooling, 10k evaluations | 0.160 ms | 0.205 ms | 0.78x, 1.3x slower |
| `radial_solver`, 3 layers (static liquid core), 300 slices | 0.44 ms | 0.51 ms | 0.87x, 1.15x slower |
| `radial_solver`, 1 layer, 200 slices | 0.45 ms | 0.50 ms | 0.90x, 1.1x slower |
| `radial_solver`, propagation matrix, 200 slices | 0.109 ms | 0.118 ms | 0.92x, 1.1x slower |

The 0.8.0 standalone radial solver wraps the world path. It builds a temporary world from the supplied arrays, solves that world's equation of state, and integrates the Love-number equations against the same dense structure the world path uses.

> [!TIP]
> We recommend the world-based approach. It is as fast as the standalone solver, and reusing a constructed `BaseWorld` is much faster than repeated calls to the standalone radial solver: its equation of state is solved once, its Love solves are cached per degree and frequency, and `calc_tides` reuses them across modes.

The two versions take identical integration steps on these problems, so the whole gap is the cost of each right-hand-side read. 0.8.0 evaluates the equation-of-state interpolant for gravity, the interpolated material for density and both static moduli, and the two supplied complex-modulus arrays. 0.7.X did four linear interpolations of its input arrays. The equation-of-state solve itself takes a tenth of the time, the temporary world a twentieth, and the Python-side handling about as long as the world build. The dense read is more accurate. Against the closed-form homogeneous sphere, the 0.8.0 degree-2 k2 is 2.6 times closer (6.8e-5 against 1.8e-4). The two versions agree to 1e-15 on the one-layer worlds and 6e-10 on the three-layer one at the tolerances timed here (`integration_rtol` 1e-8, `integration_atol` 1e-12, both versions).

These settings reduce the time:

- `integration_rtol` and `integration_atol` set the step count. At the `[radial_solver]` defaults (1e-6 and 1e-10, looser than the rows above), the 0.8.0 solver takes 0.25, 0.30, and 0.35 ms on the three shooting rows, faster than 0.7.X at its tighter setting. k2 moves by 4e-9 to 6e-9.
- The equation-of-state settings (`eos_rtol`, `eos_atol`, `eos_integration_method`) change the total by under 10 percent, and the slice count matters little once the searches are seeded.
- `RK45` for the Love integration is slower than `DOP853`, which reaches the tolerance in fewer steps.

The cost of a single scalar rheology call is the Python-to-C++ boundary, not the arithmetic. The numba path of 0.7.X crosses a cheaper boundary. Use the vectorized calls, where 0.8.0 is about 3x faster, for more than a handful of values. The radiogenic sweep evaluates one exponential per isotope and time in both versions. numba vectorizes those exponentials, and the C++ loop does not. The vectorized convective cooling gap has not been investigated.

### First Call

Steady-state timings leave out the startup cost. 0.7.X compiles its numba kernels on their first run and caches the machine code on disk. The first session after installing or upgrading pays the full compile, and every later session pays to load and dispatch the cached kernels. 0.8.0 has nothing to compile.

| First call | 0.7.X, first session after installing | 0.7.X, later sessions | 0.8.0 |
|---|---|---|---|
| Tidal heating, degrees 2 to 4, e^10 | 6.2 s | 1.0 s | 0.07 ms |
| 3D heating map (50 x 16 x 32 x 8 times) | 4.2 s | 0.88 s | 5.4 ms |
| Build a planet with its interior (Io) | not measured | 1.3 s | 1.8 ms |

### Threads for 3D Grids

The 3D grid methods (`calc_3d_tides`, `calc_3d_stress_strain`, `calc_3d_displacements`, and `get_3d_tidal_heating_array`) take `num_threads`, which spreads the radial solves and the per-point evaluation over threads. The default, 0, uses the number of logical processors minus 4 (at least 1). Every thread count returns identical values. `calc_tides` also spreads its Love-number solves over threads (see [Parallel Love Solves](RadialSolver/parallel.md)). 0.7.X has no equivalent.

The table times three grids of a homogeneous Io at degrees 2 to 3 with eccentricity, a non-synchronous spin, and obliquity, on 20 radii by 45 colatitudes by 90 longitudes. It uses the `tides_3d:*_1_thread` and `tides_3d:*_all_threads` tasks in `Benchmarks/Performance` on one machine with 16 hardware threads. Each figure is the lowest of three fresh processes, each taking the best of three batches.

| Grid | 1 thread | 16 threads | Change |
|---|---|---|---|
| Secular heating map | 149 ms | 22 ms | 6.7x faster |
| Stress and strain, 4 times | 266 ms | 39 ms | 6.8x faster |
| Displacements, 24 times | 302 ms | 53 ms | 5.7x faster |

Part of each call (building the waves and the strain coefficients from the radial solutions) runs on the calling thread, so the gain stays below the thread count.

## Learning TidalPy 0.8.0

Start with the [Getting Started](Overview/1_Getting_Started.md) page and the notebooks in the Demos section of the navigation: the `Basics` notebooks (configuration, world building, save and load), the `Physics` notebooks (orbits, tides, rheology, Love numbers, 3D heating, thermal interiors, truncation levels), and the `Systems` notebooks (multi-world systems, coupled thermal-orbital evolution, the Earth-Moon-Sun system). The Benchmarks section validates TidalPy against published results and tracks its performance.
