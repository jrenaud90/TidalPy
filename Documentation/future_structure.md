# Migrating from TidalPy 0.7.X

_Updated: 2026-10-08_

TidalPy 0.8.0 replaced the Python, Cython, and numba code of 0.7.X with a C++ backend wrapped by Cython. The modules, classes, functions, configuration file, and logging all changed, so 0.7.X scripts need updating. This page maps the 0.7.X API onto 0.8.0 and shows how to port common workflows. The <a href="code_map.html">interactive code map</a> shows the main classes and functions of 0.8.0 and the calls between them.

The 0.7.X API gets only bug fixes, until approx. the end of 2026. To keep it, pin the version:

```bash
pip install "TidalPy<0.8"
```

or, with conda,

```bash
conda install -c conda-forge "tidalpy<0.8"
```

## What Changed

- The physics runs in C++ behind thin Cython layers (no numba jit compiles).
- Classes store configuration and return results from explicit `solve_*`, `get_*`, and `calc_*` calls. Changing an attribute no longer updates the world, its layers, and its orbit.
- Orbital state no longer lives on a world. A `System` holds the orbits and passes them to each tidal calculation.
- Worlds and systems are described by TOML files that carry a `schema_version`, or by the equivalent Python dict.
- Every physics family follows one pattern: model classes, a `make_<family>(name, config)` factory, vectorized `calc_*` methods, a config dict, binary save and load, and `parameters`, `get_parameter`, `get_parameter_info`, and `with_parameters`. Most also have direct helper functions.
- A layer's interior is a `Material`: a solid and a liquid `Phase` (each an equation of state, a shear-modulus law, viscosity laws, default rheologies, and thermal constants) with melting curves, melt weakening, and latent heat. MatPack ships many named materials, from simplified rock and ice to peridotite, the ices, iron, and giant-planet envelopes.
- One configuration file, `TidalPy_Configs.toml`, in a data directory scoped to the minor version (`<Documents>/TidalPy/0.8.X/`).
- One logger using spdlog for both the Python and C++ side.
- Importing TidalPy no longer warns about the backend change. `TidalPy.exceptions.TidalPyDeprecationWarning` still exists, so code that filters is fine.

## Module Map

Imports are case-sensitive on every operating system: `import TidalPy.rheology` raises `ModuleNotFoundError` even though `TidalPy.Rheology` exists.

| 0.7.X | 0.8.0 | Documentation |
|---|---|---|
| `TidalPy.structures` (worlds, layers, `Orbit`) | `TidalPy.Structures` (worlds, layers, `System`) | [Structures](Structures/index.md), [System](Structures/system/system.md), [TOML schema](Structures/config/toml_schema.md) |
| `TidalPy.RadialSolver` | `TidalPy.RadialSolver` | [RadialSolver](RadialSolver/index.md) |
| `TidalPy.Material.eos` | `TidalPy.Material` (equation-of-state and shear-modulus laws, `Phase`, `Material`, MatPack); `TidalPy.Material.eos` holds the whole-planet solver | [Materials](Material/index.md), [MatPack](Material/matpack.md) |
| `TidalPy.rheology` | `TidalPy.Rheology` | [Rheology](Rheology/index.md) |
| `TidalPy.rheology.viscosity` | `TidalPy.Viscosity` | [Viscosity](Viscosity/index.md) |
| `TidalPy.rheology.partial_melt` | `TidalPy.PartialMelt` (melting curves, melt weakening, and bulk mixing, which a material combines) | [Partial Melting](PartialMelt/index.md) |
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
| `[layers.ice]`, `[layers.rock]`, `[layers.iron]` | Removed. Materials come from MatPack (see [Layers and Materials](#layers-and-materials)). `[layers]` holds only `material`, the MatPack material of a layer that names none (default `simple_rock`). |
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
TidalPy.save_config("my_run_config.toml")                                   # Record the settings of a run
TidalPy.reinit(provided_config="default")                                   # Back to your saved file
```

## Logging

0.7.X logged through Python's `logging` module (`TidalPy.logger.get_logger`). 0.8.0 logs through one C++ logger (spdlog), which Python code reaches through `TidalPy.Utilities.logging`:

```python
from TidalPy.Utilities.logging import log_info, set_console_level

set_console_level("debug")                            # Show debug messages in the console for this session
log_info("Starting the Io run")                       # Written to the same console and file as TidalPy's messages
```

Handlers attached to Python's `logging` (including pytest's `caplog`) do not see these messages. The log file keeps its name, `TidalPy_<YYYYMMDD-HHMMSS>.log`, and is written only when `[logging] write_log_to_disk` is set.

## Worlds

The bundled worlds are new TOML files, and the set changed:

- `earth` is replaced by `earth_prem` (built from the PREM profile) beside `earth_simple`, and `earth_prem_q` which uses PREM's quality factors.
- `jupiter_simple` is new, and so are `europa_dynamic`, `luna_dynamic`, `mercury_dynamic`, and `pluto_dynamic`, which solve a liquid layer (an ocean or a fluid outer core) as a dynamic, compressible liquid.
- The warm silicate mantles of the bundled worlds now use the peridotite melting curves, with melting and pressure melting on.
- `io_simple`, `triton_simple`, `55cnc`, `55cnce`, `55cnce_simple`, and `nereid_dev` are removed.

`TidalPy.Structures.available_worlds()` lists the bundled worlds. They are copied to `<Documents>/TidalPy/0.8.X/Worlds` for editing; `TidalPy.Structures.install_worldpack(force=True)` restores the packaged copies and discards those edits.

A world computes its interior in `solve_eos`, its Love numbers in `solve_love_numbers`, and its tides in `calc_tides`, which takes the orbital state as arguments (a world in a `System` takes any argument left out from the system's orbit). `world.copy()`, `copy.deepcopy`, and `pickle` give an independent world with the same layers, models, and settings.

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

## Layers and Materials

0.7.X built a layer's interior from a material `type` (`rock`, `ice`, `iron`) whose defaults lived in the configuration file, optionally through BurnMan. 0.8.0 has one layer class, `Layer` (`TidalPy.Structures.layers`), and gives each layer a material: a MatPack name, a MatPack preset with overrides, or a full material table of phases and melting laws (see [Materials](Material/index.md) and [Layer](Structures/layers/layer.md)). The layer's physics switches (`use_thermal_expansion`, `use_melting`, `use_pressure_melting`, `use_melt_density`, `use_heating`, and `use_tides`) say how much of the material it uses; each is off by default except `use_tides`.

| 0.7.X layer key | 0.8.0 |
|---|---|
| `PhysicsLayer`, `GasLayer`, `LayerBase` | `Layer`, one class for every layer |
| `type` (`rock`, `ice`, `iron`) and its `[layers.<type>]` defaults | `material = "<MatPack name>"` (`TidalPy.Material.available_materials()`), or a `[layers.<name>.material]` table |
| `radius` (the outer radius) | one of `radius_outer_m`, `radius_fraction`, or `volume_fraction` |
| BurnMan `material`, `material_source`, `material_fractions` | the equation-of-state law of the material's phases (`birch_murnaghan`, `vinet`, `murnaghan`, `polytrope`, `modified_polytrope`, or a tabulated `interpolate` profile) |
| `density` | a `constant` equation of state, or the `simple_*` MatPack materials |
| `is_tidal` | `use_tides` |
| `temperature_mode = "user-defined"`, `temperature_fixed` | `temperature_k` with no cooling model: one temperature throughout |
| `temperature_mode = "adiabatic"`, `temperature_top` | `temperature_k` with a `convection` cooling model, in a solve with `solve_temperature=True`: the layer's temperature applies at the top of its adiabatic interior |
| `shear_modulus`, `thermal_conductivity`, `thermal_expansion`, `heat_fusion` | the phase's `shear_modulus` law and `thermal_conductivity_w_mk`, the equation of state's `thermal_expansion_1_k`, and the material's `latent_heat_j_kg` |
| `solid_viscosity`, `liquid_viscosity` | the `shear_viscosity` of the material's `solid` and `liquid` phases |
| `partial_melting` (`model`, `solidus`, `liquidus`) | the material's `melting` table (`solidus` and `liquidus` curves and a `weakening` law), used when the layer sets `use_melting` |
| `rheology` | the layer's `shear_rheology` table, or the default `shear_rheology` of the material's phase |
| `radiogenics`, `cooling` | `[layers.<name>.radiogenics]` and `[layers.<name>.cooling]` tables, or the layer's `radiogenics` and `cooling` properties in Python |

```python
from TidalPy.Material import available_materials, load_material
from TidalPy.Structures import build_world
from TidalPy.Structures.layers import Layer

print(available_materials("rocky"))   # The MatPack names in one category

# A world from a dict, matches the format of the toml files
io_like = build_world({
    "schema_version": "0.2.0",
    "name": "Io-like",
    "type": "terrestrial",
    "radius_m": 1.8215e6,
    "mass_kg": 8.93e22,
    "layers": {
        "core": {
            "radius_fraction": 0.45,
            "material": "simple_iron_core",   # A MatPack name
            "use_tides": False,               # Was is_tidal
            "temperature_k": 1800.0},
        "mantle": {
            "radius_fraction": 1.0,
            "material": {
                "preset": "peridotite",       # A MatPack material with an override
                "solid": {"shear_rheology": {"model": "maxwell"}}},
            "temperature_k": 1600.0,
            "use_melting": True,              # Its melting curves and weakening take part
            "use_pressure_melting": True,
            "cooling": {"model": "convection"}}}})
io_like.solve_eos(
    solve_temperature=True,
    surface_temperature=110.0)                # An adiabatic mantle under a conducting boundary layer

# The same mantle built in Python
mantle = Layer(
    "mantle",
    1,
    8.2e5,
    1.8215e6,
    material=load_material(
        "peridotite",
        solid={"shear_rheology": {"model": "maxwell"}}),
    temperature=1600.0,
    use_melting=True,
    use_pressure_melting=True,
    cooling="convection")
```

## Systems

`System` replaces `Orbit` (`PhysicsOrbit`). It returns the derivative rates for the current state, and `System.evolve` integrates a world's orbit, spin, and layer temperatures about its host over time (Demo S02); for any other time evolution, integrate the rates with an integrator of your choice.

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

`TidalPy.toolbox` (`quick_tidal_dissipation` and `quick_dual_body_tidal_dissipation`) has no direct replacement. But, `calc_world_evolution` returns the same outputs for any world in a `System` (the tidal heating, the potential derivatives, and the orbit and spin rates). The world's `love_method` sets how its Love numbers are found:

- `"radial_solver"` (the shooting method, and the default) integrates the radial equations through every layer.
- `"homogeneous"` treats each tidal layer as a homogeneous sphere of its averaged material, with no radial solve, as `quick_tidal_dissipation` did.

The bundled Io has a core, a mantle, and a thin asthenosphere that does almost all of the dissipating. Continuing the [Systems](#systems) example:

```python
for love_method in ("radial_solver", "homogeneous"):
    io.set_tide_config(love_method=love_method)
    rates = system.calc_world_evolution(io)
    print(love_method, rates["tidal_heating"])       # [W]: 9.33e13 (radial solver), 1.37e14 (homogeneous)
    print(rates["dU_dM"], rates["dU_dw"], rates["dU_dO"])
    print(rates["da_dt"], rates["de_dt"], rates["dspin_dt"])   # [m s-1], [s-1], [rad s-2]
```

The radial solver reproduces, _e.g._, Io very well (calibrated using Lainey et al. 2009), in about 1.2 ms per call. The homogeneous method weights the layers by volume, which only approximates a thin weak layer inside a stiffer planet. It produces less accurate results (sometimes significantly so) but is about 10x faster. Use it for quick sweeps and first estimates. See [Love numbers](Tides/love/love_numbers.md) for the other methods.

A one-layer world built from bare numbers reproduces `quick_tidal_dissipation`. The world below is its default homogeneous Maxwell body at degree 2 and eccentricity truncation 2, and the result matches the homogeneous-sphere closed form.

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
            "radius_fraction": 1.0,
            "temperature_k": 1600.0,
            "material": {
                "solid": {
                    "eos": {
                        "model": "constant",
                        "reference_density_kg_m3": density},
                    "shear_modulus": {
                        "model": "constant",
                        "shear_modulus_pa": 6.0e10},
                    "shear_viscosity": {
                        "model": "constant",
                        "reference_viscosity_pas": 1.0e19}}},
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

For dual-body dissipation, as in `quick_dual_body_tidal_dissipation`, make each body the other's tidal host with `system.set_tidal_host(host, target)` and call `system.calc_pair_evolution(target)`. Each body dissipates through its own tide model. A gas giant built from a file whose `[tides]` table names no model takes the builder default, `fixed_dt`; the bundled gas giants (`jupiter`, `jupiter_simple`, `neptune`) set `fixed_q`. The per-degree Love numbers that `quick_tidal_dissipation` returned come from `target.get_tidal_love_k(l, m, p, q)` after a `calc_tides` call, or from the closed-form functions in `TidalPy.Tides.love`.

## Radial Solver

`radial_solver` takes the same positional arguments and nearly the same keywords as in 0.7.X. The differences:

- `use_prop_matrix=True` is now `love_method="propagation_matrix"`.
- `use_kamata=True` is now `starting_method="kamata"` (and `use_kamata=False`, the default, is `starting_method="takeuchi"`). The same key replaces `use_kamata` in the `[radial_solver]` configuration and world tables, and two new methods join it: `"power_series"` and `"unity"` (see [Starting Conditions](RadialSolver/starting_conditions.md)). Takeuchi is still the recommended starting condition for most problems.
- The solver settings (`integration_method`, `integration_rtol`, `integration_atol`, `expected_size`, the `eos_*` arguments, and the rest) default to `None`, which reads the `[radial_solver]` and `[eos_solver]` sections of the configuration. The packaged `integration_rtol` and `integration_atol` are both 3e-8, against 1e-5 and 1e-8 in 0.7.X.
- `solve_for` is a tuple of case-insensitive strings (_e.g._, `("tidal", "loading")`).
- Invalid inputs raise `ValueError` instead of `ArgumentException` or `UnknownModelError`.
- When several boundary conditions are solved for, `k`, `h`, and `l` are complex128 arrays; after a failed solve they are complex128 NaN arrays (float64 in 0.7.X).
- `moi_factor` is now the conventional $C/(M R^2)$, 0.4 for a uniform sphere. The 0.7.X value, $C/(0.4 M R^2)$, is now `moi_sphere_ratio`. Code that reads `moi_factor` gets a value 2.5 times smaller.
- The solution is evaluated at any radius through dense output (`get_radial_solution(radius)`), and `plot_ys` and `plot_interior` take `show_plot` and plotting keywords.
- The input builders `build_rs_input_homogeneous_layers` and `build_rs_input_from_data` keep their argument names. Their rheology arguments take `TidalPy.Rheology` models or model names, and one model can stand in for every layer. `perform_checks` is accepted and ignored: inputs are always validated.

```python
import numpy as np

from TidalPy.RadialSolver import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.Rheology import Elastic, Maxwell

# Build the radial solver inputs with a helper (see "Helper Functions" below)
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

The rheology, viscosity, partial-melt, cooling, radiogenics, and luminosity functions of 0.7.X are now classes in 0.8.0, built by name with a helper builder function (called factories). Config keys carry their units (`reference_viscosity_pas`, `solidus_k`), and each model refuses a key it does not read with `ValueError`, naming the closest key it does.

| 0.7.X | 0.8.0 |
|---|---|
| `find_rheology(name)` | `make_rheology(name, config)`, which returns a model instance |
| `rheology(frequency, modulus, viscosity)` | `rheology.calc_complex_modulus(modulus, viscosity, frequency)` (note the order) |
| `vectorize_frequency`, `vectorize_modulus_viscosity` | `calc_complex_modulus_vectorize_frequency`, `_vectorize_modulus`, `_vectorize_all` |
| `Voigt(args=(5.0, 0.02))`, `Andrade(args=(0.3, 1.0))`, `change_args` | Keyword parameters: `Voigt(voigt_modulus_frac=5.0, voigt_viscosity_frac=0.02)`, `Andrade(alpha=0.3, zeta=1.0)` |
| `Newton`, `SundbergCooper`; the names `voigtkelvin` and `sundbergcooper` | `Viscous`, `Sundberg`; the factory takes `newton`, `voigt-kelvin`, and `sundberg-cooper` as aliases |
| the complex compliance functions (`rheology.complex_compliance`) | removed: the models return the complex modulus, whose reciprocal is the compliance |
| `rheology.viscosity` functions (`arrhenius`, `reference`, `constant`) | `make_viscosity(name, config)` and `calc_viscosity(temperature, pressure)` |
| `rheology.partial_melt` (`spohn`, `henning`, `calculate_melt_fraction`) | `make_melting_curve(name, config)` for the solidus and liquidus and `make_melt_weakening(name, config)` (`none`, `spohn`, `henning`), combined in a `Material`; `Material.calc_state(pressure, temperature, use_melting=True)` returns the melt fraction and the weakened shear modulus and viscosity (see the note below) |
| `cooling` functions (`convection`, `conduction`, `off`) | `make_cooling(name, config)` and `calc_cooling`, or the direct functions `convective`, `conductive`, `cooling_off`; `off` reports a boundary-layer thickness of NaN, where 0.7.X gave half the layer thickness |
| `radiogenics` functions (`isotope`, `fixed`, `off`) with times in Myr | `make_radiogenics(name, config)` and `calc_heating(time, mass)` with times in seconds; isotope sets are named datasets |
| `stellar.luminosity_from_mass` | `TidalPy.Stellar.mass_to_luminosity(mass)` or `make_luminosity("mass_to_luminosity")` |

Both melt-weakening laws are continuous where 0.7.X stepped: across the breakdown band (`crit_melt_frac` to `crit_melt_frac + crit_melt_frac_width`) the aggregate blends into the liquid's values instead of jumping to them at the band's end. Spohn starts from the solid's own values at the solidus, unless given the absolute anchors `fs_visc_log10_at_solidus = 15.875` and `fs_shear_log10_at_solidus = 10.65` of the 0.7.X fit.

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

- The eccentricity functions are unsquared: `TidalPy.Tides.eccentricity_func(eccentricity, degree_l, truncation)`, with `TidalPy.Tides.eccentricity.eccentricity_squared_func` for the squares. 0.7.X tabulated the squared functions per degree and level (`eccentricity_funcs_l2_trunc10` and so on).
- As in 0.7.X, truncation level $N$ keeps every product of two eccentricity functions, and so every heating term, through $e^N$. Levels 2 to 10 and 20 keep the same terms in both versions.
- 0.8.0 tabulates levels 2, 4, 6, 8, 10, 20, and 50, plus `"exact"`. A configured level of 12 to 18 or 22 is promoted to the next tabulated level with a warning, and the direct functions raise `NotImplementedError` for it. The default level is now 10 (it was 6). `TidalPy.Tides.eccentricity.recommend_eccentricity_truncation` picks a level for a given eccentricity.
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

Timings of TidalPy 0.7.6 against 0.8.0 on one machine; treat the ratios as rough magnitudes and measure your own workload.

| Task | 0.7.X | 0.8.0 | Change |
|---|---|---|---|
| Build a planet with its interior (Io, 3 layers) | 187 ms | 0.84 ms | 220x faster |
| Orbit-averaged 3D heating map (50 x 16 x 32), 1 thread | 13.4 ms | 2.14 ms | 6.3x faster |
| Global tidal heating, degrees 2 to 4, e^10 | 0.051 ms | 0.0135 ms | 3.8x faster |
| Rheology, 10k complex moduli | 0.176 ms | 0.057 ms | 3.1x faster |
| `radial_solver`, 1 layer, 200 slices | 0.46 ms | 0.53 ms | 1.1x slower |
| Convective cooling, one evaluation | 0.17 us | 1.50 us | 9x slower |
| First call of a tidal heating calculation | 1.1 to 7.0 s (numba compile or cache load) | 0.08 ms | |

The largest gains are where 0.7.X called BurnMan or compiled numba kernels; 0.8.0 has nothing to compile on first call. The standalone `radial_solver` is about as fast as in 0.7.X, because it builds a temporary world and runs a more accurate EOS solve.

> [!TIP]
> We recommend the world-based approach. It is as fast as the standalone solver, and reusing a constructed `BaseWorld` is much faster than repeated standalone calls: its equation of state is solved once, its Love solves are cached per degree and frequency, and `calc_tides` reuses them across modes.

A single scalar call is dominated by the crossing from Python, so use the vectorized calls (about 3x faster than 0.7.X) for more than a handful of values. The direct cooling functions (`convective`, `conductive`) also build their model on every call (about a microsecond); a model object's `calc_cooling` does not.

### Threads for 3D Grids

The 3D grid methods (`calc_3d_tides`, `calc_3d_stress_strain`, `calc_3d_displacements`, and `get_3d_tidal_heating_array`) take `num_threads`, which 0.7.X lacked. The default, 0, uses the number of logical processors minus 4 (at least 1), and every thread count returns identical values. On a 16-thread machine, 16 threads run a 3D grid about 6 to 7 times faster than one. On a very small call, starting the threads costs more than they save, so pass `num_threads=1` for many small calls. `calc_tides` also spreads its Love solves over threads (see [Parallel Love Solves](RadialSolver/parallel.md)).

## Learning TidalPy 0.8.0

Start with the [Getting Started](Overview/1_Getting_Started.md) page and the notebooks in the Demos section: `Basics` (configuration, world building, save and load), `Physics` (orbits, tides, rheology, Love numbers, 3D heating, thermal interiors, truncation levels), and `Systems` (multi-world systems, coupled thermal-orbital evolution, the Earth-Moon-Sun system). The Benchmarks section validates TidalPy against published results and tracks its performance.
