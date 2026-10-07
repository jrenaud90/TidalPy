# Worlds (`Structures.worlds`)

_Updated: 2026-10-06_

The world classes are the top-level structural objects in TidalPy. A world owns its identity, orbital and thermal scalars, bulk geometry, spin model, and tide model, and an ordered stack of [layers](../layers/layer.md), which may be empty. It runs the whole-planet equation-of-state (EOS), thermal, radial (Love number), and tidal solves, and holds the heat sources that act inside its layers.

## Inheritance

```
TidalPyBaseClass
  └── StructureBase
        └── BaseWorld
              ├── TerrestrialWorld
              ├── GasGiantWorld
              └── StarWorld
```

| Class | TOML `type` | Purpose |
|-------|-------------|---------|
| `BaseWorld` | `layered` | A world: layers, identity, albedo/emissivity/obliquity/spin, bulk geometry, equilibrium temperature, spin and tide models, heat sources, and the EOS, Love, and tidal solves. |
| `TerrestrialWorld` | `terrestrial` | A rocky or icy planet or moon; a `BaseWorld` with its own type and binary id. |
| `GasGiantWorld` | `gasgiant` | A gas giant; a `BaseWorld` with its own type and binary id. |
| `StarWorld` | `star` | A star; adds the effective temperature and luminosity, linked by the Stefan-Boltzmann law. Usually has no layers. |

## `BaseWorld`

```python
from TidalPy.Structures.worlds import BaseWorld

welcome_to_earth = BaseWorld(
    name="Earth",
    radius=6.371e6,
    mass=5.972e24,
    world_type="terrestrial",
    albedo=0.3,
    emissivity=1.0,
    obliquity=0.41,
    spin_frequency=7.29e-5,
)                                                     # A world with no layers yet
```

**Read-only properties:** `name`, `world_type`, `radius`, `mass`, `albedo`, `emissivity`, `obliquity`, `spin_frequency`, and `moment_of_inertia_factor` ($C / (M R^2)$ of the solved structure: the solved moment of inertia over the solved mass times the radius squared; NaN before a successful `solve_eos`, while a spin calculation uses the spin model's factor, see [Other World Members](#other-world-members)).

**Methods**

| Method | Returns | Description |
|--------|---------|-------------|
| `calc_surface_gravity()` | float [m/s²] | $GM/R^{2}$. |
| `calc_escape_velocity()` | float [m/s] | $\sqrt{2GM/R}$. |
| `calc_mean_density()` | float [kg/m³] | $M / (\tfrac{4}{3}\pi R^{3})$. |
| `calc_equilibrium_temperature(F)` | float or array [K] | $\left[(1-A)\,F/(4\varepsilon\sigma)\right]^{1/4}$ (fast rotator, $F$ = insolation flux [W/m²], a float or an array of any shape). |
| `set_spin_frequency(ω)` | - | Set rotation rate [rad/s]. A `System` reads it when it builds the world's tidal state; `calc_tides` takes its spin rate as an argument, so a tidal result already solved is left alone. |
| `set_obliquity(θ)` | - | Set axial obliquity [rad]. The same holds as for the spin rate. |
| `copy()` | world | An independent world of the same class with the same state, built through the binary record (see [Binary Serialization](#binary-serialization)). `copy.copy`, `copy.deepcopy`, and `pickle` use it, so a world can be sent to a process pool. |
| `summary()` | str | A multi-line description: the world, then one row per layer with its name, radii [km], state, material phases, solved density range [kg/m³], temperature [K], and shear rheology. |

`repr(world)` is one line with the class, name, radius [km], mass [kg], and number of layers, as in `TerrestrialWorld('Io', radius_km=1821.49, mass_kg=8.9298e+22, num_layers=3)`.

```python
from TidalPy.Structures import build_world

io = build_world("io")                                # A bundled world
io.solve_eos()                                        # The solved structure fills the summary's densities
print(io.summary())                                   # One row per layer
print(io.moment_of_inertia_factor)                    # C / (M R^2) of the solved structure
```

`get_config_dict()` returns the world as the TOML builder's world table: `schema_version`, `name`, `type` (the builder's world type, from `get_builder_world_type()`), `radius_m`, `mass_kg`, `albedo`, `emissivity`, `obliquity_rad`, `spin_frequency_rad_s`, `moment_of_inertia_factor` (the spin model's, written when the world's source configuration gave it or it differs from the `[worlds]` default for the world type, which a rebuild takes again), a `tides` table when a tide model is attached (`global_tidal_model`, its per-degree parameters, and the settings from `get_tide_config()`), the `eos_solver` and `radial_solver` tables when the world pins any solver key, and a `layers` table when the world has layers. Each layer entry is the layer's own config dict (scalars, `material`, and model sub-tables; `mass_kg` only once the layer has a mass) without the standalone-only keys the builder derives itself, so `build_world(world.get_config_dict())` rebuilds the same structure. A world constructed directly in Python with no tide model attached writes no `tides` table, so its rebuild takes the builder's default tide model for its type (`rheology` for a terrestrial or layered world, `fixed_dt` for a gas giant, `fixed_q` for a star). `save_config` / `save_binary` / `load_binary` are inherited from `TidalPyBaseClass`; the binary methods take a `str` or `os.PathLike` path. `save_to_toml` validates the dict against the schema before writing when no build configuration is retained.

## Layers

A world holds an ordered (inner-to-outer) stack of [layers](../layers/layer.md).

```python
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds import TerrestrialWorld

world = TerrestrialWorld("Earth", 6.371e6, 5.972e24)
world.add_layer(Layer(
    "core",
    radius_outer=3.485e6,
    material="simple_iron_core",
    temperature=4000.0))                              # Innermost layer: index 0, starting at r = 0
world.add_layer(Layer(
    "mantle",
    radius_outer=6.371e6,
    material="simple_rock",
    temperature=1600.0))                              # Index 1, starting where the core ends
```

A layer constructed without a `layer_index` or a `radius_inner` takes them from its place in the stack when it is added. The positional form `Layer("core", 0, 0.0, 3.485e6, ...)` gives both, and `add_layer` then checks them.

**Layer management**

| Member | Description |
|--------|-------------|
| `add_layer(layer)` | Add a layer inner-to-outer. A layer constructed without a `layer_index` gets its place in the stack (the number of layers already added), and one without a `radius_inner` starts at the current outermost radius (0 for the first layer). Ownership of the layer (and its attached models) transfers into the world; the passed wrapper stays usable as a non-owning view of that layer, like the one `world.<layer name>` returns. Raises `ValueError` if the layer was already added, if a given index is not its place in the stack, if it does not continue the stack (the innermost starts at 0), if it reaches past the world radius, or if its name is taken; a refused layer is left as it was constructed. |
| `num_layers` | Number of layers (property). |
| `calc_total_mass()` | Sum of the layer masses [kg]; equals `planet_mass_eos` after a successful EOS solve. |
| `calc_internal_heating(time)` | Sum of the radiogenic heating [W] of every layer with a radiogenics model, at a time [s], whatever the layer's `use_heating`. Uses each layer's `mass`, so solve the EOS first when the layers were built without one. |
| `validate_layers()` | `True` if every boundary is continuous and the innermost starts at 0. |

**Accessing layers**

Python reaches a world's layers through non-owning views:

| Access | Returns |
|--------|---------|
| `world.get_layer(i)` / `world[i]` | the layer at index `i` (0 = innermost; negative indices count from the outermost), or with a name, the layer of that name (`world["mantle"]`) |
| `world[a:b]` | a list of the sliced layer views |
| `world.layers` | a list of all layer views, inner to outer |
| `world.<layer_name>` | the layer with that name (e.g. `world.mantle`) |
| `for layer in world:` | iterate the layers inner-to-outer; `len(world)` is the layer count (a world is always true, even with no layers) |

Every world method that takes a layer (`get_layer`, `world[...]`, `get_layer_tidal_heating`, `get_layer_tidal_scale`, `set_prescribed_heating`, `calc_layer_temperature_rate`, `calc_layer_thermal_capacity`, `calc_layer_latent_capacity`) takes its index or its name the same way: an index out of range raises `IndexError`, a name the world does not have raises `KeyError` listing the layer names (with the closest one), and anything else raises `TypeError`.

Each view is a `Layer` with the layer's full API (`world.mantle.temperature = 1700.0`, `world.mantle.get_tidal_heating()`, `world.core.calc_complex_shear_modulus(r, ω)`). A setting changed through a view reaches the world, which forgets its solved structure when the solve reads that setting (see [Changes That Clear a Solve](../layers/layer.md#changes-that-clear-a-solve)). The view keeps the world alive, so it is safe to hold. Views are built once and cached (rebuilt only when a layer is added), so repeated access returns the same object (`world.mantle is world.mantle`). Access by layer name runs only after normal attribute lookup, so defined members win, and it ignores names starting with `_`.

## Equation of State

Each layer's [material](../layers/layer.md#material) is its density source. Once every layer has one (`all_materials_set`), `BaseWorld.solve_eos(...)` integrates the planet's radial structure from center to surface. It populates every layer's profile (density, gravity, pressure, temperature, heat flow, and the material state) and sets each layer's `mass` (and so its `density_bulk`) to the mass the solved profile places in the layer. Before the solve, a layer's mass is the value it was constructed with. The TOML builder uses 0.0 when a file gives none, which is the usual case.

```python
result = world.solve_eos(surface_pressure=0.0)  # EOSResult: a dict of profile arrays and scalars
rho    = world.get_density(5.0e6)               # [kg/m³] at radius 5000 km
g      = world.get_gravity(world.radius)        # Surface gravity [m/s²]
p0     = world.get_pressure(0.0)                # Central pressure [Pa]
print(result["success"], world.planet_mass_eos) # True and the solved mass [kg]
```

The central pressure is found by a secant iteration on the surface-pressure mismatch. The first step assumes a unit slope (exact for an incompressible planet), and later steps use the slope measured between iterations, so a compressible planet converges in a few steps. For large, soft planets the surface pressure can at first fall as the central pressure rises. While the measured slope is not positive, the step doubles each pass. Once the mismatch has changed sign, the root is bracketed and found by secant, false-position, and bisection steps (as in Brent's method). Each step integrates the whole structure. Only the converged integration keeps the dense output the profile getters and Love solves read, so one that converges without it is run again. When every layer's density is set by its radius and temperature alone (`constant` and `interpolate` equations of state, and no melt-density mixing in a melting range that moves with pressure, as in a profile world or the standalone `radial_solver`), nothing else in the structure depends on the pressure, so it is left out of the integrator's step control. Every step then takes the same radii, the unit-slope step is exact to rounding, and the integration after it keeps its output: a solve from scratch takes two integrations instead of three or more. The pressure's accuracy then follows from the steps gravity, mass, and moment of inertia need. The integration runs in non-dimensional units (the planet radius, its bulk density, and $1/\sqrt{\pi G \rho}$ as the length, density, and time units), so the tolerances mean the same thing for every planet. Every result is returned in SI.

**`solve_eos(*, surface_pressure=0.0, slices_per_layer=None, G_to_use=None, integration_method=None, rtol=None, atol=None, pressure_tol=None, max_iters=None, nondimensionalize=None, temperature=None, solve_temperature=None, surface_temperature=None, reset_layer_masses=False, verbose=False, time=None, max_thermal_passes=None, thermal_tol=None, raise_on_fail=False) -> EOSResult`**

Every argument is keyword-only. Every solver setting left as `None` takes the `[eos_solver]` value of the TidalPy configuration (see [Configurations](../../Overview/2_TidalPy_Configurations.md)), the same defaults the standalone `radial_solver` uses, unless the world pins it (see [Other World Members](#other-world-members)). `G_to_use` left as `None` is the configured gravitational constant. With `raise_on_fail=True` a failed solve raises `SolutionFailedError` (a `RuntimeError`) with the solve's message instead of returning `success = False`. `pressure_tol` is relative to the central-pressure scale $(2/3) \pi G \rho^2 R^2$ and must stay above `rtol`, the integrator's own noise on the surface pressure. Hitting `max_iters` sets `max_iters_hit` in the result and runs one last pass. The solve succeeds only if that pass meets `pressure_tol`, and otherwise fails as described below. `temperature` gives every layer one temperature \[K\] for this solve in place of its own. `solve_temperature`, `surface_temperature`, `time`, `max_thermal_passes`, and `thermal_tol` belong to the [thermal solve](#temperature-and-heat-flow), and `reset_layer_masses` to the [layers that hold their mass](#layer-size).

### Solved State

A re-solve starts its central-pressure iteration from the world's last converged solve instead of from a uniform sphere, so after a small change (a new layer temperature, a later time) it usually converges in one or two passes. The result then depends on the call history at the level of the tolerances. A world solved twice at different surface pressures can differ from a fresh one by a few parts in $10^9$ in its central pressure and about one part in $10^9$ in $k_2$. Build a fresh world when results must reproduce bit for bit.

A solve commits its result only when it finishes. Materials and models are immutable and shared, so a solution released from it (`release_radial_solution`) keeps the materials it was solved with, and a later solve never changes a solution already provided to the user. Anything that changes what the solve reads clears every solved result:

- A change to the layer stack (`add_layer`, `load_binary`) or a layer moved with `set_radii`
- A layer's material, temperature, material switches, `use_heating`, `is_volume_fixed`, cooling model, or radiogenics model
- A layer `state` change that decides whether the layer can change state
- `set_prescribed_heating`

The profile getters then return NaN, `eos_solved` is `False`, `zones` is empty, and a Love solve or a `rheology` `calc_tides` raises errors until `solve_eos` runs again, as before the first solve. Every Love solve reads the radial-solver flags (`state`, `is_static`, `is_incompressible`) and the rheologies afresh, so changing those keeps the solved structure. A failed solve also leaves the world unsolved, keeping only its diagnostics (`success`, `message`) and NaN profile arrays.

A successful solve whose mass (`planet_mass_eos`) differs from the world's stated `mass` by more than 1 percent logs a warning, once per world, naming both. The two then describe different planets: the Love numbers, tides, and moment of inertia follow the solved structure, while the orbit in a `System`, `calc_surface_gravity`, and `calc_mean_density` use the stated mass. Adjust the layers' materials or radii, or the stated mass, until they agree.

A world runs one heavy call at a time: threads sharing a world take turns on `solve_eos`, the Love solves, `calc_tides`, the 3D calls, `release_radial_solution`, and `load_binary`. These calls release the GIL, so separate worlds can run in parallel. A read of a property called during another thread's `solve_eos` waits and then reads the new profile. A call on an array of radii takes one turn, so all its values come from one solve. Give each thread its own world when reads must run in parallel. The setters of a world's tide model and settings (`set_tide_model`, `set_tide_config`, `set_spin_model`, `set_solver_defaults`) and of its layers' models and flags take the same turns, so a change waits for a running solve instead of landing in the middle of it. The lock covers each call on its own, not a solve followed by a read of its result, so threads sharing a world can read each other's results. [Parallel Love Solves](../../RadialSolver/parallel.md) shows how to avoid that and how to run solves on thread and process pools.

A solve that reaches `max_iters` with its surface pressure still off the target by more than `pressure_tol` has found no hydrostatic structure. It returns `success = False` with the message "no hydrostatic structure", sets `max_iters_hit`, and leaves the world unsolved. A converged solve fails the same way when its enclosed mass differs from the world's stated mass by more than the factor `[numerical] maximum_eos_mass_ratio` (default 10) either way. Layers with no hydrostatic structure near that mass (a core much too dense for its radius, say) can still meet the surface pressure on a collapsed branch at an absurd central pressure. The message says how far the mass is off.

A converged solve also fails when a layer is in tension past what its material's pressure law represents (Birch-Murnaghan or Vinet). The law sees the pressure less the thermal pressure $\alpha_0 K_0 (T - T_\mathrm{ref})$ of a layer with `use_thermal_expansion`, so a hot layer whose material has a large thermal expansivity or a large $K_0'$ can fall below the law's tension limit. There the density is held at the law's smallest compression and the bulk modulus is near zero. The message names the layer and gives the pressures. A layer past the law's compression end (a Birch-Murnaghan $K_0'$ below 4 turns over) is held at the law's largest compression there, and the solve logs a warning and stands.

The result is an `EOSResult` (`TidalPy.Structures.worlds.EOSResult`), a `dict` subclass that behaves as a plain `dict` (item access, equality with a `dict` of the same entries, copying, and pickling). Its `repr` is a short summary instead of every profile array: `success`, `iterations`, `message`, the planet mass \[kg\], radius \[m\], and central pressure \[Pa\], and the list of keys, so a notebook cell that ends in a solve prints a few lines. It contains:
- `success`, `message`, `iterations` (central-pressure steps), `structure_integrations` (integrations of the whole structure over every thermal pass, the repeat of a converged one included: the measure of the solve's cost), `max_iters_hit`, and `pressure_error` \[Pa\].
- The profile arrays: `radius`, `gravity`, `pressure`, `mass`, `moi`, `density`, `temperature`, `heat_flow`.
- The scalar results: `surface_gravity`, `surface_pressure`, `central_pressure`, `planet_mass`, `planet_moi`.
- The solid and liquid zones, `zones`, as the [`zones`](#pieces-and-zones) property reports them.
- The thermal iteration report: `thermal_passes` and `thermal_converged`.
- The per-layer lists, one entry per layer, inner to outer:
  - `layer_radius_outer` \[m\]: where each layer ended (it moves for a layer that holds its mass).
  - `layer_temperature`, `layer_node_temperature` (at the interface above the layer), `layer_top_temperature` and `layer_base_temperature` (the two ends of a convecting interior, whose top is the layer's own temperature) \[K\].
  - `layer_heat_flow_in`, `layer_heat_flow_out`, and `layer_heating` \[W\], with each heat source's part of the heating in `layer_heating_radiogenic`, `layer_heating_tidal`, and `layer_heating_prescribed`.
  - `layer_boundary_thickness` \[m\], `layer_rayleigh_number`, `layer_nusselt_number`, `layer_boundary_fallback`, `layer_magma_ocean`, `layer_reference_pressure` \[Pa\], `layer_reference_viscosity` \[Pa s\], and `layer_reference_melt_fraction`: the cooling model's profile (see [Temperature and Heat Flow](#temperature-and-heat-flow)).
  - `layer_in_thermal_network`: whether the layer has a temperature of its own.

A layer's temperature rate and heat capacities are world methods rather than result entries: `calc_layer_temperature_rate`, `calc_layer_thermal_capacity`, and `calc_layer_latent_capacity` compute them from the last solve on the first call after it, since their quadratures cost more than an isothermal solve (see [Temperature and Heat Flow](#temperature-and-heat-flow)).

### Pieces and Zones

The solve integrates each layer in pieces. A piece ends at the first of three conditions:

- A radius: the top of a layer that holds its volume, or the end of one of its thermal segments.
- An enclosed mass: the top of a layer that holds its mass (see [Layer Size](#layer-size)).
- A change of state: the radius where the material of a layer that can change state crosses from solid to liquid or back.

A layer can change state when its `state` is `"auto"`, its `use_melting` is on, and its material has a solid and a liquid phase (`Layer.can_change_state`). Its state at a radius is set by its post-melt rigidity $\mu / (\bar{\rho} g R)$, with $\bar{\rho}$, $g$, and $R$ the world's stated bulk density, surface gravity, and radius: liquid where it falls to `[numerical] minimum_solid_rigidity` (10$^{-6}$ by default) or below, which a fully molten material always does, and solid elsewhere. A partially molten material therefore stays solid until its shear modulus is far too low to integrate as a solid. The integration locates each crossing as a root, to the integration tolerance, and restarts there in the other state, so the boundaries do not depend on `slices_per_layer`.

Consecutive pieces of one layer in one state form a zone. `zones` lists them, inner to outer, as dicts with `layer` (the layer's name), `radius_inner` and `radius_outer` \[m\], `mass_inner` and `mass_outer` (the enclosed mass at each end) \[kg\], and `state` (`"solid"` or `"liquid"`). A layer that cannot change state is one zone, liquid for a liquid layer (`Layer.is_liquid`) and solid otherwise. `molten_regions` lists the liquid zones of layers that are not liquid throughout, the stretches melting made liquid, as `(layer_name, radius_inner, radius_outer)`, and an info-level log message names each one. Both are empty before a solve.

```python
world = TerrestrialWorld("Hot-Earth", 6.371e6, 5.972e24)
world.add_layer(Layer(
    "core", 0, 0.0, 3.485e6,
    material="simple_iron_core",
    temperature=4000.0))
world.add_layer(Layer(
    "mantle", 1, 3.485e6, 6.371e6,
    material="peridotite",
    temperature=1950.0,
    use_melting=True,
    use_pressure_melting=True))                     # Its solidus rises with pressure
world.solve_eos()                                   # Finds where the mantle's rigidity crosses the threshold

for zone in world.zones:
    print(zone["layer"], zone["state"], zone["radius_inner"], zone["radius_outer"])
print(world.molten_regions)                         # The mantle's top 52 km is a magma ocean
```

The radial (Love number) solver takes each zone as a layer of its own. A solid zone keeps its layer's `is_static` and `is_incompressible` choices, and a liquid zone is solved with the liquid equations under the same choices, so a melting mantle (`is_static` true by default) is a static liquid where it is molten. A zone thinner than `[numerical] minimum_zone_fraction` (10$^{-7}$ of the world radius by default, about 0.6 m in an Earth-sized world) takes the state of its thicker neighbor in the layer, since the radial solver cannot integrate across it. Treating a solid of rigidity $10^{-6}$ as a liquid changes the Love numbers by about that fraction. Outside the radial solve the layer is one layer, and its 3D heating takes nothing from its liquid zones (a liquid carries no shear dissipation).

`zones` is what the solve found. A layer forced to `"solid"` or `"liquid"` after the solve is one zone in that state to the next Love solve, which reads the flags afresh, while `zones` keeps the solved states until the next `solve_eos`. A change between `"auto"` and a forced state that decides whether the layer can change state makes the world forget its solve.

> [!NOTE]
> Only a layer that can change state is split. A solid layer given a near-zero shear modulus some other way (a material without melting, or `use_melting` off) is solved as a solid, and the solver may fail on it or return $\mathrm{Im}[k]$ with the wrong sign and only a conditioning warning.

### Layer Size

A layer holds its volume by default, so a world solves with the boundaries it was built with. A layer with `is_volume_fixed = False` holds its mass instead: the integration ends it where the enclosed mass reaches its mass, inside one solve. Every layer above it moves with it and keeps its own volume, unless it also holds its mass, and the world radius becomes the top of the outermost layer. The layers and the world radius change only when the solve succeeds. The thermal segments of a layer that holds its mass end at fixed fractions of its mass, so its boundary layers move with it.

The mass a layer holds is its `mass` (`mass_kg` in a file) when that is positive. Otherwise the first solve keeps the layer's volume and the layer then holds the mass that solve finds inside its boundaries. `solve_eos(reset_layer_masses=True)` forgets the held masses, so the current boundaries set them again.

```python
world = TerrestrialWorld("Expanding", 2.0e6, 1.0e23)
world.add_layer(Layer(
    "core", 0, 0.0, 1.0e6,
    material="simple_iron_core",
    temperature=1600.0))
world.add_layer(Layer(
    "mantle", 1, 1.0e6, 2.0e6,
    material="simple_rock",
    temperature=1600.0,
    use_thermal_expansion=True,                     # Its density falls as it warms
    is_volume_fixed=False))                         # It holds its mass, not its size
world.solve_eos()                                   # The first solve sets the mass it holds
world.mantle.temperature = 2100.0
result = world.solve_eos()                          # The mantle expands to hold the same mass
print(result["layer_radius_outer"])                 # [m] where the boundaries ended
print(world.radius)                                 # Follows the outermost layer, about 2.01e6
```

A layer that holds its mass reads NaN above its top (the integration that found the top is not read past it).

### Temperature and Heat Flow

Each layer carries its own temperature (`Layer.temperature`) and its [cooling model](../../Cooling/cooling_models.md) sets how heat moves inside it. From these the solve calculates a temperature profile, the heat flow through every radius, and each layer's rate of temperature change. `surface_temperature` \[K\] is the temperature the outermost layer radiates to. When it is left out, no heat leaves the world.

Each cooling model builds its layer's profile. A layer with no cooling model is isothermal, as with `off`.

| Model | Profile inside the layer | Where its temperature applies |
|---|---|---|
| `off` | Isothermal: one temperature throughout, and no modeled gradient, so the layer conducts perfectly. | Everywhere |
| `conduction` | Two conducting halves, $T = T_0 - (L / 4 \pi k)(1/r_0 - 1/r)$. | The mid-radius |
| `convection` | A conducting boundary layer at the base and the top, sized by the model's Nusselt scaling, around an adiabatic interior. | The top of the interior |

The adiabat follows

$$\frac{dT}{dr} = -\frac{(\alpha + \alpha_L)\, g\, T}{c_p}$$

with the material's expansivity $\alpha$ and heat capacity $c_p$ at the local pressure, and $\alpha_L$ the latent heat's share of the expansivity inside a melting range whose curves follow the pressure (`latent_expansion` of `Material.calc_state`, zero outside the range or with `use_pressure_melting` off). The interior warms downward from the layer's temperature, so its base sits at $T \exp\left(\int (\alpha + \alpha_L) g / c_p \, dr\right)$. A layer whose base carries no heat (the innermost layer, or one above a layer outside the network) has no boundary layer at its base. The layer's temperature is applied to the top of the interior, under the upper boundary layer: the upper-mantle temperature of parameterized convection (Stevenson et al. 1983; Schubert et al. 2001).

The convection model evaluates its Rayleigh number's viscosity at that reference point: the top of the interior, at the layer's temperature and the solved pressure there (`layer_reference_pressure`, `layer_reference_viscosity`, `layer_reference_melt_fraction`). The gravity, density, and thermal constants are read at the mid-layer. An interior that is liquid at the reference point (fully molten, or past the rigidity threshold of [Pieces and Zones](#pieces-and-zones)) is a magma ocean and takes the liquid scaling $\mathrm{Nu} = a_\mathrm{liquid} \mathrm{Ra}^{\beta_\mathrm{liquid}}$ (Solomatov 2000), which `layer_magma_ocean` reports. A conducting layer reads its conductivity, heat capacity, and expansivity at the mid-radius pressure and its own temperature.

A layer with no temperature of its own (a temperature that is not positive, such as the 0 K default) is neither a heat sink nor a source, and a neighbor keeps its own temperature at the shared interface. The layer is isothermal at its placeholder temperature, and `layer_in_thermal_network` reports it `False`. Where that leaves its material at the cold, rigid limit of its viscosity law, the first solve logs a warning naming the layer, once until the layer's temperature is set again. The `temperature` argument overrides every layer's temperature with one number.

A convecting layer's Rayleigh number uses the temperature drop across both of its boundary layers: from the top of the layer below (the end of its adiabat, when that layer convects) to the layer's own temperature, plus from that temperature to the layer above or to `surface_temperature`. A mantle at the temperature of the layer above it therefore still convects when the core below it is hotter. The innermost layer, and a layer above one outside the network, has only the upper drop and only the upper boundary layer, since no heat crosses its base. The cooling model's boundary-layer thickness $d / \mathrm{Nu}$ carries its flux across the whole drop, so a layer with two boundary layers gives each half of it and a layer with one gives it the whole, at most 40 percent of the layer each. The flux through the top of a layer whose two drops are equal is then the cooling model's `cooling_flux`, and the boundary layers of a sub-critical layer ($\mathrm{Nu} = 1$) are close to the two conducting halves of a `conduction` layer. A cooling model that gives no usable thickness (usually from a NaN viscosity at the layer's temperature) leaves each boundary layer at 40 percent, `layer_boundary_fallback` reports it, and the solve logs a warning.

The layers form a chain of thermal resistances. A conducting spherical shell between $r_a$ and $r_b$ has

$$R = \frac{1}{4 \pi k} \left( \frac{1}{r_a} - \frac{1}{r_b} \right)$$

and the heat flow through an interface is $L = \Delta T / R$ across the two resistances facing it. Both layer temperatures are inputs, so the flow entering a layer generally differs from the flow leaving it. The difference is the heat the layer stores or releases, reported by `calc_layer_temperature_rate(layer)`:

$$\left(C + C_\mathrm{latent}\right) \frac{dT}{dt} = L_\mathrm{in} - L_\mathrm{out} + H$$

with $H$ \[W\] the heat generated inside the layer (`layer_heating`, see [Heat Sources](#heat-sources)), $C$ \[J K$^{-1}$\] the heat the layer's profile stores per kelvin of its temperature (`calc_layer_thermal_capacity(layer)`), and $C_\mathrm{latent}$ \[J K$^{-1}$\] the latent heat of its zone boundaries (`calc_layer_latent_capacity(layer)`). With the neighbors' interface temperatures and the heating held, a change $\delta T$ of the layer's temperature moves its profile by $S(r)\,\delta T$, so

$$C = \int \rho\, c_p\, S \, 4 \pi r^2 \, dr$$

over the layer, with $c_p$ the material's effective heat capacity at the solved pressure and temperature. $S$ is one through an isothermal layer and $T(r)/T$ along a convecting interior, which scales with the temperature at its top (Stevenson et al. 1983). Across a conducting stretch the change solves Laplace's equation whatever the heating, so $S$ runs linearly in $1/r$ from zero at the interface it meets to one at the end the layer holds ($T_\mathrm{base}/T$ at the base of a convecting interior); below an insulated base it is one. A convecting mantle therefore stores more than $M c_p$ per kelvin of its upper-mantle temperature (about 1.3 times for the bundled Earth). The integral is split where the profile crosses a melting curve. A material that melts over a range carries its latent heat in its effective heat capacity, spread over the range, so only the part of the layer inside the range stores it. A material that melts at one temperature (its solidus and liquidus the same curve) has no range to spread it over, so the boundary between the layer's solid and liquid zones carries it instead: as the layer's temperature changes the boundary moves and melts or freezes mass. With $G(r) = T(r) - T_m(P(r))$ zero at the boundary $r_b$,

$$C_\mathrm{latent} = \frac{L\, 4 \pi r_b^2\, \rho_\mathrm{solid}\, S}{\left| dG/dr \right|}$$

with $L$ the material's latent heat \[J kg$^{-1}$\] and $S$ the sensitivity of the temperature at the boundary to the layer's temperature (one in an isothermal layer, $T(r)/T$ on an adiabat). It assumes the neighbors' interface temperatures hold while the layer's temperature changes, and that the melt that forms takes the solid's density at the boundary. It is zero for a layer without such a boundary.

Where neither side of an interface has a resistance (an `off` layer under another `off` layer, or an `off` outermost layer under the surface), nothing holds a temperature contrast across it. The lower layer then stores no heat, instead it passes on the heat entering it plus the heat it generates, and the interface keeps its temperature. An `off` outermost layer therefore loses all of that heat through the surface, its temperature rate is zero, and its profile ends at its own temperature rather than at `surface_temperature`.

A world whose layers are all at one temperature and generate no heat has no profile to integrate. The solve then keeps its four structure variables and returns what it would with `solve_temperature=False`, at the same cost. The profile queries still report each layer's own temperature. Otherwise the solve adds temperature and heat flow as two more state variables and iterates. The first pass is isothermal. Each later pass integrates the profile and then relaxes the boundary layers, interface temperatures, and heat flows against the structure it produced. `thermal_passes` counts the passes and `thermal_converged` reports whether they settled: the largest relative change in the interface temperatures and heat flows between two passes fell below `thermal_tol`. `max_thermal_passes` caps the passes. Both are `[eos_solver]` settings that a call or the world may override. A solve that uses every pass without settling logs a warning and keeps the last pass. Layers that hold their mass move in the same passes and end on the grid of the last pass, so every reported slice and interface lies in its own layer.

Each layer's material is evaluated at the solved temperature and pressure of every radius, so an Arrhenius viscosity is stiff where the profile is cold, and a melting layer melts where its profile crosses its solidus. A layer with `use_thermal_expansion` also has its density follow the profile.

```python
from TidalPy.Radiogenics import make_radiogenics

world = TerrestrialWorld("Io-like", 1.8215e6, 8.93e22)
world.add_layer(Layer(
    "core", 0, 0.0, 8.1e5,
    material="simple_iron_core",
    temperature=1900.0))                           # Isothermal: no cooling model
world.add_layer(Layer(
    "mantle", 1, 8.1e5, 1.8215e6,
    material="simple_rock",
    temperature=1600.0,                            # [K] at the top of its adiabatic interior
    cooling="convection",
    radiogenics=make_radiogenics("isotope"),
    use_heating=True))                             # Its radiogenics heat it
result = world.solve_eos(
    surface_temperature=110.0)                     # [K] what the outermost layer radiates to

print(world.get_temperature(0.9 * world.radius))   # [K] on the solved profile
print(world.get_heat_flow(world.radius))           # [W] leaving the world
print(result["layer_nusselt_number"])              # The mantle convects: Nu of about 7
print(world.calc_layer_temperature_rate("mantle"))  # [K/s] from the mantle's heat imbalance
```

### Heat Sources

A layer with `use_heating` on is heated by the world's three heat sources during a solve that carries temperature. Each gives a volumetric heating $h$ \[W m$^{-3}$\], and the structure integrates

$$\frac{dL}{dr} = 4 \pi r^2 h$$

so the heat flow grows through a heated layer and its conducting stretches bend: a uniformly heated conducting shell follows $T = B + A/r - h r^2 / 6k$. The resistance chain includes the same heating. With $H(r)$ the heat generated between the base of a conducting stretch and $r$, the flow leaving its top is the flow entering plus $H$, and the temperature drop across it is $L_\mathrm{base} R + \int H / (4 \pi r^2 k) \, dr$. Both follow from the heating and the solved density, so the solved profile still passes through every layer temperature and reaches `surface_temperature`. A heated world is a thermal solve even when its layers share one temperature. With `solve_temperature=False` there is no heat flow, so the heating is ignored (a note in the debug log).

| Source | Heating | Set by |
|---|---|---|
| Radiogenic | The layer's [radiogenics model](../../Radiogenics/radiogenics_models.md) gives a specific rate $\epsilon$ \[W kg$^{-1}$\] at the solve's `time` \[s\], and $h = \epsilon \rho$ with the local density, so it is exact for a layer whose mass is an output of the solve. `time=None` uses each model's own reference time. | `Layer.radiogenics` |
| Tidal | The heating \[W\] the world's last `calc_tides` put in each layer. With a radial-solver Love method it follows the radial profile of the per-layer heating integral (piecewise linear in the layer's radial fraction, renormalized so the layer receives exactly its heating); otherwise it is spread by mass. | `calc_tides` |
| Prescribed | A power \[W\] spread over the layer by mass, or a specific rate \[W kg$^{-1}$\]. A power reaches the layer to the thermal tolerance whatever mass the solve gives it. | `set_prescribed_heating` |

The report gives each layer's heating from each source (`layer_heating_radiogenic`, `layer_heating_tidal`, `layer_heating_prescribed`) and their sum (`layer_heating`), all zero for a layer without `use_heating` and in a solve that carries no temperature.

| Member | Description |
|---|---|
| `set_prescribed_heating(layer, power=None, specific_rate=None)` | Prescribe a layer's heating by its name or index: a `power` \[W\] or a `specific_rate` \[W kg$^{-1}$\]. With neither it clears the layer's prescribed heating. Both, an infinite value, or a layer the world does not have raise `ValueError`. The EOS solve reads it, so the world forgets its solved structure. |
| `prescribed_heating` | The prescribed heating by layer name, as `{"power": watts}` or `{"specific_rate": watts_per_kg}`; layers without one are left out. |
| `tidal_heat_source` | The tidal source: each layer's heating \[W\] from the last `calc_tides`, by layer name; empty before one. |
| `clear_tidal_heating()` | Forget the tidal source, so later solves heat no layer tidally. |
| `get_heating(radius)` | The volumetric heating \[W m$^{-3}$\] the last solve's sources give at a radius (float or array), every source summed: zero in a layer without `use_heating`, NaN before a successful solve or outside the world. |
| `calc_layer_temperature_rate(layer)` | A layer's temperature rate \[K s$^{-1}$\] by its name or index: $(C + C_\mathrm{latent})\,dT/dt = L_\mathrm{in} - L_\mathrm{out} + H$ with the heat flows, radiogenic heat, and prescribed heat of the last solve, and the tidal heat of the latest `calc_tides`. NaN before a solve. |
| `calc_layer_thermal_capacity(layer)` | The heat $C$ \[J K$^{-1}$\] a layer's profile stores per kelvin of its temperature, by its name or index. Computed from the last solve on the first call after it and kept until the next; NaN before a solve. |
| `calc_layer_latent_capacity(layer)` | The latent heat $C_\mathrm{latent}$ \[J K$^{-1}$\] the layer's zone boundaries absorb per kelvin of its temperature; zero without such a boundary, NaN before a solve. |

The tidal source depends on the solved interior through the tides, so it is kept across solves: a later `solve_eos` keeps it, while a failed `calc_tides`, `clear_tidal_heating`, `add_layer`, and `load_binary` forget it. A layer with `use_tides` off dissipates nothing (the radial solver treats it as elastic), so it takes no tidal heating and its heat is in neither the total nor any other layer. A radial-solver tide with `layer_tidal_heating` off gives no per-layer split, so the source spreads the total over the tidal layers by mass. A step of an evolution is `solve_eos`, then `calc_tides`, then the temperature rates: `calc_layer_temperature_rate` takes the tides just computed, with no second solve, so the first step already has them and no tidal term is added by hand. Neither the tidal source nor the prescribed heating is saved with the world.

```python
import math

from TidalPy.Tides.classes import make_tide

world.set_tide_model(make_tide("rheology"))              # Love numbers from the radial solver
world.set_tide_config(
    max_degree_l=2,
    eccentricity_truncation=2,
    obliquity_truncation=0)
orbital_frequency = 2.0 * math.pi / (1.769 * 86400.0)    # [rad s-1]

world.solve_eos(surface_temperature=110.0)               # One evolution step: the structure first,
world.calc_tides(
    orbital_frequency=orbital_frequency,
    spin_frequency=orbital_frequency,
    eccentricity=0.0041,
    obliquity=0.0,
    semi_major_axis=4.217e8,
    host_mass=1.898e27)                                  # then the tides of that structure,
print(world.tidal_heat_source)                           # [W] by layer, about 8.3e11 in the mantle
print(world.calc_layer_temperature_rate("mantle"))       # then the rates [K/s], tides included

world.set_prescribed_heating(
    "mantle",
    specific_rate=1.0e-11)                               # [W/kg] beside the radiogenic and tidal heat
result = world.solve_eos(surface_temperature=110.0)
print(result["layer_heating_tidal"], result["layer_heating_prescribed"])
print(world.get_heating(1.5e6))                          # [W/m^3] every source summed
```

### Profile Queries

After a successful solve the world reads its profiles at any radius \[m\]. Each getter takes a float or an `np.ndarray` and returns a float or an array of the same shape, NaN where unsolved.

| Member | Returns | Description |
|--------|---------|-------------|
| `get_density(r)` | float [kg/m³] | Density at radius `r` [m]. |
| `get_gravity(r)` | float [m/s²] | Gravitational acceleration at `r`. |
| `get_pressure(r)` | float [Pa] | Pressure at `r`. |
| `get_temperature(r)` | float [K] | Temperature at `r` on the solved profile. |
| `get_heat_flow(r)` | float [W] | Heat flowing outward through the sphere of radius `r`. |
| `get_heating(r)` | float [W/m³] | Volumetric heating of the last solve's heat sources. |
| `get_state(r)` | dict | Density, gravity, pressure, the post-melt moduli and viscosities, and the melt fraction, from one evaluation of the solved state per radius. |
| `eos_solved` | bool | `True` once profiles are populated. |
| `all_materials_set` | bool | `True` once every layer has a material. |
| `surface_gravity_eos`, `central_pressure`, `planet_mass_eos`, `planet_moi_eos` | float | Scalar results of the last solve (NaN if unsolved). |
| `zones`, `molten_regions` | list | The solid and liquid zones of the last solve, and the liquid zones of layers not liquid throughout (see [Pieces and Zones](#pieces-and-zones)). |

The getters delegate to the layer that contains `r`, so each layer also exposes the same getters (see [Profile Getters](../layers/layer.md#profile-getters)).

### C++ API

The solve runs in C++ and can be called directly from C++:

```cpp
tidalpy::c_WorldEOSSolveConfig cfg;          // starts from the [eos_solver] section of the configuration
cfg.surface_pressure = 1.0e5;                // override only what this solve changes
world.solve_eos(cfg);                        // populates each layer's c_LayerEOSData
double rho = world.get_density(5.0e6);
```

`c_BaseWorld::solve_eos(const c_WorldEOSSolveConfig&)` lays out each layer's segments (its cooling model's profile in a thermal solve), estimates the bulk density, converts to the non-dimensional solve units, and calls `Material/eos`'s `c_solve_eos` (`solver_.hpp`: CyRK ODE integration inside the central-pressure iteration). Each layer is described by a `c_EOSLayerBounds` (its radii, the mass it holds, and the rigidity-margin event of a layer that can change state), and `c_EOSPassIntegrator` integrates one pass in pieces with CyRK terminal events: `c_eos_mass_event` (`ode_.hpp`) for an enclosed mass and `c_preeval_rigidity_margin` for a change of state. The material is read through `c_preeval_material` (`Material/eos/methods/material_preeval_.hpp`): the density alone while iterating, the thermal properties in a thermal solve, and the full state for dense output. The converged pass is committed (`c_commit_eos_pass`) into a `c_EOSSolution` (`eos_solution_.hpp`) that holds the pieces (`c_EOSPiece`), reads each piece only inside its own span, and builds the zones (`c_EOSZone`). The solve throws `std::invalid_argument` on bad input (`ValueError` in Cython via `except +`).

`c_BaseWorld` exposes `get_density/get_gravity/get_pressure(double) const` (delegating to the containing layer via `find_layer_for_radius`), the vectorized `get_eos_fields(field_indices, num_fields, radii, num_radii, values_out) const` (entries of the `eos_layout_.hpp` evaluation layout at every radius, field-major, under one hold of the call lock) and `calc_complex_moduli(is_shear, radii, num_radii, frequency, moduli_out) const`, the result accessors `get_eos_success/get_eos_message/get_eos_iterations/get_eos_pressure_error/` `get_surface_gravity_eos/get_surface_pressure_eos/get_central_pressure/` `get_planet_mass_eos/get_planet_moi_eos`, the zones (`get_zones()` as solved, `get_radial_zones()` under the layers' current flags, `get_molten_regions()`, `get_is_liquid_at(radius)`), the heat sources (`set_prescribed_heating`, `get_tidal_heat_source`, `clear_tidal_heating`, `get_heating`, `calc_layer_temperature_rate`, and the lazily computed, cached `calc_layer_thermal_capacity` and `calc_layer_latent_capacity`), the retained full solution via `get_eos_solution() -> const c_EOSSolution*`, and `get_eos_solved()` / `get_all_materials_set()`. The heat sources are in `heating_.hpp` (`c_Heating`, with `c_HeatSourceKind` `Radiogenic`, `Tidal`, and `Prescribed`), and the thermal network, with `c_zone_boundary_latent_capacity`, in `thermal_layout_.hpp`. The Cython `BaseWorld.solve_eos` wrapper converts the integration-method string to the CyRK enum, fills `c_WorldEOSSolveConfig`, calls the C++ method under `nogil`, and builds the Python result dict from the report.

## Viscoelastic Properties

After a successful `solve_eos`, the world and each layer expose the viscoelastic getters the radial solver uses.

| Member | Returns | Description |
|--------|---------|-------------|
| `get_shear_modulus(r)` | float [Pa] | Post-melt static shear modulus at `r`. |
| `get_bulk_modulus(r)` | float [Pa] | Post-melt adiabatic bulk modulus at `r`. |
| `get_shear_viscosity(r)` | float [Pa·s] | Post-melt shear viscosity at `r`. |
| `get_bulk_viscosity(r)` | float [Pa·s] | Post-melt bulk viscosity at `r`. |
| `calc_complex_shear_modulus(r, ω)` | complex [Pa] | Rheology-derived complex shear modulus at radius `r` [m], frequency `ω` [rad/s]. |
| `calc_complex_bulk_modulus(r, ω)` | complex [Pa] | Rheology-derived complex bulk modulus. |

All of the above accept a float or `np.ndarray` for `r` (and `ω`); array inputs return an `np.ndarray` of the same shape.

## Calculating Love Numbers

`BaseWorld.solve_love_numbers(...)` calculates the viscoelastic-gravitational tidal or loading Love numbers $k$, $h$, $l$ at a tidal forcing frequency.

`love_method` selects the method (names are case-insensitive; aliases in parentheses):

| `love_method` | What it does |
|---|---|
| `radial_solver` (`shooting`, `rs`; default) | Numerically integrates the radial ODEs from the center to the surface. Works for arbitrary multi-layer, solid/liquid, static/dynamic, compressible/incompressible worlds. |
| `propagation_matrix` (`prop_matrix`, `pm`, `prop`) | Quasi-analytic matrix propagation, restricted to a single solid, static, incompressible layer. An incompatible world fails the solve gracefully (`love_success` is `False`, `love_error_code` non-zero). `core_model` selects the core starting condition. |
| `homogeneous` (`homogen`) | Quasi-homogeneous: each tidal layer is treated as a homogeneous incompressible planet made of its own averaged material, and the world's Love numbers are the sum of the layers' weighted by their tidal scales (see Quasi-Homogeneous Love Numbers below). |
| `cpl` | The same, with each layer's static (unrelaxed) averaged shear modulus and a constant phase lag: $k$, $h$, $l$ are multiplied by $(1 - i/Q)$ so $-\mathrm{Im}[k] = \mathrm{Re}[k]/Q$. $Q$ is `fixed_q` (argument or `[tides]` config) or, when unset, the attached tide model's fixed Q for the degree. |
| `ctl` | Similar to `cpl` but with a constant time lag, $(1 - i\,\omega\,\Delta t)$; $\Delta t$ is `fixed_dt` or the tide model's fixed time lag. |
| `laterally_inhomogeneous` (`3d`, `lat_inhom`) | Reserved for a laterally inhomogeneous (3D) Love solver, which is not implemented; raises `NotImplementedError`. |

### Quasi-Homogeneous Love Numbers

For the `homogeneous`, `cpl`, and `ctl` methods, each layer uses the homogeneous-sphere formulas, $k_{l} = \frac{3}{2(l-1)}\,\frac{1}{1+\bar{\mu}_{l}}$, $h_{l} = \frac{2l+1}{2(l-1)}\,\frac{1}{1+\bar{\mu}_{l}}$, $l_{l} = \frac{3}{2l(l-1)}\,\frac{1}{1+\bar{\mu}_{l}}$ with $\bar{\mu}_{l} = \frac{2l^{2} + 4l + 3}{l}\,\frac{\mu}{\rho g R}$, where $\mu$ is the layer's rheology at the forcing frequency applied to its averaged moduli and viscosity, and $\rho$, $g$, $R$ are the planet's EOS bulk density, EOS surface gravity, and radius. This allows these methods to keep some of the layered structure without a radial solve (so it is much faster). Each tidal layer is averaged into a homogeneous material. Its post-melt shear and bulk moduli are volume averaged. Its post-melt viscosity is log-volume averaged, so a viscosity that spans decades across the layer is averaged in log space. A homogeneous planet made of that material, with the planet's radius, bulk density, and surface gravity, has Love numbers $k_i$. The layer's tidal scale $s_i$ weights them, and the world's Love numbers are

$$k = \sum_i s_i k_i$$

(and likewise $h$ and $l$). The tidal scale is the layer's `tidal_scale`, or its volume over the planet's when none is set. It is 0 for a layer with `use_tides` off, which then takes no part. With the default scales, a one-layer planet gets exactly the homogeneous-sphere value. A small, highly dissipative layer (an asthenosphere, for example) adds dissipation in proportion to its volume instead of making the whole planet dissipate like it. The tidal heating is linear in $-\mathrm{Im}[k]$, so each layer's heating is the heating of its own term $s_i k_i$, constant within the layer, and the layers sum to the total. `love_layer_parts` lists each layer's $s_i$, $k_i$, $h_i$, $l_i$, and complex shear modulus. `get_layer_tidal_scale(index)` returns the scale in use.

### Input Checks

Every Love solve checks its input first and raises `ValueError` for:
- Degree below 2 in a tidal or free-surface solve. Degree 1 is a translation of the body. A loading-only solve may use degree 1, whose load Love numbers depend on the reference frame.
- Frequency that is not finite and positive, or that lies outside the `[numerical]` `minimum_frequency` to `maximum_frequency` range. The Love numbers at $-\omega$ are the complex conjugates of those at $\omega$.
- `start_radius_tol` outside (0, 1).
- Negative `starting_radius` or `max_step`.

`max_step` \[m\] is converted into the integration's units. The `homogeneous`, `cpl`, and `ctl` methods take the bulk density from the solved mass (the structure the EOS surface gravity comes from), not from the declared mass.

The analytic methods have no radial functions (the radial-function getters, such as `get_love_radial_y`, return NaN) and no depth-resolved solution, so the 3D stress/strain/heating path (`calc_3d_tides`, `get_3d_tidal_heating`) raises `RuntimeError` when an analytic method is used. Free-function versions of the formulas live in `TidalPy.Tides.love` (`calc_homogeneous_love_numbers`, `calc_effective_rigidity`, `apply_fixed_q`, `apply_fixed_dt`; see [Love numbers](../../Tides/love/love_numbers.md)).

The world's default method is set with `set_tide_config(love_method=..., love_fixed_q=..., love_fixed_dt=...)` or the `[tides]` keys `love_method`, `love_fixed_q`, `love_fixed_dt_s` in a world TOML file. `solve_love_numbers` takes the method per call.

### Example

```python
from TidalPy.Material import Material, Phase

world = BaseWorld("planet", 6.0e6, 3.6e24)
world.add_layer(Layer(
    "mantle", 0, 0.0, 6.0e6,
    material=Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": 4000.0, "bulk_modulus_pa": 1.3e11},
        shear_modulus={"model": "constant", "shear_modulus_pa": 6.0e10},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e21})),
    shear_rheology="maxwell"))                   # A uniform Maxwell planet
world.solve_eos()                                # Required before any Love solve
result = world.solve_love_numbers(
    frequency=1.0e-5)                            # [rad s-1]

print(result["love_number_k"])                   # Complex k2, also world.love_number_k
print(world.love_number_h, world.love_number_l)
```

The moduli and the viscosity are properties of the layer's material, not of the rheology model. A rheology class only holds model parameters (the Andrade exponent, the Voigt fractions, etc.) and uses the modulus and viscosity as arguments. A solid layer with no shear modulus has no strength, and the solve fails.

**`solve_love_numbers( frequency=1e-5, degree_l=2, solve_for='tidal', core_model=0, use_kamata=None, nondimensionalize=None, starting_radius=0.0, start_radius_tol=None, integration_method=None, rtol=None, atol=None, scale_rtols=None, max_num_steps=None, expected_size=None, max_ram_MB=None, max_step=0.0, verbose=False, warnings=True, love_method=None, fixed_q=None, fixed_dt=None, raise_on_fail=False) -> dict`**

Every solver setting left as `None` takes the `[radial_solver]` value of the TidalPy configuration (see [Configurations](../../Overview/2_TidalPy_Configurations.md)), the same defaults the standalone `radial_solver` and the world's own tidal solves use. `love_method`, `fixed_q`, and `fixed_dt` left as `None` take the world's `[tides]` settings. `solve_eos` must be called first: raises `ValueError` if the EOS has not been solved. Returns a dict (`success`, `error_code`, `message`, `love_method`, `love_number_k/h/l`). The results are also stored on the world and read through the properties below. With `raise_on_fail=True` a failed solve raises `SolutionFailedError` (a `RuntimeError`) with its message.

### Frequency Sweeps

**`calc_love_numbers(frequencies, **solve_love_numbers_kwargs) -> dict`**

Runs `solve_love_numbers` at each frequency \[rad s-1\] with the same other arguments and returns `frequency` and, each shaped like it, `k`, `h`, `l` (complex), `success`, and `message`. A solve that fails gives NaN at its frequency and the sweep goes on, so only bad input (an unsolved EOS, a frequency outside the allowed range) raises.

```python
import numpy as np

periods = np.logspace(0, 3, 30) * 86400.0             # 1 to 1000 days [s]
sweep = world.calc_love_numbers(
    2.0 * np.pi / periods,
    degree_l=2)                                       # One Love solve per frequency
neg_imag_k2 = -sweep["k"].imag                        # NaN where a solve failed
```

The radial solver integrates each [zone](#pieces-and-zones) as a layer, and does not interpolate between EOS slices. At the exact radius the integrator requests:
- Gravity, pressure, mass, and moment of inertia come from the world's dense EOS solution.
- The density and the static moduli and viscosities come from the same solution (the layer's material evaluated them as the structure was integrated).
- The complex moduli come from the layer's rheology applied to those static values at that solve's frequency.

The Love numbers therefore converge with the integration tolerance alone and are entirely independent of `slices_per_layer`, whether or not the moduli and viscosities vary with depth. `slices_per_layer` still sets the size of the returned profile arrays. The propagation-matrix method propagates across those slices, so its Love numbers _do_ depend on it.

`release_radial_solution()` returns the world's last radial solve as a `RadialSolverSolution`, with its `result` grid filled on the solve's radius grid (`sample_radii()` gives those radii, so the two plot together). A radius above the surface has no solution and reads NaN.

### Love-Number Properties

These describe the world's last `solve_love_numbers` (or `solve_love_numbers_supplied`) call. The per-frequency Love solves of `calc_tides` and the 3D paths use their own workspaces. They leave these properties, and what `release_radial_solution` returns, untouched. After only a `calc_tides`, the properties still read as unsolved (`love_success` is `False` and `love_error_code` is -100).

| Property | Type | Description |
|----------|------|-------------|
| `love_solved` | bool | `True` while a successful solve is held. A later `solve_eos` clears it: Love numbers describe the structure they were solved with, so every Love number and radial-function getter returns NaN until the next solve. |
| `love_success` | bool | `True` if the last solve converged and still describes the structure. |
| `love_error_code` | int | Solver error code (0 = success; < 0 = failure; -100 = no Love solve has run since the world was built or last solved its EOS). |
| `love_message` | str | Human-readable solver message. |
| `love_num_ytypes` | int | Number of independent solution types (boundary-condition models requested). |
| `love_number_k`, `love_number_h`, `love_number_l` | complex | Love numbers for the first boundary condition at the solved degree. Equivalent to `get_love_number_k(0)` and friends. |
| `love_q_k`, `love_lag_k` | float | Quality factor $Q = -s\,\lvert k \rvert / \mathrm{Im}(k)$ and phase lag $\arctan_2(-s\,\mathrm{Im}(k), \lvert\mathrm{Re}(k)\rvert)$ \[rad\] of `love_number_k`, with $s$ the sign of $\mathrm{Re}(k)$: the definitions `RadialSolverSolution.Q_k` and `lag_k` use. $Q$ is infinite and the lag 0 for a purely elastic $k$; both are NaN when $k$ is. |
| `love_method` | str | Canonical name of the method the last solve used. |
| `love_surface_amplification` | float | Conditioning of the surface boundary-condition solve, recorded on every shooting solve whether or not `warnings` is on; near 1 is healthy, 0 after an analytic solve. |
| `love_surface_rcond` | float | Reciprocal condition number of the surface boundary-condition system, the rank measures if the solution constants are undetermined. |
| `love_effective_shear_modulus`, `love_tidal_volume` | complex, float | The tidal-scale-weighted mean of the layers' complex shear moduli \[Pa\] and the volume of the layers that took part \[m$^3$\] in the last quasi-homogeneous solve; NaN after a radial-solver solve. |
| `love_layer_parts` | list of dict | Each tidal layer's part of the last quasi-homogeneous solve: `layer`, `tidal_scale`, `love_number_k`, `love_number_h`, `love_number_l`, and `shear_modulus`. Empty after a radial-solver solve. |

**`get_love_number_k(ytype_idx=0) -> complex`**, **`get_love_number_h(ytype_idx=0) -> complex`**, **`get_love_number_l(ytype_idx=0) -> complex`**

Return the Love numbers for boundary-condition model index `ytype_idx` (0 = first requested, usually tidal).

**`get_love_surface_y(ytype_idx, y_idx) -> complex`**

Raw radial function y₁…y₆ at the surface for solution type `ytype_idx`, function index `y_idx` (0–5).

**`get_love_radial_y(radius, ytype_idx=0, y_idx=0) -> complex or array`**

The same radial function at any radius \[m\] (a float, or an array giving a complex array of its shape), evaluated from the solver's dense interpolants, so it is accurate between grid slices. Returns NaN if the solve failed, if an analytic Love method was used (those have no radial functions), or if the radius sits below the solver's starting radius.

### Love-Number C++ API

```cpp
tidalpy::c_LoveSolveConfig cfg = world.make_love_solve_config();   // [radial_solver] defaults + the world's [tides]
cfg.frequency = 1.0e-5;                                            // [rad/s]
cfg.degree_l  = 2;
world.solve_love_numbers(cfg);   // delegates to the cached radial solver

std::complex<double> k2 = world.get_love_number_k(0);
```

A default-constructed `c_LoveSolveConfig` (and `c_WorldEOSSolveConfig`) reads the `[radial_solver]` (`[eos_solver]`) section of the shared runtime config, so C++ callers and the tide paths start from the same defaults as Python callers.

`c_BaseWorld::solve_love_numbers(const c_LoveSolveConfig&)` delegates to a cached helper, `c_WorldRadialSolver` (held by `p_love.radial_solver`). The helper separates the frequency-independent setup (built once and reused) from the frequency-dependent work (recomputed on every call), which speeds up frequency sweeps and orbital evolution.

The non-dimensionalization is frequency-independent (the `c_NonDimensionalScales` time scale is $1/\sqrt{\pi G \bar{\rho}}$ for the bulk density $\bar{\rho}$, not $1/\omega$). Only the complex moduli and the shooting integration change between calls at different frequencies.

1. Validates `eos_solved` and `tidalpy_config_ptr`.
2. If the cache does not match the current EOS grid/assumptions, `build_cache` captures (once), for the solver's layers (the world's zones, `get_radial_zones()`; the providers map each solver layer back to its world layer): the non-dim radius/density/gravity/pressure/mass/moi arrays, per-layer metadata (solid/liquid, static, incompressible) and slice partitioning, the non-dim scalars (`G`, bulk density, surface pressure), and a reused `c_RadialSolutionStorage` whose internal `c_EOSSolution` arrays serve as the scratch buffers. The cache is invalidated automatically whenever `solve_eos` re-runs.
3. Per call: the world installs a material-state provider (`c_EOSSolution::MaterialEval`, a type-erased callable carrying that solve's frequency) that the shooting solve calls for density and the complex moduli at each integration radius, and also fills the dimensional moduli scratch at the slice radii for the propagation-matrix method and the array outputs. The helper non-dimensionalizes the scratch in place, re-applies the cached non-dimensional structure arrays, and runs the selected solver.
4. Re-dimensionalizes the y-solution, restores the SI surface gravity, and calls `c_RadialSolutionStorage::find_love()`.

A Love solve writes everything it produces (the cached radial solver and its storage, or the quasi-homogeneous results) into a `c_LoveWorkspace`. The world holds one for `solve_love_numbers`, which `get_love_number_k/h/l`, `get_love_surface_y`, and the status accessors read.

## Global (1D) Tidal Dissipation

`BaseWorld.calc_tides(...)` calculates the body's total tidal heating and three orbital potential partial derivatives by summing over the active tidal modes (the global or "1D potential" approach). A tide model (see [Global Tidal Dissipation](../../Tides/global_tides.md)) supplies the per-mode dissipation multiplier $-\mathrm{Im}[k_{l}]$. The world runs the global-potential engine for its stored `[tides]` config and the supplied orbital and spin state, collapses, and resolves each layer's share of the heating. It also records each layer's heating as the world's [tidal heat source](#heat-sources).

`calc_tides` returns a dict: `tidal_heating` \[W\], the potential derivatives `dU_dM`, `dU_dw`, and `dU_dO` \[J kg$^{-1}$ rad$^{-1}$\], `num_tidal_modes`, and `layer_tidal_heating` (each layer's heating \[W\] by name). The same results stay on the world for the getters below until the next solve.

### Orbital State

The six orbital-state arguments (`orbital_frequency`, `spin_frequency`, `eccentricity`, `obliquity`, `semi_major_axis`, `host_mass`) keep their order, and each one left as `None` comes from `get_tide_state()`, the state the system the world belongs to gives it: its orbit about its tidal host, the host's mass, and its own spin and obliquity. A world of a `System` is solved in its current state by `world.calc_tides()`, and `world.calc_tides(eccentricity=0.05)` changes one value. A world outside a system, or without a tidal host, has no state to take, so an argument left out raises `ValueError` naming it. The 3D methods (`get_3d_tidal_heating`, `get_3d_tidal_heating_array`, `calc_3d_displacements`, `calc_3d_stress_strain`, `calc_3d_tides`) take the state the same way; their point and grid arguments are then given by keyword.

A tide solve warns, once per world, when the spin is within $10^{-3}$ of the orbital mean motion but not equal to it. A slightly non-synchronous spin adds a slow forcing term at the difference frequency, which can change the heating by orders of magnitude (a bundled Io whose spin is 0.016 percent off the mean motion of its rounded semi-major axis heats about four times more than a synchronous one). A synchronous world should spin at exactly its orbital frequency (`world.set_spin_frequency(system.calc_orbital_frequency(world))`). The check sits on the C++ path every tide solve takes, so a `System` evolution warns too.

### Tidal Heating of Each Layer

A layer's heating depends on the source of the Love numbers:

| Love numbers from | Layer heating |
|-------------------|---------------|
| The radial solver (`radial_solver`, `propagation_matrix`) | The volume integral of the radial solution's orbit-averaged heating density over the layer: the same integral as `calc_3d_tides` with every axis summed, scaled so the tidal layers sum to the 1D total. A liquid zone carries no shear dissipation, and a layer with `use_tides` off, solved with purely real moduli, takes none. It costs about as much as the global solve again, so `set_tide_config(layer_tidal_heating=False)` (or `layer_tidal_heating = false` in `[tides]`) skips it and leaves every layer's heating NaN. |
| The quasi-homogeneous methods (`homogeneous`, `cpl`, `ctl`) | The heating of the layer's own term $s_i k_i$ (see Quasi-Homogeneous Love Numbers), constant within the layer; the layers sum to the total. |
| An analytic tide model (`cpl`, `ctl`, `ctl_q` tide models, which describe the whole body) | The total times the layer's tidal scale $s_i$. |

The world builder normally sets the tide model and configuration from the `[tides]` TOML table, with per-family defaults (star: `fixed_q`, gasgiant: `fixed_dt`, terrestrial: `rheology`). The table may hold the per-degree lists of every analytic model, and the builder passes the tide model only the ones it reads (`tide_config_keys`): a `fixed_dt` world ignores `fixed_q`. They can also be set directly:

```python
world.set_tide_model(make_tide(
    "cpl",
    fixed_k=[0.3],
    fixed_q=[50.0]))                     # An analytic tide model
world.set_tide_config(
    max_degree_l=2,
    eccentricity_truncation=2,
    obliquity_truncation=0)

result = world.calc_tides(
    orbital_frequency=2.05e-5,
    spin_frequency=2.05e-5,
    eccentricity=0.0041,
    obliquity=0.0,
    semi_major_axis=4.2e8,
    host_mass=1.898e27)                  # Tides of an Io-like orbit about Jupiter

result["tidal_heating"]                  # Total heating [W], also world.get_tidal_heating()
world.get_tidal_potential_derivatives()  # {"dU_dM": ..., "dU_dw": ..., "dU_dO": ...} [J kg-1 rad-1]
world.get_layer_tidal_heating("mantle")  # = world heating × the mantle's tidal scale (an analytic tide model)
```

`set_tide_model` also takes a model name (`world.set_tide_model("fixed_q")`, the model at its `[tides]` defaults) or a table like a world file's `[tides]` table, with the model under `model` or `global_tidal_model`, its per-degree lists, and any tide settings, which it applies too. `world.tide_model` returns the attached model (or `None`), and setting it calls `set_tide_model`.

`set_tide_config` takes every key `get_tide_config` returns (`eccentricity_trunc_lvl`, `obliquity_trunc_lvl`, `love_fixed_dt_s`, and the rest) as well as its own argument names, so `world.set_tide_config(**world.get_tide_config())` restores a configuration; a call gives each setting under one name only. `temporary_tide_config(**settings)` applies settings for a `with` block and restores the previous configuration on exit, also when the block raises:

```python
with world.temporary_tide_config(eccentricity_truncation=10):
    high_order = world.calc_tides(
        orbital_frequency=2.05e-5,
        spin_frequency=2.05e-5,
        eccentricity=0.0041,
        obliquity=0.0,
        semi_major_axis=4.2e8,
        host_mass=1.898e27)              # Level 10 inside the block, level 2 again after it
```

A new tide configuration clears the last tidal result, so read the result inside the block.

For a synchronous, low-eccentricity body the `cpl` result reproduces the standard CPL rate $\frac{21}{2}\,\frac{k_{2}}{Q}\,\frac{G M_{h}^{2} R^{5} n e^{2}}{a^{6}}$ (host mass $M_{h}$).

### Rheology Model

The analytic models (`cpl`/`ctl`/`ctl_q`) take $-\mathrm{Im}[k_{l}]$ from their fixed per-degree parameters and need no interior solution. The `rheology` model derives $-\mathrm{Im}[k_{l}(\omega)]$ from the world radial solver. `calc_tides` runs the global-potential engine and then, for each unique tidal frequency, solves the world's complex Love numbers (reusing the frequency-independent radial-solver cache). It feeds the per-mode `k_l` into the collapse and keeps the full `k`/`h`/`l` set per mode for inspection.

```python
world.solve_eos(temperature=1500.0)              # Required for the rheology model
world.set_tide_model(make_tide("rheology"))
world.set_tide_config(
    max_degree_l=2,
    eccentricity_truncation=2,
    obliquity_truncation=0)
world.calc_tides(
    orbital_frequency=2.0e-5,
    spin_frequency=1.0e-5,
    eccentricity=0.01,
    obliquity=0.0,
    semi_major_axis=4.0e8,
    host_mass=1.9e27)                            # A non-synchronous spin

world.get_tidal_heating()                        # Total heating [W]
world.get_tidal_love_k(2, 2, 0, 0)               # Complex k for the (l,m,p,q) = (2,2,0,0) mode
```

`calc_tides` raises `RuntimeError` if the `rheology` model is selected but the EOS has not been solved, or if a per-frequency radial solve fails.

| Member | Returns | Description |
|--------|---------|-------------|
| `set_tide_model(tide)` | - | Attach a tide model: a `TideBase`, a model name, or a `[tides]`-style table (see above); `None` detaches it. The world shares a model rather than copying it (models are not changed in place), so one model can serve several worlds. |
| `tide_model` | `TideBase` or None | The attached model; settable, through `set_tide_model`. |
| `tide_model_set` | bool | Whether a model is attached. |
| `set_tide_config(min_degree_l=None, max_degree_l=None, eccentricity_truncation=None, obliquity_truncation=None, layer_tidal_heating=None, eccentricity_exact_tolerance=None, love_method=None, love_fixed_q=None, love_fixed_dt=None, *, eccentricity_trunc_lvl=None, obliquity_trunc_lvl=None, love_fixed_dt_s=None)` | - | Change the stored `[tides]` truncation/degree settings (`eccentricity_truncation` a level or `"exact"`, with `eccentricity_exact_tolerance` its mode range; see [Eccentricity Functions](../../Tides/Eccentricity.md); `obliquity_truncation` 0 or `"off"`, 2, 4, or `"gen"`; see [Obliquity Functions](../../Tides/Obliquity.md)), whether `calc_tides` resolves each layer's heating on the radial-solver path, and the world's default Love-number method (see [Calculating Love Numbers](#calculating-love-numbers)). Only the arguments given change; the rest keep their current values (`get_tide_config()`). The three keyword-only names are the `get_tide_config` keys of the same settings; giving both names of one setting raises `ValueError`. A NaN `love_fixed_q` or `love_fixed_dt` clears it. |
| `get_tide_config()` | dict | The stored settings under the builder's `[tides]` key names (`*_trunc_lvl`). |
| `temporary_tide_config(**settings)` | context manager | Applies `set_tide_config(**settings)` for a `with` block and restores the previous configuration on exit. |
| `calc_tides(orbital_frequency=None, spin_frequency=None, eccentricity=None, obliquity=None, semi_major_axis=None, host_mass=None)` | dict | Run the global tidal solve; each argument left as `None` comes from `get_tide_state()` (see [Orbital State](#orbital-state)). |
| `tides_solved` | bool | Whether a successful solve is held. A new tide model or tide configuration clears it, and so does `solve_eos`, after which the heating and potential-derivative getters return NaN and each layer's `get_tidal_heating()` does too, until the next `calc_tides`. The tidal heat source is kept (see [Heat Sources](#heat-sources)). |
| `get_tidal_heating()` | float [W] | Total global tidal heating (NaN if unsolved). |
| `get_tidal_heat_flux()` | float [W m⁻²] | The heating over the surface area, $\dot{E} / (4 \pi R^2)$ (NaN if unsolved). |
| `get_tidal_potential_derivatives()` | dict | `dU_dM`, `dU_dw`, `dU_dO` [J kg⁻¹ rad⁻¹], by the mean anomaly, the argument of pericenter, and the longitude of the node. |
| `get_tidal_dU_dM_minus_dw()` | float | The per-mode sum of `dUdM - dUdw` [J kg⁻¹ rad⁻¹], to pass as `dU_dM_minus_dw` to `OrbitSolver.calc_de_dt`: at small eccentricity the two separate sums nearly cancel in $\dot{e}$. NaN if unsolved. |
| `get_num_tidal_modes()` | int | Active modes summed. |
| `get_layer_tidal_heating(layer=None)` | float [W] or dict | The heating the last `calc_tides` put in a layer, given by index or name (see Tidal Heating of Each Layer); NaN before one. With no layer, every layer's heating by name. |
| `get_layer_tidal_scale(layer)` | float | The tidal scale of a layer given by index or name: its `tidal_scale`, or its volume over the planet's when none is set; 0 with `use_tides` off. |
| `get_tidal_love_k(l, m, p, q)` | complex | Per-mode radial-solver `k_l` (rheology only; NaN for analytic models). |

Each layer also stores its own tidal heating: `layer.get_tidal_heating()` (C++ `c_Layer::get_tidal_heating()`) returns the same value as `get_layer_tidal_heating`.

## Other World Members

**Radial queries.** Each takes a radius [m], scalar or array, and returns the value there from the solver's dense interpolants, not by re-interpolating the gridded arrays. The `calc_*` forms solve the EOS first when it is unsolved or `force_recalc=True`.

| Member | Returns |
|---|---|
| `calc_density`, `calc_gravity`, `calc_pressure` | Density [kg m-3], gravity [m s-2], and pressure [Pa]. |
| `calc_shear_modulus`, `calc_bulk_modulus` | Static moduli [Pa], after any melt weakening. |
| `calc_shear_viscosity`, `calc_bulk_viscosity` | Viscosities [Pa s], after any melt weakening. |
| `calc_static_viscoelastics`, `get_static_viscoelastics` | The static (unrelaxed) moduli and viscosities together. |
| `get_melt_fraction` | Melt fraction of the material; `0.0` where it does not melt. |

**Geometry.** `calc_surface_area(radius)`, `calc_volume_sphere(radius)`, and `calc_volume_shell(outer, inner)` are the shared spherical helpers every structure inherits.

**Spin and orbit.** `set_spin_model(spin)` attaches a spin model. `get_moment_of_inertia()` returns the EOS-solved moment of inertia [kg m2], or the spin model's estimate `moment_of_inertia_factor` $\times\,M R^{2}$ before a solve. `calc_spin_derivative(host_mass)` returns the spin rate of change [rad s-2] from the current tidal solution. `calc_synchronous_spin(orbital_frequency)` returns the synchronous rate [rad s-1]. See [Dynamics](../../Dynamics/dynamics.md).

**State.** `get_state(radius)` returns every EOS profile at a radius as a dict (see [Profile Queries](#profile-queries)), and `calc_state(radius, force_recalc=False)` solves the EOS first when needed.

**Three-dimensional tides.** `get_3d_tidal_heating_array(...)` is the vectorized form of `get_3d_tidal_heating`, for building a heating map; one of its `radii` and `colatitudes` may be a scalar, which pairs with every value of the other. `calc_3d_displacements(...)` returns the instantaneous displacement grid, and `calc_3d_stress_strain(...)` returns the stress and strain grids. All three, like `calc_3d_tides`, take `num_threads` (default 0, the logical processors minus 4) to spread their radial solves and per-point work over threads. See the [3D heating page](../../Tides/multilayer_3d_heating.md).

**Configuration and identity.** `source_config` is the normalized configuration the world was built from, if any. `portable_config` is the configuration as given for a world built from a `data_file` (what `save_to_toml` writes while the world is unchanged since its build), and `built_config` the world's `get_config_dict()` at the end of its build. `get_save_config(destination_dir=None)` returns what `save_to_toml` writes. `family_world_type()` returns the builder's world type for this class. `get_schema_version_str()` returns the schema version the class writes. See the [TOML schema](../config/toml_schema.md).

**Pinned solver settings.** `set_solver_defaults(eos_solver=None, radial_solver=None)` pins keys of the `[eos_solver]` and `[radial_solver]` configuration sections on this world, and `get_solver_defaults()` returns them. A world file's tables of the same names are applied here. A call's own argument wins over a pinned key, and a pinned key wins over the TidalPy configuration for every solve the world runs. An unpinned key follows the configuration. `get_config_dict()` includes the tables.

## `TerrestrialWorld`

A `BaseWorld` whose `world_type` defaults to `"terrestrial"`, with its own binary class id. It has the same API as `BaseWorld`; rocky and icy planets and moons, the bundled ones included, are built as terrestrial worlds.

## `GasGiantWorld`

A `BaseWorld` whose `world_type` defaults to `"gasgiant"`, with its own binary class id. It has the same API as `BaseWorld` and its layers are usually fluids: layers of a liquid-only material (MatPack has `simple_gas`, `h2_he_molecular`, and `h_he_metallic`).

```python
from TidalPy.Structures.worlds import GasGiantWorld

jupiter = GasGiantWorld("Jupiter", 7.0e7, 1.898e27)
jupiter.add_layer(Layer(
    "envelope", 0, 0.0, 7.0e7,
    material="simple_gas"))                     # A liquid-only material: a liquid layer
print(jupiter.envelope.is_liquid)               # True
```

## `StarWorld`

A star: a `BaseWorld` that adds an effective temperature and luminosity, kept consistent through the Stefan-Boltzmann law, $L = 4\pi R^{2}\sigma T^{4}$. A star usually has no layers. It may hold them and solve its EOS like any world, and the TOML builder accepts `[layers.<name>]` tables for a star, but it needs neither.

```python
from TidalPy.Structures.worlds import StarWorld

sun = StarWorld("Sun", 6.957e8, 1.989e30, effective_temperature=5772.0)
sun.luminosity                # ~3.83e26 W (derived from T if luminosity == 0)
sun.set_luminosity(3.828e26)  # recomputes effective_temperature
sun.effective_temperature
```

**Properties:** `effective_temperature` \[K\], `luminosity` \[W\]. **Methods:** `calc_luminosity_from_temperature(T)`, `calc_temperature_from_luminosity(L)`, `set_effective_temperature(T)`, `set_luminosity(L)`, and `calc_insolation_flux(distance, eccentricity=0.0)`, the orbit-averaged flux \[W m$^{-2}$\] $L / (4 \pi a^2 \sqrt{1 - e^2})$ at a body with semi-major axis `distance` \[m\] (floats or arrays), the flux a `System` gives its worlds. A luminosity-model hierarchy (fixed, mass-to-luminosity, power law) can be attached via `set_luminosity_model` (see [Luminosity](../../Stellar/luminosity.md)).

**Tides.** A star without layers uses the analytic tide models (`cpl`, `ctl`, `ctl_q`); the `rheology` model needs layers and a solved EOS, so `calc_tides` raises if it is selected on a star that has neither. See [Global (1D) Tidal Dissipation](#global-1d-tidal-dissipation).

**Spin.** A star carries a spin model like every world, so a `System` evolves its spin rate. Without a solved EOS its moment of inertia is `moment_of_inertia_factor` $\times\,M R^{2}$. The builder takes the factor from `[worlds.star]` in the TidalPy configuration, 0.0754 (an $n = 3$ polytrope, a Sun-like star, as the `[tides.star]` Love numbers assume); a fully convective M dwarf is closer to 0.205 ($n = 1.5$). Set `moment_of_inertia_factor` in the star's TOML, or call `set_spin_model`, to change it.

## Binary Serialization

A world's binary file ([Binary Serialization](../../Utilities/binary.md)) holds the whole world structure, every layer with its material and models, and every setting that affects calculations.

- Every world type: the tide model and its configuration (`get_tide_config()`: degrees, truncations, Love method, `love_fixed_q`, `love_fixed_dt`, `layer_tidal_heating`), the spin model's `moment_of_inertia_factor`, the pinned solver settings (`get_solver_defaults()`), and every layer.
- `StarWorld`: the luminosity model.

Solved state is not saved. The EOS profile (`c_LayerEOSData`), the zones, Love numbers, and tide results must be recomputed with `solve_eos` and `calc_tides` after a load. The tidal heat source is a runtime input, like the orbit, and is not saved either: a load clears it. The prescribed heating is saved: the binary record holds each layer's, and `get_config_dict` writes it in a `[prescribed_heating]` table. Loading into an existing world replaces all of these, so a record saved without a tide model leaves the world without one.

The layers' materials and models are restored, so the EOS solve needs nothing re-attached. Layer views taken before the load (`world.<name>`, `get_layer`) refer to replaced layers, so take new ones. A layer that belongs to a world cannot be loaded in place. Load the world, or a standalone layer. A file of another class is refused with both classes named ("it is a GasGiantWorld file, not a TerrestrialWorld one"). Paths may be a `str` or an `os.PathLike` such as a `pathlib.Path`.

`TidalPy.Structures.load_world(path, force=False)` reads a file into a new world of the class that saved it, with no placeholder object, and `build_world(path)` does the same for a binary file. `load_world` raises `IOError` naming what a file holds when it is not a world's.

`copy()` (and `copy.copy`, `copy.deepcopy`, and `pickle`) goes through the same record held in memory, so a copy holds exactly what a saved and loaded world holds, and it is unsolved and belongs to no system. The configurations the world was built from (`source_config`, `portable_config`, `built_config`) are copied with it. `TidalPy.Structures.worlds.world_from_bytes(world_class, record, configs=None, force=False)` rebuilds a world from such a record.

```python
import pickle
from TidalPy.Structures import load_world

world.save_binary("world.tpyb")                       # A str or an os.PathLike path
loaded = load_world("world.tpyb")                     # A new world of the saved class
twin = world.copy()                                   # Independent: same class, layers, models, and settings
twin.solve_eos()                                      # Copies are unsolved
restored = pickle.loads(pickle.dumps(world))          # A pickle round trip, as a process pool makes
```

## References

- Schubert, G., Turcotte, D. L., and Olson, P. (2001). *Mantle Convection in the Earth and Planets*. Cambridge University Press. The upper-mantle temperature of parameterized convection.
- Solomatov, V. S. (2000). Fluid dynamics of a terrestrial magma ocean. In *Origin of the Earth and Moon*, 323-338. University of Arizona Press. The liquid convection scaling.
- Stevenson, D. J., Spohn, T., and Schubert, G. (1983). Magnetism and thermal evolution of the terrestrial planets. *Icarus*, 54(3), 466-489. Parameterized convection.
