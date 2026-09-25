# Worlds (`structures_x.worlds`)

_Updated: 2026-09-24_

The world classes are the top-level structural objects in TidalPy. A world owns its identity, orbital and thermal scalars, and bulk geometry; a layered world also owns an ordered stack of [layers](../layers/base_layer.md) and runs the whole-planet equation-of-state and radial (Love number) solves.

## Inheritance

```
TidalPyBaseClass
  └── StructureBase
        └── BaseWorld
              ├── LayeredWorld
              │     └── GasGiantWorld
              └── StarWorld
```

| Class | Layers? | EOS? | Purpose |
|-------|---------|------|---------|
| `BaseWorld` | no | no | Identity, albedo/emissivity/obliquity/spin, bulk geometry, equilibrium temperature. |
| `StarWorld` | no | no | A star; effective temperature and luminosity linked by the Stefan-Boltzmann law. |
| `GasGiantWorld` | yes | yes | Gas giant; a `LayeredWorld` with its own type/binary id. |
| `LayeredWorld` | yes | yes | Terrestrial/layered body; owns layers, aggregates mass and heating, runs the whole-planet EOS solve. |

## `BaseWorld`

```python
from TidalPy.structures_x.worlds import BaseWorld

welcome_to_earth = BaseWorld(
    name="Earth", radius=6.371e6, mass=5.972e24,
    world_type="terrestrial", albedo=0.3, emissivity=1.0,
    obliquity=0.41, spin_frequency=7.29e-5,
)
```

**Read-only properties:** `name`, `world_type`, `radius`, `mass`, `albedo`, `emissivity`, `obliquity`, `spin_frequency`.

**Methods**

| Method | Returns | Description |
|--------|---------|-------------|
| `calc_surface_gravity()` | float [m/s²] | $GM/R^{2}$. |
| `calc_escape_velocity()` | float [m/s] | $\sqrt{2GM/R}$. |
| `calc_mean_density()` | float [kg/m³] | $M / (\tfrac{4}{3}\pi R^{3})$. |
| `calc_equilibrium_temperature(F)` | float [K] | $\left[(1-A)\,F/(4\varepsilon\sigma)\right]^{1/4}$ (fast rotator, $F$ = insolation flux [W/m²]). |
| `set_spin_frequency(ω)` | - | Set rotation rate [rad/s]. A `System` reads it when it builds the world's tidal state; `calc_tides` takes its spin rate as an argument, so a tidal result already solved is left alone. |
| `set_obliquity(θ)` | - | Set axial obliquity [rad]. The same holds as for the spin rate. |

`get_config_dict()` returns the world as the TOML builder's world table: `schema_version`, `name`, `type` (the builder's world type, from `get_builder_world_type()`), `radius`, `mass`, `albedo`, `emissivity`, `obliquity`, `spin_frequency`, and a `tides` table when a tide model is attached (`global_tidal_model`, its per-degree parameters, and the settings from `get_tide_config()`). `save_config` / `save_binary` / `load_binary` are inherited from `TidalPyBaseClass`; `save_to_toml` validates the dict against the schema before writing when no build configuration is retained.

Binary class id 200 (`BinaryClassID::BaseWorld`).

## `LayeredWorld`

A world built from an ordered (inner-to-outer) stack of layers.

```python
from TidalPy.structures_x.worlds import LayeredWorld
from TidalPy.structures_x.layers import SolidLiquidLayer

world = LayeredWorld("Earth", 6.371e6, 5.972e24, world_type="terrestrial")
world.add_layer(SolidLiquidLayer("core",   0, 0.0,     3.485e6, 1.932e24))
world.add_layer(SolidLiquidLayer("mantle", 1, 3.485e6, 6.371e6, 4.040e24))
```

**Layer management**

| Member | Description |
|--------|-------------|
| `add_layer(layer)` | Add a layer inner-to-outer. Ownership of the layer (and its attached physics models) transfers into the world; the passed wrapper stays usable as a non-owning view of that layer, like the one `world.<layer name>` returns. Raises `ValueError` if the layer was already added or if its inner radius is not continuous with the current outermost radius (innermost must start at 0). A rejected layer is not consumed. |
| `num_layers` | Number of layers (property). |
| `calc_total_mass()` | Sum of the layer masses [kg]; equals `planet_mass_eos` after a successful EOS solve. |
| `calc_internal_heating(time)` | Sum of the layer radiogenic heating [W]; only `SolidLiquidLayer`s with an attached radiogenics model contribute. Uses each layer's `mass`, so solve the EOS first when the layers were built without one. |
| `validate_layers()` | `True` if every boundary is continuous and the innermost starts at 0. |

**Accessing layers**

A built world owns its layers; Python reaches them through non-owning views:

| Access | Returns |
|--------|---------|
| `world.get_layer(i)` / `world[i]` | the layer at index `i` (0 = innermost; negative indices allowed) |
| `world[a:b]` | a list of the sliced layer wrappers |
| `world.layers` | a list of all layer wrappers, inner to outer |
| `world.<layer_name>` | the layer with that name (e.g. `world.mantle`) |
| `for layer in world:` | iterate the layers inner-to-outer; `len(world)` is the layer count |

Each returns a non-owning view dispatched to the matching subclass (`PhysicsLayer`/`SolidLiquidLayer`/`GasLayer`/`BaseLayer`), so the layer's full API is available (`world.mantle.shear_modulus_static`, `world.mantle.get_tidal_heating()`, `world.core.calc_complex_shear_modulus(r, ω)`, ...). The world still owns the C++ layer; the view keeps the world alive, so it is safe to hold but must not be mutated through. The views are built once and cached (rebuilt only when a layer is added), so repeated access returns the same object (`world.mantle is world.mantle`). Attribute access by name is only consulted after normal attribute lookup (defined members win) and ignores names starting with `_`.

`get_config_dict()` adds a `layers` table keyed by layer name to the `BaseWorld` keys. Each entry is the layer's own config dict (`class`, scalars, attached-model sub-tables) without the standalone-only keys the builder derives itself, so `build_world(world.get_config_dict())` rebuilds the same structure.

Binary class id 201 (`BinaryClassID::LayeredWorld`). See [Binary serialization](#binary-serialization) below.

### Equation of State

Each layer carries a [material EOS model](../../material_x/material_eos.md) (its density source), attached with `BaseLayer.set_eos(model)`. Once every layer has one, `LayeredWorld.solve_eos(...)` integrates the planet's radial structure from center to surface, populates every layer's density, gravity, and pressure profile, and sets each layer's `mass` (and so its `density_bulk`) to the mass the solved profile places between the layer's radii. Until then a layer's mass is whatever it was constructed with; the TOML builder uses 0.0 when a file gives none, which is the usual case.

```python
from TidalPy.structures_x.worlds import LayeredWorld
from TidalPy.structures_x.layers import BaseLayer
from TidalPy.Material_x.eos import make_material_eos

world  = LayeredWorld("Earth", 6.371e6, 5.972e24)
core   = BaseLayer("core", 0, 0.0, 3.485e6, 0.0)
mantle = BaseLayer("mantle", 1, 3.485e6, 6.371e6, 0.0)
core.set_eos(make_material_eos("constant", {"reference_density_kg_m3": 11000.0}))
mantle.set_eos(make_material_eos("constant", {"reference_density_kg_m3": 4500.0}))
world.add_layer(core)
world.add_layer(mantle)

result = world.solve_eos(surface_pressure=0.0)  # dict of profile arrays + scalars
rho    = world.get_density(5.0e6)               # [kg/m³] at radius 5000 km
g      = world.get_gravity(world.radius)        # surface gravity [m/s²]
p0     = world.get_pressure(0.0)                # central pressure [Pa]
```

The solver carries pressure as a radial state variable, so analytic density-from-pressure models (Birch-Murnaghan, Vinet) are evaluated inline; the constant and interpolated models ignore pressure. The central pressure is found by a secant iteration on the surface-pressure mismatch: the first step assumes a unit slope (exact for an incompressible planet) and later steps use the slope measured between iterations, so a compressible planet converges in a few steps. The iteration is safeguarded for large, soft planets, whose surface pressure can at first fall as the central pressure rises: while the measured slope is not positive the step doubles each pass, and once the mismatch has changed sign the root is bracketed and found by secant, false-position, and bisection steps (as in Brent's method). A body whose materials cannot hold it up has no root; the solve then reports "no hydrostatic structure" and leaves the world unsolved. A world that has been solved before starts the iteration from its last central pressure instead of from a uniform sphere, so a re-solve after a small change (a new layer temperature, a later time) usually converges on its first pass. The integration runs in non-dimensional units (the planet radius, its bulk density, and $1/\sqrt{\pi G \rho}$ as the length, density, and time units) so the tolerances mean the same thing for every planet; every result is returned in SI.

**`solve_eos(surface_pressure=0.0, slices_per_layer=None, G_to_use=-1.0, integration_method=None, rtol=None, atol=None, pressure_tol=None, max_iters=None, nondimensionalize=None, temperature=None, solve_temperature=None, surface_temperature=None, reset_layer_masses=False, verbose=False, time=None) -> dict`**

Every solver setting left as `None` takes the `[eos_solver]` value of the TidalPy configuration (see [Configurations](../../Overview/2_TidalPy_Configurations.md)), the same defaults the standalone `radial_solver` uses. `pressure_tol` is relative to the central-pressure scale $(2/3) \pi G \rho^2 R^2$ and must stay above `rtol`, the integrator's own noise on the surface pressure. Hitting `max_iters` logs a warning, sets `max_iters_hit` in the result, and keeps the last iteration's profile.

#### Solved State

A re-solve starts its central-pressure iteration from the world's last converged solve, so a small change settles in a pass or two. The answer then depends on the call history at the level of the tolerances: a world solved twice at different surface pressures and back can differ from a fresh one by a few parts in $10^9$ in its central pressure and about one part in $10^9$ in k2. Build a fresh world when results must reproduce bit for bit.

A solve works on a copy of every layer's material model and commits its result to the world only when it finishes. The solved profile, and any solution released from it (`release_radial_solution`), keeps reading exactly the materials it was solved with: a model swapped or edited on a layer afterwards (`set_eos`, `set_partial_melt`, a new viscosity model) takes effect at the next `solve_eos`, and a later solve never changes a solution already handed out. Anything that changes the layer stack leaves the world unsolved: `add_layer` and `load_binary` clear every solved result (the profile getters return NaN and a Love solve raises until `solve_eos` runs again). A failed solve does the same, keeping only its diagnostics (`success`, `message`) and returning NaN profile arrays, and a solve that raises puts floating layers back where they were.

A world takes one of its heavy calls at a time: two threads sharing a world take turns on `solve_eos`, the Love solves, `calc_tides`, the 3D calls, `release_radial_solution`, and `load_binary`, while separate worlds run in parallel (these calls release the GIL). The reads of the solved state take the same turns: the radius getters of the world (`get_density`, `get_state`, `calc_complex_shear_modulus`, and the rest of the profile queries), the same getters on its layer views (`world.mantle.get_density(r)`), and the Love-number and radial-function getters. A getter called while another thread runs `solve_eos` waits for the solve and then reads the new profile, so it never sees a profile being replaced. One call reads a whole array of radii under a single turn, so every value it returns comes from one solve and a loop over an array pays for one lock. Two threads reading one world also take turns; give each thread its own world when the reads themselves must run in parallel. Do not change a world's layers or models from one thread while another thread uses it.

The turns come from one recursive mutex per world (`c_WorldCallLock` in C++). A world hands it to each layer it takes in (`add_layer`, `load_binary`), and the layer's profile getters take it too. It is recursive because the locked calls read the profile through the same getters on the same thread (`calc_tides` evaluates the complex moduli at every radius it integrates). A shared (reader-writer) lock cannot be taken again by the thread that holds it, so it would deadlock on those paths. The lock is only held inside a C++ call and never while waiting for the GIL, so a scalar getter keeps the GIL on its fast path without risking a deadlock, and an array getter releases it.

Raises `ValueError` if the world has no layers, any layer lacks an EOS model, `slices_per_layer < 2`, the outermost layer does not end at the world radius, or the integration method is unknown. A solve that reaches `max_iters` with its surface pressure still off the target by more than `pressure_tol` has found no hydrostatic structure: it returns `success = False`, sets `max_iters_hit`, and leaves the world unsolved. So does a converged solve whose enclosed mass differs from the world's stated mass by more than the factor `[numerical] maximum_eos_mass_ratio` (default 10) either way: layers with no hydrostatic structure near that mass (a core much too dense for its radius, say) can still meet the surface pressure on a collapsed branch at an absurd central pressure, and the message says how far the mass is off. The returned dict contains `success`, `message`, `iterations`, `max_iters_hit`, `pressure_error` \[Pa\], the profile arrays (`radius`, `gravity`, `pressure`, `mass`, `moi`, `density`, `temperature`, `heat_flow`), the per-layer lists (`layer_radius_outer`, `layer_temperature`, `layer_heat_flow_in`, `layer_heat_flow_out`, `layer_heating`, `layer_temperature_rate`, and the network detail `layer_node_temperature` (at the interface above the layer), `layer_top_temperature` (where a convecting layer's adiabat ends), `layer_boundary_thickness` \[m\], `layer_rayleigh_number`, `layer_nusselt_number`, and `layer_in_thermal_network`), the iteration report (`thermal_passes`, `thermal_converged`, `geometry_converged`), and the scalar results (`surface_gravity`, `surface_pressure`, `central_pressure`, `planet_mass`, `planet_moi`).

#### Layer Size

A layer holds its volume by default, so the boundaries a world was built with are the boundaries it solves with. Setting `is_volume_fixed = false` on a layer makes it hold its mass instead: the solve moves its outer radius using its EOS-derived density and constant mass. Every layer above it moves with it, each keeping its own volume (unless they too are not volume fixed). The world radius follows the outermost layer.

The mass a floating layer holds is its `mass_kg`, or, when its configuration gives none, the mass its first solve finds inside the boundaries it was built with. `solve_eos(reset_layer_masses=True)` takes the current geometry as the new reference.

Each pass measures how far the layer is from that mass and steps its outer radius by the mass it is short of over the slope of the enclosed mass, $dm/dr = 4 \pi r^2 \rho$, which is exact to first order, so a few passes settle it. `geometry_converged` says whether they did, and `layer_radius_outer` reports where the boundaries ended up. A layer whose mass cannot fit inside its own base raises `RuntimeError`.

```python
world.core.is_volume_fixed = False     # the core holds its mass, not its size
result = world.solve_eos()
result["layer_radius_outer"]           # [m] where the boundaries settled
world.radius                           # follows the outermost layer
```

#### Temperature and Heat Flow

Each layer carries its own temperature (`temperature_k`, see [PhysicsLayer](../layers/physics_layer.md)) and its [cooling model](../../cooling_x/cooling_models.md) says how heat moves inside it. The solve turns those into a temperature profile, the heat flowing through every radius, and the rate each layer's temperature changes at. `temperature` overrides every layer's value with one number, and `surface_temperature` \[K\] is what the outermost layer radiates to; left out, no heat leaves the world.

What each cooling model makes of its layer:

| Model | Profile inside the layer | Where its temperature applies |
|---|---|---|
| `off` | Isothermal: one temperature throughout, and no modeled gradient, so the layer conducts perfectly. | Everywhere |
| `conduction` | Two conducting halves, $T = T_0 - (L / 4 \pi k)(1/r_0 - 1/r)$. | The mid-radius |
| `convection` | A conducting boundary layer at the base and the top, sized by the model's Nusselt scaling, around an adiabatic interior, $dT/dr = -\alpha g T / c_p$. | The base of the interior |

A layer with no temperature of its own, a geometry-only `BaseLayer` or a layer whose temperature is not a positive number (the 0 K default, say), takes no part in the network: its interfaces carry no heat, so it is neither a heat sink nor a source for its neighbors, and a neighbor keeps its own temperature at the shared interface. It is isothermal at its placeholder temperature, and `layer_in_thermal_network` reports it `False`. A world-wide `temperature` gives every layer that temperature, geometry-only layers included.

A convecting layer's Rayleigh number takes the temperature drop across both of its boundary layers: from the top of the layer below (the end of its adiabat, when that layer convects) to the layer's own temperature, plus from that temperature to the layer above or to `surface_temperature`. The innermost layer has only the upper drop. Both boundary layers take the one thickness the resulting Nusselt number gives, at most 40 percent of the layer each, so a mantle at the temperature of the layer above it still convects when the core below it is hotter. A cooling model that gives no usable thickness (a NaN viscosity at the layer's temperature is the usual cause) leaves each boundary layer at that 40 percent, and the solve logs a warning naming the layer.

The layers form a chain of thermal resistances. A conducting spherical shell between $r_a$ and $r_b$ has

$$R = \frac{1}{4 \pi k} \left( \frac{1}{r_a} - \frac{1}{r_b} \right)$$

and the heat flow through an interface is $L = \Delta T / R$ across the two resistances facing it. Both layer temperatures are inputs, so that flow is generally not the same entering a layer as leaving it. The difference is the heat the layer stores or releases, which is what `layer_temperature_rate` reports:

$$M c_p \frac{dT}{dt} = L_\mathrm{in} - L_\mathrm{out} + H$$

with $H$ \[W\] the heat generated inside the layer (`layer_heating`), zero unless the layer is heated.

**Internal heating**

A layer with `use_heating` set is heated by the world's heat sources. The radiogenic source takes the layer's [radiogenics model](../../radiogenics_x/radiogenics_models.md) as a specific rate $\epsilon$ \[W kg$^{-1}$\] at the solve's `time` \[s\] and heats the layer at $h = \epsilon \rho$ \[W m$^{-3}$\] with the local density, so it is exact for a layer whose mass is an output of the solve. `time=None` takes each model's own reference time. The structure solve integrates

$$\frac{dL}{dr} = 4 \pi r^2 h$$

so the heat flow grows through a heated layer and its conducting stretches bend: a uniformly heated conducting shell follows $T = B + A/r - h r^2 / 6k$. The resistance chain accounts for the same heating. With $H(r)$ the heat generated between the base of a conducting stretch and $r$, the flow leaving its top is the flow entering plus $H$, and the temperature drop across it is $L_\mathrm{base} R + \int H / (4 \pi r^2 k) \, dr$. Both follow from the heating and the solved density, so the solved profile still passes through every layer temperature and arrives at `surface_temperature`. A heated world is a thermal solve even when its layers share one temperature. With `solve_temperature=False` there is no heat flow to act through, so the heating is ignored and a warning is logged.

A world whose layers are all at one temperature has no profile to integrate. The solve then keeps its four structure variables and returns exactly what it returns with `solve_temperature=False`, at the same cost; the profile queries still report each layer's own temperature. Otherwise the solve adds temperature and heat flow as two more state variables and repeats: the first pass is isothermal, and each later pass integrates the profile and then relaxes the boundary layers, interface temperatures, and heat flows against the structure it produced. `thermal_passes` counts them and `thermal_converged` says whether they settled; a solve that uses every pass without settling logs a warning and keeps the last pass. Layers that hold their mass move in the same passes (`geometry_converged`), and they end on the grid the last pass solved on, so every reported slice and interface lies in its own layer.

The viscosity and partial-melt models of every layer are evaluated at the solved temperature of each slice, so an Arrhenius layer is stiff where the profile is cold. A layer whose `use_thermal_eos` is set also passes that temperature to its EOS model, so its density follows the profile.

```python
world.mantle.temperature = 1600.0            # [K] the layer's own temperature
world.mantle.set_cooling(make_cooling("convection"))
result = world.solve_eos(
    surface_temperature=250.0)               # [K] what the outermost layer radiates to

world.get_temperature(0.9 * world.radius)    # [K] on the solved profile
world.get_heat_flow(world.radius)            # [W] leaving the world
result["layer_temperature_rate"]             # [K/s] per layer, from its heat imbalance

world.mantle.use_heating = True              # The mantle's radiogenics model now heats it
result = world.solve_eos(
    surface_temperature=250.0,
    time=1.0e17)                             # [s] when the radiogenics models are evaluated
result["layer_heating"]                      # [W] generated inside each layer
```

**Profile queries (after a successful solve)**

`get_state(radius)` (the density, gravity, pressure, static moduli, viscosities, and melt fraction as a dict) and `get_static_viscoelastics(radius)` read one evaluation of the solved state per radius, so they cost about what one single-quantity getter does and return exactly what the single getters return.

| Member | Returns | Description |
|--------|---------|-------------|
| `get_density(r)` | float [kg/m³] | Density at radius `r` [m] (NaN if unsolved). |
| `get_gravity(r)` | float [m/s²] | Gravitational acceleration at `r`. |
| `get_pressure(r)` | float [Pa] | Pressure at `r`. |
| `get_temperature(r)` | float [K] | Temperature at `r` on the solved profile. |
| `get_heat_flow(r)` | float [W] | Heat flowing outward through the sphere of radius `r`. |
| `eos_solved` | bool | `True` once profiles are populated. |
| `all_eos_set` | bool | `True` once every layer has an EOS model. |
| `surface_gravity_eos`, `central_pressure`, `planet_mass_eos`, `planet_moi_eos` | float | Scalar results of the last solve (NaN if unsolved). |
| `molten_regions` | list | Molten stretches of solid layers in the last solve, as `(layer_name, radius_inner, radius_outer)` with radii in \[m\]. The radial solver treats each as a static liquid (see Molten Stretches in a Solid Layer under Calculating Love Numbers). |

`get_density` / `get_gravity` / `get_pressure` delegate to the layer that contains `r` (radii beyond the surface clamp to the outermost layer), so the individual layers expose the same getters independently.

#### C++ API

The entire solve runs in C++ and can be driven directly from other C++ code:

```cpp
tidalpy::c_WorldEOSSolveConfig cfg;          // starts from the [eos_solver] section of the configuration
cfg.surface_pressure = 1.0e5;                // override only what this solve changes
world.solve_eos(cfg);                        // populates each layer's c_LayerEOSData
double rho = world.get_density(5.0e6);
```

`c_LayeredWorld::solve_eos(const c_WorldEOSSolveConfig&)` generates the per-layer radius grids, estimates the bulk density, converts to the non-dimensional solve units, calls `Material_x/eos`'s `c_solve_eos` (CyRK ODE integration with the secant surface-pressure iteration), returns the solution to SI, and slices it into each layer's `c_LayerEOSData`. It throws `std::invalid_argument` on bad input (surfaced as `ValueError` in Cython via `except +`). The per-layer density source is `tidalpy::c_MaterialEOSBase`, attached via `c_BaseLayer::set_eos(std::unique_ptr<c_MaterialEOSBase>)` (non-owning observer through `get_eos()`; `get_eos_set()`); during integration the `c_preeval_material_eos` pre-eval (in `Material_x/eos/methods/material_.hpp`, paired with the `c_MaterialEOSInput` struct holding the model pointer and the unit scales of the solve) dispatches to `model->calc_density(pressure, temperature, radius)` in SI. The retained integrators carry only the four structure variables; the density, moduli, and viscosities at any radius are evaluated on demand from the interpolated state by the same EOS function, so a dense evaluation is a plain polynomial evaluation plus one model call.

`c_LayeredWorld` exposes `get_density/get_gravity/get_pressure(double) const` (delegating to the containing layer via `find_layer_for_radius`), the vectorized `get_eos_fields(field_indices, num_fields, radii, num_radii, values_out) const` (entries of the `eos_layout_.hpp` evaluation layout at every radius, field-major, under one hold of the call lock) and `calc_complex_moduli(is_shear, radii, num_radii, frequency, moduli_out) const`, the result accessors `get_eos_success/get_eos_message/get_eos_iterations/get_eos_pressure_error/` `get_surface_gravity_eos/get_surface_pressure_eos/get_central_pressure/` `get_planet_mass_eos/get_planet_moi_eos`, the retained full solution via `get_eos_solution() -> const c_EOSSolution*`, and `get_eos_solved()` / `get_all_eos_set()`. The Cython `LayeredWorld.solve_eos` wrapper only converts the integration-method string to the CyRK enum, fills `c_WorldEOSSolveConfig`, calls the C++ method under `nogil`, and builds the Python result dict from the retained solution.

### Viscoelastic Properties (after EOS solve)

Once `solve_eos` has succeeded, the world and each layer expose the viscoelastic getters the radial solver needs.

| Member | Returns | Description |
|--------|---------|-------------|
| `get_shear_modulus(r)` | float [Pa] | Post-melt static shear modulus at `r`. |
| `get_bulk_modulus(r)` | float [Pa] | Post-melt static bulk modulus at `r`. |
| `get_shear_viscosity(r)` | float [Pa·s] | Post-melt shear viscosity at `r`. |
| `get_bulk_viscosity(r)` | float [Pa·s] | Post-melt bulk viscosity at `r`. |
| `calc_complex_shear_modulus(r, ω)` | complex [Pa] | Rheology-derived complex shear modulus at radius `r` [m], frequency `ω` [rad/s]. |
| `calc_complex_bulk_modulus(r, ω)` | complex [Pa] | Rheology-derived complex bulk modulus. |

All of the above accept a float or `np.ndarray` for `r` (and `ω`); array inputs return an `np.ndarray` of the same shape.

### Calculating Love Numbers

`LayeredWorld.solve_love_numbers(...)` computes the viscoelastic-gravitational Love tidal or loading numbers $k$, $h$, $l$ for a given tidal forcing frequency. `solve_eos` must be called first.

The `love_method` argument selects how the Love numbers are obtained (names are case-insensitive; aliases in parentheses):

| `love_method` | What it does |
|---|---|
| `radial_solver` (`shooting`, `rs`; default) | Numerically integrates the radial ODEs from the center to the surface. Works for arbitrary multi-layer, solid/liquid, static/dynamic, compressible/incompressible worlds. |
| `propagation_matrix` (`prop_matrix`, `pm`, `prop`) | Quasi-analytic matrix propagation, restricted to a single solid, static, incompressible layer. An incompatible world fails the solve gracefully (`love_success` is `False`, `love_error_code` non-zero). `core_model` selects the core starting condition. |
| `homogeneous` (`homogen`) | Quasi-homogeneous: each tidal layer is treated as a homogeneous incompressible planet made of its own averaged material, and the world's Love numbers are the sum of the layers' weighted by their tidal scales (see Quasi-Homogeneous Love Numbers below). Each layer uses the homogeneous-sphere formulas, $k_{l} = \frac{3}{2(l-1)}\,\frac{1}{1+\bar{\mu}_{l}}$, $h_{l} = \frac{2l+1}{2(l-1)}\,\frac{1}{1+\bar{\mu}_{l}}$, $l_{l} = \frac{3}{2l(l-1)}\,\frac{1}{1+\bar{\mu}_{l}}$ with $\bar{\mu}_{l} = \frac{2l^{2} + 4l + 3}{l}\,\frac{\mu}{\rho g R}$, where $\mu$ is the layer's rheology at the forcing frequency applied to its averaged moduli and viscosity, and $\rho$, $g$, $R$ are the planet's EOS bulk density, EOS surface gravity, and radius. Fast; no radial functions. |
| `cpl` | The same, with each layer's static (unrelaxed) averaged shear modulus and a constant phase lag: $k$, $h$, $l$ are multiplied by $(1 - i/Q)$ so $-\mathrm{Im}[k] = \mathrm{Re}[k]/Q$. $Q$ is `fixed_q` (argument or `[tides]` config) or, when unset, the attached tide model's fixed Q for the degree. |
| `ctl` | Similar to `cpl` but with a constant time lag, $(1 - i\,\omega\,\Delta t)$; $\Delta t$ is `fixed_dt` or the tide model's fixed time lag. |
| `laterally_inhomogeneous` (`3d`, `lat_inhom`) | Reserved for a future 3D Love solver; raises `NotImplementedError`. |

#### Quasi-Homogeneous Love Numbers

The `homogeneous`, `cpl`, and `ctl` methods keep some of the layered structure without a radial solve. Each tidal layer is averaged into a homogeneous material: a volume average of its post-melt shear and bulk moduli and a log-volume average of its post-melt viscosity, so a viscosity that spans decades across the layer is averaged in the space it varies in. A homogeneous planet of the planet's radius, bulk density, and surface gravity, made of that material, has Love numbers $k_i$; the layer's tidal scale $s_i$ then weighs them, and the world's Love numbers are

$$k = \sum_i s_i k_i$$

(and likewise $h$ and $l$). The tidal scale is the layer's `tidal_scale`, or its volume over the planet's when none is set, and 0 for a layer that is not tidal. With the default scales a one-layer planet gets exactly the homogeneous-sphere value, and a small, highly dissipative layer, an asthenosphere for example, adds dissipation in proportion to its volume instead of making the whole planet dissipate like it. Because the tidal heating is linear in $-\mathrm{Im}[k]$, each layer's heating is the heating of its own term $s_i k_i$, constant within the layer, and the layers sum to the total. Gas layers carry no shear modulus and take no part. `love_layer_parts` lists each layer's $s_i$, $k_i$, $h_i$, $l_i$, and complex shear modulus, and `get_layer_tidal_scale(index)` gives the scale in use.

Every Love solve checks its input first and raises `ValueError` for a degree below 2 in a tidal or free-surface solve (degree 1 is a translation of the body; a solve for loading alone may ask for degree 1, whose load Love numbers depend on the reference frame), a frequency that is not finite and positive or lies outside the `[numerical]` `minimum_frequency` to `maximum_frequency` range (the Love numbers at $-\omega$ are the complex conjugates of those at $\omega$), a `start_radius_tol` outside (0, 1), or a negative `starting_radius` or `max_step`. `max_step` [m] is converted into the units the integration runs in. The `homogeneous`, `cpl`, and `ctl` methods take the bulk density from the solved mass, the same structure the EOS surface gravity comes from, not from the declared mass.

The analytic methods report the tidal-scale-weighted mean of the layers' complex shear moduli and the volume of the layers that took part through `love_effective_shear_modulus` and `love_tidal_volume`, return `love_surface_amplification = 0`, and give NaN for the radial-function getters (`get_radial_solution_y`, ...). They have no depth-resolved solution, so the 3D stress/strain/heating path (`calc_3d_tides`, `get_3d_tidal_heating`) raises `RuntimeError` while an analytic method is the world's configured method. Free-function versions of the formulas live in `TidalPy.Tides_x.love` (`calc_homogeneous_love_numbers`, `calc_effective_rigidity`, `apply_fixed_q`, `apply_fixed_dt`; see [Love numbers](../../Tides_x/love/love_numbers.md)).

The world's default method, used whenever its tide model asks for Love numbers inside `calc_tides`, is set with `set_tide_config(love_method=..., love_fixed_q=..., love_fixed_dt=...)` or the matching `[tides]` keys `love_method`, `love_fixed_q`, `love_fixed_dt_s` in a world TOML file; `solve_love_numbers` takes the method per call.

```python
from TidalPy.structures_x.worlds import LayeredWorld
from TidalPy.structures_x.layers import SolidLiquidLayer
from TidalPy.Material_x.eos import make_material_eos
from TidalPy.rheology_x import make_rheology
from TidalPy.viscosity_x import make_viscosity

world = LayeredWorld("planet", 6.0e6, 4.2e24)
layer = SolidLiquidLayer("mantle", 0, 0.0, 6.0e6, 4.2e24)
layer.set_eos(make_material_eos(
    "constant", {"reference_density_kg_m3": 4000.0, "shear_modulus_static_pa": 6.0e10, "bulk_modulus_static_pa": 1.3e11}))
layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e21}))
layer.set_shear_rheology(make_rheology("maxwell"))
world.add_layer(layer)

world.solve_eos()
result = world.solve_love_numbers(frequency=1.0e-5)

print(result["love_number_k"])   # complex k2, also world.love_number_k
print(world.love_number_h, world.love_number_l)
```

The moduli and the viscosity are properties of the layer, not of the rheology model: a rheology model holds only its own shape parameters (the Andrade exponent, the Voigt fractions), and reads the modulus and viscosity it is handed. A layer with no shear modulus and no viscosity model deforms as if it had no strength, and the solve fails rather than guessing.

**`solve_love_numbers( frequency=1e-5, degree_l=2, solve_for='tidal', core_model=0, use_kamata=None, nondimensionalize=None, starting_radius=0.0, start_radius_tol=None, integration_method=None, rtol=None, atol=None, scale_rtols=None, max_num_steps=None, expected_size=None, max_ram_MB=None, max_step=0.0, verbose=False, warnings=True, love_method=None, fixed_q=None, fixed_dt=None) -> dict`**

Every solver setting left as `None` takes the `[radial_solver]` value of the TidalPy configuration (see [Configurations](../../Overview/2_TidalPy_Configurations.md)), the same defaults the standalone `radial_solver` and the world's own tidal solves use; `love_method`, `fixed_q`, and `fixed_dt` left as `None` take the world's `[tides]` settings. Raises `ValueError` if the EOS has not yet been solved. Returns a dict (`success`, `error_code`, `message`, `love_method`, `love_number_k/h/l`); the results are also stored internally and accessed through the properties below.

The solver interpolates nothing between EOS slices. Gravity, pressure, mass, and moment of inertia come from the world's dense EOS solution, the density and the static moduli and viscosities from the same solution (the layer's material evaluated them as the structure was integrated), and the complex moduli from the layer's rheology applied to those static values at that solve's frequency, all at the exact radius the integrator asks for. The Love numbers therefore converge with the integration tolerance alone and are independent of `slices_per_layer` exactly, not just to within an interpolation error, whether or not the moduli and viscosities vary with depth.

`slices_per_layer` still sizes the profile arrays the solve returns and the `[layers.*]` array properties, and the propagation-matrix method still propagates across those slices, so it remains a real knob for those. It no longer affects a shooting-method Love number.

`release_radial_solution()` hands the world's last radial solve out as a `RadialSolverSolution`, with its `result` grid filled on the solve's radius grid (`sample_radii()` gives those radii, so the two plot together). A radius above the surface has no solution and reads NaN.

#### Molten Stretches in a Solid Layer

A solid layer's partial-melt model can weaken part of the layer past use as a solid, for example the base of a mantle over a hot core. Past the critical melt fraction the post-melt shear modulus falls steeply to the model's `liquid_shear` floor (10$^{-5}$ Pa by default), and the solid equations, which divide by the shear modulus, cannot be integrated through it. After every EOS solve the world marks as molten each stretch of a layer with a partial-melt model where the modulus sits at that floor or where its rigidity $\mu / (\bar{\rho} g R)$ (planet bulk density, surface gravity, and radius) is below `minimum_solid_rigidity` in the `[numerical]` section of the TidalPy configuration (10$^{-6}$ by default). The stretches are found on the EOS slices and their edges refined by bisection on the dense profile. `molten_regions` lists them as `(layer_name, radius_inner, radius_outer)`, and a message at the info log level names each one.

The radial solver splits the layer at those edges and solves each molten stretch as a static liquid, which reads only the density and gravity. The partial-melt model does not set a liquid's bulk modulus there (it changes the bulk modulus only when its `bulk_melt_weakening` switch is on), which a compressible dynamic liquid would read. The solid parts keep the layer's own flags, and each stretch takes its share of the layer's slices, at least five, for the solution's output grid. Treating a solid of rigidity $10^{-6}$ as a liquid changes the Love numbers by about that fraction. The split changes only the radial solve: the layer is one layer everywhere else, its 3D heating takes nothing from the molten stretch (a liquid carries no shear dissipation there), and a layer declared liquid is never split.

> [!NOTE]
> Only a layer with a partial-melt model is split. A solid layer given a near-zero shear modulus some other way is solved as a solid, and the solver may fail on it or return $\mathrm{Im}[k]$ with the wrong sign and only a conditioning warning.

**Love-number properties (after a successful solve)**

These describe the world's own last `solve_love_numbers` (or `solve_love_numbers_supplied`) call. The Love solves `calc_tides` and the 3D paths run for each tidal frequency use workspaces of their own and leave these untouched, so after only a `calc_tides` they still read as unsolved (`love_success` is `False` and `love_error_code` is -100). A failed Love solve inside `calc_tides` raises `RuntimeError` instead.

| Property | Type | Description |
|----------|------|-------------|
| `love_solved` | bool | `True` while a successful solve is held. A later `solve_eos` clears it: Love numbers describe the structure they were solved with, so every Love number and radial-function getter returns NaN until the next solve. |
| `love_success` | bool | `True` if the last solve converged and still describes the structure. |
| `love_error_code` | int | Solver error code (0 = success; < 0 = failure; -100 = no Love solve has run since the world was built or last solved its EOS). |
| `love_message` | str | Human-readable solver message. |
| `love_num_ytypes` | int | Number of independent solution types (boundary-condition models requested). |
| `love_number_k`, `love_number_h`, `love_number_l` | complex | Love numbers for the first boundary condition at the solved degree. Equivalent to `get_love_number_k(0)` and friends. |
| `love_method` | str | Canonical name of the method the last solve used. |
| `love_surface_amplification` | float | Conditioning of the surface boundary-condition solve, recorded on every shooting solve whether or not `warnings` is on; near 1 is healthy, 0 after an analytic solve. |
| `love_surface_rcond` | float | Reciprocal condition number of the surface boundary-condition system, the rank measure the amplification cannot give: near machine epsilon the solution constants are undetermined and the solve fails with error code -13 (below `[numerical] minimum_surface_rcond`); below the integration rtol it draws a conditioning warning. NaN before a shooting solve and for the other methods. The standalone solver reports it as `surface_solve_rcond`. |
| `love_effective_shear_modulus`, `love_tidal_volume` | complex, float | The tidal-scale-weighted mean of the layers' complex shear moduli [Pa] and the volume of the layers that took part [m3] in the last quasi-homogeneous solve; NaN after a radial-solver solve. |
| `love_layer_parts` | list of dict | Each tidal layer's part of the last quasi-homogeneous solve: `layer`, `tidal_scale`, `love_number_k`, `love_number_h`, `love_number_l`, and `shear_modulus`. Empty after a radial-solver solve. |

**`get_love_number_k(ytype_idx=0) -> complex`**, **`get_love_number_h(ytype_idx=0) -> complex`**, **`get_love_number_l(ytype_idx=0) -> complex`**

Return the Love numbers for boundary-condition model index `ytype_idx` (0 = first requested, usually tidal).

**`get_love_surface_y(ytype_idx, y_idx) -> complex`**

Raw radial function y₁…y₆ at the surface for solution type `ytype_idx`, function index `y_idx` (0–5).

**`get_love_radial_y(radius, ytype_idx=0, y_idx=0) -> complex`**

The same radial function at any radius [m], evaluated from the solver's dense interpolants, so it is accurate between grid slices too. Returns NaN if the solve failed, if an analytic Love method was used (those have no radial functions), or if the radius sits below the solver's starting radius.

#### C++ API

```cpp
tidalpy::c_LoveSolveConfig cfg = world.make_love_solve_config();   // [radial_solver] defaults + the world's [tides]
cfg.frequency = 1.0e-5;                                            // [rad/s]
cfg.degree_l  = 2;
world.solve_love_numbers(cfg);   // delegates to the cached radial solver

std::complex<double> k2 = world.get_love_number_k(0);
```

A default-constructed `c_LoveSolveConfig` (and `c_WorldEOSSolveConfig`) reads the `[radial_solver]` (`[eos_solver]`) section of the shared runtime config, so C++ callers and the tide paths start from the same defaults as Python callers.

`c_LayeredWorld::solve_love_numbers(const c_LoveSolveConfig&)` delegates to a cached helper, `c_WorldRadialSolver` (held by `p_radial_solver`), that separates the frequency-independent setup (built once and reused) from the frequency-dependent work (recomputed on every call), since the Love-number solve is the hot loop for frequency sweeps and orbital evolution.

The non-dimensionalization is itself frequency-independent (the `c_NonDimensionalScales` time scale is $1/\sqrt{\pi G \bar{\rho}}$ for the bulk density $\bar{\rho}$, not $1/\omega$), so the only quantities that change between calls at different frequencies are the complex moduli and the shooting integration.

1. Validates `eos_solved` and `tidalpy_config_ptr`.
2. If the cache does not match the current EOS grid/assumptions, `build_cache` captures (once), for the solver's layers (the world's layers, with a solid layer split at the edges of its molten stretches; `get_radial_segments()` lists the stretches and the providers map each solver layer back to its world layer): the non-dim radius/density/gravity/pressure/mass/moi arrays, per-layer metadata (solid/liquid, static, incompressible) and slice partitioning, the non-dim scalars (`G`, bulk density, surface pressure), and a reused `c_RadialSolutionStorage` whose internal `c_EOSSolution` arrays serve as the scratch buffers. The cache is invalidated automatically whenever `solve_eos` re-runs.
3. Per call: the world installs a material-state provider (`c_EOSSolution::MaterialEval`, a type-erased callable carrying that solve's frequency) that the shooting solve calls for density and the complex moduli at each integration radius, and also fills the dimensional moduli scratch at the slice radii for the propagation-matrix method and the array outputs. The helper non-dimensionalizes the scratch in place, re-applies the cached non-dim structure arrays via `inject_from_world_eos`, and runs the selected solver.
4. Re-dimensionalizes the y-solution, restores the SI surface gravity, and calls `c_RadialSolutionStorage::find_love()`.

A Love solve writes everything it produces (the cached radial solver and its storage, or the quasi-homogeneous results) into a `c_LoveWorkspace`. The world holds one for `solve_love_numbers`, which `get_love_number_k/h/l`, `get_love_surface_y`, and the status accessors read; `calc_tides` and the 3D paths solve into workspaces of their own.

### Global (1D) Tidal Dissipation

`LayeredWorld.calc_tides(...)` computes the body's total tidal heating and three orbital potential partial derivatives by summing over the active tidal modes (the global / "1D potential" approach). A tide model (see [Global Tidal Dissipation](../../Tides_x/global_tides.md)) supplies the per-mode dissipation multiplier $-\mathrm{Im}[k_{l}]$; the world runs the global-potential engine for its stored `[tides]` config + the supplied orbital/spin state, collapses, and resolves each layer's share of the heating. The Love solves `calc_tides` runs are its own: they leave the world's last `solve_love_numbers` result, and what `release_radial_solution` hands out, untouched.

#### Tidal Heating of Each Layer

How a layer's heating is found depends on where the Love numbers come from:

| Love numbers from | Layer heating |
|-------------------|---------------|
| The radial solver (`radial_solver`, `propagation_matrix`) | The volume integral of the radial solution's orbit-averaged heating density over the layer: the same integral as `calc_3d_tides` with every axis summed, scaled so the layers sum to the 1D total. A liquid layer carries no shear dissipation and takes 0. It costs about as much as the global solve again, so `set_tide_config(layer_tidal_heating=False)` (or `layer_tidal_heating = false` in `[tides]`) skips it and leaves every layer's heating NaN. |
| The quasi-homogeneous methods (`homogeneous`, `cpl`, `ctl`) | The heating of the layer's own term $s_i k_i$ (see Quasi-Homogeneous Love Numbers), constant within the layer; the layers sum to the total. |
| An analytic tide model (`cpl`, `ctl`, `ctl_q` tide models, which describe the whole body) | The total times the layer's tidal scale $s_i$. |

A layer that is not tidal (`is_tidal = false`) has a tidal scale of 0 and takes no part in the quasi-homogeneous methods.

The tide model + config are normally wired by the world builder from the `[tides]` TOML table (with per-family defaults: star = `fixed_q`, gasgiant = `fixed_dt`, terrestrial = `rheology`), but can also be set directly:

```python
from TidalPy.Tides_x.classes import make_tide

world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]}))
world.set_tide_config(max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)

world.calc_tides(orbital_frequency=2.05e-5, spin_frequency=2.05e-5, eccentricity=0.0041,
                 obliquity=0.0, semi_major_axis=4.2e8, host_mass=1.898e27)

world.get_tidal_heating()                # total heating [W]
world.get_tidal_potential_derivatives()  # (dUdM, dUdw, dUdO) [J kg-1 rad-1]
world.get_layer_tidal_heating(0)         # = world heating × layer 0's tidal scale (an analytic tide model)
```

For a synchronous, low-eccentricity body the `cpl` result reproduces the standard CPL rate $\frac{21}{2}\,\frac{k_{2}}{Q}\,\frac{G M_{h}^{2} R^{5} n e^{2}}{a^{6}}$ (host mass $M_{h}$).

### Rheology Model

The analytic models (`cpl`/`ctl`/`ctl_q`) take $-\mathrm{Im}[k_{l}]$ from their fixed per-degree parameters and need no interior solution. The `rheology` model instead derives $-\mathrm{Im}[k_{l}(\omega)]$ from the world radial solver: `calc_tides` runs the global-potential engine, then for each unique tidal frequency it solves the world's complex Love numbers (reusing the frequency-independent radial-solver cache), feeds the per-mode `k_l` into the collapse, and retains the full `k`/`h`/`l` suite per mode for inspection. Because it runs the radial solver, the EOS must be solved first:

```python
world.solve_eos(G_to_use=G, temperature=1500.0)   # required for the rheology model
world.set_tide_model(make_tide("rheology"))
world.set_tide_config(max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
world.calc_tides(orbital_frequency=2.0e-5, spin_frequency=1.0e-5, eccentricity=0.01,
                 obliquity=0.0, semi_major_axis=4.0e8, host_mass=1.9e27)

world.get_tidal_heating()             # total heating [W]
world.get_tidal_love_k(2, 2, 0, 0)    # complex k₂ for the (l,m,p,q) = (2,2,0,0) mode
```

`calc_tides` raises `RuntimeError` if the `rheology` model is selected but the EOS has not been solved, or if a per-frequency radial solve fails.

| Member | Returns | Description |
|--------|---------|-------------|
| `set_tide_model(tide)` | - | Attach a tide model (transfers ownership). |
| `tide_model_set` | bool | Whether a model is attached. |
| `set_tide_config(min_degree_l=None, max_degree_l=None, eccentricity_truncation=None, obliquity_truncation=None, layer_tidal_heating=None, eccentricity_exact_tolerance=None, love_method=None, love_fixed_q=None, love_fixed_dt=None)` | - | Change the stored `[tides]` truncation/degree settings (`eccentricity_truncation` a level or `"exact"`, with `eccentricity_exact_tolerance` its mode range; see [Eccentricity Functions](../../Tides_x/eccentricity.md); `obliquity_truncation` 0 or `"off"`, 2, 4, or `"gen"`; see [Obliquity Functions](../../Tides_x/obliquity.md)), whether `calc_tides` resolves each layer's heating on the radial-solver path, and the world's default Love-number method (see the RadialSolver section). Only the arguments given change; the rest keep their current values (`get_tide_config()`). A NaN `love_fixed_q` or `love_fixed_dt` clears it. |
| `get_tide_config()` | dict | The stored settings under the builder's `[tides]` key names (`*_trunc_lvl`). |
| `calc_tides(orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass)` | - | Run the global tidal solve. |
| `tides_solved` | bool | Whether a successful solve is held. A new tide model or tide configuration clears it, and so does a layered world's `solve_eos`, after which the heating and potential-derivative getters return NaN and each layer's `get_tidal_heating()` does too, until the next `calc_tides`. |
| `get_tidal_heating()` | float [W] | Total global tidal heating (NaN if unsolved). |
| `get_tidal_potential_derivatives()` | tuple | `(dUdM, dUdw, dUdO)` [J kg⁻¹ rad⁻¹]. |
| `get_tidal_dU_dM_minus_dw()` | float | The per-mode sum of `dUdM - dUdw` [J kg⁻¹ rad⁻¹], to pass as `dU_dM_minus_dw` to `OrbitSolver.calc_de_dt`: at small eccentricity the two separate sums nearly cancel in $\dot{e}$. NaN if unsolved. |
| `get_num_tidal_modes()` | int | Active modes summed. |
| `get_layer_tidal_heating(index)` | float [W] | The heating the last `calc_tides` put in the layer (see Tidal Heating of Each Layer); NaN before one. |
| `get_layer_tidal_scale(index)` | float | The layer's tidal scale: its `tidal_scale`, or its volume over the planet's when none is set; 0 when not tidal. |
| `get_tidal_love_k(l, m, p, q)` | complex | Per-mode radial-solver `k_l` (rheology only; NaN for analytic models). |

Each layer also stores its own tidal heating: `layer.get_tidal_heating()` (C++ `c_BaseLayer::get_tidal_heating()`) returns the same value as `get_layer_tidal_heating`.

### Other World Members

The remaining public surface, grouped by what it is for.

**Radial queries after an EOS solve.** Each takes a radius [m] and returns the value there, read from the solver's dense interpolants rather than by re-interpolating the gridded arrays. They accept a scalar or an array.

| Member | Returns |
|---|---|
| `calc_gravity`, `calc_pressure` | Gravity [m s-2] and pressure [Pa]. |
| `calc_shear_modulus`, `calc_bulk_modulus` | Complex moduli [Pa] at a forcing frequency. |
| `calc_shear_viscosity`, `calc_bulk_viscosity` | Viscosities [Pa s], after any melt weakening. |
| `calc_static_viscoelastics`, `get_static_viscoelastics` | The static (unrelaxed) moduli and viscosities together. |
| `get_melt_fraction` | Melt fraction from the material's partial-melt model; `0.0` where it has none. |

**Geometry.** `calc_surface_area(radius)`, `calc_volume_sphere(radius)`, and `calc_volume_shell(outer, inner)` are the shared spherical helpers every structure inherits.

**Spin and orbit.** `set_spin_model(spin)` attaches a spin model; `get_moment_of_inertia()` returns the EOS-solved moment of inertia [kg m2] (or the spin model's estimate `moment_of_inertia_factor` $\times\,M R^{2}$ before a solve); `calc_spin_derivative(host_mass)` gives the spin rate of change [rad s-2] from the current tidal solution; `calc_synchronous_spin(orbital_frequency)` returns the synchronous rate [rad s-1]. See [Dynamics](../../dynamics_x/dynamics.md).

**State.** `get_state()` returns the world's current scalar state as a dict, and `calc_state()` recomputes it. `get_state()` is the cheap read.

**Three-dimensional tides.** `get_3d_tidal_heating_array(...)` is the vectorized form of `get_3d_tidal_heating` and the efficient way to build a heating map; `calc_3d_displacements(...)` and `calc_3d_stress_strain(...)` return the instantaneous displacement grid and the stress and strain grids. All three, like `calc_3d_tides`, take `num_threads` (default 1) to spread their per-point work over threads, and all are described on the [3D heating page](../../Tides_x/multilayer_3d_heating.md).

**Configuration and identity.** `source_config` is the normalized configuration the world was built from, when it was built from one, and `portable_config` the configuration as given for a world built from a `data_file` (what `save_to_toml` writes); `family_world_type()` gives the builder's world type for this class; `get_schema_version_str()` reports the schema version the class writes. See the [TOML schema](../config/toml_schema.md).

**Pinned solver settings.** `set_solver_defaults(eos_solver=None, radial_solver=None)` pins keys of the `[eos_solver]` and `[radial_solver]` sections of the TidalPy configuration on this world (a world file's tables of the same names arrive here), and `get_solver_defaults()` returns them. A pinned key wins over the configuration for every solve the world runs, a call's own argument wins over the pinned key, and an unpinned key keeps following the configuration; `get_config_dict()` carries the tables.

## `GasGiantWorld`

A `LayeredWorld` whose `world_type` defaults to `"gasgiant"` and which uses a dedicated binary class id 202 (`BinaryClassID::GasGiantWorld`). Same API as `LayeredWorld`; typically populated with `GasLayer`s.

```python
from TidalPy.structures_x.worlds import GasGiantWorld
from TidalPy.structures_x.layers import GasLayer

jupiter = GasGiantWorld("Jupiter", 7.0e7, 1.898e27)
jupiter.add_layer(GasLayer("envelope", 0, 0.0, 7.0e7, 1.898e27))
```

## `StarWorld`

A star: no layers, no EOS. Effective temperature and luminosity are kept consistent through the Stefan-Boltzmann law, $L = 4\pi R^{2}\sigma T^{4}$.

```python
from TidalPy.structures_x.worlds import StarWorld

sun = StarWorld("Sun", 6.957e8, 1.989e30, effective_temperature=5772.0)
sun.luminosity                # ~3.83e26 W (derived from T if luminosity == 0)
sun.set_luminosity(3.828e26)  # recomputes effective_temperature
```

**Properties:** `effective_temperature` [K], `luminosity` [W]. **Methods:** `calc_luminosity_from_temperature(T)`, `calc_temperature_from_luminosity(L)`, `set_effective_temperature(T)`, `set_luminosity(L)`. A luminosity-model hierarchy (fixed, mass-to-luminosity, power law) can be attached via `set_luminosity_model` (see `stellar_x/luminosity.md`).

**Tides.** The analytic tide pipeline (`set_tide_model`/`set_tide_config`/`calc_tides` and the `get_tidal_*` accessors) lives on `BaseWorld`, so a star dissipates tidally too, with the analytic models only (`cpl`, `ctl`, `ctl_q`). The `rheology` model needs the radial solver and a layered interior, so `calc_tides` raises if it is selected on a star. A star has no layers, so there is no per-layer heating distribution. See [Global tidal dissipation](#global-1d-tidal-dissipation).

Binary class id 203 (`BinaryClassID::StarWorld`).

## Binary Serialization

A `LayeredWorld` (and `GasGiantWorld`) serializes its `BaseWorld` fields and a layer count, then each layer's own complete binary record in index order. Because each layer recursively serializes its attached material EOS, rheology, viscosity, partial-melt, cooling, and radiogenics models (see [Binary serialization](../../utilities_x/binary_x.md)), a single `save_binary` / `load_binary` round-trips the entire world graph: no Python reconstruction step is needed. On load, each layer is rebuilt as the correct concrete subclass via the layer binary-dispatch factory (`c_layer_from_binary`).

The record also carries every setting that changes a result, so a loaded world computes what the saved one did:

- Every world type: the tide model and its configuration (`get_tide_config()`: degrees, truncations, Love method, `love_fixed_q`, `love_fixed_dt`, `layer_tidal_heating`).
- `LayeredWorld` and `GasGiantWorld`: the spin model's `moment_of_inertia_factor` and the pinned solver settings (`get_solver_defaults()`).
- `StarWorld`: the luminosity model.

Solved state is not saved: the EOS profile, Love numbers, and tide results are recomputed with `solve_eos` and `calc_tides` after a load. Loading into an existing world replaces all of these, so a record saved without a tide model leaves the world without one.

```python
world.save_binary("earth.tpyb")
reloaded = LayeredWorld("placeholder", 1.0, 1.0)
reloaded.load_binary("earth.tpyb")
assert reloaded.num_layers == world.num_layers
assert reloaded.calc_internal_heating(0.0) == world.calc_internal_heating(0.0)
```

EOS profile data (`c_LayerEOSData`) is never serialized; it is repopulated by running the whole-planet EOS solve after loading. The layers' material EOS models are restored, so that solve needs nothing re-attached. A load into a world that was solved clears everything solved, and layer views taken from the world before the load (`world.<name>`, `get_layer`) refer to layers the load replaced, so take new ones. A file must hold a record of the class it is loaded into: a `LayeredWorld` file into a `LayeredWorld`, a `PhysicsLayer` file into a `PhysicsLayer`, a model's file into the same model; anything else raises `IOError` before the object is touched, as does a corrupt or truncated file. A layer that belongs to a world cannot be loaded in place; load the world, or a standalone layer.
