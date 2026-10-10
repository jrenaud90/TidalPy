# Worlds (`Structures.worlds`)

_Updated: 2026-10-09_

A world is TidalPy's top-level structure: a planet, moon, or star. It holds its identity, orbital and thermal scalars, bulk geometry, spin model, tide model, and an ordered stack of [layers](../layers/layer.md) (possibly empty). It runs the whole-planet equation-of-state (EOS), thermal, Love-number, and tidal solves, and holds the heat sources that act inside its layers. Build one by name or from a file with `build_world` (see the [TOML Schema](../config/toml_schema.md)), or in Python.

## Quick Start

```python
import math

from TidalPy.Structures import build_world

io = build_world("io")                                 # A bundled world, or a path to a world file
result = io.solve_eos()                                # The interior: density, gravity, pressure, ...
print(io.summary())                                    # One row per layer
print(result["success"], io.moment_of_inertia_factor)  # C / (M R^2) of the solved structure

orbital_frequency = 2.0 * math.pi / (1.769 * 86400.0)  # [rad s-1]
love = io.solve_love_numbers(
    frequency=orbital_frequency,
    degree_l=2)                                        # Complex k, h, l at one frequency
print(love["love_number_k"], io.love_q_k)

tides = io.calc_tides(
    orbital_frequency=orbital_frequency,
    spin_frequency=orbital_frequency,
    eccentricity=0.0041,
    obliquity=0.0,
    semi_major_axis=4.217e8,
    host_mass=1.898e27)                                # Tides about Jupiter
print(tides["tidal_heating"])                          # Total [W]
print(tides["layer_tidal_heating"])                    # [W] by layer name
```

## Building a World

```python
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds import TerrestrialWorld

world = TerrestrialWorld(
    name="Earth",
    radius=6.371e6,
    mass=5.972e24,
    albedo=0.3,
    emissivity=1.0,
    obliquity=0.41,
    spin_frequency=7.29e-5)
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

The classes (TOML `type` in parentheses) are `TerrestrialWorld` (`terrestrial`), `GasGiantWorld` (`gasgiant`), `StarWorld` (`star`, usually without layers), and their base `BaseWorld` (`layered`, with a `world_type` argument). All share the `BaseWorld` API (see [World Classes](#world-classes)).

**Read-only properties:** `name`, `world_type`, `radius`, `mass`, `albedo`, `emissivity`, `obliquity`, `spin_frequency`, and `moment_of_inertia_factor`, the solved structure's $C / (M R^2)$, NaN before a successful `solve_eos` (a spin calculation uses the spin model's factor; see [Other World Members](#other-world-members)).

| Method | Returns | Description |
|--------|---------|-------------|
| `calc_surface_gravity()` | float [m/s²] | $GM/R^{2}$. |
| `calc_escape_velocity()` | float [m/s] | $\sqrt{2GM/R}$. |
| `calc_mean_density()` | float [kg/m³] | $M / (\tfrac{4}{3}\pi R^{3})$. |
| `calc_equilibrium_temperature(F)` | float or array [K] | $\left[(1-A)\,F/(4\varepsilon\sigma)\right]^{1/4}$ (fast rotator; $F$ = insolation flux [W/m²], float or array). |
| `set_spin_frequency(ω)`, `set_obliquity(θ)` | - | Set the rotation rate [rad/s] or obliquity [rad], which a `System` reads to build the tidal state. A solved tidal result is left alone. |
| `copy()` | world | An independent, unsolved world of the same class and state (see [Binary Serialization](#binary-serialization)); `copy.copy`, `copy.deepcopy`, and `pickle` use it, so a world can go to a process pool. |
| `summary()` | str | The world, then one row per layer: name, radii [km], state, material phases, solved density range, temperature, and shear rheology. |

`repr(world)` is one line: `TerrestrialWorld('Io', radius_km=1821.49, mass_kg=8.9298e+22, num_layers=3)`.

### Layers

A world holds an inner-to-outer stack of [layers](../layers/layer.md). A layer built without a `layer_index` or `radius_inner` takes them from its place in the stack; the positional form `Layer("core", 0, 0.0, 3.485e6, ...)` gives both, and `add_layer` checks them.

| Member | Description |
|--------|-------------|
| `add_layer(layer)` | Add a layer, inner to outer; the world then owns it and its models, and the passed object becomes a view of it. Raises `ValueError` (leaving the layer as built) if it was already added, if a given index is not its place, if it does not continue the stack (the innermost starts at 0), if it reaches past the world radius, or if its name is taken. |
| `num_layers` | Number of layers. |
| `calc_total_mass()` | Sum of the layer masses [kg]; equals `planet_mass_eos` after a successful EOS solve. |
| `calc_internal_heating(time)` | Radiogenic heating [W] of every layer with a radiogenics model at a time [s], whatever its `use_heating`, from each layer's `mass` (so solve the EOS first if they have none). |
| `validate_layers()` | `True` if every boundary is continuous and the innermost starts at 0. |
| `world[i]`, `world.get_layer(i)`, `world.<layer_name>` | The layer at index `i` (0 = innermost, negative from the outermost) or by name (`world["mantle"]`, `world.mantle`). |
| `world[a:b]`, `world.layers`, `for layer in world` | Lists and iteration, inner to outer; `len(world)` is the layer count (a world is always true). |

Every world method that takes a layer takes its index or name: an index out of range raises `IndexError`, an unknown name `KeyError` (listing the names, with the closest), and anything else `TypeError`.

Each is a non-owning view with the full `Layer` API (`world.mantle.temperature = 1700.0`). A change through a view reaches the world, which forgets its solved structure when the solve reads that setting (see [Changes That Clear a Solve](../layers/layer.md#changes-that-clear-a-solve)). A view keeps the world alive, and repeated access returns the same object until a layer is added. Name access runs after normal attribute lookup (defined members win) and ignores names starting with `_`.

## Equation of State

Once every layer has a [material](../layers/layer.md#material) (`all_materials_set`), `solve_eos(...)` integrates the radial structure from center to surface. It fills every layer's profile (density, gravity, pressure, temperature, heat flow, material state) and sets each layer's `mass` (and `density_bulk`) to the mass the solved profile places in it, replacing the value it was built with (0.0 from a file that gives none, the usual case).

```python
result = world.solve_eos(surface_pressure=0.0)  # EOSResult: a dict of profile arrays and scalars
rho    = world.get_density(5.0e6)               # [kg/m³] at radius 5000 km
g      = world.get_gravity(world.radius)        # Surface gravity [m/s²]
p0     = world.get_pressure(0.0)                # Central pressure [Pa]
print(result["success"], world.planet_mass_eos) # True and the solved mass [kg]
```

A secant iteration on the surface-pressure mismatch finds the central pressure in a few whole-structure integrations (where a large, soft planet's surface pressure first falls, the step doubles until the root is bracketed, then closes as in Brent's method). When every density depends only on radius and temperature (`constant` and `interpolate` equations of state, no melt-density mixing in a pressure-dependent melting range), a solve from scratch takes two integrations instead of three or more. The integration is non-dimensional (planet radius, bulk density, and $1/\sqrt{\pi G \rho}$), so the tolerances mean the same for every planet; results are in SI.

**`solve_eos(*, surface_pressure=0.0, slices_per_layer=None, G_to_use=None, integration_method=None, rtol=None, atol=None, pressure_tol=None, max_iters=None, nondimensionalize=None, temperature=None, solve_temperature=None, surface_temperature=None, reset_layer_masses=False, verbose=False, time=None, max_thermal_passes=None, thermal_tol=None, raise_on_fail=False) -> EOSResult`**

A setting left as `None` takes the world's [pinned value](#other-world-members), else the `[eos_solver]` [configuration](../../Overview/2_TidalPy_Configurations.md) value (the standalone `radial_solver`'s defaults); `G_to_use` defaults to the configured constant. `pressure_tol` is relative to the central-pressure scale $(2/3) \pi G \rho^2 R^2$ and must stay above `rtol`, the integrator's noise on the surface pressure. `temperature` gives every layer one temperature \[K\]. `solve_temperature`, `surface_temperature`, `time`, `max_thermal_passes`, and `thermal_tol` belong to the [thermal solve](#temperature-and-heat-flow), `reset_layer_masses` to [Layer Size](#layer-size). `raise_on_fail=True` raises `SolutionFailedError` (a `RuntimeError`) instead of returning `success = False`. See [`solve_eos` Result](worlds.md#solve_eos-result) for the keys.

### Solved State

A re-solve starts from the last converged solve, so after a small change (a layer temperature, a later time) it usually converges in one or two passes. Results then depend on the call history at the tolerance level (a world solved at two surface pressures can differ from a fresh one by a few parts in $10^9$ in central pressure and about one in $10^9$ in $k_2$), so build a fresh world when results must reproduce bit for bit.

A solve commits only when it finishes, and a solution released by `release_radial_solution` keeps the (immutable) materials it was solved with, so a later solve never changes it. Changing anything the solve reads clears every solved result: the layer stack (`add_layer`, `load_binary`), `set_prescribed_heating`, and the layer settings in [Changes That Clear a Solve](../layers/layer.md#changes-that-clear-a-solve). The profile getters then return NaN, `eos_solved` is `False`, `zones` is empty, and a Love solve or `rheology` `calc_tides` raises until `solve_eos` runs again. Love solves read `state`, `is_static`, `is_incompressible`, and the rheologies afresh, so changing those keeps the structure. A failed solve leaves the world unsolved, with only `success`, `message`, and NaN profiles.

A world runs one heavy call at a time (its solves, `release_radial_solution`, `load_binary`, and the setters of its models, settings, and layers), so a change waits for a running solve and a call on an array of radii reads one solve. These calls release the GIL, so separate worlds run in parallel. A solve and a later read are separate turns, so give each thread its own world (see [World-Attached Solves](../../RadialSolver/parallel.md#world-attached-solves)).

### Pieces and Zones

The solve integrates each layer in pieces, each ending at a radius (the top of a layer that holds its volume, or the end of a thermal segment), an enclosed mass (the top of a layer that [holds its mass](#layer-size)), or a change of state.

A layer can change state when its `state` is `"auto"`, `use_melting` is on, and its material has both phases (`Layer.can_change_state`). It is liquid where its post-melt rigidity $\mu / (\bar{\rho} g R)$ (with the world's stated bulk density, surface gravity, and radius) is at or below `[numerical] minimum_solid_rigidity` (10$^{-6}$ by default), as a fully molten material always is, and solid elsewhere. A partially molten material therefore stays solid until its shear modulus is far too low to integrate as a solid. Each crossing is located as a root to the integration tolerance and the integration restarts there, so the boundaries do not depend on `slices_per_layer`.

Consecutive pieces of one layer in one state form a zone. `zones` lists them inner to outer as dicts with `layer`, `radius_inner` and `radius_outer` \[m\], `mass_inner` and `mass_outer` (enclosed mass) \[kg\], and `state` (`"solid"` or `"liquid"`). A layer that cannot change state is one zone, liquid if `Layer.is_liquid`. `molten_regions` lists the liquid zones melting made in layers not liquid throughout, as `(layer_name, radius_inner, radius_outer)`, each also named in an info-level log message. Both are empty before a solve.

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

The Love solver takes each zone as a layer with its layer's `is_static` and `is_incompressible`, a liquid zone with the liquid equations, so a melting mantle (`is_static` by default) is a static liquid where molten. A zone thinner than `[numerical] minimum_zone_fraction` (10$^{-7}$ of the world radius by default, about 0.6 m for Earth) takes its thicker neighbor's state, since the radial solver cannot integrate across it. Treating a solid of rigidity $10^{-6}$ as a liquid changes the Love numbers by about that fraction. Elsewhere the layer is one layer, and its 3D heating takes nothing from its liquid zones. Each zone reads the material on its own side of a zone boundary, so a static liquid zone's interface conditions take the liquid's density there, as they would for a liquid layer. The profile getters give a radius on the boundary the lower zone's values.

A layer forced to `"solid"` or `"liquid"` after the solve is one zone in that state to the next Love solve, while `zones` keeps the solved states until the next `solve_eos`. A change between `"auto"` and a forced state that decides whether the layer can change state makes the world forget its solve.

> [!NOTE]
> Only a layer that can change state is split. A solid layer given a near-zero shear modulus some other way (a material without melting, or `use_melting` off) is solved as a solid, and the solver may fail on it or return $\mathrm{Im}[k]$ with the wrong sign and only a conditioning warning.

### Layer Size

A layer holds its volume by default. One with `is_volume_fixed = False` holds its mass: the integration ends it where the enclosed mass reaches its mass, inside one solve. The layers above move with it, keeping their volumes (unless they also hold their mass), and the world radius becomes the top of the outermost layer, but only if the solve succeeds. Its thermal segments end at fixed fractions of its mass, so its boundary layers move with it.

The mass held is the layer's `mass` (`mass_kg` in a file) when positive; otherwise the first solve keeps its volume and it holds the mass found there. `solve_eos(reset_layer_masses=True)` lets the current boundaries set the masses again. A layer that holds its mass reads NaN above its top.

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

### When the EOS Solve Fails

- **No hydrostatic structure:** `max_iters` reached (plus one last pass) with the surface pressure off by more than `pressure_tol`. Returns `success = False` with "no hydrostatic structure", sets `max_iters_hit`, and leaves the world unsolved.
- **Mass far off:** a converged solve fails when its mass differs from the stated mass by more than the factor `[numerical] maximum_eos_mass_ratio` (default 10) either way, as when layers with no hydrostatic structure near that mass (a core much too dense for its radius) meet the surface pressure on a collapsed branch at an absurd central pressure. The message says how far off.
- **Mass slightly off:** more than 1 percent between `planet_mass_eos` and `mass` warns once per world. Love numbers, tides, and moment of inertia follow the solved structure; a `System` orbit, `calc_surface_gravity`, and `calc_mean_density` use the stated mass. Adjust materials, radii, or the stated mass until they agree.
- **Tension:** a layer in tension past what its pressure law (Birch-Murnaghan or Vinet) represents fails the solve; the message names the layer and pressures. The law sees the pressure less the thermal pressure $\alpha_0 K_0 (T - T_\mathrm{ref})$ of a layer with `use_thermal_expansion`, so a hot layer with a large expansivity or $K_0'$ can get there. Its density is then held at the smallest compression, with a near-zero bulk modulus.
- **Over-compression:** past the law's compression end (a Birch-Murnaghan $K_0'$ below 4 turns over) the density is held at the largest compression; the solve warns and stands.
- **Thermal passes not settled:** warns and keeps the last pass (see [Thermal Passes](#thermal-passes)).

## Calculating Love Numbers

`solve_love_numbers(...)` calculates the tidal or loading Love numbers $k$, $h$, $l$ at one forcing frequency; `calc_love_numbers` sweeps many. `solve_eos` must run first. `love_method` (case-insensitive; aliases in parentheses) selects the method:

| `love_method` | What it does |
|---|---|
| `radial_solver` (`shooting`, `rs`; default) | Integrates the radial ODEs from center to surface, for any multi-layer, solid/liquid, static/dynamic, compressible/incompressible world. |
| `propagation_matrix` (`prop_matrix`, `pm`, `prop`) | Quasi-analytic matrix propagation for one solid, static, incompressible layer; other worlds fail gracefully (`love_success` `False`, `love_error_code` non-zero). `core_model` selects the core starting condition. |
| `homogeneous` (`homogen`) | Each tidal layer as a homogeneous incompressible planet, summed with tidal-scale weights (see [Quasi-Homogeneous Love Numbers](#quasi-homogeneous-love-numbers)). |
| `cpl` | The same on the static (unrelaxed) averaged shear modulus, times $(1 - i/Q)$, so $-\mathrm{Im}[k] = \mathrm{Re}[k]/Q$; $Q$ is `fixed_q` (argument or `[tides]`), else the tide model's for the degree. |
| `ctl` | Like `cpl` with a constant time lag, $(1 - i\,\omega\,\Delta t)$; $\Delta t$ is `fixed_dt` or the tide model's. |
| `laterally_inhomogeneous` (`3d`, `lat_inhom`) | Reserved for a 3D solver; raises `NotImplementedError`. |

The default method comes from `set_tide_config(love_method=..., love_fixed_q=..., love_fixed_dt=...)` or the `[tides]` keys `love_method`, `love_fixed_q`, `love_fixed_dt_s`.

```python
from TidalPy.Material import Material, Phase
from TidalPy.Structures.worlds import BaseWorld

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

The moduli and viscosity belong to the material; a rheology holds only model parameters (the Andrade exponent, the Voigt fractions). A solid layer with no shear modulus has no strength, and the solve fails.

**`solve_love_numbers(frequency=1e-5, degree_l=2, solve_for='tidal', core_model=0, starting_method=None, nondimensionalize=None, starting_radius=0.0, start_radius_tol=None, integration_method=None, rtol=None, atol=None, scale_rtols=None, max_num_steps=None, expected_size=None, max_ram_MB=None, max_step=0.0, verbose=False, warnings=True, love_method=None, fixed_q=None, fixed_dt=None, raise_on_fail=False, love_only=False) -> dict`**

A setting left as `None` takes the `[radial_solver]` [configuration](../../Overview/2_TidalPy_Configurations.md) value (as the standalone `radial_solver` and the world's tidal solves do); `love_method`, `fixed_q`, and `fixed_dt` take the world's `[tides]` settings. Raises `ValueError` if the EOS is unsolved. Returns `success`, `error_code`, `message`, `love_method`, and `love_number_k/h/l`, which also stay on the world ([Love-Number Properties](#love-number-properties)). `raise_on_fail=True` raises `SolutionFailedError` (a `RuntimeError`) on failure. `max_step` is in meters.

`love_only=True` makes the shooting method faster by reading each layer's solution at its top, with no dense output. The Love numbers, surface values, and step counts are unchanged, but `get_love_radial_y` and a released solution's `result` and `get_radial_solution` raise `ValueError` until a solve without it. `calc_tides` and the 3D paths always keep the radial functions; the propagation matrix ignores the argument.

The radial solver integrates each [zone](#pieces-and-zones) as a layer, reading the structure and static moduli from the dense EOS solution at each radius (no interpolation between slices) and the complex moduli from the rheology, so its Love numbers converge with the integration tolerance alone and do not depend on `slices_per_layer`. The propagation matrix propagates across those slices, so its Love numbers _do_ depend on it.

`release_radial_solution()` returns the last radial solve as a `RadialSolverSolution`, its `result` on the solve's radius grid (`sample_radii()`); above the surface it reads NaN.

### Frequency Sweeps

**`calc_love_numbers(frequencies, **solve_love_numbers_kwargs) -> dict`** runs `solve_love_numbers` at each frequency \[rad s-1\] and returns `frequency` and, shaped like it, `k`, `h`, `l` (complex), `success`, and `message`. A failed solve gives NaN at its frequency and the sweep goes on; only bad input (an unsolved EOS, a frequency outside the allowed range) raises. The frequency-independent setup is reused, so only the complex moduli and the integration are redone per frequency.

```python
import numpy as np

periods = np.logspace(0, 3, 30) * 86400.0             # 1 to 1000 days [s]
sweep = world.calc_love_numbers(
    2.0 * np.pi / periods,
    degree_l=2)                                       # One Love solve per frequency
neg_imag_k2 = -sweep["k"].imag                        # NaN where a solve failed
```

### Quasi-Homogeneous Love Numbers

The `homogeneous`, `cpl`, and `ctl` methods give each layer the homogeneous-sphere formulas, $k_{l} = \frac{3}{2(l-1)}\,\frac{1}{1+\bar{\mu}_{l}}$, $h_{l} = \frac{2l+1}{2(l-1)}\,\frac{1}{1+\bar{\mu}_{l}}$, $l_{l} = \frac{3}{2l(l-1)}\,\frac{1}{1+\bar{\mu}_{l}}$ with $\bar{\mu}_{l} = \frac{2l^{2} + 4l + 3}{l}\,\frac{\mu}{\rho g R}$, where $\mu$ is the layer's rheology at the forcing frequency applied to its averaged moduli and viscosity, and $\rho$, $g$, $R$ are the planet's EOS bulk density (from the solved mass, not the declared one), EOS surface gravity, and radius. They keep some layered structure without a radial solve, so they are much faster.

Each tidal layer's post-melt moduli are volume averaged and its post-melt viscosity log-volume averaged (so a viscosity spanning decades is averaged in log space). A homogeneous planet of that material, with the planet's radius, bulk density, and surface gravity, has Love numbers $k_i$; with the layer's tidal scale $s_i$,

$$k = \sum_i s_i k_i$$

(likewise $h$ and $l$). The tidal scale is the layer's `tidal_scale`, else its volume over the planet's, and 0 with `use_tides` off. With default scales a one-layer planet gets exactly the homogeneous-sphere value, and a small, highly dissipative layer (an asthenosphere) adds dissipation in proportion to its volume instead of making the whole planet dissipate like it. Heating is linear in $-\mathrm{Im}[k]$, so each layer's heating is that of its own term $s_i k_i$, constant within the layer, and the layers sum to the total.

`love_layer_parts` lists each layer's $s_i$, $k_i$, $h_i$, $l_i$, and complex shear modulus. These methods have no radial functions (the radial-function getters return NaN), so the 3D path (`calc_3d_tides`, `get_3d_tidal_heating`) raises `RuntimeError` after one. Free-function versions are in `TidalPy.Tides.love` ([Love numbers](../../Tides/love/love_numbers.md)).

### Input Checks

Every Love solve raises `ValueError` for:

- Degree below 2 in a tidal or free-surface solve (degree 1 is a translation of the body). A loading-only solve may use degree 1, in the reference frame `degree1_frame` (see [Degree-1 Load Love Numbers](../../RadialSolver/calculating_love_numbers.md#degree-1-load-love-numbers)).
- A loading or free-surface solve with `homogeneous`, `cpl`, or `ctl`, which give tidal Love numbers only.
- A frequency that is not finite and positive, or outside `[numerical]` `minimum_frequency` to `maximum_frequency` (the tide paths instead solve a mode at $\sqrt{\omega^2 + \omega_c^2}$, with $\omega_c$ the world's continuation frequency, `calc_continuation_frequency()`, and scale its dissipation by $|\omega|$ over that; see [Numerical Settings](../../Overview/2_TidalPy_Configurations.md#numerical-settings)). The Love numbers at $-\omega$ are the complex conjugates of those at $\omega$.
- `start_radius_tol` outside (0, 1), or a negative `starting_radius` or `max_step`.

## Global (1D) Tidal Dissipation

`calc_tides(...)` sums the active tidal modes (the global or "1D potential" approach) for the total tidal heating and three orbital potential derivatives. The tide model ([Global Tidal Dissipation](../../Tides/global_tides.md)) supplies each mode's $-\mathrm{Im}[k_{l}]$: the analytic models (`cpl`, `ctl`, `ctl_q`) from fixed per-degree parameters, with no interior, and the `rheology` model from the world's complex Love numbers at each unique tidal frequency, keeping every mode's `k`, `h`, `l` (`get_tidal_love_k(l, m, p, q)`). It raises `RuntimeError` if `rheology` has no solved EOS or a per-frequency radial solve fails.

It returns `tidal_heating` \[W\], `dU_dM`, `dU_dw`, `dU_dO` \[J kg$^{-1}$ rad$^{-1}$\], `num_tidal_modes`, and `layer_tidal_heating` (\[W\] by layer name). These stay on the world ([Tide Members](#tide-members)), and the layer heating becomes the [tidal heat source](#heat-sources).

```python
from TidalPy.Tides.classes import make_tide

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

For a synchronous, low-eccentricity body the `cpl` result reproduces the standard CPL rate $\frac{21}{2}\,\frac{k_{2}}{Q}\,\frac{G M_{h}^{2} R^{5} n e^{2}}{a^{6}}$ (host mass $M_{h}$).

The builder sets the tide model from the `[tides]` table, by default `fixed_q` for a star, `fixed_dt` for a gas giant, and `rheology` for a terrestrial world. The table may hold every analytic model's per-degree lists; the model gets only those it reads (`tide_config_keys`), so a `fixed_dt` world ignores `fixed_q`. `set_tide_model` also takes a name (`"fixed_q"`, at its `[tides]` defaults) or a `[tides]`-style table (the model under `model` or `global_tidal_model`, plus lists and settings, which it applies). `world.set_tide_config(**world.get_tide_config())` restores a configuration, and `temporary_tide_config` restores it after a `with` block, also when the block raises:

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

A new tide configuration clears the last tidal result, so read it inside the block.

### Orbital State

Each of the six orbital-state arguments (`orbital_frequency`, `spin_frequency`, `eccentricity`, `obliquity`, `semi_major_axis`, `host_mass`, in that order) left as `None` comes from `get_tide_state()`, the state the world's `System` gives it (its orbit about its tidal host, the host's mass, its spin and obliquity). So `world.calc_tides()` solves a system's world as it is, and `world.calc_tides(eccentricity=0.05)` changes one value. Without a system or tidal host, a missing argument raises `ValueError` naming it. The 3D methods take the state the same way, with their point and grid arguments by keyword.

A tide solve (in a `System` evolution too) warns once per world when the spin is within $10^{-3}$ of the orbital mean motion but not equal to it. The slow forcing at the difference frequency can change the heating by orders of magnitude: a bundled Io spinning 0.016 percent off the mean motion of its rounded semi-major axis heats about four times more than a synchronous one. Set a synchronous spin exactly (`world.set_spin_frequency(system.calc_orbital_frequency(world))`).

### Tidal Heating of Each Layer

| Love numbers from | Layer heating |
|-------------------|---------------|
| The radial solver (`radial_solver`, `propagation_matrix`) | The volume integral of the orbit-averaged heating density over the layer (the `calc_3d_tides` integral, every axis summed), scaled so the tidal layers sum to the 1D total. Liquid zones and layers with `use_tides` off (purely real moduli) take none. It costs about another global solve; `layer_tidal_heating=False` (in `set_tide_config` or `[tides]`) skips it, leaving every layer NaN. |
| The quasi-homogeneous methods (`homogeneous`, `cpl`, `ctl`) | The heating of the layer's own term $s_i k_i$, constant within the layer; the layers sum to the total. |
| An analytic tide model (`cpl`, `ctl`, `ctl_q` tide models, which describe the whole body) | The total times the layer's tidal scale $s_i$. |

`layer.get_tidal_heating()` gives the same value as `get_layer_tidal_heating`.

## Temperature and Heat Flow

Each layer has its own temperature (`Layer.temperature`), and its [cooling model](../../Cooling/cooling_models.md) sets how heat moves inside it. From these the solve calculates a temperature profile, the heat flow through every radius, and each layer's temperature rate. The outermost layer radiates to `surface_temperature` \[K\]; left out, no heat leaves the world.

| Model | Profile inside the layer | Where its temperature applies |
|---|---|---|
| `off`, or no model | Isothermal, with no modeled gradient, so the layer conducts perfectly. | Everywhere |
| `conduction` | Two conducting halves, $T = T_0 - (L / 4 \pi k)(1/r_0 - 1/r)$. | The mid-radius |
| `convection` | A conducting boundary layer at the base and top, sized by the model's Nusselt scaling, around an adiabatic interior. | The top of the interior |

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

Materials are evaluated at the solved temperature and pressure, so an Arrhenius viscosity stiffens where the profile is cold, a melting layer melts where it crosses the solidus, and with `use_thermal_expansion` the density follows it.

### Convecting Layers

A convecting interior follows the adiabat

$$\frac{dT}{dr} = -\frac{(\alpha + \alpha_L)\, g\, T}{c_p}$$

with the material's expansivity $\alpha$ and heat capacity $c_p$ at the local pressure, and $\alpha_L$ the latent heat's share of the expansivity inside a melting range whose curves follow the pressure (`latent_expansion` of `Material.calc_state`; zero outside the range or with `use_pressure_melting` off). The layer's temperature applies at the top of the interior, under the upper boundary layer (the upper-mantle temperature of parameterized convection; Stevenson et al. 1983; Schubert et al. 2001), and the interior warms downward to $T \exp\left(\int (\alpha + \alpha_L) g / c_p \, dr\right)$ at its base. A layer whose base carries no heat (the innermost, or one above a layer outside the network) has no lower boundary layer.

The Rayleigh number's viscosity is taken at the layer's temperature and the solved pressure near the top of the interior (`layer_reference_pressure`, `layer_reference_viscosity`, `layer_reference_melt_fraction`; see [Convective Reference State](../../Cooling/cooling_models.md#convective-reference-state)), and gravity, density, and the thermal constants at mid-layer. An interior liquid there (fully molten, or past the rigidity threshold of [Pieces and Zones](#pieces-and-zones)) is a magma ocean (`layer_magma_ocean`) with the liquid scaling $\mathrm{Nu} = a_\mathrm{liquid} \mathrm{Ra}^{\beta_\mathrm{liquid}}$ (Solomatov 2000). A conducting layer reads its thermal properties at the mid-radius pressure and its own temperature.

The Rayleigh number uses the drop across both boundary layers: from the top of the layer below (the end of its adiabat, if it convects) to the layer's temperature, plus from that to the layer above or `surface_temperature`. A mantle as warm as the layer above therefore still convects over a hotter core. A layer whose base carries no heat has only the upper drop and boundary layer. The thickness $d / \mathrm{Nu}$ carries the flux across the whole drop, so two boundary layers take half each and a single one the whole, at most 40 percent of the layer each. With equal drops the flux through the top is the model's `cooling_flux`, and a sub-critical layer ($\mathrm{Nu} = 1$) is close to a `conduction` layer. A model that gives no usable thickness (usually a NaN viscosity) leaves each boundary layer at 40 percent, sets `layer_boundary_fallback`, and warns. [Thermal Profiles](../../Cooling/cooling_models.md#thermal-profiles) shows how each model builds its profile.

### Heat Balance and Capacity

The layers form a chain of thermal resistances. A conducting spherical shell between $r_a$ and $r_b$ has

$$R = \frac{1}{4 \pi k} \left( \frac{1}{r_a} - \frac{1}{r_b} \right)$$

and the heat flow through an interface is $L = \Delta T / R$ across the two resistances facing it. Both layer temperatures are inputs, so the flow into a layer generally differs from the flow out; the difference is stored or released, giving `calc_layer_temperature_rate(layer)`:

$$\left(C + C_\mathrm{latent}\right) \frac{dT}{dt} = L_\mathrm{in} - L_\mathrm{out} + H$$

with $H$ \[W\] the heat generated inside (`layer_heating`, see [Heat Sources](#heat-sources)), $C$ \[J K$^{-1}$\] the heat its profile stores per kelvin of its temperature (`calc_layer_thermal_capacity`), and $C_\mathrm{latent}$ \[J K$^{-1}$\] the latent heat of its zone boundaries (`calc_layer_latent_capacity`). With the neighbors' interface temperatures and the heating held, a change $\delta T$ moves the profile by $S(r)\,\delta T$, so

$$C = \int \rho\, c_p\, S \, 4 \pi r^2 \, dr$$

over the layer, with $c_p$ the effective heat capacity at the solved pressure and temperature. $S$ is one through an isothermal layer and $T(r)/T$ along a convecting interior, which scales with the temperature at its top (Stevenson et al. 1983). Across a conducting stretch the change solves Laplace's equation whatever the heating, so $S$ runs linearly in $1/r$ from zero at the interface it meets to one at the end the layer holds ($T_\mathrm{base}/T$ at the base of a convecting interior); below an insulated base it is one. A convecting mantle therefore stores more than $M c_p$ per kelvin of its upper-mantle temperature (about 1.3 times for the bundled Earth). The integral is split where the profile crosses a melting curve.

A material that melts over a range carries its latent heat in its effective heat capacity, so only the part inside the range stores it. One that melts at a single temperature (solidus and liquidus the same curve) has no range, so the boundary between its solid and liquid zones carries it: as the layer's temperature changes, the boundary moves and melts or freezes mass. With $G(r) = T(r) - T_m(P(r))$ zero at the boundary $r_b$,

$$C_\mathrm{latent} = \frac{L\, 4 \pi r_b^2\, \rho_\mathrm{solid}\, S}{\left| dG/dr \right|}$$

with $L$ the latent heat \[J kg$^{-1}$\] and $S$ the sensitivity of the boundary temperature to the layer's (one if isothermal, $T(r)/T$ on an adiabat). It assumes the neighbors' interface temperatures hold and the new melt takes the solid's density at the boundary; it is zero without such a boundary. Both capacities are computed on the first call after a solve (their quadratures cost more than an isothermal solve).

Where neither side of an interface has a resistance (an `off` layer under another, or an `off` outermost layer under the surface), nothing holds a temperature contrast: the lower layer stores no heat, passing on what enters plus what it generates, and the interface keeps its temperature. An `off` outermost layer thus loses all of it through the surface, its rate is zero, and its profile ends at its own temperature, not `surface_temperature`.

A layer whose temperature is not positive (such as the 0 K default) is outside the thermal network (`layer_in_thermal_network` is `False`): neither sink nor source, isothermal at its placeholder, and its neighbors keep their own temperatures at the shared interfaces. If that leaves its material at the cold, rigid limit of its viscosity law, the first solve warns, naming it, once until its temperature is set again.

### Thermal Passes

A world whose layers share one temperature and generate no heat has no profile to integrate, so the solve costs and returns what it would with `solve_temperature=False`, and the profile queries report each layer's temperature. Otherwise temperature and heat flow become two more state variables and the solve iterates: an isothermal first pass, then passes that integrate the profile and relax the boundary layers, interface temperatures, and heat flows against it. `thermal_passes` counts them, and `thermal_converged` says whether the largest relative change in interface temperatures and heat flows between two passes fell below `thermal_tol`. `max_thermal_passes` caps them; both are `[eos_solver]` settings a call or the world may override. A solve that uses every pass without settling warns and keeps the last pass. Layers that hold their mass move in the same passes and end on the last pass's grid, so every reported slice and interface lies in its own layer.

### Heat Sources

In a solve that carries temperature, the world's three sources heat each layer with `use_heating` on. Each gives a volumetric heating $h$ \[W m$^{-3}$\], and the structure integrates

$$\frac{dL}{dr} = 4 \pi r^2 h$$

so heat flow grows through a heated layer and its conducting stretches bend: a uniformly heated conducting shell follows $T = B + A/r - h r^2 / 6k$. With $H(r)$ the heat generated between the base of a conducting stretch and $r$, the flow out of its top is the flow in plus $H$, and the drop across it is $L_\mathrm{base} R + \int H / (4 \pi r^2 k) \, dr$, so the solved profile still passes through every layer temperature and reaches `surface_temperature`. A heated world is a thermal solve even when its layers share one temperature. With `solve_temperature=False` the heating is ignored (a debug-log note).

| Source | Heating | Set by |
|---|---|---|
| Radiogenic | The [radiogenics model](../../Radiogenics/radiogenics_models.md)'s $\epsilon$ \[W kg$^{-1}$\] at the solve's `time` \[s\] (`None`: each model's reference time), with $h = \epsilon \rho$, exact even when the layer's mass is a solve output. | `Layer.radiogenics` |
| Tidal | What the last `calc_tides` put in each layer \[W\]: with a radial-solver Love method, along the radial profile of the per-layer heating integral (piecewise linear in radial fraction, renormalized to the layer's total); otherwise by mass. | `calc_tides` |
| Prescribed | A power \[W\] spread by mass (reached to the thermal tolerance whatever mass the solve gives), or a specific rate \[W kg$^{-1}$\]. | `set_prescribed_heating` |

The result reports each layer's heating by source (`layer_heating_radiogenic`, `layer_heating_tidal`, `layer_heating_prescribed`) and summed (`layer_heating`), all zero without `use_heating` or temperature.

| Member | Description |
|---|---|
| `set_prescribed_heating(layer, power=None, specific_rate=None)` | A `power` \[W\] or `specific_rate` \[W kg$^{-1}$\]; neither clears it. Both, an infinite value, or an unknown layer raise `ValueError`. Clears the solved structure. |
| `prescribed_heating` | By layer name, `{"power": watts}` or `{"specific_rate": watts_per_kg}`. |
| `tidal_heat_source`, `clear_tidal_heating()` | Each layer's tidal heating \[W\] from the last `calc_tides` (empty before one), and forgetting it. |
| `get_heating(radius)` | Volumetric heating \[W m$^{-3}$\] of every source: zero without `use_heating`, NaN before a successful solve or outside the world. |
| `calc_layer_temperature_rate(layer)` | $dT/dt$ \[K s$^{-1}$\] from the last solve and the latest `calc_tides`; NaN before a solve. |
| `calc_layer_thermal_capacity(layer)`, `calc_layer_latent_capacity(layer)` | $C$ and $C_\mathrm{latent}$ \[J K$^{-1}$\]; NaN before a solve. |

The tidal source depends on the interior through the tides, so it survives a later `solve_eos`; a failed `calc_tides`, `clear_tidal_heating`, `add_layer`, and `load_binary` forget it. A layer with `use_tides` off is elastic, so its heat is in neither the total nor any other layer. With `layer_tidal_heating` off a radial-solver tide spreads the total over the tidal layers by mass. The tidal source is not saved with the world; the prescribed heating is (see [Binary Serialization](#binary-serialization)).

An evolution step is `solve_eos`, `calc_tides`, then the temperature rates, which take the tides just computed with no second solve:

```python
import math

from TidalPy.Tides.classes import make_tide

world.set_tide_model(make_tide("rheology"))              # Love numbers from the radial solver
world.set_tide_config(
    max_degree_l=2,
    eccentricity_truncation=2,
    obliquity_truncation=0)
orbital_frequency = 2.0 * math.pi / (1.769 * 86400.0)    # [rad s-1]

world.solve_eos(surface_temperature=110.0)               # The structure first,
world.calc_tides(
    orbital_frequency=orbital_frequency,
    spin_frequency=orbital_frequency,
    eccentricity=0.0041,
    obliquity=0.0,
    semi_major_axis=4.217e8,
    host_mass=1.898e27)                                  # then its tides,
print(world.tidal_heat_source)                           # [W] by layer, about 8.3e11 in the mantle
print(world.calc_layer_temperature_rate("mantle"))       # then the rates [K/s], tides included

world.set_prescribed_heating(
    "mantle",
    specific_rate=1.0e-11)                               # [W/kg] beside the radiogenic and tidal heat
result = world.solve_eos(surface_temperature=110.0)
print(result["layer_heating_tidal"], result["layer_heating_prescribed"])
print(world.get_heating(1.5e6))                          # [W/m^3] every source summed
```

## Reference

### `solve_eos` Result

An `EOSResult` (`TidalPy.Structures.worlds.EOSResult`) is a `dict` subclass that behaves as a plain `dict` (item access, equality, copying, pickling). Its `repr` is a short summary (`success`, `iterations`, `message`, planet mass, radius, central pressure, and the keys). It contains:

- `success`, `message`, `iterations` (central-pressure steps), `structure_integrations` (whole-structure integrations over every thermal pass, the measure of cost), `max_iters_hit`, `pressure_error` \[Pa\].
- Profile arrays: `radius`, `gravity`, `pressure`, `mass`, `moi`, `density`, `temperature`, `heat_flow`.
- Scalars: `surface_gravity`, `surface_pressure`, `central_pressure`, `planet_mass`, `planet_moi`.
- [`zones`](#pieces-and-zones), and the thermal report `thermal_passes` and `thermal_converged`.
- Per-layer lists, inner to outer:
  - `layer_radius_outer` \[m\]: where each layer ended (it moves for a layer that holds its mass).
  - `layer_temperature`, `layer_node_temperature` (at the interface above), `layer_top_temperature` and `layer_base_temperature` (the ends of a convecting interior, whose top is the layer's temperature) \[K\].
  - `layer_heat_flow_in`, `layer_heat_flow_out`, `layer_heating`, `layer_heating_radiogenic`, `layer_heating_tidal`, `layer_heating_prescribed` \[W\].
  - The cooling model's profile: `layer_boundary_thickness` \[m\], `layer_rayleigh_number`, `layer_nusselt_number`, `layer_boundary_fallback`, `layer_magma_ocean`, `layer_reference_pressure` \[Pa\], `layer_reference_viscosity` \[Pa s\], `layer_reference_melt_fraction`.
  - `layer_in_thermal_network`: whether the layer has a temperature of its own.

### Profile Queries

After a successful solve each getter takes a radius \[m\] (float or `np.ndarray`) and returns the same shape, NaN where unsolved, from the solver's dense interpolants (not by re-interpolating the gridded arrays).

| Member | Returns |
|--------|---------|
| `get_density`, `get_gravity`, `get_pressure`, `get_shear_modulus`, `get_bulk_modulus`, `get_shear_viscosity`, `get_bulk_viscosity`, `get_melt_fraction`, `get_static_viscoelastics`, `get_state`, `calc_complex_shear_modulus(r, ω)`, `calc_complex_bulk_modulus(r, ω)` | The layer [Profile Getters](../layers/layer.md#profile-getters), read from the layer that contains `r`. |
| `get_temperature`, `get_heat_flow`, `get_heating` | Temperature [K] on the solved profile, heat [W] outward through the sphere of radius `r`, and volumetric heating [W/m³]. |
| `eos_solved`, `all_materials_set` | Whether profiles are populated, and whether every layer has a material. |
| `surface_gravity_eos`, `central_pressure`, `planet_mass_eos`, `planet_moi_eos` | Scalars of the last solve (NaN if unsolved). |
| `zones`, `molten_regions` | See [Pieces and Zones](#pieces-and-zones). |

`calc_density`, `calc_gravity`, `calc_pressure`, `calc_shear_modulus`, `calc_bulk_modulus`, `calc_shear_viscosity`, `calc_bulk_viscosity`, `calc_static_viscoelastics`, and `calc_state(radius, force_recalc=False)` solve the EOS first when it is unsolved or `force_recalc=True`.

### Love-Number Properties

These describe the last `solve_love_numbers` (or `solve_love_numbers_supplied`) call. The Love solves of `calc_tides` and the 3D paths leave them (and `release_radial_solution`) untouched: after only a `calc_tides` they read as unsolved.

| Property | Type | Description |
|----------|------|-------------|
| `love_solved` | bool | `True` while a successful solve is held. A later `solve_eos` clears it, and every Love number and radial-function getter returns NaN until the next solve. |
| `love_success` | bool | Whether the last solve converged and still describes the structure. |
| `love_error_code` | int | 0 = success, < 0 = failure, -100 = no Love solve since the world was built or last solved its EOS. |
| `love_message`, `love_method` | str | The solver message, and the canonical name of the method used. |
| `love_num_ytypes` | int | Number of solution types (boundary-condition models requested). |
| `love_number_k`, `love_number_h`, `love_number_l` | complex | Love numbers of the first boundary condition (`get_love_number_k(0)` and friends). |
| `love_q_k`, `love_lag_k` | float | $Q = -s\,\lvert k \rvert / \mathrm{Im}(k)$ and lag $\arctan_2(-s\,\mathrm{Im}(k), \lvert\mathrm{Re}(k)\rvert)$ \[rad\] of `love_number_k`, $s$ the sign of $\mathrm{Re}(k)$ (as `RadialSolverSolution.Q_k` and `lag_k`). $Q$ is infinite and the lag 0 for an elastic $k$; both NaN when $k$ is. |
| `love_surface_amplification` | float | Conditioning of the surface boundary-condition solve, recorded on every shooting solve; near 1 is healthy, 0 after an analytic solve. |
| `love_surface_rcond` | float | Reciprocal condition number of the surface boundary-condition system: near machine precision, the solution constants are undetermined (the standalone `surface_solve_rcond`). |
| `love_surface_frame_residual` | float | How far a degree-1 loading solve leaves the surface condition its reference frame replaced (the standalone `surface_frame_residual`); NaN otherwise. |
| `love_effective_shear_modulus`, `love_tidal_volume` | complex, float | Tidal-scale-weighted mean complex shear modulus \[Pa\] and the participating volume \[m$^3$\] of the last quasi-homogeneous solve; NaN after a radial-solver solve. |
| `love_layer_parts` | list of dict | Per tidal layer of the last quasi-homogeneous solve: `layer`, `tidal_scale`, `love_number_k`, `love_number_h`, `love_number_l`, `shear_modulus`; empty after a radial-solver solve. |
| `get_love_number_k/h/l(ytype_idx=0)` | complex | Love numbers for boundary-condition model `ytype_idx` (0 = first requested, usually tidal). |
| `get_love_surface_y(ytype_idx, y_idx)` | complex | Radial function y₁…y₆ (`y_idx` 0–5) at the surface. |
| `get_love_radial_y(radius, ytype_idx=0, y_idx=0)` | complex or array | The radial function at any radius \[m\], accurate between grid slices. NaN after a failed or analytic solve or below the starting radius; `ValueError` after a `love_only` solve. |

### Tide Members

| Member | Returns | Description |
|--------|---------|-------------|
| `set_tide_model(tide)`, `tide_model` | -, `TideBase` or None | Attach a model (a `TideBase`, name, or `[tides]`-style table; `None` detaches), or read it (setting `tide_model` calls `set_tide_model`). Worlds share models (models are not changed in place). |
| `tide_model_set` | bool | Whether a model is attached. |
| `set_tide_config(min_degree_l=None, max_degree_l=None, eccentricity_truncation=None, obliquity_truncation=None, layer_tidal_heating=None, eccentricity_exact_tolerance=None, love_method=None, love_fixed_q=None, love_fixed_dt=None, *, eccentricity_trunc_lvl=None, obliquity_trunc_lvl=None, love_fixed_dt_s=None)` | - | Change the given `[tides]` settings: `eccentricity_truncation` (a level or `"exact"`, with `eccentricity_exact_tolerance` its mode range; [Eccentricity Functions](../../Tides/Eccentricity.md)), `obliquity_truncation` (0 or `"off"`, 2, 4, `"gen"`; [Obliquity Functions](../../Tides/Obliquity.md)), the per-layer split, and the default [Love method](#calculating-love-numbers). The keyword-only names are the `get_tide_config` keys; both names of one setting raise `ValueError`. A NaN `love_fixed_q` or `love_fixed_dt` clears it. |
| `get_tide_config()` | dict | The settings under the builder's `[tides]` keys (`*_trunc_lvl`). |
| `calc_tides(orbital_frequency=None, spin_frequency=None, eccentricity=None, obliquity=None, semi_major_axis=None, host_mass=None)` | dict | The global tidal solve ([Orbital State](#orbital-state)). |
| `tides_solved` | bool | Whether a successful solve is held. A new tide model or configuration, or `solve_eos`, clears it, and the heating and potential getters (and each layer's `get_tidal_heating()`) return NaN until the next `calc_tides`; the [tidal heat source](#heat-sources) is kept. |
| `get_tidal_heating()`, `get_tidal_heat_flux()` | float | Total heating [W] and $\dot{E} / (4 \pi R^2)$ [W m⁻²]; NaN if unsolved. |
| `get_tidal_potential_derivatives()` | dict | `dU_dM`, `dU_dw`, `dU_dO` [J kg⁻¹ rad⁻¹]: by mean anomaly, argument of pericenter, and longitude of the node. |
| `get_tidal_dU_dM_minus_dw()` | float | The per-mode sum of `dUdM - dUdw`, to pass as `dU_dM_minus_dw` to `OrbitSolver.calc_de_dt`: at small eccentricity the separate sums nearly cancel in $\dot{e}$. NaN if unsolved. |
| `get_num_tidal_modes()` | int | Active modes summed. |
| `get_layer_tidal_heating(layer=None)` | float [W] or dict | A layer's heating ([Tidal Heating of Each Layer](#tidal-heating-of-each-layer)), NaN before a solve; with no layer, all by name. |
| `get_layer_tidal_scale(layer)` | float | The layer's `tidal_scale`, else its volume over the planet's; 0 with `use_tides` off. |
| `get_tidal_love_k(l, m, p, q)` | complex | Per-mode radial-solver `k_l` (rheology only; NaN for analytic models). |

### Other World Members

**Geometry.** `calc_surface_area(radius)`, `calc_volume_sphere(radius)`, and `calc_volume_shell(outer, inner)`, shared by every structure.

**Spin.** `set_spin_model(spin)` attaches a spin model. `get_moment_of_inertia()` returns the EOS-solved moment of inertia [kg m2], or the spin model's `moment_of_inertia_factor` $\times\,M R^{2}$ before a solve. `calc_spin_derivative(host_mass)` returns the spin rate of change [rad s-2] from the current tidal solution, and `calc_synchronous_spin(orbital_frequency)` the synchronous rate [rad s-1] ([Dynamics](../../Dynamics/dynamics.md)).

**Three-dimensional tides.** `get_3d_tidal_heating_array(...)` vectorizes `get_3d_tidal_heating` for a heating map (one of `radii` and `colatitudes` may be a scalar), `calc_3d_displacements(...)` returns the displacement grid, and `calc_3d_stress_strain(...)` the stress and strain grids. Like `calc_3d_tides`, they take `num_threads` (default 0, the logical processors minus 4). See [3D heating](../../Tides/multilayer_3d_heating.md).

**Pinned solver settings.** `set_solver_defaults(eos_solver=None, radial_solver=None)` pins `[eos_solver]` and `[radial_solver]` keys on this world (a world file's tables of those names land here), and `get_solver_defaults()` returns them. A call's argument beats a pinned key, which beats the TidalPy configuration.

### Configuration Files

`get_config_dict()` returns the world as the builder's world table, so `build_world(world.get_config_dict())` rebuilds the same structure (see [Round Trip](../config/toml_schema.md#round-trip)). Beside the scalars (`schema_version`, `name`, `type` from `get_builder_world_type()`, `radius_m`, `mass_kg`, `albedo`, `emissivity`, `obliquity_rad`, `spin_frequency_rad_s`) it holds:

- `moment_of_inertia_factor` (the spin model's), when the source configuration gave it or it differs from the `[worlds]` default for the world type, which a rebuild takes again.
- A `tides` table when a tide model is attached (`global_tidal_model`, its per-degree parameters, and `get_tide_config()`). Without one, a rebuild takes the type's default (`rheology` for terrestrial or layered, `fixed_dt` for a gas giant, `fixed_q` for a star).
- `eos_solver` and `radial_solver` tables when a solver key is pinned, and a `prescribed_heating` table when a layer has one.
- A `layers` table of each layer's config dict ([Layer Serialization](../layers/layer.md#serialization)).

`save_to_toml` writes `get_save_config(destination_dir=None)`, validated against the schema when no build configuration is retained. `source_config`, `portable_config` (for a world built from a `data_file`, what `save_to_toml` writes while the world is unchanged), and `built_config` (`get_config_dict()` after the build) keep the configurations from the build. `family_world_type()` gives the builder's world type for the class and `get_schema_version_str()` its schema version. `save_config`, `save_binary`, and `load_binary` (a `str` or `os.PathLike` path) come from `TidalPyBaseClass`.

### Binary Serialization

A world's [binary file](../../Utilities/binary.md) holds every layer with its material and models, and every setting that affects calculations: the tide model and `get_tide_config()`, the spin model's `moment_of_inertia_factor`, the pinned solver settings, each layer's prescribed heating, and a `StarWorld`'s luminosity model. Solved state (profiles, zones, Love numbers, tide results) and the tidal heat source are not saved: rerun `solve_eos` and `calc_tides` after a load. Loading into an existing world replaces all of these (a record without a tide model leaves none), and layer views taken before it refer to replaced layers. A layer that belongs to a world cannot be loaded in place. A file of another class is refused with both classes named ("it is a GasGiantWorld file, not a TerrestrialWorld one").

`TidalPy.Structures.load_world(path, force=False)` reads a file into a new world of the saved class (as does `build_world(path)` for a binary file), raising `IOError` naming what a non-world file holds. `copy()`, `copy.copy`, `copy.deepcopy`, and `pickle` use the same record in memory, plus `source_config`, `portable_config`, and `built_config`; a copy is unsolved and belongs to no system. `TidalPy.Structures.worlds.world_from_bytes(world_class, record, configs=None, force=False)` rebuilds a world from a record.

```python
import pickle
from TidalPy.Structures import load_world

world.save_binary("world.tpyb")                       # A str or an os.PathLike path
loaded = load_world("world.tpyb")                     # A new world of the saved class
twin = world.copy()                                   # Same class, layers, models, and settings
twin.solve_eos()                                      # Copies are unsolved
restored = pickle.loads(pickle.dumps(world))          # A pickle round trip, as a process pool makes
```

### World Classes

`TerrestrialWorld` (rocky and icy planets and moons, the bundled ones included), `GasGiantWorld`, and `StarWorld` derive from `BaseWorld` (a `StructureBase` and `TidalPyBaseClass`), each with its own default `world_type` and binary class id. A gas giant's layers are usually of a liquid-only material (MatPack has `simple_gas`, `h2_he_molecular`, and `h_he_metallic`):

```python
from TidalPy.Structures.worlds import GasGiantWorld

jupiter = GasGiantWorld("Jupiter", 7.0e7, 1.898e27)
jupiter.add_layer(Layer(
    "envelope", 0, 0.0, 7.0e7,
    material="simple_gas"))                     # A liquid-only material: a liquid layer
print(jupiter.envelope.is_liquid)               # True
```

`StarWorld` adds an effective temperature and luminosity linked by the Stefan-Boltzmann law, $L = 4\pi R^{2}\sigma T^{4}$. It may hold layers and solve its EOS (the builder accepts `[layers.<name>]` tables for a star) but needs neither.

```python
from TidalPy.Structures.worlds import StarWorld

sun = StarWorld("Sun", 6.957e8, 1.989e30, effective_temperature=5772.0)
sun.luminosity                # ~3.83e26 W (derived from T if luminosity == 0)
sun.set_luminosity(3.828e26)  # recomputes effective_temperature
sun.effective_temperature
```

- Also `calc_luminosity_from_temperature(T)`, `calc_temperature_from_luminosity(L)`, `set_effective_temperature(T)`, and `calc_insolation_flux(distance, eccentricity=0.0)`: the orbit-averaged flux \[W m$^{-2}$\] $L / (4 \pi a^2 \sqrt{1 - e^2})$ at semi-major axis `distance` \[m\] (floats or arrays), which a `System` gives its worlds. `set_luminosity_model` attaches a luminosity model (fixed, mass-to-luminosity, power law; [Luminosity](../../Stellar/luminosity.md)).
- **Tides:** a star without layers uses `cpl`, `ctl`, or `ctl_q`; `calc_tides` raises if `rheology` is selected without layers and a solved EOS.
- **Spin:** a `System` evolves a star's spin. Without a solved EOS its moment of inertia is `moment_of_inertia_factor` $\times\,M R^{2}$, from `[worlds.star]`: 0.0754 (an $n = 3$ polytrope, a Sun-like star, as the `[tides.star]` Love numbers assume). A fully convective M dwarf is closer to 0.205 ($n = 1.5$); set `moment_of_inertia_factor` in the star's TOML, or call `set_spin_model`.

### C++ API

```cpp
tidalpy::c_WorldEOSSolveConfig cfg;          // starts from the [eos_solver] section of the configuration
cfg.surface_pressure = 1.0e5;                // override only what this solve changes
world.solve_eos(cfg);
double rho = world.get_density(5.0e6);

tidalpy::c_LoveSolveConfig love_cfg = world.make_love_solve_config();   // [radial_solver] defaults + the world's [tides]
love_cfg.frequency = 1.0e-5;                                            // [rad/s]
love_cfg.degree_l  = 2;
world.solve_love_numbers(love_cfg);
std::complex<double> k2 = world.get_love_number_k(0);
```

A default-constructed `c_WorldEOSSolveConfig` (`c_LoveSolveConfig`) reads the `[eos_solver]` (`[radial_solver]`) section of the shared runtime config, the same defaults Python uses. `solve_eos` throws `std::invalid_argument` on bad input (`ValueError` in Python). `c_BaseWorld` also has:

- Profile reads: `get_density`, `get_gravity`, `get_pressure(double) const`; `get_eos_fields(field_indices, num_fields, radii, num_radii, values_out) const` (fields of the `eos_layout_.hpp` layout, field-major, from one solve); `calc_complex_moduli(is_shear, radii, num_radii, frequency, moduli_out) const`.
- EOS results: `get_eos_success`, `get_eos_message`, `get_eos_iterations`, `get_eos_pressure_error`, `get_surface_gravity_eos`, `get_surface_pressure_eos`, `get_central_pressure`, `get_planet_mass_eos`, `get_planet_moi_eos`, `get_eos_solved()`, `get_all_materials_set()`, `get_eos_solution() -> const c_EOSSolution*`.
- Zones: `get_zones()` (as solved), `get_radial_zones()` (under the layers' current flags, as the Love solver uses them), `get_molten_regions()`, `get_is_liquid_at(radius)`.
- Heat sources: `set_prescribed_heating`, `get_tidal_heat_source`, `clear_tidal_heating`, `get_heating`, `calc_layer_temperature_rate`, `calc_layer_thermal_capacity`, `calc_layer_latent_capacity` (kinds `c_HeatSourceKind::Radiogenic`, `Tidal`, `Prescribed`).
- Love results: `get_love_number_k/h/l`, `get_love_surface_y`, and the status accessors. The Love solve's non-dimensional time scale is $1/\sqrt{\pi G \bar{\rho}}$ (bulk density $\bar{\rho}$), not $1/\omega$, so its setup is frequency-independent: a new frequency reuses it, and a new `solve_eos` or a change of degree, layer flags, or `solve_for` rebuilds it.

## References

- Schubert, G., Turcotte, D. L., and Olson, P. (2001). *Mantle Convection in the Earth and Planets*. Cambridge University Press. The upper-mantle temperature of parameterized convection.
- Solomatov, V. S. (2000). Fluid dynamics of a terrestrial magma ocean. In *Origin of the Earth and Moon*, 323-338. University of Arizona Press. The liquid convection scaling.
- Stevenson, D. J., Spohn, T., and Schubert, G. (1983). Magnetism and thermal evolution of the terrestrial planets. *Icarus*, 54(3), 466-489. Parameterized convection.
