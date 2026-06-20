# Worlds (`structures_x.worlds`)

The world classes are the top-level structural objects in TidalPy.
A world owns its identity, orbital/thermal scalars, and bulk geometry; a *layered*
world additionally owns an ordered stack of [layers](../layers/) and is the object
on which the whole-planet equation-of-state and radial (Love number) solves run.

All values are **MKS** (radius [m], mass [kg], angles [rad], frequency [rad/s],
temperature [K], luminosity [W]).

---

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
| `StarWorld` | no | no | A star; effective temperature ↔ luminosity via Stefan-Boltzmann. |
| `GasGiantWorld` | yes | yes | Gas giant; a `LayeredWorld` with its own type/binary id. |
| `LayeredWorld` | yes | yes | Terrestrial/layered body; owns layers, aggregates mass and heating, runs the whole-planet EOS solve. |

---

## `BaseWorld`

```python
from TidalPy.structures_x.worlds import BaseWorld

w = BaseWorld(
    name="Earth", radius_m=6.371e6, mass_kg=5.972e24,
    world_type="terrestrial", albedo=0.3, emissivity=1.0,
    obliquity_rad=0.41, spin_frequency_rad_s=7.29e-5,
)
```

**Read-only properties:** `name`, `world_type`, `radius`, `mass`, `albedo`,
`emissivity`, `obliquity`, `spin_frequency`.

**Methods**

| Method | Returns | Description |
|--------|---------|-------------|
| `calc_surface_gravity()` | float [m/s²] | G·M/R². |
| `calc_escape_velocity()` | float [m/s] | √(2·G·M/R). |
| `calc_mean_density()` | float [kg/m³] | M / V_sphere(R). |
| `calc_equilibrium_temperature(F)` | float [K] | `[ (1−A)·F / (4·ε·σ) ]^(1/4)` (fast rotator, `F` = insolation flux [W/m²]). |
| `set_spin_frequency(ω)` | — | Set rotation rate [rad/s]. |
| `set_obliquity(θ)` | — | Set axial obliquity [rad]. |

`get_config_dict()` returns `name`, `world_type`, `radius_m`, `mass_kg`, `albedo`,
`emissivity`, `obliquity_rad`, `spin_frequency_rad_s`. `save_config` /
`save_binary` / `load_binary` are inherited from `TidalPyBaseClass`.

Binary class id: **200** (`BinaryClassID::BaseWorld`).

---

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
| `add_layer(layer)` | Add a layer inner-to-outer. **Ownership of the layer (and its attached physics models) transfers into the world**; the passed wrapper becomes an empty shell. Raises `ValueError` if the layer was already added or if its inner radius is not continuous with the current outermost radius (innermost must start at 0). A rejected layer is *not* consumed. |
| `num_layers` | Number of layers (property). |
| `calc_total_mass()` | Σ layer masses [kg]. |
| `calc_internal_heating(time_s)` | Σ radiogenic heating [W]; only `SolidLiquidLayer`s with an attached radiogenics model contribute. |
| `validate_layers()` | `True` if every boundary is continuous and the innermost starts at 0. |

`get_config_dict()` adds `num_layers` and a `layers` list (each a geometry-level
dict, in index order) to the `BaseWorld` keys.

Binary class id: **201** (`BinaryClassID::LayeredWorld`). See
[Binary serialization](#binary-serialization) below.

### Equation of state

Each layer carries a [material EOS model](../../material_x/material_eos.md) (its
density source), attached with `BaseLayer.set_eos(model)`. Once every layer has
one, `LayeredWorld.solve_eos(...)` integrates the planet's radial structure from
center to surface and populates every layer's density/gravity/pressure profile.

```python
from TidalPy.structures_x.worlds import LayeredWorld
from TidalPy.structures_x.layers import BaseLayer
from TidalPy.Material_x.eos import make_material_eos

world  = LayeredWorld("Earth", 6.371e6, 5.972e24)
core   = BaseLayer("core",   0, 0.0,     3.485e6, 0.0)
mantle = BaseLayer("mantle", 1, 3.485e6, 6.371e6, 0.0)
core.set_eos(make_material_eos("constant", {"reference_density_kg_m3": 11000.0}))
mantle.set_eos(make_material_eos("constant", {"reference_density_kg_m3": 4500.0}))
world.add_layer(core)
world.add_layer(mantle)

result = world.solve_eos(surface_pressure=0.0)   # dict of profile arrays + scalars
rho = world.get_density(5.0e6)                    # [kg/m³] at radius 5000 km
g   = world.get_gravity(world.radius)             # surface gravity [m/s²]
p0  = world.get_pressure(0.0)                      # central pressure [Pa]
```

The solver carries pressure as a radial state variable, so analytic
density-from-pressure models (Birch-Murnaghan, Vinet) are evaluated inline; the
constant and interpolated models ignore pressure. A convergence loop on the
surface pressure fixes the central pressure.

**`solve_eos(surface_pressure=0.0, slices_per_layer=100, G_to_use=-1.0,
integration_method='DOP853', rtol=1e-6, atol=1e-10, pressure_tol=1e-3,
max_iters=100, temperature=0.0, verbose=False) -> dict`**

Raises `ValueError` if the world has no layers, any layer lacks an EOS model,
`slices_per_layer < 2`, or the integration method is unknown. The returned dict
contains `success`, `message`, `iterations`, `pressure_error`, the profile arrays
(`radius`, `gravity`, `pressure`, `mass`, `moi`, `density`), and the scalar
results (`surface_gravity`, `surface_pressure`, `central_pressure`,
`planet_mass`, `planet_moi`).

**Profile queries (after a successful solve)**

| Member | Returns | Description |
|--------|---------|-------------|
| `get_density(r)` | float [kg/m³] | Density at radius `r` [m] (NaN if unsolved). |
| `get_gravity(r)` | float [m/s²] | Gravitational acceleration at `r`. |
| `get_pressure(r)` | float [Pa] | Pressure at `r`. |
| `eos_solved` | bool | `True` once profiles are populated. |
| `all_eos_set` | bool | `True` once every layer has an EOS model. |
| `surface_gravity_eos`, `central_pressure`, `planet_mass_eos`, `planet_moi_eos` | float | Scalar results of the last solve (NaN if unsolved). |

`get_density` / `get_gravity` / `get_pressure` delegate to the layer that contains
`r` (radii beyond the surface clamp to the outermost layer), so the individual
layers expose the same getters independently.

#### C++ API

The entire solve runs in C++ and can be driven directly from other C++ code:

```cpp
tidalpy::c_WorldEOSSolveConfig cfg;          // surface_pressure, slices_per_layer,
cfg.integration_method = ODEMethod::DOP853;  // G_to_use, rtol/atol, ...
world.solve_eos(cfg);                        // populates each layer's c_LayerEOSData
double rho = world.get_density(5.0e6);
```

`c_LayeredWorld::solve_eos(const c_WorldEOSSolveConfig&)` generates the per-layer
radius grids, estimates the bulk density, calls `Material_x/eos`'s `c_solve_eos`
(CyRK ODE integration with the surface-pressure loop), and slices the result into
each layer's `c_LayerEOSData`. It throws `std::invalid_argument` on bad input
(surfaced as `ValueError` in Cython via `except +`). The per-layer density source
is `tidalpy::c_MaterialEOSBase`, attached via
`c_BaseLayer::set_eos(std::unique_ptr<c_MaterialEOSBase>)` (non-owning observer
through `get_eos()`; `get_eos_set()`); during integration the
`c_preeval_material_eos` pre-eval (in `Material_x/eos/methods/material_.hpp`,
paired with the `c_MaterialEOSInput` struct holding the model pointer) dispatches
to `model->calc_density(pressure, temperature, radius)`.

`c_LayeredWorld` exposes `get_density/get_gravity/get_pressure(double) const`
(delegating to the containing layer via `find_layer_for_radius`), the result
accessors `get_eos_success/get_eos_message/get_eos_iterations/get_eos_pressure_error/`
`get_surface_gravity_eos/get_surface_pressure_eos/get_central_pressure/`
`get_planet_mass_eos/get_planet_moi_eos`, the retained full solution via
`get_eos_solution() -> const c_EOSSolution*`, and `get_eos_solved()` /
`get_all_eos_set()`. The Cython `LayeredWorld.solve_eos` wrapper only converts the
integration-method string to the CyRK enum, fills `c_WorldEOSSolveConfig`, calls
the C++ method under `nogil`, and builds the Python result dict from the retained
solution.

### Viscoelastic Properties (after EOS solve)

Once `solve_eos` has succeeded, the world and each layer expose the full
complex-moduli / viscoelastic getter surface, which is needed by the Love number / radial solver.

| Member | Returns | Description |
|--------|---------|-------------|
| `get_shear_modulus(r)` | float [Pa] | Post-melt static shear modulus at `r`. |
| `get_bulk_modulus(r)` | float [Pa] | Post-melt static bulk modulus at `r`. |
| `get_shear_viscosity(r)` | float [Pa·s] | Post-melt shear viscosity at `r`. |
| `get_bulk_viscosity(r)` | float [Pa·s] | Post-melt bulk viscosity at `r`. |
| `calc_complex_shear_modulus(r, ω)` | complex [Pa] | Rheology-derived complex shear modulus at radius `r` [m], frequency `ω` [rad/s]. |
| `calc_complex_bulk_modulus(r, ω)` | complex [Pa] | Rheology-derived complex bulk modulus. |

All of the above accept a float or `np.ndarray` for `r` (and `ω`);
array inputs return an `np.ndarray` of the same shape.

---

### RadialSolver (Calculating Love Numbers)

`LayeredWorld.solve_love_numbers(...)` computes the viscoelastic-gravitational Love
tidal or loading numbers k, h, l for a given tidal forcing frequency. `solve_eos`
must be called first.

Two radial-solve methods are available behind the same call, selected by
`use_prop_matrix`:

* **Shooting method** (`use_prop_matrix=False`, default) — numerically integrates
  the radial ODEs from the center to the surface. Works for arbitrary multi-layer,
  solid/liquid, static/dynamic, compressible/incompressible worlds.
* **Propagation-matrix method** (`use_prop_matrix=True`) — a quasi-analytic
  matrix-propagation technique restricted to a **single solid, static,
  incompressible layer** (the classic homogeneous-sphere case). An incompatible
  world fails the solve gracefully (`love_success` is `False`, `love_error_code`
  non-zero) rather than raising. `core_model` selects the core starting condition.

```python
from TidalPy.structures_x.worlds import LayeredWorld
from TidalPy.structures_x.layers import SolidLiquidLayer
from TidalPy.rheology_x import make_rheology
from TidalPy.Material_x.eos import make_material_eos

world = LayeredWorld("planet", 6e6, 4.2e24)
layer = SolidLiquidLayer("mantle", 0, 0.0, 6e6, 4.2e24)
layer.set_eos(make_material_eos("constant", {"reference_density_kg_m3": 4000.0}))
layer.set_shear_rheology(make_rheology("maxwell", {
    "shear_modulus_pa": 6e10, "viscosity_pa_s": 1e21
}))
world.add_layer(layer)

world.solve_eos()
world.solve_love_numbers(frequency_rad_s=1e-5)

print(world.love_k2)   # complex k₂
print(world.love_h2)   # complex h₂
print(world.love_l2)   # complex l₂
```

**`solve_love_numbers(
   frequency_rad_s=1e-5,
   degree_l=2,
   solve_tidal=True,
   use_kamata=True,
   nondimensionalize=True,
   starting_radius=0.0,
   start_radius_tol=1e-4,
   integration_method='DOP853',
   rtol=1e-6,
   atol=1e-10,
   scale_rtols=True,
   max_num_steps=500000,
   expected_size=500,
   max_ram_MB=500,
   max_step=0.0,
   verbose=False,
   warnings=True,
   eos_rtol=1e-6,
   eos_atol=1e-10,
   eos_pressure_tol=1e-3,
   eos_max_iters=100) -> None`**

Raises `RuntimeError` if the EOS has not yet been solved. Results are stored
internally and accessed through the properties below.

**Love-number properties (after a successful solve)**

| Property | Type | Description |
|----------|------|-------------|
| `love_solved` | bool | `True` if a solve has been attempted. |
| `love_success` | bool | `True` if the last solve converged. |
| `love_error_code` | int | Solver error code (0 = success; < 0 = failure). |
| `love_message` | str | Human-readable solver message. |
| `love_num_ytypes` | int | Number of independent solution types (boundary-condition models requested). |
| `love_k2`, `love_h2`, `love_l2` | complex | Degree-2 Love numbers for the default (tidal) boundary condition. |

**`get_love_number_k(ytype_idx=0) -> complex`**,
**`get_love_number_h(ytype_idx=0) -> complex`**,
**`get_love_number_l(ytype_idx=0) -> complex`**

Return the Love numbers for boundary-condition model index `ytype_idx` (0 = first
requested, usually tidal).

**`get_love_surface_y(ytype_idx, y_idx) -> complex`**

Raw radial function y₁…y₆ at the surface for solution type `ytype_idx`, function
index `y_idx` (0–5).

#### C++ API

```cpp
tidalpy::c_LoveSolveConfig cfg;
cfg.frequency_rad_s   = 1.0e-5;          // [rad/s]
cfg.degree_l          = 2;
cfg.nondimensionalize = true;
cfg.integration_method = ODEMethod::DOP853;
world.solve_love_numbers(cfg);           // delegates to the cached radial solver

std::complex<double> k2 = world.get_love_number_k(0);
```

`c_LayeredWorld::solve_love_numbers(const c_LoveSolveConfig&)` delegates to a
cached helper, `c_WorldRadialSolver` (held by `p_radial_solver`), that separates the
**frequency-independent** setup (built once and reused) from the
**frequency-dependent** work (recomputed cheaply on every call). This matters
because the Love-number solve is the hot loop for frequency sweeps and orbital
evolution.

The non-dimensionalization is itself frequency-independent (the
`c_NonDimensionalScales` time scale is `1/(π·G·ρ_bulk)`, not `1/ω`), so the only
quantities that change between calls at different frequencies are the complex
moduli and the shooting integration.

1. Validates `eos_solved` and `tidalpy_config_ptr`.
2. If the cache does not match the current EOS grid/assumptions, `build_cache`
   captures (once): the non-dim radius/density/gravity/pressure/mass/moi arrays,
   per-layer metadata (solid/liquid, static, incompressible) and slice
   partitioning, the non-dim scalars (`G`, bulk density, surface pressure), and a
   **reused** `c_RadialSolutionStorage` whose internal `c_EOSSolution` arrays serve
   as the scratch buffers. The cache is invalidated automatically whenever
   `solve_eos` re-runs.
3. Per call: `calc_complex_shear/bulk_modulus` fills the dimensional moduli scratch
   at the requested frequency; the helper non-dimensionalizes them in place,
   re-applies the cached non-dim structure arrays via `inject_from_world_eos`
   (which sets `p_use_array_interp` so the EOS interpolates from the arrays — no
   CyRK dense output needed), and runs the selected solver.
4. Re-dimensionalizes the y-solution, restores the SI surface gravity, and calls
   `c_RadialSolutionStorage::find_love()`.

`get_love_number_k/h/l`, `get_love_surface_y`, and the status accessors read
through `p_radial_solver->get_storage()`.

### Global (1D) Tidal Dissipation

`LayeredWorld.calc_tides(...)` computes the body's total tidal heating and three
orbital potential partial derivatives by summing over the active tidal modes
(the global / "1D potential" approach). A tide model (see
[Global Tidal Dissipation](../../Tides_x/global_tides.md)) supplies the per-mode
dissipation multiplier `−Im[k_l]`; the world runs the global-potential engine for its
stored `[tides]` config + the supplied orbital/spin state, collapses, and distributes the
heat to layers by each layer's `tidal_scale`.

The tide model + config are normally wired by the world builder from the `[tides]` TOML
table (with per-family defaults: star = `fixed_q`, gasgiant = `fixed_dt`, terrestrial = `rheology`),
but can also be set directly:

```python
from TidalPy.Tides_x.classes import make_tide

world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]}))
world.set_tide_config(max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)

world.calc_tides(orbital_frequency=2.05e-5, spin_frequency=2.05e-5, eccentricity=0.0041,
                 obliquity=0.0, semi_major_axis=4.2e8, host_mass=1.898e27)

world.get_tidal_heating()                # total heating [W]
world.get_tidal_potential_derivatives()  # (dUdM, dUdw, dUdO) [J kg-1 rad-1]
world.get_layer_tidal_heating(0)         # = world heating × layer 0's tidal_scale
```

For a synchronous, low-eccentricity body the `cpl` result reproduces the classic CPL rate
`(21/2)(k₂/Q)·G·M_host²·R⁵·n·e²/a⁶`.

### Rheology model

The analytic models (`cpl`/`ctl`/`ctl_q`) take `−Im[k_l]` from their fixed per-degree
parameters and need no interior solution. The `rheology` model instead derives `−Im[k_l(ω)]`
from the world radial solver: `calc_tides` runs the global-potential engine, then for each
unique tidal frequency it solves the world's complex Love numbers (reusing the frequency-
independent radial-solver cache), feeds the per-mode `k_l` into the collapse, and retains the
full `k`/`h`/`l` suite per mode for inspection. Because it runs the radial solver, the EOS must
be solved first:

```python
world.solve_eos(G_to_use=G, temperature=1500.0)   # required for the rheology model
world.set_tide_model(make_tide("rheology"))
world.set_tide_config(max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
world.calc_tides(orbital_frequency=2.0e-5, spin_frequency=1.0e-5, eccentricity=0.01,
                 obliquity=0.0, semi_major_axis=4.0e8, host_mass=1.9e27)

world.get_tidal_heating()             # total heating [W]
world.get_tidal_love_k(2, 2, 0, 0)    # complex k₂ for the (l,m,p,q) = (2,2,0,0) mode
```

`calc_tides` raises `RuntimeError` if the `rheology` model is selected but the EOS has not been
solved, or if a per-frequency radial solve fails.

| Member | Returns | Description |
|--------|---------|-------------|
| `set_tide_model(tide)` | — | Attach a tide model (transfers ownership). |
| `tide_model_set` | bool | Whether a model is attached. |
| `set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=10)` | — | Set the stored `[tides]` truncation/degree. |
| `calc_tides(orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass)` | — | Run the global tidal solve. |
| `tides_solved` | bool | Whether a solve has succeeded. |
| `get_tidal_heating()` | float [W] | Total global tidal heating (NaN if unsolved). |
| `get_tidal_potential_derivatives()` | tuple | `(dUdM, dUdw, dUdO)` [J kg⁻¹ rad⁻¹]. |
| `get_num_tidal_modes()` | int | Active modes summed. |
| `get_layer_tidal_heating(index)` | float [W] | Per-layer heating = world heating × `tidal_scale`. |
| `get_tidal_love_k(l, m, p, q)` | complex | Per-mode radial-solver `k_l` (rheology only; NaN for analytic models). |

Each layer also stores its own tidal heating: the C++ `c_BaseLayer::get_tidal_heating()` returns
the world heating scaled by that layer's contribution.

---

## `GasGiantWorld`

A `LayeredWorld` whose `world_type` defaults to `"gasgiant"` and which uses a
dedicated binary class id **202** (`BinaryClassID::GasGiantWorld`). Same API as
`LayeredWorld`; typically populated with `GasLayer`s.

```python
from TidalPy.structures_x.worlds import GasGiantWorld
from TidalPy.structures_x.layers import GasLayer

jupiter = GasGiantWorld("Jupiter", 7.0e7, 1.898e27)
jupiter.add_layer(GasLayer("envelope", 0, 0.0, 7.0e7, 1.898e27))
```

---

## `StarWorld`

A star: no layers, no EOS. Effective temperature and luminosity are kept
consistent through the Stefan-Boltzmann law `L = 4·π·R²·σ·T⁴`.

```python
from TidalPy.structures_x.worlds import StarWorld

sun = StarWorld("Sun", 6.957e8, 1.989e30, effective_temperature_k=5772.0)
sun.luminosity            # ~3.83e26 W (derived from T if luminosity_w == 0)
sun.set_luminosity(3.828e26)  # recomputes effective_temperature
```

**Properties:** `effective_temperature` [K], `luminosity` [W].
**Methods:** `calc_luminosity_from_temperature(T)`,
`calc_temperature_from_luminosity(L)`, `set_effective_temperature(T)`,
`set_luminosity(L)`. A richer luminosity-model hierarchy (mass→luminosity) is a
later phase.

Binary class id: **203** (`BinaryClassID::StarWorld`).

---

## Binary serialization

A `LayeredWorld` (and `GasGiantWorld`) serializes its `BaseWorld` fields and a
layer count, then **each layer's own complete binary record in index order**.
Because each layer recursively serializes its attached rheology / cooling /
radiogenics models (see [Binary serialization](../../utilities_x/binary_x.md)), a
single `save_binary` / `load_binary` round-trips the entire world graph — no
Python reconstruction step is needed. On load, each layer is rebuilt as the
correct concrete subclass via the layer binary-dispatch factory
(`c_layer_from_binary`).

```python
world.save_binary("earth.tpyb")
reloaded = LayeredWorld("placeholder", 1.0, 1.0)
reloaded.load_binary("earth.tpyb")
assert reloaded.num_layers == world.num_layers
assert reloaded.calc_internal_heating(0.0) == world.calc_internal_heating(0.0)
```

EOS profile data (`c_LayerEOSData`) is never serialized; it is repopulated by
running the whole-planet EOS solve after loading.
