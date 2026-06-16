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

`LayeredWorld.solve_love_numbers(...)` runs the shooting-method radial integration
to compute the viscoelastic-gravitational Love tidal or loading numbers k, h, l for a given tidal
forcing frequency. `solve_eos` must be called first.

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
world.solve_love_numbers(cfg);           // populates p_love_solution

std::complex<double> k2 = world.get_love_number_k(0);
```

`c_LayeredWorld::solve_love_numbers(const c_LoveSolveConfig&)` is fully
self-contained:

1. Validates `eos_solved` and `tidalpy_config_ptr`.
2. Copies the EOS radius/density/gravity/pressure/mass/moi arrays from the stored
   `c_EOSSolution` and calls `calc_complex_shear/bulk_modulus` at the requested
   frequency to build per-slice complex-moduli arrays.
3. Builds a `c_NonDimensionalScales` object and non-dimensionalizes all arrays
   (if `cfg.nondimensionalize`).
4. Creates a `c_RadialSolutionStorage` (owned by `p_love_solution`) and populates
   its `c_EOSSolution` via `inject_from_world_eos`, which sets the
   `p_use_array_interp` flag so the EOS solution interpolates from the copied
   arrays (no CyRK dense-output required in the radial-solve context).
5. Calls `c_shooting_solver` directly (no intermediate `c_radial_solver`).
6. Re-dimensionalizes and calls `c_RadialSolutionStorage::find_love()`.

The `c_RadialSolutionStorage` is owned by `p_love_solution
(unique_ptr<c_RadialSolutionStorage>)` on `c_LayeredWorld`; `get_love_number_k/h/l`,
`get_love_surface_y`, and the status accessors delegate to it.

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
