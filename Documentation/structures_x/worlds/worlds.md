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
| `GasGiantWorld` | yes | yes (later phase) | Gas giant; a `LayeredWorld` with its own type/binary id. |
| `LayeredWorld` | yes | yes (later phase) | Terrestrial/layered body; owns layers, aggregates mass and heating. |

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
