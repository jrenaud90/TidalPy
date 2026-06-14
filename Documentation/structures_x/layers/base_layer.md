# BaseLayer

`TidalPy.structures_x.layers.BaseLayer`

## Overview

`BaseLayer` is the geometry-only base for all TidalPy layer types. It stores the
inner and outer radii, total mass, and an optional material identifier for one
spherically symmetric shell inside a planetary body.

All spatial data is stored and returned in **MKS units** (meters, kilograms,
seconds). Derived geometry (thickness, volume, surface areas) is computed at
construction and accessible via read-only properties.

A **material EOS model** (the layer's density source) is attached with `set_eos`.
An **EOS profile** (density, gravity, and pressure as a function of radius) is then
populated by the world-level EOS solve
([`LayeredWorld.solve_eos`](../worlds/worlds.md#equation-of-state)), or directly via
`update_eos_data`. Until populated, all EOS getters return `NaN`.

---

## Inheritance

```
TidalPyBaseClass
  └── StructureBase
        └── BaseLayer
```

`BaseLayer` inherits binary serialization (`save_binary`, `load_binary`) and TOML
config saving (`save_config`, `get_config_dict`) from `TidalPyBaseClass` via
`StructureBase`.

---

## Constructor

```python
BaseLayer(
    name:           str,
    layer_index:    int,
    radius_inner_m: float,
    radius_outer_m: float,
    mass_kg:        float,
    material_name:  str   = "",
    is_tidal:       bool  = True,
    tidal_scale:    float = 1.0,
)
```

### Parameters

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| `name` | `str` | — | Human-readable identifier (e.g. `"mantle"`). |
| `layer_index` | `int` | — | Zero-based index; innermost layer = 0. |
| `radius_inner_m` | `float` | m | Inner boundary radius. |
| `radius_outer_m` | `float` | m | Outer boundary radius. |
| `mass_kg` | `float` | kg | Total layer mass. |
| `material_name` | `str` | — | Material identifier (e.g. `"perovskite"`). Optional. |
| `is_tidal` | `bool` | — | Whether this layer dissipates tidal energy. Default `True`. |
| `tidal_scale` | `float` | — | Dimensionless scale on tidal heating. Default `1.0`. |

---

## Properties

### Geometry (read-only, MKS)

| Property | Units | Description |
|----------|-------|-------------|
| `name` | — | Layer name. |
| `layer_index` | — | Zero-based index. |
| `radius` / `radius_outer` | m | Outer boundary radius. |
| `radius_inner` | m | Inner boundary radius. |
| `thickness` | m | `radius_outer - radius_inner`. |
| `mass` | kg | Total layer mass. |
| `volume` | m³ | Spherical shell volume. |
| `surface_area_outer` | m² | Outer surface area. |
| `surface_area_inner` | m² | Inner surface area. |
| `material_name` | — | Material identifier. |
| `is_tidal` | — | Tidal dissipation flag. |
| `tidal_scale` | — | Tidal heating scale factor. |

### EOS profile

| Property | Description |
|----------|-------------|
| `eos_data_populated` | `True` after the EOS profile has been populated (by the world EOS solve or `update_eos_data`). |
| `eos_set` | `True` after a material EOS model has been attached via `set_eos`. |

---

## Methods

### `set_eos(model)`

Attach a [material EOS model](../../material_x/material_eos.md) (the per-layer
density source). Ownership of the C++ model transfers into the layer; the passed
wrapper becomes an empty shell. The model is consumed by the world-level
[`solve_eos`](../worlds/worlds.md#equation-of-state), which integrates the planet
structure and populates this layer's EOS profile.

```python
from TidalPy.Material_x.eos import make_material_eos

layer.set_eos(make_material_eos("birch_murnaghan", {
    "reference_density_kg_m3": 4500.0,
    "reference_bulk_modulus_pa": 2.5e11,
    "bulk_modulus_derivative": 4.0,
}))
```

Raises `ValueError` if the model has already been attached or moved.

### `update_eos_data(radius_m, density_kgm3, gravity_ms2, pressure_pa)`

Populate the EOS profile directly from sorted radius arrays (normally done for you
by the world EOS solve; useful for tests or manual construction).

```python
import numpy as np

r   = np.linspace(r_inner, r_outer, 200)
rho = ...  # density profile [kg/m^3]
g   = ...  # gravity profile [m/s^2]
p   = ...  # pressure profile [Pa]

layer.update_eos_data(r, rho, g, p)
```

**Notes:**
- In normal workflow this is called automatically by the world EOS solve
  ([`LayeredWorld.solve_eos`](../worlds/worlds.md#equation-of-state)).
- All sequences must be the same length and `radius_m` must be sorted ascending.
- Linear interpolation is used; values are clamped at the layer boundaries.

### `get_density(radius_m)` → float

Density at `radius_m` [kg/m³]. Returns `NaN` if EOS data not populated.

### `get_gravity(radius_m)` → float

Gravitational acceleration at `radius_m` [m/s²]. Returns `NaN` if not populated.

### `get_pressure(radius_m)` → float

Pressure at `radius_m` [Pa]. Returns `NaN` if not populated.

### Inherited geometry calculations (from `StructureBase`)

Pure-function methods that do not depend on stored state:

```python
layer.calc_surface_area(r)          # 4πr² [m²]
layer.calc_volume_sphere(r)         # (4/3)πr³ [m³]
layer.calc_volume_shell(r_out, r_in) # shell volume [m³]
layer.calc_surface_gravity(m, r)    # G·m/r² [m/s²]
layer.calc_mean_density(m, v)       # m/v [kg/m³]
layer.calc_escape_velocity(m, r)    # √(2Gm/r) [m/s]
```

### Binary I/O (inherited)

```python
layer.save_binary("layer.tpyb")

restored = BaseLayer("placeholder", 0, 0.0, 1.0, 1.0)
restored.load_binary("layer.tpyb")
```

**Note:** EOS profile data is NOT included in the binary file. It must be
repopulated after loading.

### TOML config (inherited)

```python
layer.save_config("layer.toml")
cfg = layer.get_config_dict()  # → dict with all construction parameters
```

---

## Example

```python
from TidalPy.structures_x.layers import BaseLayer

mantle = BaseLayer(
    name           = "mantle",
    layer_index    = 1,
    radius_inner_m = 3.485e6,   # CMB radius [m]
    radius_outer_m = 6.371e6,   # Earth surface [m]
    mass_kg        = 4.043e24,  # mantle mass [kg]
    material_name  = "perovskite",
)

print(f"Thickness:     {mantle.thickness / 1e3:.0f} km")
print(f"Volume:        {mantle.volume:.3e} m³")
print(f"EOS populated: {mantle.eos_data_populated}")
# EOS populated: False

# After EOSHandler runs (or for testing):
mantle.update_eos_data(
    radius_m     = [3.485e6, 6.371e6],
    density_kgm3 = [5560.0,  3300.0],
    gravity_ms2  = [10.68,    9.81],
    pressure_pa  = [1.36e11,  0.0],
)
print(f"EOS populated: {mantle.eos_data_populated}")
# EOS populated: True
print(f"Density at CMB: {mantle.get_density(3.485e6):.0f} kg/m³")
```
