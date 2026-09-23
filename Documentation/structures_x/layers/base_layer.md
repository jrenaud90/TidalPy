# BaseLayer

_Updated: 2026-09-23_

`TidalPy.structures_x.layers.BaseLayer` is the geometry-only base for all TidalPy layer types. It stores the inner and outer radii \[m\], total mass \[kg\], and an optional material identifier for one spherically symmetric shell inside a planetary body. Derived geometry (thickness, volume, surface areas) is computed at construction and accessible through read-only properties.

A material EOS model (the layer's density source) is attached with `set_eos`. An EOS profile (density, gravity, and pressure as a function of radius) is then populated by the world-level EOS solve ([`LayeredWorld.solve_eos`](../worlds/worlds.md#equation-of-state)), or directly with `update_eos_data`. Until populated, all EOS getters return `NaN`.

## Inheritance

```
TidalPyBaseClass
  └── StructureBase
        └── BaseLayer
```

`BaseLayer` inherits binary serialization (`save_binary`, `load_binary`) and TOML config saving (`save_config`, `get_config_dict`) from `TidalPyBaseClass` via `StructureBase`.

## Constructor

```python
BaseLayer(
    name:               str,
    layer_index:        int,
    radius_inner:       float,
    radius_outer:       float,
    mass:               float,
    material_name:      str   = "",
    is_tidal:           bool  = True,
    tidal_scale:        float = None,
)
```

### Parameters

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| `name` | `str` | - | Human-readable identifier (e.g. `"mantle"`). |
| `layer_index` | `int` | - | Zero-based index; innermost layer = 0. |
| `radius_inner` | `float` | m | Inner boundary radius. |
| `radius_outer` | `float` | m | Outer boundary radius. |
| `mass` | `float` | kg | Total layer mass. Overwritten by each successful world EOS solve. |
| `material_name` | `str` | - | Material identifier (e.g. `"perovskite"`). Optional. |
| `is_tidal` | `bool` | - | Whether this layer dissipates tidal energy. Default `True`. |
| `tidal_scale` | `float` | - | The layer's share of the planet in the quasi-homogeneous Love methods (`homogeneous`, `cpl`, `ctl`) and of an analytic tide model's heating; `None` (default) takes its volume fraction. See [Worlds](../worlds/worlds.md). |

## Properties

### Geometry

Read-only properties.

| Property | Units | Description |
|----------|-------|-------------|
| `name` | - | Layer name. |
| `layer_index` | - | Zero-based index. |
| `radius` / `radius_outer` | m | Outer boundary radius. |
| `radius_inner` | m | Inner boundary radius. |
| `thickness` | m | `radius_outer - radius_inner`. |
| `mass` | kg | Total layer mass. Set at construction, then overwritten by each successful world EOS solve with the mass the solved profile places between the layer's radii. |
| `volume` | m³ | Spherical shell volume. |
| `density_bulk` | kg/m³ | Bulk density, `mass / volume`; follows the EOS-set mass. |
| `surface_area_outer` | m² | Outer surface area. |
| `surface_area_inner` | m² | Inner surface area. |
| `material_name` | - | Material identifier. |
| `is_tidal` | - | Tidal dissipation flag. |
| `tidal_scale` | - | The configured tidal scale, or `None` when the layer takes its volume fraction. Settable. |

### EOS Profile

| Property | Description |
|----------|-------------|
| `eos_data_populated` | `True` after the EOS profile has been populated (by the world EOS solve or `update_eos_data`). |
| `eos_set` | `True` after a material EOS model has been attached via `set_eos`. |

## Methods

### `set_eos(model)`

Attach a [material EOS model](../../material_x/material_eos.md) (the per-layer density source). Ownership of the C++ model transfers into the layer; the passed wrapper becomes an empty shell. The model is consumed by the world-level [`solve_eos`](../worlds/worlds.md#equation-of-state), which integrates the planet structure and populates this layer's EOS profile.

```python
from TidalPy.Material_x.eos import make_material_eos

layer.set_eos(make_material_eos("birch_murnaghan", {
    "reference_density_kg_m3": 4500.0,
    "reference_bulk_modulus_pa": 2.5e11,
    "bulk_modulus_derivative": 4.0,
}))
```

Raises `ValueError` if the model has already been attached or moved.

### `update_eos_data(...)`

Populate the EOS profile directly from sorted radius arrays (normally done for you by the world EOS solve; useful for tests or manual construction).

```python
import numpy as np

r   = np.linspace(r_inner, r_outer, 200)
rho = ...  # density profile [kg/m^3]
g   = ...  # gravity profile [m/s^2]
p   = ...  # pressure profile [Pa]

layer.update_eos_data(r, rho, g, p)
```

In the normal workflow the world EOS solve ([`LayeredWorld.solve_eos`](../worlds/worlds.md#equation-of-state)) calls this. All sequences must be the same length and `radius` must be sorted ascending. Linear interpolation is used; values are clamped at the layer boundaries.

### `get_density(radius)` -> float or ndarray

Density at `radius` [kg/m³]. Returns `NaN` if EOS data not populated.

### `get_gravity(radius)` -> float or ndarray

Gravitational acceleration at `radius` [m/s²]. Returns `NaN` if not populated.

### `get_pressure(radius)` -> float or ndarray

Pressure at `radius` [Pa]. Returns `NaN` if not populated.

### Viscoelastic profile getters -> float or ndarray

After the world EOS solve populates the layer, the radius-resolved viscoelastic state is readable through the same getter names the world exposes: `get_shear_modulus`, `get_bulk_modulus`, `get_shear_viscosity`, `get_bulk_viscosity` (all after the partial-melt step), `get_melt_fraction`, and the shorthand bundles `get_static_viscoelastics(radius)` (the post-melt 4-tuple) and `get_state(radius)` (all profiles as a dict). All return `NaN` before the profile is populated.

Nothing here is stored on a grid, and nothing is calculated by the layer. The layer's material (its EOS model) evaluated these properties while the structure was integrated, and each getter reads them back from that solved EOS at the radius asked for, so a getter and the solve always agree. Every layer class has them, because the material belongs to the EOS model rather than to the layer class. A profile supplied by hand through `update_eos_data` carries density, gravity, and pressure alone, so these report `NaN` for it.

| Getter | Returns |
|---|---|
| `get_shear_modulus`, `get_bulk_modulus` | Static moduli [Pa] after melt weakening. |
| `get_shear_viscosity`, `get_bulk_viscosity` | Viscosities [Pa s] after melt weakening. |
| `get_melt_fraction` | Melt fraction from the attached partial-melt model; `0.0` without one. |
| `get_static_viscoelastics(radius)` | The post-melt four-tuple in one call. |
| `get_state(radius)` | Every profile at that radius as a dict. |

`viscoelastic_populated` says whether these are meaningful yet: it is `False` until the world's EOS solve gives the layer a structure profile to read, and every getter returns NaN before then.

Every profile getter on this page accepts a float or an `np.ndarray` of radii and returns a matching scalar or same-shape array, evaluated in a C loop:

```python
import numpy as np
radii = np.linspace(3.5e6, 6.3e6, 100)
rho = mantle.get_density(radii)        # ndarray, shape (100,)
mu, eta_mu, kk, eta_k = mantle.get_static_viscoelastics(radii)
```

### Inherited Geometry Calculations

Pure-function methods that do not depend on stored state (Inherited from `StructureBase`):

```python
layer.calc_surface_area(r)           # 4πr² [m²]
layer.calc_volume_sphere(r)          # (4/3)πr³ [m³]
layer.calc_volume_shell(r_out, r_in) # shell volume [m³]
layer.calc_surface_gravity(m, r)     # G·m/r² [m/s²]
layer.calc_mean_density(m, v)        # m/v [kg/m³]
layer.calc_escape_velocity(m, r)     # √(2Gm/r) [m/s]
```

### Binary I/O

```python
layer.save_binary("layer.tpyb")

restored = BaseLayer("placeholder", 0, 0.0, 1.0, 1.0)
restored.load_binary("layer.tpyb")
```

> [!NOTE]
> An attached material EOS model is saved and restored with the layer, but the EOS profile data it produces is not; re-run the world's `solve_eos` after loading.

### TOML Config

```python
layer.save_config("layer.toml")
cfg = layer.get_config_dict()  # -> dict with all construction parameters
```

The dict follows the world builder's layer schema: `class` names the layer class (`base`, `physics`, `solidliquid`, or `gas`), the scalar keys are the constructor parameters, and each attached physics model is a sub-table keyed by `model` (`material` here, the layer's material EOS model; subclasses add their own). `name` and `radius_inner` belong to a standalone layer only; a world drops them when it nests the layer under its name (`LAYER_STANDALONE_CONFIG_KEYS`).

```python
layer.set_eos(ConstantDensityEOS(reference_density=4400.0))
layer.get_config_dict()["material"]  # {'model': 'constant', 'reference_density_kg_m3': 4400.0}
```

### Tidal Bookkeeping

| Member | Description |
|---|---|
| `get_tidal_heating()` | Tidal heating deposited in this layer [W], set by the world's tidal solve. |
| `tidal_scale` | The layer's share of the planet in the quasi-homogeneous Love methods and of an analytic tide model's heating (its volume fraction when unset). Settable; see [Worlds](../worlds/worlds.md). |
| `is_tidal` | Whether the layer takes any tidal heating at all. A non-tidal layer always gets zero. |
| `get_schema_version_str()` | The schema version this class reads and writes, for configuration and binary compatibility. |

## Example

```python
from TidalPy.structures_x.layers import BaseLayer

mantle = BaseLayer(
    name        = "mantle",
    layer_index = 1,
    radius_inner = 3.485e6,   # CMB radius [m]
    radius_outer = 6.371e6,   # Earth surface [m]
    mass        = 4.043e24,  # mantle mass [kg]
    material_name  = "perovskite",
)

print(f"Thickness:     {mantle.thickness / 1e3:.0f} km")
print(f"Volume:        {mantle.volume:.3e} m³")
print(f"EOS populated: {mantle.eos_data_populated}")
# EOS populated: False

# After the world's EOS solve runs (or for testing):
mantle.update_eos_data(
    radius       = [3.485e6, 6.371e6],
    density_kgm3 = [5560.0,  3300.0],
    gravity_ms2  = [10.68,    9.81],
    pressure     = [1.36e11,  0.0],
)
print(f"EOS populated: {mantle.eos_data_populated}")
# EOS populated: True
print(f"Density at CMB: {mantle.get_density(3.485e6):.0f} kg/m³")
```
