# BaseLayer

_Updated: 2026-09-30_

`TidalPy.Structures.layers.BaseLayer` holds one spherically symmetric shell inside a planetary body: its inner and outer radii \[m\], total mass \[kg\], and an optional material identifier, the radial-solver assumptions, the layer temperature, the material, and the shear and bulk rheology. Derived geometry (thickness, volume, surface areas) is computed at construction and read through properties.

The material is defined via the layer's EOS model (see [Material EOS Models](../../Material/material_eos.md)), attached with `set_eos`. It holds the density law, the static moduli, the shear law, the viscosities, and the viscosity and partial-melt models, and the world-level EOS solve ([`BaseWorld.solve_eos`](../worlds/worlds.md#equation-of-state)) evaluates it as it integrates. The solve populates the layer's EOS profile (density, gravity, pressure, and the viscoelastic state as a function of radius).

When a rheology model (a `RheologyBase` subclass) is attached with `set_shear_rheology` or `set_bulk_rheology`, `calc_complex_shear_modulus` and `calc_complex_bulk_modulus` apply it to the static modulus and viscosity the solved EOS reports. Until then they return the static modulus as a purely real complex number, which is perfectly elastic behavior.

## Inheritance

```
TidalPyBaseClass
  └── StructureBase
        └── BaseLayer
              ├── SolidLiquidLayer
              └── GasLayer
```

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
    is_volume_fixed:    bool  = True,
    tidal_scale:        float = None,
    is_solid:           bool  = True,
    is_static:          bool  = True,
    is_incompressible:  bool  = False,
    temperature:        float = 0.0,
    use_thermal_eos:    bool  = False,
    use_heating:        bool  = False,
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
| `is_volume_fixed` | `bool` | - | `False` lets the layer grow or shrink to hold its mass during an EOS solve. Default `True`. |
| `tidal_scale` | `float` | - | The layer's share of the planet in the quasi-homogeneous Love methods (`homogeneous`, `cpl`, `ctl`) and of an analytic tide model's heating; `None` (default) takes its volume fraction. See [Worlds](../worlds/worlds.md). |
| `is_solid`, `is_static`, `is_incompressible` | `bool` | - | Radial-solver assumptions; see [Layer Assumptions](#layer-assumptions). Defaults `True`, `True`, `False`. |
| `temperature` | `float` | K | Layer temperature at which the material's viscosity and melt models are evaluated. Default `0.0`, the cold rigid limit of the viscosity laws. |
| `use_thermal_eos` | `bool` | - | Let the material's density law see the temperature, so the density and bulk modulus depend on it. Default `False`. |
| `use_heating` | `bool` | - | Let the world's heat sources act inside the layer during a thermal EOS solve. Default `False`. |

The static moduli and viscosities and the shear law are arguments of the EOS model, not of the layer:

```python
from TidalPy.Material.eos import ConstantDensityEOS
from TidalPy.Structures.layers import BaseLayer

mantle = BaseLayer("mantle", 1, 3.485e6, 6.371e6, 4.043e24)
mantle.set_eos(ConstantDensityEOS(
    reference_density=4500.0,
    shear_modulus_static=1.67e11,
    bulk_modulus_static=3.57e11))
```

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
| `is_volume_fixed` | - | `False` if the layer holds its mass rather than its volume. Writable. |
| `tidal_scale` | - | The configured tidal scale, or `None` when the layer takes its volume fraction. Writable. |

### Material and Models

Read-only properties.

| Property | Units | Description |
|----------|-------|-------------|
| `shear_modulus_static` | Pa | The material's unrelaxed shear modulus, read from the attached EOS model. NaN when none is attached. |
| `bulk_modulus_static` | Pa | The material's unrelaxed bulk modulus, as above. |
| `shear_viscosity_static` | Pa·s | The material's static shear viscosity, as above. |
| `bulk_viscosity_static` | Pa·s | The material's static bulk viscosity, as above. |
| `eos_set` | - | `True` after a material EOS model is attached. |
| `shear_rheology_set`, `bulk_rheology_set` | - | `True` after the corresponding rheology model is attached. |
| `shear_viscosity_set`, `bulk_viscosity_set` | - | `True` when the material holds the corresponding viscosity model. |
| `partial_melt_set` | - | `True` when the material holds a partial-melt model. |
| `eos_data_populated` | - | `True` after the EOS profile is populated (by the world EOS solve or `update_eos_data`). |
| `viscoelastic_populated` | - | `True` after the world EOS solve gives the layer its viscoelastic profile. |

### Layer State

| Property | Units | Description |
|----------|-------|-------------|
| `temperature` | K | Layer temperature. |
| `use_thermal_eos` | - | `True` if the material's density law receives the temperature. |
| `use_heating` | - | `True` if the world's heat sources act inside the layer during a thermal EOS solve (see [Worlds](../worlds/worlds.md)). |

The EOS solve reads all three, so writing one on a layer of a solved world leaves the world unsolved until its next `solve_eos` (see [Solved State](../worlds/worlds.md#solved-state)).

### Layer Assumptions

These three flags decide which equations the radial solver uses inside this layer. They are constructor arguments, layer keys in a world TOML (see [TOML schema](../config/toml_schema.md)), and writable after construction. A liquid layer is static unless `is_static` is set `False`.

| Property | Meaning |
|---|---|
| `is_solid` | Solid rather than liquid. A liquid layer has no shear strength and contributes fewer independent solutions to the radial solve. |
| `is_static` | The quasi-static assumption, dropping the inertial terms. See Beuthe (2015). |
| `is_incompressible` | The incompressible assumption. |

```python
mantle.is_incompressible = True    # e.g. to use the propagation-matrix method
mantle.is_incompressible = False
```

## Methods

### `set_eos(model)`

Attach a [material EOS model](../../Material/material_eos.md), the layer's material and density source. Ownership of the C++ model transfers into the layer; the passed wrapper becomes an empty shell, and attaching it again raises `ValueError`. The layer's viscosity and partial-melt models are held by its material, so replacing the material keeps the ones attached before unless the new model carries its own. On a layer of a solved world, a new material leaves the world unsolved until its next `solve_eos` (see [Solved State](../worlds/worlds.md#solved-state)).

```python
from TidalPy.Material.eos import make_material_eos

crust = BaseLayer("crust", 2, 6.371e6, 6.400e6, 1.0e22)
crust.set_eos(make_material_eos("birch_murnaghan", {
    "reference_density_kg_m3": 4500.0,
    "reference_bulk_modulus_pa": 2.5e11,
    "bulk_modulus_derivative": 4.0,
}))
```

### `set_shear_rheology(rheology)` / `set_bulk_rheology(rheology)`

Attach a rheology model (a `RheologyBase` subclass such as `Maxwell()` or `make_rheology("andrade")`) used to compute the complex shear or bulk modulus. Ownership of the underlying C++ model is transferred into the layer; the passed wrapper becomes an empty shell and must not be reused (attaching it again raises `ValueError`).

```python
from TidalPy.Rheology import Maxwell, make_rheology

mantle.set_shear_rheology(Maxwell())
mantle.set_bulk_rheology(make_rheology("andrade", {"alpha": 0.3}))
```

### `set_shear_viscosity(model)` / `set_bulk_viscosity(model)` / `set_partial_melt(model)`

These helpers hand the models to the layer's EOS (material) rather than store it on the layer.

```python
from TidalPy.Viscosity import make_viscosity
from TidalPy.PartialMelt import make_partial_melt

mantle.set_shear_viscosity(make_viscosity("reference", {
    "reference_viscosity_pas": 1.0e19,
    "reference_temperature_k": 1400.0}))
mantle.set_partial_melt(make_partial_melt("henning"))
```

A viscosity model from [`Viscosity`](../../Viscosity/viscosity_models.md) turns the temperature and pressure into the viscosity the rheology then uses. Without one the material falls back to its static viscosity, which is NaN unless specifically set. A partial-melt model from [`PartialMelt`](../../PartialMelt/partial_melt_models.md) weakens the modulus and the viscosity between the solidus and the liquidus.

### `calc_complex_shear_modulus(frequency)` -> complex

Complex shear modulus \[Pa\] at the given tidal forcing frequency \[rad s$^{-1}$\], from the material's static constants. With a shear rheology attached it is the complex modulus $\mu^*(\omega)$ that model returns for the static shear modulus and viscosity; without one it is `shear_modulus_static + 0j`.

### `calc_complex_shear_modulus(radius, frequency)` -> complex or ndarray

Radius-resolved form: applies the shear rheology to the static modulus and viscosity the solved EOS reports at `radius`, exactly like the world-level [`BaseWorld.calc_complex_shear_modulus`](../worlds/worlds.md). `radius` may be a float (returns `complex`) or an `np.ndarray` of radii (returns a same-shape complex array). Returns `NaN` before the world EOS solve populates the layer.

### `calc_complex_bulk_modulus(...)` -> complex or ndarray

Complex bulk modulus \[Pa\], in both the material-constant `(frequency)` and the radius-resolved `(radius, frequency)` forms, following the same rules as `calc_complex_shear_modulus`.

### `update_eos_data(radius, density, gravity, pressure)`

Populate the density, gravity, and pressure profile directly from sorted arrays. The world EOS solve does this in the normal workflow. This direct form is for tests and manual construction. All sequences must be the same length and `radius` sorted ascending. Values are interpolated linearly and clamped at the layer boundaries.

### Profile Getters -> float or ndarray

| Getter | Returns |
|---|---|
| `get_density`, `get_gravity`, `get_pressure` | Density \[kg m$^{-3}$\], gravity \[m s$^{-2}$\], and pressure \[Pa\]. |
| `get_shear_modulus`, `get_bulk_modulus` | Static moduli \[Pa\] after melt weakening. |
| `get_shear_viscosity`, `get_bulk_viscosity` | Viscosities \[Pa s\] after melt weakening. |
| `get_melt_fraction` | Melt fraction from the attached partial-melt model; `0.0` without one. |
| `get_static_viscoelastics(radius)` | The post-melt four-tuple in one call. |
| `get_state(radius)` | Every profile at that radius as a dict. |

Nothing here is stored on a grid, and nothing is calculated by the layer. The layer's material evaluated these properties while the structure was integrated, and each getter reads them back from that solved EOS at the radius asked for, so a getter and the solve always agree. A profile supplied through `update_eos_data` carries density, gravity, and pressure alone, so the viscoelastic getters report `NaN` for it.

Every getter accepts a float or an `np.ndarray` of radii and returns a matching scalar or same-shape array, evaluated in a C loop:

```python
import numpy as np

mantle.update_eos_data([3.485e6, 6.371e6], [5560.0, 3300.0], [10.68, 9.81], [1.36e11, 0.0])
radii = np.linspace(3.5e6, 6.3e6, 100)
rho = mantle.get_density(radii)        # ndarray, shape (100,)
```

A layer that belongs to a world reads its profile under that world's call lock, so a getter on a layer view takes turns with the world's `solve_eos` on another thread, as the world's own getters do (see the threading note under [Equation of State](../worlds/worlds.md#equation-of-state)). One call reads a whole array under a single turn; in C++ that is `c_BaseLayer::get_eos_fields(field_indices, num_fields, radii, num_radii, values_out)`, and `c_BaseLayer::calc_complex_moduli` is the matching complex-modulus form. The setters of what the world's solves read (the material, the models, the flags, and the layer state) take the same turns. A standalone layer has no lock to take.

### Inherited Geometry Calculations

Pure-function methods from `StructureBase` that do not depend on stored state:

```python
mantle.calc_surface_area(6.371e6)                # 4πr² [m²]
mantle.calc_volume_shell(6.371e6, 3.485e6)       # shell volume [m³]
mantle.calc_escape_velocity(5.972e24, 6.371e6)   # √(2Gm/r) [m/s]
```

`calc_volume_sphere`, `calc_surface_gravity`, and `calc_mean_density` follow the same pattern.

### Configuration Dict

```python
cfg = mantle.get_config_dict()   # every construction parameter and attached model
mantle.save_config("mantle.toml")
```

The dict follows the world builder's layer schema. `class` names the layer class (`base`, `solidliquid`, or `gas`), the scalar keys are the constructor parameters (`temperature_k` for the temperature), and each attached model is a sub-table keyed by `model`. The `material` table is the EOS model with its static constants, its shear law, and its own `shear_viscosity`, `bulk_viscosity`, and `partial_melt` tables; `shear_rheology` and `bulk_rheology` are the rheologies. `name` and `radius_inner_m` belong to a standalone layer only; a world drops them when it nests the layer under its name (`LAYER_STANDALONE_CONFIG_KEYS`).

### Tidal Bookkeeping

| Member | Description |
|---|---|
| `get_tidal_heating()` | Tidal heating deposited in this layer \[W\], set by the world's tidal solve. |
| `tidal_scale` | The layer's share of the planet in the quasi-homogeneous Love methods and of an analytic tide model's heating (its volume fraction when unset). See [Worlds](../worlds/worlds.md). |
| `is_tidal` | Whether the layer takes any tidal heating at all. |

## Example

```python
import math
from TidalPy.Material.eos import ConstantDensityEOS
from TidalPy.Rheology import Maxwell
from TidalPy.Structures.layers import BaseLayer

mantle = BaseLayer(
    name          = "mantle",
    layer_index   = 1,
    radius_inner  = 3.485e6,
    radius_outer  = 6.371e6,
    mass          = 4.043e24,
    material_name = "perovskite",
)
# The material: density law, static moduli, static viscosities.
mantle.set_eos(ConstantDensityEOS(
    reference_density      = 4500.0,
    shear_modulus_static   = 1.67e11,
    bulk_modulus_static    = 3.57e11,
    shear_viscosity_static = 1.0e21,
    bulk_viscosity_static  = 2.0e21,
))

freq = 2.0 * math.pi / (1.77 * 86400.0)   # Io's orbital frequency [rad/s]
mu   = mantle.calc_complex_shear_modulus(freq)
print(f"Thickness:              {mantle.thickness / 1e3:.0f} km")
print(f"Complex shear modulus:  {mu.real:.3e} + {mu.imag:.3e}j Pa")   # no rheology has been set yet: imaginary part 0.0

mantle.set_shear_rheology(Maxwell())
mu = mantle.calc_complex_shear_modulus(freq)
print(f"With a Maxwell rheology: {mu.real:.3e} + {mu.imag:.3e}j Pa")
```
