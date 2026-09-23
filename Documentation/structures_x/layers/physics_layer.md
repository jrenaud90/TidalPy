# PhysicsLayer

_Updated: 2026-09-23_

`TidalPy.structures_x.layers.PhysicsLayer` extends `BaseLayer` with what a tidal calculation needs from the layer itself: the radial-solver assumptions, the layer temperature, the Love numbers, and the shear and bulk rheology.

The material does not live here. The static moduli, the shear law, the viscosities, and the viscosity and partial-melt models all belong to the layer's EOS model (see [Material EOS Models](../../material_x/material_eos.md)), which the whole-planet EOS solve evaluates as it integrates. There is one path for those properties, so a layer and the solve that used it cannot disagree. After a solve they are read back with getters; nothing is recalculated.

The rheology is the one thing that knows a forcing frequency. When a rheology model (a `RheologyBase` subclass) is attached with `set_shear_rheology` or `set_bulk_rheology`, `calc_complex_shear_modulus` and `calc_complex_bulk_modulus` apply it to the static modulus and viscosity the solved EOS reports. Until then they return the static modulus as a purely real complex number, which is perfectly elastic behavior.

## Inheritance

```
TidalPyBaseClass
  └── StructureBase
        └── BaseLayer
              └── PhysicsLayer
```

## Constructor

```python
PhysicsLayer(
    name:               str,
    layer_index:        int,
    radius_inner:       float,
    radius_outer:       float,
    mass:               float,
    material_name:      str     = "",
    is_tidal:           bool    = True,
    is_volume_fixed:    bool    = True,
    tidal_scale:        float   = None,
    love_number_k:      complex = 0+0j,
    love_number_h:      complex = 0+0j,
    love_number_l:      complex = 0+0j,
    is_solid:           bool    = True,
    is_static:          bool    = True,
    is_incompressible:  bool    = False,
    temperature:        float   = 0.0,
    use_thermal_eos:    bool    = False,
    use_heating:        bool    = False,
)
```

### Parameters

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| `name` | `str` | - | Human-readable layer name. |
| `layer_index` | `int` | - | Zero-based index; innermost layer = 0. |
| `radius_inner` | `float` | m | Inner boundary radius. |
| `radius_outer` | `float` | m | Outer boundary radius. |
| `mass` | `float` | kg | Total layer mass. Overwritten by each successful world EOS solve. |
| `material_name` | `str` | - | Material identifier. Optional. |
| `is_tidal` | `bool` | - | Tidal dissipation flag. Default `True`. |
| `is_volume_fixed` | `bool` | - | `False` lets the layer grow or shrink to hold its mass during an EOS solve. Default `True`. |
| `tidal_scale` | `float` | - | The layer's share of the planet in the quasi-homogeneous Love methods (`homogeneous`, `cpl`, `ctl`) and of an analytic tide model's heating; `None` (default) takes its volume fraction. See [Worlds](../worlds/worlds.md). |
| `love_number_k` | `complex` | - | Potential Love number k (placeholder). Default `0+0j`. |
| `love_number_h` | `complex` | - | Radial displacement Love number h (placeholder). Default `0+0j`. |
| `love_number_l` | `complex` | - | Tangential displacement Love number l (placeholder). Default `0+0j`. |
| `is_solid`, `is_static`, `is_incompressible` | `bool` | - | Radial-solver assumptions; see Layer Assumptions below. Defaults `True`, `True`, `False`. |
| `temperature` | `float` | K | Layer temperature at which the material's viscosity and melt models are evaluated. Default `0.0`, the cold rigid limit of the viscosity laws. |
| `use_thermal_eos` | `bool` | - | Let the material's density law see the temperature, so the density and bulk modulus depend on it. Default `False`. |
| `use_heating` | `bool` | - | Let the world's heat sources act inside the layer during a thermal EOS solve. Default `False`. |

The static moduli and viscosities and the shear law are arguments of the EOS model, not of the layer:

```python
from TidalPy.Material_x.eos import ConstantDensityEOS

mantle = PhysicsLayer("mantle", 1, 3.485e6, 6.371e6, 4.043e24)
mantle.set_eos(ConstantDensityEOS(
    reference_density=4500.0, shear_modulus_static=1.67e11, bulk_modulus_static=3.57e11))
```

## Properties

### Inherited from BaseLayer

_Read-only properties._

See [BaseLayer](base_layer.md) for the full list: `name`, `layer_index`, `radius`, `radius_inner`, `radius_outer`, `thickness`, `mass`, `volume`, `density_bulk`, `surface_area_inner`, `surface_area_outer`, `material_name`, `is_tidal`, `tidal_scale`, `eos_data_populated`.

### Mechanical

_Read-only properties._

| Property | Units | Description |
|----------|-------|-------------|
| `shear_modulus_static` | Pa | The material's unrelaxed shear modulus, read from the attached EOS model. NaN when none is attached. |
| `bulk_modulus_static` | Pa | The material's unrelaxed bulk modulus, as above. |
| `shear_viscosity_static` | Pa·s | The material's static shear viscosity, as above. |
| `bulk_viscosity_static` | Pa·s | The material's static bulk viscosity, as above. |
| `love_numbers` | - | All three Love numbers as a `LoveNumbers` object. |
| `love_number_k` | - | Potential Love number k. Returns `complex`. |
| `love_number_h` | - | Radial displacement Love number h. Returns `complex`. |
| `love_number_l` | - | Tangential displacement Love number l. Returns `complex`. |
| `shear_rheology_set`, `bulk_rheology_set` | - | `True` after the corresponding rheology model is attached. |
| `shear_viscosity_set`, `bulk_viscosity_set` | - | `True` when the material holds the corresponding viscosity model. |
| `partial_melt_set` | - | `True` when the material holds a partial-melt model. |

### Layer State

| Property | Units | Description |
|----------|-------|-------------|
| `temperature` | K | Layer temperature. Writable. |
| `use_thermal_eos` | - | `True` if the material's density law receives the temperature. Writable. |
| `use_heating` | - | `True` if the world's heat sources act inside the layer during a thermal EOS solve (see [Worlds](../worlds/worlds.md)). Writable. |

### Layer Assumptions

These three flags decide which equations the radial solver uses inside this layer. They are constructor arguments, layer keys in a world TOML (see [TOML schema](../config/toml_schema.md)), and writable after construction. A liquid layer is static unless `is_static` is set `False`.

| Property | Meaning |
|---|---|
| `is_solid` | Solid rather than liquid. A liquid layer has no shear strength and contributes fewer independent solutions to the radial solve. |
| `is_static` | The quasi-static assumption, dropping the inertial terms. See Beuthe (2015). |
| `is_incompressible` | The incompressible assumption. |

```python
layer.is_incompressible = True    # e.g. to use the propagation-matrix method
```

## Methods

### `set_shear_rheology(rheology)` / `set_bulk_rheology(rheology)`

Attach a rheology model (a `RheologyBase` subclass such as `Maxwell()` or `make_rheology("andrade")`) used to compute the complex shear / bulk modulus. Ownership of the underlying C++ model is transferred into the layer; the passed Python wrapper becomes an empty, non-owning shell and must not be reused (attempting to attach it again raises `ValueError`).

```python
from TidalPy.rheology_x import Maxwell, make_rheology
mantle.set_shear_rheology(Maxwell())
mantle.set_bulk_rheology(make_rheology("andrade", {"alpha": 0.3}))
```

### `set_shear_viscosity(model)` / `set_bulk_viscosity(model)` / `set_partial_melt(model)`

Helpers. The material owns these models, so each call hands the model to the layer's EOS model rather than storing it on the layer; they exist so a layer can be configured in one place. Attach the EOS first (`set_eos`): with none there is no material to give the model to, and the call raises `ValueError`. The same methods are on the EOS model itself, for a material configured before it is attached.

```python
from TidalPy.viscosity_x import make_viscosity
from TidalPy.partial_melt_x import make_partial_melt

mantle.set_shear_viscosity(make_viscosity("reference", {
    "reference_viscosity_pas": 1.0e19,
    "reference_temperature_k": 1400.0}))
mantle.set_partial_melt(make_partial_melt("henning"))
```

A viscosity model from [`viscosity_x`](../../viscosity_x/viscosity_models.md) turns the temperature and pressure into the viscosity the rheology then uses; without one the material falls back to its static viscosity, which is NaN unless you set it. A partial-melt model from [`partial_melt_x`](../../partial_melt_x/partial_melt_models.md) weakens the modulus and the viscosity between the solidus and the liquidus.

### `calc_complex_shear_modulus(frequency)` -> complex

Complex shear modulus \[Pa\] at the given tidal forcing frequency, from the material's static constants.

When a shear rheology model is attached the result is the complex modulus $\mu^*(\omega)$ returned by that model (evaluated from the static shear modulus, static shear viscosity, and frequency). Without a rheology model the return value is `shear_modulus_static + 0j`.

```python
mu = mantle.calc_complex_shear_modulus(2.0 * math.pi / 86400.0)
print(f"Re(μ) = {mu.real:.3e} Pa,  Im(μ) = {mu.imag:.3e} Pa")
```

### `calc_complex_shear_modulus(radius, frequency)` -> complex or ndarray

Radius-resolved form: applies the shear rheology to the static modulus and viscosity the solved EOS reports at `radius`, exactly like the world-level [`LayeredWorld.calc_complex_shear_modulus`](../worlds/worlds.md). This is the only step of the chain that depends on frequency, and it is what the radial Love-number solve does at every radius it visits. `radius` may be a float (returns `complex`) or an `np.ndarray` of radii (returns a same-shape complex array). Returns `NaN` before the world EOS solve populates the layer.

```python
import numpy as np
radii = np.linspace(3.5e6, 6.3e6, 100)
mu_of_r = mantle.calc_complex_shear_modulus(radii, 2.0 * math.pi / 86400.0)
```

### `calc_complex_bulk_modulus(...)` -> complex or ndarray

Complex bulk modulus [Pa]; both the material-constant `(frequency)` and the radius-resolved `(radius, frequency)` forms, with the same delegation logic as `calc_complex_shear_modulus`.

### Inherited from `BaseLayer`

`update_eos_data`, `get_density`, `get_gravity`, `get_pressure`, the static material getters (`get_shear_modulus`, `get_bulk_modulus`, `get_shear_viscosity`, `get_bulk_viscosity`, `get_melt_fraction`, `get_state`), `calc_surface_area`, `calc_volume_sphere`, `calc_volume_shell`, `calc_surface_gravity`, `calc_mean_density`, `calc_escape_velocity`, `save_binary`, `load_binary`, `save_config`, `get_config_dict`.

`get_config_dict()` adds the three layer-assumption flags, `temperature_k`, `use_thermal_eos`, `use_heating`, the Love-number components, and a sub-table for each attached rheology (`shear_rheology`, `bulk_rheology`), keyed by `model` exactly as the world builder reads it. The `material` table (from `BaseLayer`) carries the EOS model with its static constants, its shear law, and its own `shear_viscosity`, `bulk_viscosity`, and `partial_melt` tables.

## Binary Serialization

`save_binary` / `load_binary` serialize all `BaseLayer` fields (see [BaseLayer](base_layer.md)) followed by six doubles for the Love numbers (`love_number_k` re+im, `love_number_h` re+im, `love_number_l` re+im), one byte each for `is_solid`, `is_static`, and `is_incompressible`, then one double for `temperature` and one byte each for `use_thermal_eos` and `use_heating`.

Following the scalar payload, an optional sub-model section is written: one-byte presence flags for the material EOS model and the shear and bulk rheology, each followed (when set) by that model's own binary record. The EOS record carries the material with it: the static constants, the shear law, and its viscosity and partial-melt models. On load, attached models are reconstructed recursively via each module's binary-dispatch factory, so a saved layer round-trips with its models intact (verify with `eos_set`, `shear_rheology_set`, `shear_viscosity_set`, and `partial_melt_set`). See [Binary serialization](../../utilities_x/binary_x.md) for the encoding.

Binary class id 101 (`BinaryClassID::PhysicsLayer`).

The EOS profile data is not serialized; re-run the world's `solve_eos` after loading.

## Example

```python
import math
from TidalPy.Material_x.eos import ConstantDensityEOS
from TidalPy.rheology_x import Maxwell
from TidalPy.structures_x.layers import PhysicsLayer

mantle = PhysicsLayer(
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
print(f"Complex shear modulus:  {mu.real:.3e} + {mu.imag:.3e}j Pa")   # no rheology yet: imaginary part 0.0

mantle.set_shear_rheology(Maxwell())
mu = mantle.calc_complex_shear_modulus(freq)
print(f"With a Maxwell rheology: {mu.real:.3e} + {mu.imag:.3e}j Pa")
```
