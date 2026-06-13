# PhysicsLayer

`TidalPy.structures_x.layers.PhysicsLayer`

## Overview

`PhysicsLayer` extends `BaseLayer` with the static mechanical properties needed
for tidal calculations: (shear/bulk) modulus and dynamic (shaer/bulk) viscosity.

When a rheology model (`c_RheologyBase` subclass, Phase 5) is attached via the
C++ `set_shear/bulk_rheology` methods, `calc_complex_shear/bulk_modulus` returns
the frequency-dependent complex modulus from the constitutive law. Until then
(or when no rheology is set), the methods return the static modulus as a purely
real complex number, equivalent to perfectly elastic behavior.

All values are in **MKS units** (meters, kilograms, seconds, pascals).

---

## Inheritance

```
TidalPyBaseClass
  └── StructureBase
        └── BaseLayer
              └── PhysicsLayer
```

---

## Constructor

```python
PhysicsLayer(
    name:                       str,
    layer_index:                int,
    radius_inner_m:             float,
    radius_outer_m:             float,
    mass_kg:                    float,
    material_name:              str     = "",
    is_tidal:                   bool    = True,
    tidal_scale:                float   = 1.0,
    shear_modulus_static_pa:    float   = 0.0,
    bulk_modulus_static_pa:     float   = 0.0,
    shear_viscosity_static_pas: float   = 0.0,
    bulk_viscosity_static_pas:  float   = 0.0,
    love_number_k:              complex = 0+0j,
    love_number_h:              complex = 0+0j,
    love_number_l:              complex = 0+0j,
)
```

### Parameters

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| `name` | `str` | — | Human-readable layer name. |
| `layer_index` | `int` | — | Zero-based index; innermost layer = 0. |
| `radius_inner_m` | `float` | m | Inner boundary radius. |
| `radius_outer_m` | `float` | m | Outer boundary radius. |
| `mass_kg` | `float` | kg | Total layer mass. |
| `material_name` | `str` | — | Material identifier. Optional. |
| `is_tidal` | `bool` | — | Tidal dissipation flag. Default `True`. |
| `tidal_scale` | `float` | — | Dimensionless tidal heating scale. Default `1.0`. |
| `shear_modulus_static_pa` | `float` | Pa | Unrelaxed shear modulus. Default `0.0`. |
| `bulk_modulus_static_pa` | `float` | Pa | Unrelaxed bulk modulus. Default `0.0`. |
| `shear_viscosity_static_pas` | `float` | Pa·s | Reference shear viscosity. Default `0.0`. |
| `bulk_viscosity_static_pas` | `float` | Pa·s | Reference bulk viscosity. Default `0.0`. |
| `love_number_k` | `complex` | — | Potential Love number k (placeholder). Default `0+0j`. |
| `love_number_h` | `complex` | — | Radial displacement Love number h (placeholder). Default `0+0j`. |
| `love_number_l` | `complex` | — | Tangential displacement Love number l (placeholder). Default `0+0j`. |

---

## Properties

### Inherited from BaseLayer (read-only, MKS)

See [BaseLayer](base_layer.md) for the full list: `name`, `layer_index`,
`radius`, `radius_inner`, `radius_outer`, `thickness`, `mass`, `volume`,
`surface_area_inner`, `surface_area_outer`, `material_name`, `is_tidal`,
`tidal_scale`, `eos_data_populated`.

### Mechanical (read-only, MKS)

| Property | Units | Description |
|----------|-------|-------------|
| `shear_modulus_static` | Pa | Unrelaxed shear modulus. |
| `bulk_modulus_static` | Pa | Unrelaxed bulk modulus. |
| `shear_viscosity_static` | Pa·s | Reference shear viscosity. |
| `bulk_viscosity_static` | Pa·s | Reference bulk viscosity. |
| `love_numbers` | — | All three Love numbers as a `LoveNumbers` object. |
| `love_number_k` | — | Potential Love number k. Returns `complex`. |
| `love_number_h` | — | Radial displacement Love number h. Returns `complex`. |
| `love_number_l` | — | Tangential displacement Love number l. Returns `complex`. |
| `shear_rheology_set` | — | `True` after a shear rheology model is attached. |
| `bulk_rheology_set` | — | `True` after a bulk rheology model is attached. |

---

## Methods

### `calc_tidal_susceptibility()` → float

Geometrical tidal susceptibility [m³]: (3/2) · r⁵ / (G · m²).

Newton's G is sourced from the TidalPy global configuration. Returns `0.0`
if the config is uninitialised or mass is zero.

```python
chi = mantle.calc_tidal_susceptibility()
print(f"Tidal susceptibility: {chi:.3e} m³")
```

**References:** Kaula (1964); Eggleton et al. (1998).

### `calc_complex_shear_modulus(frequency_rad_s)` → complex

Complex shear modulus [Pa] at the given tidal forcing frequency.

When a shear rheology model is attached the result is 1 / J(ω), where J(ω)
is the complex compliance returned by the rheology model. Without a rheology
model the return value is `shear_modulus_static + 0j`.

```python
mu = mantle.calc_complex_shear_modulus(2.0 * math.pi / 86400.0)
print(f"Re(μ) = {mu.real:.3e} Pa,  Im(μ) = {mu.imag:.3e} Pa")
```

### `calc_complex_bulk_modulus(frequency_rad_s)` → complex

Complex bulk modulus [Pa] at the given tidal forcing frequency.  Same
delegation logic as `calc_complex_shear_modulus`.

### Inherited from BaseLayer

`update_eos_data`, `get_density`, `get_gravity`, `get_pressure`,
`calc_surface_area`, `calc_volume_sphere`, `calc_volume_shell`,
`calc_surface_gravity`, `calc_mean_density`, `calc_escape_velocity`,
`save_binary`, `load_binary`, `save_config`, `get_config_dict`.

---

## Binary serialization

`save_binary` / `load_binary` serialize all `BaseLayer` fields (see
[BaseLayer](base_layer.md)) followed by ten doubles in order:
`shear_modulus_static_pa`, `bulk_modulus_static_pa`, `shear_viscosity_static_pas`,
`bulk_viscosity_static_pas`, then `love_number_k` re+im, `love_number_h` re+im,
`love_number_l` re+im (6 doubles total for the Love numbers).

Binary class ID: **101** (`BinaryClassID::PhysicsLayer`).

**Rheology objects and EOS profile data are NOT serialized.**

---

## Example

```python
import math
from TidalPy.structures_x.layers import PhysicsLayer

mantle = PhysicsLayer(
    name                        = "mantle",
    layer_index                 = 1,
    radius_inner_m              = 3.485e6,
    radius_outer_m              = 6.371e6,
    mass_kg                     = 4.043e24,
    material_name               = "perovskite",
    shear_modulus_static_pa     = 1.67e11,
    bulk_modulus_static_pa      = 3.57e11,
    shear_viscosity_static_pas  = 1.0e21,
    bulk_viscosity_static_pas   = 2.0e21,
)

freq = 2.0 * math.pi / (1.77 * 86400.0)   # Io's orbital frequency [rad/s]
mu   = mantle.calc_complex_shear_modulus(freq)
chi  = mantle.calc_tidal_susceptibility()

print(f"Thickness:              {mantle.thickness / 1e3:.0f} km")
print(f"Tidal susceptibility:   {chi:.3e} m³")
print(f"Complex shear modulus:  {mu.real:.3e} + {mu.imag:.3e}j Pa")
# Rheology not yet set, so imaginary part is 0.0
```
