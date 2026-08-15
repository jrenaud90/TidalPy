# SolidLiquidLayer

`TidalPy.structures_x.layers.SolidLiquidLayer`

## Overview

`SolidLiquidLayer` extends `PhysicsLayer` with thermomechanical behavior:
phase-change tracking, Arrhenius viscosity, dynamic shear modulus, thermal
transport, and optional sub-model hooks for radiogenic heating and
convective/conductive cooling.

Key physics:

- **Melt fraction** — power-law interpolation between solidus and liquidus:
  φ = clamp((T − T_s) / (T_l − T_s), 0, 1)^n

- **Arrhenius viscosity** with pressure correction and partial-melt reduction:
  η = η_ref · exp((E_a + P·V_a)/(R·T) − E_a/(R·T_ref)) · exp(−C · φ)

- **Dynamic shear modulus** — G_eff = G_static · (1 − φ)

- **Thermal transport** — constant-k conductivity; diffusivity κ = k/(ρ_ref·c_p);
  adiabatic gradient α·T·g/c_p (requires EOS data for gravity).

- **Conductive heat flux** — F = k · (T_base − T_top) / h.

The solidus/liquidus temperatures are constant: the melt curve carries no
pressure dependence.

All values are in **MKS units** (meters, kilograms, seconds, pascals, kelvin).

---

## Inheritance

```
TidalPyBaseClass
  └── StructureBase
        └── BaseLayer
              └── PhysicsLayer
                    └── SolidLiquidLayer
```

---

## Constructor

```python
SolidLiquidLayer(
    name:                          str,
    layer_index:                   int,
    radius_inner_m:                float,
    radius_outer_m:                float,
    mass_kg:                       float,
    material_name:                 str   = "",
    is_tidal:                      bool  = True,
    tidal_scale:                   float = 1.0,
    shear_modulus_static_pa:       float = 0.0,
    bulk_modulus_static_pa:        float = 0.0,
    viscosity_static_pas:          float = 0.0,
    love_number_re:                float = 0.0,
    love_number_im:                float = 0.0,
    thermal_conductivity_ref_w_mk: float = 4.0,
    thermal_expansion_ref_1_k:     float = 3e-5,
    heat_capacity_ref_j_kgk:       float = 1200.0,
    activation_energy_j_mol:       float = 300e3,
    activation_volume_m3_mol:      float = 5e-6,
    solidus_temperature_k:         float = 1600.0,
    liquidus_temperature_k:        float = 2000.0,
    melt_fraction_exponent:        float = 1.0,
    reference_density_kg_m3:       float = 3500.0,
    reference_temperature_k:       float = 1600.0,
    melt_viscosity_reduction:      float = 25.0,
)
```

### Parameters

| Parameter | Units | Description |
|-----------|-------|-------------|
| `name` | — | Human-readable layer name. |
| `layer_index` | — | Zero-based index; innermost layer = 0. |
| `radius_inner_m` | m | Inner boundary radius. |
| `radius_outer_m` | m | Outer boundary radius. |
| `mass_kg` | kg | Total layer mass. |
| `material_name` | — | Material identifier. Default `""`. |
| `is_tidal` | — | Tidal dissipation flag. Default `True`. |
| `tidal_scale` | — | Dimensionless tidal heating scale. Default `1.0`. |
| `shear_modulus_static_pa` | Pa | Unrelaxed shear modulus. Default `0.0`. |
| `bulk_modulus_static_pa` | Pa | Unrelaxed bulk modulus. Default `0.0`. |
| `viscosity_static_pas` | Pa·s | Reference dynamic viscosity (at `reference_temperature_k`, P=0). Default `0.0`. |
| `love_number_re` | — | Real part of complex Love number (placeholder). Default `0.0`. |
| `love_number_im` | — | Imaginary part of complex Love number (placeholder). Default `0.0`. |
| `thermal_conductivity_ref_w_mk` | W/(m·K) | Reference thermal conductivity. Default `4.0`. |
| `thermal_expansion_ref_1_k` | 1/K | Reference thermal expansion coefficient α. Default `3e-5`. |
| `heat_capacity_ref_j_kgk` | J/(kg·K) | Reference specific heat capacity c_p. Default `1200.0`. |
| `activation_energy_j_mol` | J/mol | Arrhenius activation energy E_a. Default `300e3`. |
| `activation_volume_m3_mol` | m³/mol | Arrhenius activation volume V_a. Default `5e-6`. |
| `solidus_temperature_k` | K | Solidus temperature T_s. Default `1600.0`. |
| `liquidus_temperature_k` | K | Liquidus temperature T_l. Default `2000.0`. |
| `melt_fraction_exponent` | — | Exponent n in the melt-fraction formula. Default `1.0` (linear). |
| `reference_density_kg_m3` | kg/m³ | Reference density ρ_ref for thermal diffusivity. Default `3500.0`. |
| `reference_temperature_k` | K | Reference temperature T_ref for Arrhenius viscosity. Default `1600.0`. |
| `melt_viscosity_reduction` | — | Coefficient C in exp(−C·φ) melt-viscosity reduction. Default `25.0`. |

---

## Properties

### Inherited from PhysicsLayer

See [PhysicsLayer](physics_layer.md): `shear_modulus_static`, `bulk_modulus_static`,
`viscosity_static`, `love_number`, `shear_rheology_set`, `bulk_rheology_set`.

### Inherited from BaseLayer

See [BaseLayer](base_layer.md): `name`, `layer_index`, `radius`, `radius_inner`,
`radius_outer`, `thickness`, `mass`, `volume`, `surface_area_inner`,
`surface_area_outer`, `material_name`, `is_tidal`, `tidal_scale`,
`eos_data_populated`.

### Thermal / melt (read-only)

| Property | Units | Description |
|----------|-------|-------------|
| `thermal_conductivity_ref` | W/(m·K) | Reference thermal conductivity. |
| `thermal_expansion_ref` | 1/K | Reference thermal expansion coefficient α. |
| `heat_capacity_ref` | J/(kg·K) | Reference specific heat capacity c_p. |
| `activation_energy` | J/mol | Arrhenius activation energy E_a. |
| `activation_volume` | m³/mol | Arrhenius activation volume V_a. |
| `solidus_temperature` | K | Solidus temperature T_s. |
| `liquidus_temperature` | K | Liquidus temperature T_l. |
| `melt_fraction_exponent` | — | Melt-fraction exponent n. |
| `reference_density` | kg/m³ | Reference density for thermal diffusivity. |
| `reference_temperature` | K | Reference temperature for Arrhenius viscosity. |
| `melt_viscosity_reduction` | — | Exponential melt-viscosity reduction coefficient C. |
| `cooling_set` | — | `True` after a cooling sub-model is attached. |
| `radiogenics_set` | — | `True` after a radiogenics sub-model is attached. |

---

## Methods

### `calc_melt_fraction(temperature_k, pressure_pa=0.0)` → float

Volumetric melt fraction φ ∈ [0, 1]:

```
τ = clamp((T − T_s) / (T_l − T_s), 0, 1)
φ = τ^n
```

`pressure_pa` is accepted for interface uniformity and is unused: the melt
curve carries no pressure dependence.

```python
phi = layer.calc_melt_fraction(3200.0)        # T = 3200 K, P = 0
```

### `calc_viscosity(temperature_k, pressure_pa=0.0)` → float

Effective dynamic viscosity [Pa·s]:

```
η = η_ref · exp(clamp((E_a + P·V_a)/(R·T) − E_a/(R·T_ref), −100, 100))
          · exp(clamp(−C · φ, −100, 0))
```

The Arrhenius exponent is clamped to [−100, 100] to prevent overflow.
Returns `η_ref` when T = 0 K.

```python
eta = layer.calc_viscosity(3000.0, 1e11)      # T = 3000 K, P = 100 GPa
```

### `calc_shear_modulus(temperature_k, pressure_pa=0.0)` → float

Effective shear modulus [Pa] accounting for partial melt:

```
G_eff = G_static · (1 − φ)
```

```python
G = layer.calc_shear_modulus(3000.0)
```

### `calc_thermal_conductivity(temperature_k)` → float

Returns the reference thermal conductivity k [W/(m·K)].
Temperature dependence is not modeled.

### `calc_thermal_diffusivity(temperature_k)` → float

Thermal diffusivity [m²/s] = k / (ρ_ref · c_p).

### `calc_adiabatic_temperature_gradient(temperature_k, pressure_pa=0.0)` → float

Adiabatic temperature gradient [K/m] = α · T · g / c_p.

Gravity g is read from the EOS profile at the outer radius. Returns `0.0`
when EOS data has not been populated via `update_eos_data`.

```python
layer.update_eos_data(radii, densities, gravities, pressures)
grad = layer.calc_adiabatic_temperature_gradient(3000.0)
```

### `calc_heat_flux_conductive(temperature_base_k, temperature_top_k)` → float

Conductive heat flux [W/m²]:

```
F = k · (T_base − T_top) / h
```

where h is the layer thickness. Returns `0.0` for zero-thickness layers.

```python
flux = layer.calc_heat_flux_conductive(T_base=3500.0, T_top=1500.0)
```


### `set_cooling(cooling)` / `set_radiogenics(radiogenics)`

Attach a cooling (`CoolingBase`) or radiogenics (`RadiogenicsBase`) sub-model.
Ownership of the underlying C++ model is **transferred** into the layer; the passed
Python wrapper becomes an empty, non-owning shell and must not be reused (raises
`ValueError` if re-attached). Shear/bulk rheology are attached via the inherited
`set_shear_rheology` / `set_bulk_rheology` (see [PhysicsLayer](physics_layer.md)).

```python
from TidalPy.cooling_x import make_cooling
from TidalPy.radiogenics_x import IsotopeRadiogenics
layer.set_cooling(make_cooling("convection"))
layer.set_radiogenics(IsotopeRadiogenics.from_dataset("modern_day_chondritic"))
```


### `calc_radiogenic_heating(time_s, mass_kg)` → float

Radiogenic heating power [W] from the attached sub-model.
Returns `0.0` when no radiogenics sub-model has been attached.

### Inherited from PhysicsLayer / BaseLayer

`calc_tidal_susceptibility`, `calc_complex_shear_modulus`,
`calc_complex_bulk_modulus`, `update_eos_data`, `get_density`, `get_gravity`,
`get_pressure`, `calc_surface_area`, `calc_volume_sphere`, `calc_volume_shell`,
`calc_surface_gravity`, `calc_mean_density`, `calc_escape_velocity`,
`save_binary`, `load_binary`, `save_config`, `get_config_dict`.

---

## `get_config_dict()` → dict

Returns all configuration values as a Python dictionary (MKS).  Includes all
`BaseLayer` + `PhysicsLayer` keys plus the 11 SolidLiquidLayer parameters.

| Key | Units | Description |
|-----|-------|-------------|
| `thermal_conductivity_ref_w_mk` | W/(m·K) | Reference thermal conductivity. |
| `thermal_expansion_ref_1_k` | 1/K | Reference thermal expansion coefficient. |
| `heat_capacity_ref_j_kgk` | J/(kg·K) | Reference specific heat capacity. |
| `activation_energy_j_mol` | J/mol | Arrhenius activation energy. |
| `activation_volume_m3_mol` | m³/mol | Arrhenius activation volume. |
| `solidus_temperature_k` | K | Solidus temperature. |
| `liquidus_temperature_k` | K | Liquidus temperature. |
| `melt_fraction_exponent` | — | Melt-fraction exponent. |
| `reference_density_kg_m3` | kg/m³ | Reference density. |
| `reference_temperature_k` | K | Reference temperature. |
| `melt_viscosity_reduction` | — | Melt-viscosity reduction coefficient. |

---

## Binary serialization

`save_binary` / `load_binary` serialize fields in this order:

1. All `BaseLayer` fields (name, geometry, mass).
2. All `PhysicsLayer` mechanical fields (G, K, η, love number re+im).
3. The 11 SolidLiquidLayer doubles in constructor order.
4. An optional sub-model section: presence flags + recursive binary records for the
   shear rheology, bulk rheology, cooling, and radiogenics models (in that order).

On load, every attached sub-model is reconstructed recursively via the
rheology / cooling / radiogenics binary-dispatch factories, so a saved layer
round-trips with all of its physics intact (verify with `shear_rheology_set`,
`cooling_set`, `radiogenics_set`, `calc_complex_shear_modulus`, and
`calc_radiogenic_heating`). See [Binary serialization](../../utilities_x/binary_x.md)
for the encoding.

Binary class ID: **102** (`BinaryClassID::SolidLiquidLayer`).

**EOS profile data is NOT serialized** (repopulate it after loading by running the
EOS handler).

---

## Example

```python
import math
from TidalPy.structures_x.layers import SolidLiquidLayer

mantle = SolidLiquidLayer(
    name                          = "mantle",
    layer_index                   = 1,
    radius_inner_m                = 3.485e6,
    radius_outer_m                = 6.371e6,
    mass_kg                       = 4.043e24,
    material_name                 = "perovskite",
    shear_modulus_static_pa       = 1.67e11,
    bulk_modulus_static_pa        = 3.57e11,
    viscosity_static_pas          = 1.0e21,
    thermal_conductivity_ref_w_mk = 4.5,
    heat_capacity_ref_j_kgk       = 1200.0,
    thermal_expansion_ref_1_k     = 2.0e-5,
    activation_energy_j_mol       = 300.0e3,
    activation_volume_m3_mol      = 5.0e-6,
    solidus_temperature_k         = 3000.0,
    liquidus_temperature_k        = 4000.0,
    reference_density_kg_m3       = 4000.0,
    reference_temperature_k       = 3000.0,
    melt_viscosity_reduction      = 25.0,
)

T = 3200.0    # K — slightly above solidus
P = 1e11      # Pa — ~100 GPa

phi = mantle.calc_melt_fraction(T, P)
eta = mantle.calc_viscosity(T, P)
G   = mantle.calc_shear_modulus(T, P)
k   = mantle.calc_thermal_conductivity(T)
kap = mantle.calc_thermal_diffusivity(T)
F   = mantle.calc_heat_flux_conductive(T_base=3500.0, T_top=1500.0)

print(f"Thickness:             {mantle.thickness / 1e3:.0f} km")
print(f"Melt fraction:         {phi:.3f}")
print(f"Effective viscosity:   {eta:.3e} Pa·s")
print(f"Effective shear mod:   {G:.3e} Pa")
print(f"Thermal conductivity:  {k:.2f} W/(m·K)")
print(f"Thermal diffusivity:   {kap:.3e} m²/s")
print(f"Conductive heat flux:  {F:.3f} W/m²")

# Binary save/load
mantle.save_binary("mantle.tpyb")
restored = SolidLiquidLayer("placeholder", 0, 0.0, 1.0, 1.0)
restored.load_binary("mantle.tpyb")
assert restored.solidus_temperature == mantle.solidus_temperature
```
