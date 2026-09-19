# SolidLiquidLayer

_Updated: 2026-09-19_

`TidalPy.structures_x.layers.SolidLiquidLayer` extends `PhysicsLayer` with thermomechanical behavior: melt-fraction tracking, an Arrhenius viscosity, a melt-weakened shear modulus, thermal transport, and optional sub-models for radiogenic heating and convective or conductive cooling.

The physics:

- Melt fraction: power-law interpolation between the solidus $T_{s}$ and the liquidus $T_{l}$, $\phi = \left[\operatorname{clamp}\!\left(\frac{T - T_{s}}{T_{l} - T_{s}}, 0, 1\right)\right]^{n}$.
- Arrhenius viscosity with a pressure correction and a partial-melt reduction, $\eta = \eta_\mathrm{ref}\exp\!\left(\frac{E_{a} + P V_{a}}{R T} - \frac{E_{a}}{R T_\mathrm{ref}}\right)\exp(-C\phi)$.
- Melt-weakened shear modulus, $G_\mathrm{eff} = G_\mathrm{static}(1 - \phi)$.
- Thermal transport: a constant conductivity $k$; the diffusivity $\kappa = k/(\rho_\mathrm{ref}\,c_{p})$; the adiabatic gradient $\alpha T g / c_{p}$ (requires EOS data for gravity).
- Conductive heat flux, $F = k\,(T_\mathrm{base} - T_\mathrm{top})/h$.

The solidus and liquidus temperatures are constant: the melt curve carries no pressure dependence.

## Inheritance

```
TidalPyBaseClass
  └── StructureBase
        └── BaseLayer
              └── PhysicsLayer
                    └── SolidLiquidLayer
```

## Constructor

```python
SolidLiquidLayer(
    name:                     str,
    layer_index:              int,
    radius_inner:             float,
    radius_outer:             float,
    mass:                     float,
    material_name:            str     = "",
    is_tidal:                 bool    = True,
    tidal_scale:              float   = 1.0,
    shear_modulus_static:     float   = 0.0,
    bulk_modulus_static:      float   = 0.0,
    shear_viscosity_static:   float   = nan,
    bulk_viscosity_static:    float   = nan,
    love_number_k:            complex = 0+0j,
    love_number_h:            complex = 0+0j,
    love_number_l:            complex = 0+0j,
    thermal_conductivity_ref: float   = 4.0,
    thermal_expansion_ref:    float   = 3.0e-5,
    heat_capacity_ref:        float   = 1200.0,
    activation_energy:        float   = 300.0e3,
    activation_volume:        float   = 5.0e-6,
    solidus_temperature:      float   = 1600.0,
    liquidus_temperature:     float   = 2000.0,
    melt_fraction_exponent:   float   = 1.0,
    reference_density:        float   = 3500.0,
    reference_temperature:    float   = 1600.0,
    melt_viscosity_reduction: float   = 25.0,
    tidal_scale_method:       str     = "user_provided",
)
```

### Parameters

| Parameter | Units | Description |
|-----------|-------|-------------|
| `name` | — | Human-readable layer name. |
| `layer_index` | — | Zero-based index; innermost layer = 0. |
| `radius_inner` | m | Inner boundary radius. |
| `radius_outer` | m | Outer boundary radius. |
| `mass` | kg | Total layer mass. Overwritten by each successful world EOS solve. |
| `material_name` | — | Material identifier. Default `""`. |
| `is_tidal` | — | Tidal dissipation flag. Default `True`. |
| `tidal_scale` | — | Dimensionless tidal heating scale. Default `1.0`. |
| `shear_modulus_static` | Pa | Unrelaxed shear modulus. Default `0.0`. |
| `bulk_modulus_static` | Pa | Unrelaxed bulk modulus. Default `0.0`. |
| `shear_viscosity_static` | Pa·s | Reference shear viscosity (at `reference_temperature`, P=0). Default `nan`, meaning unset: attach a viscosity model instead, or a viscous rheology returns NaN. |
| `bulk_viscosity_static` | Pa·s | Reference bulk viscosity. Default `nan`, as above. |
| `love_number_k`, `love_number_h`, `love_number_l` | — | Per-layer complex Love numbers, if you want to carry them on the layer. Default `0+0j`. |
| `thermal_conductivity_ref` | W/(m·K) | Reference thermal conductivity. Default `4.0`. |
| `thermal_expansion_ref` | 1/K | Reference thermal expansion coefficient α. Default `3e-5`. |
| `heat_capacity_ref` | J/(kg·K) | Reference specific heat capacity c_p. Default `1200.0`. |
| `activation_energy` | J/mol | Arrhenius activation energy E_a. Default `300e3`. |
| `activation_volume` | m³/mol | Arrhenius activation volume V_a. Default `5e-6`. |
| `solidus_temperature` | K | Solidus temperature T_s. Default `1600.0`. |
| `liquidus_temperature` | K | Liquidus temperature T_l. Default `2000.0`. |
| `melt_fraction_exponent` | — | Exponent n in the melt-fraction formula. Default `1.0` (linear). |
| `reference_density` | kg/m³ | Reference density ρ_ref for thermal diffusivity. Default `3500.0`. |
| `reference_temperature` | K | Reference temperature T_ref for Arrhenius viscosity. Default `1600.0`. |
| `melt_viscosity_reduction` | — | Coefficient C in exp(−C·φ) melt-viscosity reduction. Default `25.0`. |

## Properties

### Inherited from `PhysicsLayer`

See [PhysicsLayer](physics_layer.md): `shear_modulus_static`, `bulk_modulus_static`, `shear_viscosity_static`, `bulk_viscosity_static`, `love_numbers`, `love_number_k`, `love_number_h`, `love_number_l`, `shear_rheology_set`, `bulk_rheology_set`.

### Inherited from `BaseLayer`

See [BaseLayer](base_layer.md): `name`, `layer_index`, `radius`, `radius_inner`, `radius_outer`, `thickness`, `mass`, `volume`, `density_bulk`, `surface_area_inner`, `surface_area_outer`, `material_name`, `is_tidal`, `tidal_scale`, `eos_data_populated`.

### Thermal / Melt

_Read-only properties._

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

## Methods

### `calc_melt_fraction(temperature, pressure=0.0)` → float

Volumetric melt fraction $\phi \in [0, 1]$:

$$\phi = \left[\operatorname{clamp}\!\left(\frac{T - T_{s}}{T_{l} - T_{s}}, 0, 1\right)\right]^{n}$$

`pressure` is accepted for interface uniformity and is unused: the melt curve carries no pressure dependence.

```python
phi = layer.calc_melt_fraction(3200.0)  # T = 3200 K, P = 0
```

### `calc_viscosity(temperature, pressure=0.0)` → float

Effective dynamic viscosity [Pa·s]:

$$\eta = \eta_\mathrm{ref}\,\exp\!\left[\operatorname{clamp}\!\left(\frac{E_{a} + P V_{a}}{R T} - \frac{E_{a}}{R T_\mathrm{ref}}, -100, 100\right)\right]\exp\!\left[\operatorname{clamp}(-C\phi, -100, 0)\right]$$

The Arrhenius exponent is clamped to $[-100, 100]$ to prevent overflow. Returns `η_ref` when $T = 0$ K.

```python
eta = layer.calc_viscosity(3000.0, 1e11)  # T = 3000 K, P = 100 GPa
```

### `calc_shear_modulus(temperature, pressure=0.0)` → float

Effective shear modulus [Pa] accounting for partial melt:

$$G_\mathrm{eff} = G_\mathrm{static}\,(1 - \phi)$$

```python
G = layer.calc_shear_modulus(3000.0)
```

### `calc_thermal_conductivity(temperature)` → float

Returns the reference thermal conductivity k [W/(m·K)]. Temperature dependence is not modeled.

### `calc_thermal_diffusivity(temperature)` → float

Thermal diffusivity [m²/s], $\kappa = k/(\rho_\mathrm{ref}\,c_{p})$.

### `calc_adiabatic_temperature_gradient(temperature, pressure=0.0)` → float

Adiabatic temperature gradient [K/m], $\alpha T g / c_{p}$.

Gravity g is read from the EOS profile at the outer radius. Returns `0.0` when EOS data has not been populated via `update_eos_data`.

```python
layer.update_eos_data(radii, densities, gravities, pressures)
grad = layer.calc_adiabatic_temperature_gradient(3000.0)
```

### `calc_heat_flux_conductive(temperature_base, temperature_top)` → float

Conductive heat flux [W/m²]:

$$F = \frac{k\,(T_\mathrm{base} - T_\mathrm{top})}{h}$$

where $h$ is the layer thickness. Returns `0.0` for zero-thickness layers.

```python
flux = layer.calc_heat_flux_conductive(temperature_base=3500.0, temperature_top=1500.0)
```


### `set_cooling(cooling)` / `set_radiogenics(radiogenics)`

Attach a cooling (`CoolingBase`) or radiogenics (`RadiogenicsBase`) sub-model. Ownership of the underlying C++ model is transferred into the layer; the passed Python wrapper becomes an empty, non-owning shell and must not be reused (raises `ValueError` if re-attached). Shear/bulk rheology are attached via the inherited `set_shear_rheology` / `set_bulk_rheology` (see [PhysicsLayer](physics_layer.md)).

```python
from TidalPy.cooling_x import make_cooling
from TidalPy.radiogenics_x import IsotopeRadiogenics
layer.set_cooling(make_cooling("convection"))
layer.set_radiogenics(IsotopeRadiogenics.from_dataset("modern_day_chondritic"))
```


### `calc_radiogenic_heating(time, mass)` → float

Radiogenic heating power [W] from the attached sub-model. Returns `0.0` when no radiogenics sub-model has been attached.

### Inherited from PhysicsLayer / BaseLayer

`calc_complex_shear_modulus`, `calc_complex_bulk_modulus`, `update_eos_data`, `get_density`, `get_gravity`, `get_pressure`, `calc_surface_area`, `calc_volume_sphere`, `calc_volume_shell`, `calc_surface_gravity`, `calc_mean_density`, `calc_escape_velocity`, `save_binary`, `load_binary`, `save_config`, `get_config_dict`.

## `get_config_dict()` → dict

Returns all configuration values as a Python dictionary (MKS). Includes all `BaseLayer` + `PhysicsLayer` keys (with `class = "solidliquid"` and the attached model sub-tables) plus the 11 SolidLiquidLayer parameters, and the `cooling` and `radiogenics` sub-tables when those models are attached.

| Key | Units | Description |
|-----|-------|-------------|
| `thermal_conductivity_ref` | W/(m·K) | Reference thermal conductivity. |
| `thermal_expansion_ref` | 1/K | Reference thermal expansion coefficient. |
| `heat_capacity_ref` | J/(kg·K) | Reference specific heat capacity. |
| `activation_energy` | J/mol | Arrhenius activation energy. |
| `activation_volume` | m³/mol | Arrhenius activation volume. |
| `solidus_temperature` | K | Solidus temperature. |
| `liquidus_temperature` | K | Liquidus temperature. |
| `melt_fraction_exponent` | — | Melt-fraction exponent. |
| `reference_density` | kg/m³ | Reference density. |
| `reference_temperature` | K | Reference temperature. |
| `melt_viscosity_reduction` | — | Melt-viscosity reduction coefficient. |

## Binary Serialization

`save_binary` / `load_binary` serialize fields in this order:

1. All `BaseLayer` fields (name, geometry, mass).
2. All `PhysicsLayer` mechanical fields (G, K, η, love number re+im).
3. The 11 SolidLiquidLayer doubles in constructor order.
4. An optional sub-model section: presence flags + recursive binary records for the material EOS model, shear rheology, bulk rheology, shear viscosity, bulk viscosity, partial melt, cooling, and radiogenics models (in that order).

On load, every attached sub-model is reconstructed recursively via each module's binary-dispatch factory, so a saved layer round-trips with all of its physics intact (verify with `eos_set`, `shear_rheology_set`, `cooling_set`, `radiogenics_set`, `calc_complex_shear_modulus`, and `calc_radiogenic_heating`). See [Binary serialization](../../utilities_x/binary_x.md) for the encoding.

Binary class id 102 (`BinaryClassID::SolidLiquidLayer`).

The EOS profile data is not serialized; re-run the world's `solve_eos` after loading.

## Example

```python
import math
from TidalPy.structures_x.layers import SolidLiquidLayer

mantle = SolidLiquidLayer(
    name                     = "mantle",
    layer_index              = 1,
    radius_inner             = 3.485e6,
    radius_outer             = 6.371e6,
    mass                     = 4.043e24,
    material_name            = "perovskite",
    shear_modulus_static     = 1.67e11,
    bulk_modulus_static      = 3.57e11,
    shear_viscosity_static   = 1.0e21,
    thermal_conductivity_ref = 4.5,
    heat_capacity_ref        = 1200.0,
    thermal_expansion_ref    = 2.0e-5,
    activation_energy        = 300.0e3,
    activation_volume        = 5.0e-6,
    solidus_temperature      = 3000.0,
    liquidus_temperature     = 4000.0,
    reference_density        = 4000.0,
    reference_temperature    = 3000.0,
    melt_viscosity_reduction = 25.0,
)

T = 3200.0    # K — slightly above solidus
P = 1e11      # Pa — ~100 GPa

phi = mantle.calc_melt_fraction(T, P)
eta = mantle.calc_viscosity(T, P)
G   = mantle.calc_shear_modulus(T, P)
k   = mantle.calc_thermal_conductivity(T)
kap = mantle.calc_thermal_diffusivity(T)
F   = mantle.calc_heat_flux_conductive(temperature_base=3500.0, temperature_top=1500.0)

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
