# SolidLiquidLayer

_Updated: 2026-09-20_

`TidalPy.structures_x.layers.SolidLiquidLayer` extends `PhysicsLayer` with optional sub-models for radiogenic heating and convective or conductive cooling, and with the thermal-transport calculations that need the layer's geometry or solved profile. It adds no parameters of its own.

The physics:

- Thermal transport: a constant conductivity $k$; the diffusivity $\kappa = k/(\rho\,c_{p})$ at the layer's bulk density; the adiabatic gradient $\alpha T g / c_{p}$ (requires EOS data for gravity).
- Conductive heat flux, $F = k\,(T_\mathrm{base} - T_\mathrm{top})/h$.

The thermal constants these read ($k$, $\alpha$, $c_{p}$) are not layer parameters, and neither are viscosity, melt fraction, the static moduli, or the melt-weakened shear modulus. They belong to the layer's material, which is its EOS model (see [Material EOS Models](../../material_x/material_eos.md#the-material)): the whole-planet EOS solve evaluates them as it integrates, and the layer's radius getters read them back, so a layer and the solve it feeds cannot disagree about them. `set_shear_viscosity`, `set_bulk_viscosity`, and `set_partial_melt` are helpers that hand the model to the attached EOS. See [viscosity_models.md](../../viscosity_x/viscosity_models.md) for the Arrhenius law and its `molar_activation_energy_j_mol` and `molar_activation_volume_m3_mol` parameters.

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
    love_number_k:            complex = 0+0j,
    love_number_h:            complex = 0+0j,
    love_number_l:            complex = 0+0j,
    tidal_scale_method:       str     = "user_provided",
    is_solid:                 bool    = True,
    is_static:                bool    = True,
    is_incompressible:        bool    = False,
    temperature:                          float = 0.0,
    use_thermal_eos:                      bool  = False,
    use_heating:                          bool  = False,
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
| `love_number_k`, `love_number_h`, `love_number_l` | — | Per-layer complex Love numbers, if you want to carry them on the layer. Default `0+0j`. |
| `tidal_scale_method` | - | How the layer's share of the world's tidal heating is set. Default `"user_provided"`. |
| `is_solid`, `is_static`, `is_incompressible` | - | Radial-solver assumptions; see [PhysicsLayer](physics_layer.md). Defaults `True`, `True`, `False`. |
| `temperature`, `use_thermal_eos`, `use_heating` | | Layer-state parameters; see [PhysicsLayer](physics_layer.md). |

## Properties

### Inherited from `PhysicsLayer`

See [PhysicsLayer](physics_layer.md): `shear_modulus_static`, `bulk_modulus_static`, `shear_viscosity_static`, `bulk_viscosity_static`, `love_numbers`, `love_number_k`, `love_number_h`, `love_number_l`, `shear_rheology_set`, `bulk_rheology_set`, `is_solid`, `is_static`, `is_incompressible`, `temperature`, and `use_thermal_eos`. The four static constants read the layer's material (its EOS model).

### Inherited from `BaseLayer`

See [BaseLayer](base_layer.md): `name`, `layer_index`, `radius`, `radius_inner`, `radius_outer`, `thickness`, `mass`, `volume`, `density_bulk`, `surface_area_inner`, `surface_area_outer`, `material_name`, `is_tidal`, `tidal_scale`, `eos_data_populated`.

### Thermal

_Read-only properties._

| Property | Units | Description |
|----------|-------|-------------|
| `thermal_conductivity` | W/(m·K) | Thermal conductivity k of the layer's material (its EOS model). NaN when none is attached. |
| `thermal_expansion` | 1/K | Thermal expansivity α of the material, the same α its density law uses. |
| `heat_capacity` | J/(kg·K) | Specific heat capacity c_p of the material. |
| `cooling_set` | — | `True` after a cooling sub-model is attached. |
| `radiogenics_set` | — | `True` after a radiogenics sub-model is attached. |

## Methods

### `calc_thermal_conductivity(temperature)` -> float

Returns the material's thermal conductivity k [W/(m·K)]. Temperature dependence is not modeled.

### `calc_thermal_diffusivity(temperature)` -> float

Thermal diffusivity [m²/s], $\kappa = k/(\rho\,c_{p})$, with $\rho$ the layer's bulk density (its mass over its volume). The material itself answers for any density: `MaterialEOSBase.calc_thermal_diffusivity(density)`, which is what the thermal solve uses with the local density.

### `calc_adiabatic_temperature_gradient(temperature, pressure=0.0)` -> float

Adiabatic temperature gradient [K/m], $\alpha T g / c_{p}$.

Gravity g is read from the EOS profile at the outer radius. Returns `0.0` when EOS data has not been populated via `update_eos_data`.

```python
layer.update_eos_data(radii, densities, gravities, pressures)
grad = layer.calc_adiabatic_temperature_gradient(3000.0)
```

### `calc_heat_flux_conductive(temperature_base, temperature_top)` -> float

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


### `calc_radiogenic_heating(time, mass)` -> float

Radiogenic heating power [W] from the attached sub-model. Returns `0.0` when no radiogenics sub-model has been attached.

### Inherited from PhysicsLayer / BaseLayer

`calc_complex_shear_modulus`, `calc_complex_bulk_modulus`, `update_eos_data`, `get_density`, `get_gravity`, `get_pressure`, `calc_surface_area`, `calc_volume_sphere`, `calc_volume_shell`, `calc_surface_gravity`, `calc_mean_density`, `calc_escape_velocity`, `save_binary`, `load_binary`, `save_config`, `get_config_dict`.

## `get_config_dict()` -> dict

Returns all configuration values as a Python dictionary (MKS): all `BaseLayer` + `PhysicsLayer` keys (with `class = "solidliquid"` and the attached model sub-tables), plus the `cooling` and `radiogenics` sub-tables when those models are attached. The layer adds no scalar keys: its thermal constants are `thermal_conductivity_w_mk`, `thermal_expansion_1_k`, and `heat_capacity_j_kgk` in the `material` table.

## Binary Serialization

`save_binary` / `load_binary` serialize fields in this order:

1. All `BaseLayer` fields (name, geometry, mass).
2. All `PhysicsLayer` fields (Love numbers re+im, the three layer-assumption flags, `temperature`, `use_thermal_eos`, `use_heating`).
3. An optional sub-model section: presence flags + recursive binary records for the material EOS model, shear rheology, bulk rheology, cooling, and radiogenics models (in that order). The EOS record carries the whole material with it: the static and thermal constants, the shear law, and its viscosity and partial-melt models.

On load, every attached sub-model is reconstructed recursively via each module's binary-dispatch factory, so a saved layer round-trips with all of its physics intact (verify with `eos_set`, `shear_rheology_set`, `cooling_set`, `radiogenics_set`, `calc_complex_shear_modulus`, and `calc_radiogenic_heating`). See [Binary serialization](../../utilities_x/binary_x.md) for the encoding.

Binary class id 102 (`BinaryClassID::SolidLiquidLayer`).

The EOS profile data is not serialized; re-run the world's `solve_eos` after loading.

## Example

```python
import math
from TidalPy.Material_x.eos import ConstantDensityEOS
from TidalPy.structures_x.layers import SolidLiquidLayer

mantle = SolidLiquidLayer(
    name                     = "mantle",
    layer_index              = 1,
    radius_inner             = 3.485e6,
    radius_outer             = 6.371e6,
    mass                     = 4.043e24,
    material_name            = "perovskite",
)

# The material: density law, static moduli, static viscosity, thermal constants. Attach a viscosity and a
# partial-melt model to it
# (mantle.set_shear_viscosity(...), mantle.set_partial_melt(...)) and the same query returns their answers.
mantle.set_eos(ConstantDensityEOS(
    reference_density      = 4000.0,
    shear_modulus_static   = 1.67e11,
    bulk_modulus_static    = 3.57e11,
    shear_viscosity_static = 1.0e21,
    thermal_conductivity   = 4.5,
    heat_capacity          = 1200.0,
    thermal_expansion      = 2.0e-5,
))

T = 3200.0    # K — hot lower mantle
P = 1e11      # Pa — ~100 GPa

# Before attaching it, ask the material about a pressure and temperature directly. Once the layer sits in a solved
# world, read mantle.get_shear_modulus(radius) and friends instead: they report what the solve used.
material = ConstantDensityEOS(reference_density=4000.0, shear_modulus_static=1.67e11, shear_viscosity_static=1.0e21)
state = material.calc_material_state(pressure=P, temperature=T)
phi = state["melt_fraction"]
eta = state["shear_viscosity"]
G   = state["shear_modulus"]
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
assert restored.thermal_conductivity == mantle.thermal_conductivity
```
