# GasLayer

`TidalPy.structures_x.layers.GasLayer` (`c_GasLayer` in C++) is the ideal-gas
fluid layer class. It inherits `PhysicsLayer` and adds thermodynamic
calculations for gas and fluid envelopes such as planetary atmospheres or
gaseous mantles.  No phase changes, cooling, or radiogenics sub-models are
available, use `SolidLiquidLayer` for those features.

## Inheritance

```
c_TidalPyBaseClass
  └── c_StructureBase
        └── c_BaseLayer
              └── c_PhysicsLayer
                    └── c_GasLayer        ← this class
```

## Constructor

```python
from TidalPy.structures_x.layers.gas import GasLayer

layer = GasLayer(
    name                         = "atmosphere",
    layer_index                  = 0,
    radius_inner               = 0.0,
    radius_outer               = 7.0e7,
    mass                      = 1.0e27,
    # optional:
    material_name                = "hydrogen",
    is_tidal                     = False,
    tidal_scale                  = 1.0,
    shear_modulus_static      = 0.0,
    bulk_modulus_static       = 0.0,
    shear_viscosity_static   = nan,
    bulk_viscosity_static    = nan,
    love_number_k                = 0+0j,
    love_number_h                = 0+0j,
    love_number_l                = 0+0j,
    mean_molecular_weight = 2.0e-3,   # H₂
    adiabatic_index              = 1.4,
    reference_temperature      = 300.0,
    reference_density      = 1.0,
)
```

### Parameters

All parameters from `PhysicsLayer` are accepted plus:

| Parameter | Unit | Default | Description |
|---|---|---|---|
| `mean_molecular_weight` | kg/mol | `2e-3` | Mean molar mass of the gas |
| `adiabatic_index` | — | `1.4` | γ = c_p/c_v (ratio of specific heats) |
| `reference_temperature` | K | `300.0` | Reference temperature |
| `reference_density` | kg/m³ | `1.0` | Reference number density |

## Properties

Inherits all `BaseLayer` and `PhysicsLayer` properties, plus:

| Property | Unit | Description |
|---|---|---|
| `mean_molecular_weight` | kg/mol | Molar mass of the gas |
| `adiabatic_index` | — | γ (ratio of specific heats) |
| `reference_temperature` | K | Reference temperature |
| `reference_density` | kg/m³ | Reference density |

## Methods

### `calc_adiabatic_lapse_rate(gravity)`

Dry adiabatic lapse rate for an ideal gas [K/m]:

$$\Gamma = \frac{g \, (\gamma - 1) \, M}{\gamma \, R}$$

where $g$ is gravitational acceleration, $\gamma$ is the adiabatic index,
$M$ is the mean molar mass, and $R$ is the universal gas constant.

**Returns** `float` [K/m].

### `calc_scale_height(temperature, gravity)`

Barometric (pressure) scale height [m]:

$$H = \frac{R \, T}{g \, M}$$

**Returns** `float` [m].

### `calc_pressure_ideal_gas(temperature, density)`

Ideal gas pressure [Pa]:

$$P = \frac{\rho \, R \, T}{M}$$

**Returns** `float` [Pa].

### `calc_sound_speed(temperature)`

Adiabatic sound speed for an ideal gas [m/s]:

$$c_s = \sqrt{\frac{\gamma \, R \, T}{M}}$$

**Returns** `float` [m/s].

## Binary serialization

`save_binary(path)` / `load_binary(path, force=False)` round-trip all
configuration fields, followed by an optional sub-model section holding the
inherited shear and bulk rheology models (presence flag + recursive binary record
each; attach them with the inherited `set_shear_rheology` / `set_bulk_rheology`).
EOS profile data is never serialized (must be re-populated after loading).

Binary class ID: `103` (`BinaryClassID::GasLayer`).

## Config I/O

```python
layer.save_config("gas_layer.toml")
cfg = layer.get_config_dict()   # dict of all fields (MKS); class = "gas" plus attached-model sub-tables
```

## Literature

- Ideal gas law and scale height: standard atmospheric physics textbooks
  (e.g. Wallace & Hobbs, *Atmospheric Science*, 2006).
- Adiabatic lapse rate: Holton, *An Introduction to Dynamic Meteorology*,
  5th ed., 2004.
