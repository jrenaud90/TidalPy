# GasLayer

_Updated: 2026-09-19_

`TidalPy.structures_x.layers.GasLayer` (`c_GasLayer` in C++) is the ideal-gas fluid layer class. It inherits `PhysicsLayer` and adds thermodynamic calculations for gas and fluid envelopes such as planetary atmospheres or gaseous mantles. No phase-change, cooling, or radiogenics sub-models are available; use `SolidLiquidLayer` for those.

## Inheritance

```
c_TidalPyBaseClass
  └── c_StructureBase
        └── c_BaseLayer
              └── c_PhysicsLayer
                    └── c_GasLayer
```

## Constructor

```python
from TidalPy.structures_x.layers.gas import GasLayer

layer = GasLayer(
    name                   = "atmosphere",
    layer_index            = 0,
    radius_inner           = 0.0,
    radius_outer           = 7.0e7,
    mass                   = 1.0e27,
    # optional:
    material_name          = "hydrogen",
    is_tidal               = False,
    tidal_scale            = 1.0,
    love_number_k          = 0+0j,
    love_number_h          = 0+0j,
    love_number_l          = 0+0j,
    mean_molecular_weight  = 2.0e-3,   # H₂
    adiabatic_index        = 1.4,
    reference_temperature  = 300.0,
    reference_density      = 1.0,
)
```

### Parameters

All parameters from `PhysicsLayer` are accepted, including the material-state parameters (`temperature`, the shear law, `use_thermal_eos`), except that `is_solid` defaults to `False` (a gas carries no shear stress, so the radial solver treats it as a static liquid), plus:

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


## Binary Serialization

`save_binary(path)` / `load_binary(path, force=False)` round-trip all configuration fields, followed by an optional sub-model section holding the material EOS model and the inherited rheology, viscosity, and partial-melt models (presence flag + recursive binary record each). The EOS profile data is never serialized; re-run the world's `solve_eos` after loading.

Binary class id 103 (`BinaryClassID::GasLayer`).

## Config I/O

```python
layer.save_config("gas_layer.toml")
cfg = layer.get_config_dict()   # dict of all fields (MKS); class = "gas" plus attached-model sub-tables
```

## References

- Wallace, J. M., and Hobbs, P. V. (2006). *Atmospheric Science*, second edition. Ideal gas law and scale height.
- Holton, J. R. (2004). *An Introduction to Dynamic Meteorology*, fifth edition. Adiabatic lapse rate.
