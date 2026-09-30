# GasLayer

_Updated: 2026-09-29_

`TidalPy.Structures.layers.GasLayer` (`c_GasLayer` in C++) is the layer class for gas and fluid envelopes such as planetary atmospheres or gaseous mantles. It inherits `BaseLayer`, so its density, moduli, and viscosities come from its material exactly as for any layer, and it adds four ideal-gas parameters (mean molecular weight, adiabatic index, and a reference temperature and density). Those are stored and serialized, but no calculation reads them: they change no result. No phase-change, cooling, or radiogenics sub-models are available; use `SolidLiquidLayer` for those.

## Inheritance

```
c_TidalPyBaseClass
  └── c_StructureBase
        └── c_BaseLayer
              └── c_GasLayer
```

## Constructor

```python
from TidalPy.Structures.layers.gas import GasLayer

layer = GasLayer(
    name                   = "atmosphere",
    layer_index            = 0,
    radius_inner           = 0.0,
    radius_outer           = 7.0e7,
    mass                   = 1.0e27,
    # optional:
    material_name          = "hydrogen",
    is_tidal               = False,
    tidal_scale            = None,
    mean_molecular_weight  = 2.0e-3,   # H₂
    adiabatic_index        = 1.4,
    reference_temperature  = 300.0,
    reference_density      = 1.0,
)
```

### Parameters

All parameters from `BaseLayer` are accepted, including the material-state parameters (`temperature`, the shear law, `use_thermal_eos`, `use_heating`), except that `is_solid` defaults to `False` (a gas carries no shear stress, so the radial solver treats it as a static liquid), plus:

| Parameter | Unit | Default | Description |
|---|---|---|---|
| `mean_molecular_weight` | kg/mol | `2e-3` | Mean molar mass of the gas |
| `adiabatic_index` | - | `1.4` | γ = c_p/c_v (ratio of specific heats) |
| `reference_temperature` | K | `300.0` | Reference temperature |
| `reference_density` | kg/m³ | `1.0` | Reference density (not the layer's density, which is its material's) |

## Properties

Inherits all `BaseLayer` properties, plus:

| Property | Unit | Description |
|---|---|---|
| `mean_molecular_weight` | kg/mol | Molar mass of the gas |
| `adiabatic_index` | - | γ (ratio of specific heats) |
| `reference_temperature` | K | Reference temperature |
| `reference_density` | kg/m³ | Reference density |


## Config I/O

```python
layer.save_config("gas_layer.toml")
cfg = layer.get_config_dict()   # dict of all fields (MKS); class = "gas" plus attached-model sub-tables
```

The config dict leaves out `reference_density`: a layer file cannot carry it (see the [TOML schema](../config/toml_schema.md)), since the layer's density is its material's.

## References

- Wallace, J. M., and Hobbs, P. V. (2006). *Atmospheric Science*, second edition. Ideal gas law and scale height, for the gas description these parameters are kept for.
- Holton, J. R. (2004). *An Introduction to Dynamic Meteorology*, fifth edition. Adiabatic lapse rate.
