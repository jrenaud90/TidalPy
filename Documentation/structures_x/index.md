# Structures (`structures_x`)

_Updated: 2026-09-16_

`TidalPy.structures_x` contains the world, layer, and system classes. Worlds own layers and run the whole-planet equation-of-state and Love-number solves; a `System` links worlds together for insolation and orbital and spin evolution.

| Page | Covers |
|---|---|
| [Worlds](worlds/worlds.md) | The world classes, the equation-of-state and Love-number solves, global tidal dissipation, and binary serialization. |
| [Base Layer](layers/base_layer.md) | Layer geometry, the material EOS model, and the EOS profile getters. |
| [Physics Layer](layers/physics_layer.md) | Static moduli and viscosities, the rheology, viscosity, and partial-melt models, and the complex moduli. |
| [Solid/Liquid Layer](layers/solidliquid_layer.md) | Thermal and melt parameters, and the cooling and radiogenics models. |
| [Gas Layer](layers/gas_layer.md) | Ideal-gas thermodynamics for gas envelopes. |
| [System](system/system.md) | Linking worlds, orbital elements, insolation, and orbital and spin evolution. |
| [TOML Schema](config/toml_schema.md) | The world and system configuration files and the world builder. |
| [WorldPack](config/worldpack.md) | The bundled worlds and how they are installed and resolved by name. |

```{toctree}
:maxdepth: 1

Worlds <worlds/worlds.md>
Base Layer <layers/base_layer.md>
Physics Layer <layers/physics_layer.md>
Solid/Liquid Layer <layers/solidliquid_layer.md>
Gas Layer <layers/gas_layer.md>
System <system/system.md>
TOML Schema <config/toml_schema.md>
WorldPack <config/worldpack.md>
```
