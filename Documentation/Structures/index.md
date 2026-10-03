# Structures (`Structures`)

_Updated: 2026-10-01_

`TidalPy.Structures` contains the world, layer, and system classes. Worlds own layers and run the whole-planet equation-of-state, thermal, and Love-number solves; each layer is a shell of one material with its own physics switches; a `System` links worlds together for insolation and orbital and spin evolution.

| Page | Covers |
|---|---|
| [Worlds](worlds/worlds.md) | The world classes, the equation-of-state solve with its solid and liquid zones, the thermal solve and heat sources, the Love-number solves, global tidal dissipation, and binary serialization. |
| [Layer](layers/layer.md) | The layer class: geometry, the material, the physics switches and flags, the state, cooling and radiogenics, rheology overrides, the profile getters, and serialization. |
| [System](system/system.md) | Linking worlds, orbital elements, insolation, and orbital and spin evolution. |
| [TOML Schema](config/toml_schema.md) | The world and system configuration files and the world builder. |
| [WorldPack](config/worldpack.md) | The bundled worlds and how they are installed and resolved by name. |
| [Schema Examples](config/schema_examples.md) | Six buildable files that use every key the schema accepts, each key commented. |

```{toctree}
:maxdepth: 1

Worlds <worlds/worlds.md>
Layer <layers/layer.md>
System <system/system.md>
TOML Schema <config/toml_schema.md>
WorldPack <config/worldpack.md>
Schema Examples <config/schema_examples.md>
```
