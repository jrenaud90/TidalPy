# Structures (`Structures`)

_Updated: 2026-10-07_

`TidalPy.Structures` contains the world, layer, and system classes. A world owns its layers and runs the whole-planet equation-of-state, thermal, Love-number, and tidal solves; each layer is a shell of one material with its own physics switches; a `System` links worlds for insolation and orbital and spin evolution. Start with [Worlds](worlds/worlds.md#quick-start) to build a world and get its Love numbers and tidal heating.

| Page | Covers |
|---|---|
| [Worlds](worlds/worlds.md) | Building a world, the equation-of-state solve and its solid and liquid zones, Love numbers, global tidal dissipation, the thermal solve and heat sources, and serialization. |
| [Layer](layers/layer.md) | The layer class: geometry, material, physics switches and flags, state, cooling and radiogenics, rheology overrides, profile getters, and serialization. |
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
