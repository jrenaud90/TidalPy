# Schema Examples

_Updated: 2026-10-07_

Six files under `Documentation/Structures/config/examples/` use, between them, every key of the world and system schema, each with a comment on what it does and what else it accepts. None describes a real body. The values are chosen so that every file builds, the worlds solve their EOS and Love numbers and run their tides, and the files round-trip through `get_config_dict`; the test suite checks this and that no accepted key is missing. The [TOML schema](toml_schema.md) is the reference these files walk through, and the bundled worlds ([WorldPack](worldpack.md)) are their terse, fitted counterparts.

| File | Shows |
|---|---|
| `example_world.toml` | A layered world with every world-level key, `[eos_solver]`, `[radial_solver]`, a `[tides]` table for the rheology model, and five layers that between them use every layer key and geometry specifier, the three forms of a material (a MatPack name, a preset with overrides, and a full table), a material with both phases and its melting curves, weakening, and mixing laws, a radius-tabulated material, rheology overrides, every cooling model, the `off` and `isotope` radiogenics models with an inline isotope table, and a layer that holds its mass. |
| `example_gasgiant.toml` | A gas giant: an envelope of a liquid-only material, the per-degree lists of an analytic tide model (`ctl_q`), the long spellings of the truncation keys, and the `ctl` Love method with its stored lag. |
| `example_star.toml` | A star: the two stellar scalars, the `[luminosity]` model table, and the `cpl` tide model. |
| `example_profile_world.toml` | A world built from a radial data file: `data_file`, the ways a table refines detected layers (by index, by name, and by a radius range), and a file overriding the RK45 default such worlds pin. |
| `example_profile_q_world.toml` | A world whose radial data file's quality factors set the loss: `q_provided`, its two settings, and a radius range of layers overriding them with its own `seismic_q` table. |
| `example_system.toml` | A system: the star flag, tidal hosts, both orbits, and a member from a bundled name, a member from a file path, and a member given inline with the `fixed` radiogenics model. |

Build any of them with `build_world(path)` or `build_system(path)`. The system's file-path member resolves against the system file's folder first, so the system builds from any directory.

## A Layered World

```{literalinclude} examples/example_world.toml
:language: toml
```

## A Gas Giant

```{literalinclude} examples/example_gasgiant.toml
:language: toml
```

## A Star

```{literalinclude} examples/example_star.toml
:language: toml
```

## A World From a Radial Profile

```{literalinclude} examples/example_profile_world.toml
:language: toml
```

## A World From a Radial Profile's Quality Factors

```{literalinclude} examples/example_profile_q_world.toml
:language: toml
```

## A System

```{literalinclude} examples/example_system.toml
:language: toml
```
