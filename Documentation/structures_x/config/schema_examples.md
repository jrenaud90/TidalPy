# Schema Examples

_Updated: 2026-09-22_

Five files under `Documentation/structures_x/config/examples/` use, between them, every key the world and system schema accepts, each with a comment saying what it does and what else it takes. None describes a real body: the values are chosen so that every file builds, the worlds solve their equation of state and Love numbers and run their tides, and the files round-trip through `get_config_dict`. `Tests/Test_Structures_x/Test_Config/test_schema_examples_01.py` checks all of that and that no accepted key is missing from the set, so a key added to the schema without an example fails the suite. The bundled worlds ([WorldPack](worldpack.md)) are the terse, fitted counterparts; the [TOML schema](toml_schema.md) is the reference these files walk through.

| File | Shows |
|---|---|
| `example_world.toml` | A layered world with every world-level key, `[eos_solver]`, `[radial_solver]`, a `[tides]` table for the rheology model, and five layers that between them use every class but `base`, every geometry specifier, every scalar flag, every material law (`constant`, `bm`, `vinet`, `interpolate`), every viscosity, melt, rheology, cooling, and radiogenics model with its keys, an inline isotope table, a floating layer, and a thermal-EOS layer. |
| `example_gasgiant.toml` | A gas giant: the gas-layer scalars, the per-degree lists of an analytic tide model (`ctl_q`), the long spellings of the truncation keys, and the `ctl` Love method with its stored lag. |
| `example_star.toml` | A star: the two stellar scalars, the `[luminosity]` model table, and the `cpl` tide model. |
| `example_profile_world.toml` | A world built from a radial data file: `data_file`, the two ways a table refines a detected layer, and a file overriding the RK45 default such worlds pin. |
| `example_system.toml` | A system: the star flag, tidal hosts, both orbits, and a member from a bundled name, a member from a file path, and a member given inline. |

Build any of them with `build_world(path)` or `build_system(path)`; the system's file-path member resolves against the working directory, so build it from the `examples` directory or edit the path.

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

## A System

```{literalinclude} examples/example_system.toml
:language: toml
```
