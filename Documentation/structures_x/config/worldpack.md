# WorldPack_x and World TOML Files (`structures_x.configs.worldpack`)

_Updated: 2026-09-16_

A world file is a TOML description of a single world (a star, gas giant, or terrestrial/layered body) for the `structures_x` class system. WorldPack_x ships a small set of example world files with TidalPy, installs them into a user-editable data directory, and resolves them by name when you call `build_world("<name>")`.

The same directory also holds system files, which describe several worlds and their orbits instead of one body. The two kinds are told apart by content: a system names its members in a `[worlds.<name>]` table and a world never does, which is the test `config_kind()` applies. `available_worlds()` and `available_systems()` therefore list disjoint sets, and handing a config to the wrong builder raises a `ValueError` naming the other one.

The TOML schema for worlds and layers is described in [`toml_schema.md`](toml_schema.md).

## How it Works

The example worlds live in the package directory `TidalPy/WorldPack_x/`. They are not loaded directly from the package; instead they are copied into a per-user, per-version data directory the first time they are needed, and that data-directory copy is preferred thereafter. This lets a user edit the installed TOML to change the world they get from `build_world`, without touching the installed package.

### Locations

| Role | Path |
|------|------|
| Packaged (read-only source) | `<TidalPy package>/WorldPack_x/*` |
| User data directory (editable) | `<user documents>/TidalPy/<major>.<minor>.X/Worlds_x/*` |

The data directory is given by `TidalPy.paths.get_worlds_x_dir()`. It is scoped to the package's major.minor version (with a literal `X` patch placeholder, e.g. `0.8.X`), so every patch release of a given major.minor shares the same directory (configs and downloaded data are not duplicated on each bugfix release).

World TOMLs and their companion data files (PREM-like profiles: `.csv`, `.txt`, `.dat`) are both installed. A world's `data_file` reference is resolved by `resolve_data_file` in this order: the world TOML's own directory, the data directory, the packaged `WorldPack_x`, then the working directory.

### Install

`install_worldpack_x(force=False)` copies each packaged `*.toml` into the data directory, but only when a file of that name is not already present there. This means:

- A user's edits (or renamed/added files) in the data directory are never overwritten.
- A world newly added to the package shows up on the next run (it is absent in the data dir, so it is copied).
- Passing `force=True` re-copies every packaged world, discarding local edits.

It returns the data directory path. It runs automatically inside `resolve_world_path` and `available_worlds`, so you rarely need to call it directly.

### Name Resolution

`build_world("<name>")` (equivalently `BaseWorld.build("<name>")`) resolves a bare name through `resolve_world_path(name)`:

1. Run `install_worldpack_x()` (copy-if-absent).
2. If `Worlds_x/<name>.toml` exists, use it (the user-editable copy wins).
3. Otherwise fall back to the packaged `WorldPack_x/<name>.toml`.
4. Otherwise raise `FileNotFoundError`.

A `source` that ends in `.toml` or names an existing file is treated as a direct path; a `dict` is used as-is. So the same `build_world` entry point accepts a bundled name, a file path, or an in-memory config.

`available_worlds()` returns the sorted union of the data-directory names and the packaged names (so newly installed and packaged worlds both appear), restricted to single worlds. `available_systems()` does the same for the system files.

```python
from TidalPy.structures_x import (
    build_world, build_system, available_worlds, available_systems, install_worldpack_x)

install_worldpack_x()              # optional; copies packaged worlds into Worlds_x/
print(available_worlds())          # ['earth_prem', 'earth_simple', 'jupiter_simple', 'sol']
print(available_systems())         # ['sol_system']

earth = build_world("earth_simple")   # data-dir copy preferred, else packaged
earth.solve_eos()

system = build_system("sol_system")   # the same resolution, for a multi-world config
```

### Versioning Caveat

Because installs are copy-if-absent and the data directory is version-scoped, a within-version schema change to a bundled world does not propagate to a user who already has the old copy in `Worlds_x/`. Across versions the new version's `Worlds_x` folder starts empty, so fresh copies install. During development, delete `Worlds_x/*.toml` (or call `install_worldpack_x(force=True)`) to pick up edits.

## Adding a Bundled World

1. Write a schema-`0.2.0` world TOML and drop it in `TidalPy/WorldPack_x/<name>.toml`. `MANIFEST.in` already globs `TidalPy/WorldPack_x/*.toml`, so it is packaged on the next `uv pip install`.
2. It installs to `Worlds_x/` and becomes available as `build_world("<name>")` and in `available_worlds()` on the next run. A system file follows the same two steps and appears in `available_systems()` instead.
3. A world of wide interest can be added to the WorldPack through a pull request on TidalPy's GitHub.

The bundled worlds favor the per-material defaults: keep them small by specifying `class`, `type`, and geometry and letting `TidalPy_Configs_x.toml` supply the EOS and physics models. Override anything inline as shown in the [TOML schema](toml_schema.md).

## API Summary (`TidalPy.structures_x.configs.worldpack`)

| Function | Description |
|----------|-------------|
| `install_worldpack_x(force=False) -> str` | Copy packaged worlds and their data files into the data dir (copy-if-absent; `force` re-copies). Returns the data dir. |
| `resolve_world_path(name) -> str` | Resolve a bundled name to a TOML path (data dir preferred, then packaged). Raises `FileNotFoundError` if unknown. |
| `resolve_data_file(data_file, base_dir=None) -> str` | Resolve a world's companion `data_file` (toml dir -> data dir -> packaged -> cwd). Raises `FileNotFoundError` if unknown. |
| `available_worlds() -> list[str]` | Sorted union of data-dir and packaged world names, systems excluded. |
| `available_systems() -> list[str]` | The same for the bundled system names. |
| `config_kind(source) -> str` | `"system"` if the config has a `worlds` table, else `"world"`. Takes a path or a dict. |
| `get_worlds_x_dir() -> str` | The user data directory for `_x` worlds (`.../TidalPy/<major>.<minor>.X/Worlds_x`). |
| `PACKAGED_WORLDPACK_DIR` | Path to the packaged `WorldPack_x` directory. |

`install_worldpack_x`, `available_worlds`, and `available_systems` are also re-exported from `TidalPy.structures_x`; `build_world` and `build_system` (and `BaseWorld.build` / `System.build`) consume the resolver transparently.
