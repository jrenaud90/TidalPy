# WorldPack_x & World TOML Files (`structures_x.configs.worldpack`)

A **world file** is a TOML description of a single world (a star, gas giant, or
terrestrial/layered body) for the `structures_x` class system. **WorldPack_x** is the
mechanism that ships a small set of example world files with TidalPy, installs them
into a user-editable data directory, and resolves them by name when you call
`build_world("<name>")`.

For a description of the toml file schema used for worlds and layers please see
[`toml_schema.md`](toml_schema.md).

---

## How WorldPack_x works

The example worlds live in the package directory `TidalPy/WorldPack_x/`. They are not
loaded directly from the package; instead they are copied into a per-user, per-version
data directory the first time they are needed, and that data-directory copy is
preferred thereafter. This lets a user edit the installed TOML to change the world
they get from `build_world`, without touching the installed package.

### Locations

| Role | Path |
|------|------|
| Packaged (read-only source) | `<TidalPy package>/WorldPack_x/*` |
| User data directory (editable) | `<user documents>/TidalPy/<major>.<minor>.X/Worlds_x/*` |

The data directory is given by `TidalPy.paths.get_worlds_x_dir()`. It is scoped to
the package's **major.minor** version (with a literal `X` patch placeholder, e.g.
`0.8.X`), so every patch release of a given major.minor shares the same directory
(configs and downloaded data are not duplicated on each bugfix release).

World TOMLs and their **companion data files** (PREM-like profiles: `.csv`, `.txt`,
`.dat`) are both installed. A world's `data_file` reference is resolved by
`resolve_data_file` in this order: the world TOML's own directory, the data
directory, the packaged `WorldPack_x`, then the working directory.

### Install (copy-if-absent)

`install_worldpack_x(force=False)` copies each packaged `*.toml` into the data
directory, but only when a file of that name is **not already present** there. This
means:

- A user's edits (or renamed/added files) in the data directory are never
  overwritten.
- A world newly added to the package shows up on the next run (it is absent in the
  data dir, so it is copied).
- Passing `force=True` re-copies every packaged world, discarding local edits.

It returns the data directory path. It runs automatically inside `resolve_world_path`
and `available_worlds`, so you rarely need to call it directly.

### Name resolution (data directory first)

`build_world("<name>")` (equivalently `BaseWorld.build("<name>")`) resolves a bare
name through `resolve_world_path(name)`:

1. Run `install_worldpack_x()` (copy-if-absent).
2. If `Worlds_x/<name>.toml` exists, use it (the user-editable copy wins).
3. Otherwise fall back to the packaged `WorldPack_x/<name>.toml`.
4. Otherwise raise `FileNotFoundError`.

A `source` that ends in `.toml` or names an existing file is treated as a direct
path; a `dict` is used as-is. So the same `build_world` entry point accepts a bundled
name, a file path, or an in-memory config.

`available_worlds()` returns the sorted union of the data-directory names and the
packaged names (so newly installed and packaged worlds both appear).

```python
from TidalPy.structures_x import build_world, available_worlds, install_worldpack_x

install_worldpack_x()              # optional; copies packaged worlds into Worlds_x/
print(available_worlds())          # ['earth_simple', 'jupiter_simple', 'sol', ...]

earth = build_world("earth_simple")   # data-dir copy preferred, else packaged
earth.solve_eos()
```

### Versioning caveat

Because installs are copy-if-absent and the data directory is version-scoped, a
*within-version* schema change to a bundled world does **not** propagate to a user
who already has the old copy in `Worlds_x/`. Across versions the new version's
`Worlds_x` folder starts empty, so fresh copies install. During development, delete
`Worlds_x/*.toml` (or call `install_worldpack_x(force=True)`) to pick up edits.

---

## Adding a bundled world

1. Write a schema-`0.2.0` world TOML and drop it in `TidalPy/WorldPack_x/<name>.toml`.
   `MANIFEST.in` already globs `TidalPy/WorldPack_x/*.toml`, so it is packaged on the
   next `uv pip install`.
2. It installs to `Worlds_x/` and becomes available as `build_world("<name>")` and in
   `available_worlds()` on the next run.

The bundled worlds favor the per-material defaults: keep them small by specifying
`class` + `type` + geometry and letting `TidalPy_Configs_x.toml` supply the EOS and
physics models. Override anything inline as shown above and [`here`](toml_schema.md).

---

## Relationship to the legacy WorldPack

`TidalPy/WorldPack/` (zipped, extracted into `Worlds/`) serves the legacy
(non-`_x`) world builder and uses the older world-config schema. `WorldPack_x/`
(loose TOMLs, copied into `Worlds_x/`) serves the new `structures_x` builder and the
`0.2.0` schema. They are independent; do not mix their files.

## API summary (`TidalPy.structures_x.configs.worldpack`)

| Function | Description |
|----------|-------------|
| `install_worldpack_x(force=False) -> str` | Copy packaged worlds and their data files into the data dir (copy-if-absent; `force` re-copies). Returns the data dir. |
| `resolve_world_path(name) -> str` | Resolve a bundled name to a TOML path (data dir preferred, then packaged). Raises `FileNotFoundError` if unknown. |
| `resolve_data_file(data_file, base_dir=None) -> str` | Resolve a world's companion `data_file` (toml dir -> data dir -> packaged -> cwd). Raises `FileNotFoundError` if unknown. |
| `available_worlds() -> list[str]` | Sorted union of data-dir and packaged world names. |
| `get_worlds_x_dir() -> str` | The user data directory for `_x` worlds (`.../TidalPy/<major>.<minor>.X/Worlds_x`). |
| `PACKAGED_WORLDPACK_DIR` | Path to the packaged `WorldPack_x` directory. |

`install_worldpack_x` and `available_worlds` are also re-exported from
`TidalPy.structures_x`; `build_world` (and `BaseWorld.build`) consume the resolver
transparently.
