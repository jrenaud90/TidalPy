# WorldPack and World TOML Files (`Structures.configs.worldpack`)

_Updated: 2026-10-10_

WorldPack is TidalPy's set of example world files (stars, gas giants, terrestrial bodies, and systems), installed into a user-editable data directory and found by name. The [TOML schema](toml_schema.md) describes their keys.

```python
from TidalPy.Structures import (
    build_world, build_system, available_worlds, available_systems, install_worldpack)

install_worldpack()                   # optional; runs on import
print(available_worlds())             # 26 names (the table below)
print(available_systems())            # ['pluto_charon_system', 'sol_system']

earth = build_world("earth_simple")   # the user's copy, else the packaged file
earth.solve_eos()
system = build_system("sol_system")
```

A file with a `[worlds.<name>]` table is a system, any other a world (`config_kind()`). Giving a config to the wrong builder raises a `ValueError` naming the other one.

## Bundled Worlds

| Name | Type | Contents |
|---|---|---|
| `sol` | star | The Sun, from its effective temperature. |
| `trappist1` | star | The M8V host of the seven-planet system, with the measured luminosity stated rather than derived. |
| `jupiter_simple` | gas giant | One uniform gas layer: the mass, but not C/MR2 or k2. |
| `jupiter` | gas giant | Heavy-element core, metallic hydrogen, and molecular hydrogen, all fluid, fitted to the mass and C/MR2. |
| `neptune` | gas giant | Rock core, hot dense ices, and a hydrogen and helium envelope, all fluid, fitted to the mass and C/MR2. |
| `earth_simple` | terrestrial | Three-layer Earth with a static liquid outer core, reproducing the mass, C/MR2, and k2. |
| `earth_thermal` | terrestrial | `earth_simple` with a pressure- and temperature-dependent mantle viscosity (dry-olivine creep) for thermal evolution. It solves its temperature profile (pass `solve_eos` the surface temperature); its M2 Q is then 315, near the solid Earth's 280. |
| `earth_prem` | terrestrial | The PREM seismic profile, read from the companion `PREM.csv`. |
| `earth_prem_q` | terrestrial | `earth_prem` with PREM's own quality factors setting the loss (`q_provided`, `seismic_q`). |
| `io` | terrestrial | The Segatz et al. (1988) asthenosphere end-member, with its viscosity fitted to Io's measured heat output. |
| `europa` | terrestrial | Iron core, silicate mantle, and a solid ice Ih shell (no ocean). |
| `luna` | terrestrial | Five-layer Moon reproducing C/MR2, k2, and Q at the month and the year. |
| `mercury` | terrestrial | Fluid outer core, with the mantle rigidity fitted to the measured k2. |
| `pluto` | terrestrial | Rock core under one water layer that is an ocean below the observed 100 km ice-shell thickness. |
| `charon` | terrestrial | Rock core under one water layer, frozen through; Pluto's partner in a mutually synchronous binary. |
| `triton` | terrestrial | The same recipe with a 58 km ocean, on a retrograde synchronous orbit about Neptune. |
| `luna_dynamic` | terrestrial | `luna` with a dynamic, compressible Fe-S liquid outer core (Birch-Murnaghan), refitted to the same observables. |
| `mercury_dynamic` | terrestrial | `mercury` with a dynamic, compressible Fe-Si liquid outer core (Birch-Murnaghan). |
| `pluto_dynamic` | terrestrial | `pluto` with its ocean solved as a dynamic, compressible liquid. |
| `europa_dynamic` | terrestrial | `europa` with a 115 km dynamic, compressible ocean under a 25 km ice shell. |
| `trappist1b` to `trappist1h` | terrestrial | Liquid iron cores under silicate mantles, built from the Agol et al. (2021) masses and radii. |
| `sol_system` | system | The Sun with Earth and Jupiter. |
| `pluto_charon_system` | system | Pluto and Charon as a mutually synchronous pair (Brozovic and Jacobson 2024 orbit), lit by the Sun at Pluto's mean heliocentric orbit. |

Each file's comments say which numbers were fitted to which observable (Io's heat output, the k2 of the three liquid-core bodies, Luna's Q at the month and the year), and the test suite checks those claims. The layered bodies share these choices:

- Each layer's material is a table of fitted values ([The Material Table](toml_schema.md#the-material-table)). Each layer states a temperature and each solid layer a viscosity at it, so every solid layer dissipates. `tidal_scale` values come from the 3D heating integral.
- The liquid cores of Luna, Mercury, Earth-Simple, and the TRAPPIST-1 planets are static and of constant density. The `_dynamic` worlds solve their liquid dynamically on a pressure-dependent law, which keeps it neutrally stratified. A constant-density compressible liquid solved dynamically is unstably stratified, and its solve takes steps in proportion to the forcing period ([Dense Radial Solutions](../../RadialSolver/dense_radial_solution.md#dynamic-liquid-layers-at-long-forcing-periods)). The `_dynamic` worlds' $k_2$ at the default tolerances agrees with `rtol = 1e-12` solves to 2e-9 or better from a one-day period down to forcing frequencies of 1e-16 rad/s, where it is the static liquid's, and to 7e-7 or better with the temperature solved.
- The hydrospheres of Pluto, Charon, and Triton are one layer of MatPack's `ice_ih` with `state = "auto"`, `use_melting`, `use_pressure_melting`, and `use_melt_density` on. The layer is ice where it is colder than the pressure-dependent melting point and water where it is warmer, so its ocean grows or freezes with the layer's temperature. Pluto's temperature, 268.197 K, is fitted to a 100 km shell.
- Layers carry thermal constants, cooling models, and radiogenics (with `use_heating` on). The warm silicate mantles can melt through a `peridotite` liquid phase, the Monteux et al. (2016) melting curves, Henning weakening, and `use_melting` and `use_pressure_melting` on.
- Temperatures are prescribed. Every layered file but `earth_thermal` (which pins `true`) and the PREM worlds pins `solve_temperature = false` in `[eos_solver]`, so the fits hold at the stated temperatures. A call may still pass `solve_temperature=True`, and `System.evolve` does for a thermal run.
- The laws stay physical over the temperatures a thermal solve reaches. Rock cores and crusts use silicate creep referenced to 1e21 Pa s at 1600 K with $E$ = 300 kJ/mol, so cold rock is elastic, and ice uses the diffusion creep of `ice_ih` (about 1e14 Pa s at its melting point). Andrade creep keeps cold ice dissipating down to about 150 K. The mantles of `earth_simple` and the TRAPPIST-1 planets are the exception. Their viscosity has no pressure dependence, so when solved thermally it softens toward the hot base of the mantle, to 1e13 Pa s in `earth_simple` (which `earth_thermal` corrects) and to about 1e15 Pa s above a TRAPPIST-1 planet's 4000 K core, which lowers its Q 4 to 11 times.
- Radiogenic heating is evaluated at the solve's `time`, the time since formation on the isotope model's clock. Left out, it is today's (4600 Myr for `modern_day_chondritic`); a run that starts at `time = 0` starts at formation, with about nine times today's chondritic heating.
- The cooling models are simple, which shows in long thermal runs. A thermal Pluto has a conducting lid about 1 km thick and loses heat about 700 times faster than its radiogenic heating replaces it, so over 1 Gyr of `System.evolve` its hydrosphere cools from about 268 K to about 151 K and the ocean freezes. Its core, which has no melting law, warms from 600 K to about 1700 K. Charon behaves similarly.
- The TRAPPIST-1 cores are static liquids, though at 4000 K iron would likely be solid at the central pressures of most of them (up to 405 GPa in TRAPPIST-1b). A solid Maxwell iron core (1e20 Pa s) also solves, at every forcing frequency down to 1e-16 rad/s within 6e-9 of an `rtol = 1e-12` solve; it lowers $k_2$ at tidal periods (0.372 against 0.599 for TRAPPIST-1b at 1e-5 rad/s) and relaxes like a fluid below about 1e-12 rad/s, where its $k_2$ approaches the liquid core's.

Each interior is fitted to the stated mass, which `solve_eos()` returns. A solar-system body fits its least-known layer density, except Pluto and Charon, which fix a shared core density and fit the core radius, as a TRAPPIST-1 planet does with both of its densities fixed. A measured moment of inertia is then an independent check (Pluto, Charon, Triton, and the TRAPPIST-1 planets have none, so their files record the model value). The gas giants fit two densities to the mass and C/MR2 together and are fluid throughout (a rigid rock core of shear modulus 1e11 Pa lowers $k_2$ by a relative 7e-5 in Jupiter and 8e-4 in Neptune). Their k2 is then a prediction: Jupiter 0.534 against the measured 0.590, Neptune 0.427 against a published 0.41, and `jupiter_simple` the fluid-sphere 1.5.

## Adding a Bundled World

1. Save a schema-`0.2.0` world TOML as `TidalPy/WorldPack/<name>.toml` (`MANIFEST.in` already includes it).
2. On the next run it installs to `Worlds/` and appears in `available_worlds()` (a system file, in `available_systems()`).
3. Worlds of wide interest are welcome as pull requests on TidalPy's GitHub.

## How it Works

On import, TidalPy copies the packaged worlds and their data files into a per-user, per-version data directory with `install_worldpack()`, copying a file only when no file of that name is there. Edits and added files are never overwritten, and new packaged worlds appear on the next run. `install_worldpack(force=True)` re-copies everything, discarding edits. Files edited, added, or deleted during a session are seen by the next build; `TidalPy.reinit()` re-reads every file.

### Locations

| Role | Path |
|------|------|
| Packaged (read-only source) | `<TidalPy package>/WorldPack/*` |
| User data directory (editable) | `<user documents>/TidalPy/<major>.<minor>.X/Worlds/*` |

Every patch release of one major.minor version (e.g. `0.8.X`) shares the data directory. Bundled materials install the same way into `.../Materials/` ([MatPack](../../Material/matpack.md)). A world's `data_file` (a `.csv`, `.txt`, or `.dat` profile, format in [Building a World From a Radial Profile](toml_schema.md#building-a-world-from-a-radial-profile)) is looked for in the world file's directory, the data directory, the packaged `WorldPack`, then the working directory.

### Name Resolution

`build_world("<name>")` (or `BaseWorld.build`) uses the data directory's `Worlds/<name>.toml` if it exists (the user's copy wins), else the packaged `WorldPack/<name>.toml`, else raises `FileNotFoundError`. A `source` ending in `.toml` or naming an existing file is a path, and a `dict` is used as is.

### Stale Copies

A bundled file updated within a version does not replace the copy in `Worlds/` (a new major.minor version starts empty). TidalPy cannot tell an old copy from an edited one, so it replaces neither, but when a world, system, or data file resolves to a copy that differs from the packaged file (line endings aside), it warns once per file per session, naming both files and the two fixes: delete the copy, or `install_worldpack(force=True)`, which discards every edit. To silence it, set in `TidalPy_Configs.toml`:

```toml
[warnings]
    stale_worldpack_copy = false
```

## API Summary (`TidalPy.Structures.configs.worldpack`)

| Function | Description |
|----------|-------------|
| `install_worldpack(force=False) -> str` | Copy packaged files into the data dir (copy-if-absent; `force` re-copies). Returns the data dir. |
| `resolve_world_path(name) -> str` | A bundled name's TOML path (data dir, then packaged), or `FileNotFoundError`. |
| `resolve_data_file(data_file, base_dir=None) -> str` | A world's `data_file` path (toml dir, data dir, packaged, cwd), or `FileNotFoundError`. |
| `available_worlds() -> list[str]` | Sorted union of data-dir and packaged world names, systems excluded. |
| `available_systems() -> list[str]` | The same for the bundled system names. |
| `config_kind(source) -> str` | `"system"` if the config has a `worlds` table, else `"world"`. Takes a path or a dict. |
| `warn_if_stale_copy(data_path) -> bool` | Whether a data-directory file differs from the packaged one, with the [stale-copy](#stale-copies) warning. |
| `get_worlds_dir() -> str` | The user data directory for worlds (also `TidalPy.paths.get_worlds_dir`). |
| `PACKAGED_WORLDPACK_DIR` | Path to the packaged `WorldPack` directory. |

`install_worldpack`, `available_worlds`, and `available_systems` are also in `TidalPy.Structures`.
