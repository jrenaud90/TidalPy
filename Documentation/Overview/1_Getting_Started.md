# Getting Started with TidalPy

_Updated: 2026-10-02_

## Installation

```{include} Readme_raw.md
:start-after: How to Install
:end-before: Using TidalPy
```

## First Calculation

TidalPy builds worlds from TOML files as stacks of layers, each made of a material (an equation of state, moduli, viscosities, and melting laws) and carrying its own physics models. Several worlds are bundled (`TidalPy.Structures.available_worlds()` lists them). The example builds Io, solves its interior, places it in orbit about Jupiter, and finds its Love numbers, tidal heating, and orbital rates.

```python
from TidalPy.Structures import System, build_world

# Build a bundled world: its layers, materials, and physics models come from its TOML file
io = build_world("io")
io.solve_eos()  # Solve the interior structure (density, gravity, pressure)

# Link Io to Jupiter; orbital state lives on the system, not on the world
jupiter = build_world("jupiter_simple")
system = System("jovian")
system.add_world(jupiter)
system.add_world(
    io,
    tidal_host=jupiter,
    semi_major_axis=4.217e8,
    eccentricity=0.0041,
    synchronous=True)  # Spin Io at its mean motion

# Love numbers at Io's mean motion, its tidal forcing frequency since it rotates synchronously
orbital_frequency = system.calc_orbital_frequency(io)  # [rad s-1], from Kepler's third law
io.solve_love_numbers(
    frequency=orbital_frequency)  # Radial solve with the world's rheology
print(io.love_number_k)  # Complex k2, about 0.036 - 0.015j

rates = system.calc_world_evolution(io)  # Tidal solve plus orbital and spin rates
print(rates["tidal_heating"])  # [W], about 9e13
print(rates["da_dt"], rates["de_dt"])  # [m s-1], [s-1]
```

TidalPy computes rates only; to evolve a system in time, integrate these rates with an integrator of your choice (the demos use [CyRK](https://github.com/jrenaud90/CyRK)).

## Logging to a File

TidalPy's messages (warnings from a radial solve, for example) go to the console. To also write them to a file of your choosing for the rest of the session, call `init_logger` after importing TidalPy:

```python
from pathlib import Path

from TidalPy.Utilities.logging import init_logger, flush_logger

log_path = Path("my_runs/io_run.log").resolve()   # Any location; missing folders are created

init_logger({                                     # Replaces the current console and file outputs
    "console_level": "warning",
    "file_level": "debug",
    "log_to_file": True,
    "log_file_path": str(log_path),
})

# ... run your calculations ...

flush_logger()                                    # Write buffered info and debug lines before reading the file
```

Warnings and errors reach the file immediately. Lower levels are buffered until `flush_logger` runs or the interpreter exits, so call it only to read a log while the interpreter is still running (in a Jupyter notebook, say). The file is appended to, not overwritten. `TidalPy.reinit()` returns the logger to your configuration file's settings. To write a timestamped log file in every session, set `write_log_to_disk = true` in the `[logging]` section of the configuration file instead (see [Configurations](2_TidalPy_Configurations.md)). The [Logging](../Utilities/logging.md) page covers the levels and the other functions.

## Learning More

The module pages are the reference for each part of the package. The notebooks in the `Demos` [folder](https://github.com/jrenaud90/TidalPy/tree/main/Demos) walk through it in order, from configuration and world building to Love numbers, 3D tidal heating, and coupled thermal-orbital evolution. The `Benchmarks` [folder](https://github.com/jrenaud90/TidalPy/tree/main/Benchmarks) compares TidalPy against published results. Users coming from TidalPy 0.7.X should read the [migration guide](../future_structure.md).

For an issue, a question, or a feature idea, please open a [GitHub issue](https://github.com/jrenaud90/TidalPy/issues). TidalPy also has a Slack channel for developers and users; email [TidalPy@gmail.com](mailto:TidalPy@gmail.com) for an invitation.

## Package Structure

TidalPy is divided into these modules:

- `TidalPy.Structures`: layers, worlds (layered, gas giant, and star), and the `System` that links them; the TOML world builder and the bundled worlds. See [Structures](../Structures/index.md).
- `TidalPy.Material`: equation-of-state and shear-modulus laws, the phases and materials built from them, MatPack (the bundled materials), and the whole-planet EOS solver. See [Materials](../Material/index.md).
- `TidalPy.Rheology`: viscoelastic rheology models that return a complex modulus. See [Rheology](../Rheology/index.md).
- `TidalPy.Viscosity`: temperature- and pressure-dependent viscosity models. See [Viscosity](../Viscosity/index.md).
- `TidalPy.PartialMelt`: melting curves, melt weakening, and bulk-mixing laws that a material uses once it begins to melt. See [Partial Melting](../PartialMelt/index.md).
- `TidalPy.Cooling`: conductive and convective cooling models. See [Cooling](../Cooling/index.md).
- `TidalPy.Radiogenics`: radiogenic heating models and isotope datasets. See [Radiogenics](../Radiogenics/index.md).
- `TidalPy.RadialSolver`: the radial solver for tidal, loading, and free Love numbers and the radial functions. See [RadialSolver](../RadialSolver/index.md).
- `TidalPy.Tides`: eccentricity and obliquity functions, tide models, global (1D) tidal dissipation, and 3D tidal stress, strain, and heating. See [Tides](../Tides/index.md).
- `TidalPy.Dynamics`: spin and orbital rates. See [Dynamics](../Dynamics/index.md).
- `TidalPy.Stellar`: stellar luminosity models. See [Stellar](../Stellar/index.md).
- `TidalPy.Utilities`: logging, binary serialization, base classes, graphics, and numerical helpers. See [Utilities](../Utilities/index.md).
- `TidalPy.WorldPack` and `TidalPy.MatPack`: not modules: the bundled world, system, and material TOML files used by `build_world`, `build_system`, and `load_material`. See [WorldPack](../Structures/config/worldpack.md) and [MatPack](../Material/matpack.md).
- `TidalPy.constants`, `TidalPy.configurations`, `TidalPy.paths`: physical constants, the configuration system, and the data directory. See [Constants](../Utilities/constants.md) and [Configurations](2_TidalPy_Configurations.md).
