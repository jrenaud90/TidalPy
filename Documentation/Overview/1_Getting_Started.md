# Getting Started with TidalPy

_Updated: 2026-09-30_

## Installation

```{include} Readme_raw.md
:start-after: How to Install
:end-before: Using TidalPy
```

## First Calculation

TidalPy builds worlds out of layers, each with an equation of state and physics models, from TOML files. Several worlds are bundled with the package (`TidalPy.Structures.available_worlds()` lists them). The example below builds Io, solves its interior and Love numbers, and places it in orbit about Jupiter to find its tidal heating and orbital rates.

```python
import numpy as np

from TidalPy.Structures import build_world
from TidalPy.Structures.system import System

# Build a bundled world: its layers, materials, and physics models come from its TOML file
io = build_world("io")
io.solve_eos()                                        # Solve the interior structure (density, gravity, pressure)

# Love numbers at Io's orbital period
orbital_frequency = 2.0 * np.pi / (1.769 * 86400.0)   # [rad s-1]
io.solve_love_numbers(
    frequency=orbital_frequency)                      # Radial solve with the world's rheology
print(io.love_number_k)                               # Complex k2, about 0.036 - 0.015j

# Link Io to Jupiter; orbital state lives on the system, not on the world
jupiter = build_world("jupiter_simple")
system = System("jovian")
system.add_world(jupiter)
system.add_world(
    io,
    tidal_host=jupiter,
    semi_major_axis=4.217e8,
    eccentricity=0.0041)
io.set_spin_frequency(system.calc_orbital_frequency(io))   # Keep Io synchronous

rates = system.calc_world_evolution(io)               # Tidal solve plus orbital and spin rates
print(rates["tidal_heating"])                         # [W], about 9e13
print(rates["da_dt"], rates["de_dt"])                 # [m s-1], [s-1]
```

TidalPy computes rates only; to evolve a system in time, integrate these rates with an integrator of your choice (the demos use [CyRK](https://github.com/jrenaud90/CyRK)).

## Logging to a File

TidalPy's messages (warnings from a radial solve, for example) go to the console by default. To also write them to a file of your choosing for the rest of the session, call `init_logger` after importing TidalPy:

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

Warnings and errors are written to the file immediately. Lower levels are buffered until `flush_logger` runs or the interpreter exits (so calling `flush_logger` is not required unless you want to look at a log while an interpreter is still running, e.g., while using a Jupyter notebook). The file is appended to, not overwritten. `TidalPy.reinit()` returns the logger to the settings in your configuration file. To write a timestamped log file to the TidalPy data directory in every session, set `write_log_to_disk = true` in the `[logging]` section of the configuration file instead (see [Configurations](2_TidalPy_Configurations.md)). The levels and the other logging functions are described on the [Logging](../Utilities/logging.md) page.

## Learning More

The module pages in this documentation are the reference for each part of the package. The notebooks in the `Demos` [folder](https://github.com/jrenaud90/TidalPy/tree/main/Demos) walk through the package in order, from configuration and world building to Love numbers, 3D tidal heating, and coupled thermal-orbital evolution. The `Benchmarks` [folder](https://github.com/jrenaud90/TidalPy/tree/main/Benchmarks) compares TidalPy against published results. Users coming from TidalPy 0.7.X should read the [migration guide](../future_structure.md).

If you find an issue, have a question, or want to share an idea for a new feature, please open a GitHub issue [here](https://github.com/jrenaud90/TidalPy/issues). TidalPy also has a slack channel for developers and users. Contact us at [TidalPy@gmail.com](mailto:TidalPy@gmail.com) if you would like to be invited.

## Package Structure

TidalPy is divided into several modules, some of which rely on each other.

- `TidalPy.Structures`: layers, worlds (layered, gas giant, and star), and the `System` that links them; the TOML world builder and the bundled worlds. See [Structures](../Structures/index.md).
- `TidalPy.Material`: material equations of state and the whole-planet EOS solver. See [Material and EOS](../Material/index.md).
- `TidalPy.Rheology`: viscoelastic rheology models that return a complex modulus. See [Rheology](../Rheology/index.md).
- `TidalPy.Viscosity`: temperature- and pressure-dependent viscosity models. See [Viscosity](../Viscosity/index.md).
- `TidalPy.PartialMelt`: partial-melt models that weaken the viscosity and shear modulus. See [Partial Melting](../PartialMelt/index.md).
- `TidalPy.Cooling`: conductive and convective cooling models. See [Cooling](../Cooling/index.md).
- `TidalPy.Radiogenics`: radiogenic heating models and isotope datasets. See [Radiogenics](../Radiogenics/index.md).
- `TidalPy.RadialSolver`: the radial solver for tidal, loading, and free Love numbers and the radial functions. See [RadialSolver](../RadialSolver/index.md).
- `TidalPy.Tides`: eccentricity and obliquity functions, tide models, global (1D) tidal dissipation, and 3D tidal stress, strain, and heating. See [Tides](../Tides/index.md).
- `TidalPy.Dynamics`: spin and orbital rates. See [Dynamics](../Dynamics/index.md).
- `TidalPy.Stellar`: stellar luminosity models. See [Stellar](../Stellar/index.md).
- `TidalPy.Utilities`: logging, binary serialization, base classes, graphics, and numerical helpers. See [Utilities](../Utilities/index.md).
- `TidalPy.WorldPack`: not a module: the bundled world and system TOML files used by `build_world` and `build_system`. See [WorldPack](../Structures/config/worldpack.md).
- `TidalPy.constants`, `TidalPy.configurations`, `TidalPy.paths`: physical constants, the configuration system, and the data directory. See [Constants](../Utilities/constants.md) and [Configurations](2_TidalPy_Configurations.md).
