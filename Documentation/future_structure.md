# The Future TidalPy Structure

TidalPy's internals are being rewritten in C++. The new implementation lives in modules that carry a `_x` suffix (`structures_x`, `Tides_x`, `RadialSolver_x`, `rheology_x`, and so on) and ships side by side with the classic modules today. In a future major release the `_x` modules will become the only TidalPy: the classic modules will be removed and the new ones will drop their suffix. Nothing about the classic API changes until then, but any new development on TidalPy will happen in the `_x` modules, other than bug fixes. We highly encourage new projects to start using the `_x` modules or make plans to switch.

The 0.8.X series is the last to include the classic modules. It will continue to receive bug fixes, but no new features, until the end of 2026, and support for 0.8.X after 2026 is not guaranteed. Plan to finish porting before then.

This page explains what is different, maps the classic modules to their replacements, and shows how to port common workflows.

TidalPy announces this transition once per session when the package is imported. The notice can be silenced with:

```python
import warnings
from TidalPy.exceptions import TidalPyDeprecationWarning
warnings.filterwarnings("ignore", category=TidalPyDeprecationWarning)
```

## Why the rework?

* **Performance.** All core physics now runs in C++ (with the Eigen linear algebra library and CyRK integrators), wrapped by thin Cython layers. There is no numba JIT warmup, and hot paths avoid Python entirely.
* **Predictability.** The classic system stored planet state on Python objects and propagated changes through cascading updates, which was hard to reason about and easy to break. The new classes store configuration and return results from `calc_*` methods without mutating states.
* **Consistency.** Every physics module (rheology, cooling, radiogenics, viscosity, partial melting, equations of state, tides) follows the pattern: a C++ class hierarchy, a name-based factory (`make_<module>`), direct callable functions, vectorized variants, TOML configuration, and binary files which can be saved and loaded from disk for fast and accurate reproducibility.

## How much faster is it?

Performance tests were run with the classic and new backend. Both are timed after warm up as the best of seven batches, and each figure is the lowest of three independent runs in fresh processes, with console logging limited to errors so terminal output is not timed. The machine is an 8-core AMD desktop running Windows 11, Python 3.13, numpy 2.4, numba 0.67, scipy 1.18, and BurnMan 2.1. Ratios move with the machine and the problem size, so read them as rough magnitudes, and measure your own workload before relying on any of them.

The new backend is dramatically faster where the classic path called out to BurnMan or paid a numba compile, two to three times faster on array work and global tidal heating, about even on the radial solver, and **slower** on a few paths, which are listed too.

### Where it is faster

| Task | Classic | New | Change |
|---|---|---|---|
| Build a planet with its interior (Io, 3 layers) | 178 ms | 3.4 ms | **52x faster** |
| Radiogenic heating, one evaluation | 0.29 us | 0.031 us | **9.5x faster** |
| Radiogenic heating, 10k times | 0.087 ms | 0.019 ms | **4.5x faster** |
| Rheology, 10k complex moduli | 0.175 ms | 0.055 ms | **3.2x faster** |
| Orbit-averaged 3D heating map (50 x 16 x 32) | 13.3 ms | 4.2 ms | **3.1x faster** |
| Global tidal heating, e^2 truncation | 0.020 ms | 0.0075 ms | **2.7x faster** |
| Global tidal heating, e^4 truncation | 0.022 ms | 0.0094 ms | **2.3x faster** |
| Global tidal heating, degrees 2 to 4, e^10 | 0.047 ms | 0.021 ms | **2.2x faster** |
| Global tidal heating, e^10 truncation | 0.025 ms | 0.015 ms | **1.7x faster** |
| Build a world from config (2 layers, no interior solve) | 1.38 ms | 0.61 ms | **2.2x faster** |
| Homogeneous Love numbers (closed form) | 0.15 us | 0.069 us | **2.2x faster** |

Building a planet is the one that changes how the package feels to use. The classic path handed the interior to BurnMan, which does mineral-physics lookups and its own root finding; the new path integrates the equation of state in C++. A fresh Io went from 1.4 seconds to 4 milliseconds.

The global tidal heating rows use the homogeneous Love method, which solves the same problem as the classic `quick_tidal_dissipation`, at eccentricity truncations both backends tabulate. The new backend accepts e^1 through e^5, e^10, e^15, and e^20 and promotes any other requested level to the next tabulated one, so a request for e^6 or e^8 runs at e^10. Its cost follows the number of distinct forcing frequencies rather than the number of modes: Love numbers are solved once per frequency and degree, and the layer-averaged shear modulus the homogeneous methods need is formed once per frequency and shared by every degree, so adding degrees adds little. With the `radial_solver` Love method each frequency and degree is a full radial solve instead, and that solve dominates.

### Where it is about even

The standalone radial solver was already Cython calling CyRK, so there was little left to win, and the rewrite bought 15 to 35 percent on realistic problems while losing about 15 percent on a tiny one where call overhead dominates.

| Task | Classic | New | Change |
|---|---|---|---|
| `radial_solver`, 1 layer, 10 slices | 0.32 ms | 0.37 ms | 0.85x, slower |
| `radial_solver`, 1 layer, 200 slices | 0.47 ms | 0.41 ms | 1.15x faster |
| `radial_solver`, 3 layers, 300 slices | 0.83 ms | 0.71 ms | 1.17x faster |
| `radial_solver`, propagation matrix, 200 slices | 0.109 ms | 0.081 ms | 1.35x faster |
| Convective cooling, one evaluation | 0.16 us | 0.16 us | 1.0x, even |

### Where it is slower

| Task | Classic | New | Change |
|---|---|---|---|
| Instantaneous 3D heating map (50 x 16 x 32 x 8 times) | 13.1 ms | 28.4 ms | **0.46x, 2.2x slower** |
| Convective cooling, 10k evaluations | 0.158 ms | 0.202 ms | 0.78x, 1.3x slower |
| Rheology, one complex modulus | 0.056 us | 0.076 us | 0.74x, 1.4x slower |

One of these has a known cause. A single scalar rheology call is dominated by the Python-to-C++ boundary rather than by the arithmetic, and the numba path crosses a cheaper one; use the vectorized calls, where the new backend wins by 3x, whenever there is more than a handful of values. The vectorized convective cooling gap has not been investigated.

The instantaneous 3D map is the largest gap and does not have a tidy explanation; both paths evaluate the same grid, and the new one takes a little over twice as long. The orbit-averaged case reverses it, because the new backend has a secular form that never builds a time axis while the classic route has to average one.

### The first call

Steady-state timings hide something users feel immediately. The classic backend compiles its numba kernels the first time they run and caches the machine code on disk, so the first session after installing or upgrading pays the full compile and every later session still pays to load and dispatch the cached kernels. The new backend has nothing to compile. Each figure below is the median first call in a fresh process:

| First call | Classic, first session after installing | Classic, later sessions | New |
|---|---|---|---|
| Tidal heating, degrees 2 to 4, e^10 | 6.9 s | 1.2 s | 0.12 ms |
| 3D heating map (50 x 16 x 32 x 8 times) | 4.7 s | 1.0 s | 29 ms |
| Build a planet with its interior (Io) | not measured | 1.4 s | 4.0 ms |
| Dual-body dissipation rates | not measured | 1.1 s | no single-call equivalent |
| Build a world from config | not measured | 0.22 s | 1.2 ms |

A script that computes one 3D map and exits spends about a second in the classic backend once its cache is warm, nearly five seconds the first time after installing, and under thirty milliseconds in the new one, whatever the steady-state ratio says.

## Module map

| Classic module | Replacement | Documentation |
|----------------|-------------|---------------|
| `TidalPy.structures` (worlds, layers, burnman builds) | `TidalPy.structures_x` | [Worlds](structures_x/worlds/worlds.md), [Layers](structures_x/layers/base_layer.md), [System](structures_x/system/system.md), [TOML schema](structures_x/config/toml_schema.md) |
| `TidalPy.RadialSolver` | `TidalPy.RadialSolver_x` | [Love number solvers](RadialSolver_x/dense_radial_solution.md) |
| `TidalPy.rheology` | `TidalPy.rheology_x` | [Rheology](rheology_x/index.md) |
| `TidalPy.tides` | `TidalPy.Tides_x` (mostly via world methods) | [Global tides](Tides_x/global_tides.md), [3D heating](Tides_x/multilayer_3d_heating.md), [Love numbers](Tides_x/love/love_numbers.md) |
| `TidalPy.cooling` | `TidalPy.cooling_x` | [Cooling](cooling_x/index.md) |
| `TidalPy.radiogenics` | `TidalPy.radiogenics_x` | [Radiogenics](radiogenics_x/index.md) |
| viscosity functions in `TidalPy.rheology` | `TidalPy.viscosity_x` | [Viscosity](viscosity_x/index.md) |
| partial melting in `TidalPy.rheology` | `TidalPy.partial_melt_x` | [Partial melting](partial_melt_x/index.md) |
| burnman equations of state | `TidalPy.Material_x` | [Material and EOS](material_x/index.md) |
| `TidalPy.stellar` | `TidalPy.stellar_x` (world-attached luminosity) | [Stellar](stellar_x/index.md) |
| orbit/spin evolution in `TidalPy.orbit` | `TidalPy.dynamics_x` + the `System` class | [Dynamics](dynamics_x/index.md), [System](structures_x/system/system.md) |
| `TidalPy.utilities.graphics` (`yplot`, `planet_plot`) | `TidalPy.Utilities_x.graphics_x` (`plot_ys`, `plot_interior`) | [Graphics](utilities_x/graphics_x.md) |
| assorted helpers | `TidalPy.Utilities_x` | [Utilities](utilities_x/index.md) |

## Porting examples

### Building a world and getting Love numbers

Classic (burnman-backed OOP world):

```python
from TidalPy.structures import build_world

world = build_world("earth")
# State-based updates followed: world.orbit = ..., world.update_orbit_spin(), etc.
```

New system (TOML-configured world; solves are explicit method calls):

```python
from TidalPy.structures_x import build_world

world = build_world("earth_simple")
world.solve_eos()

result = world.solve_love_numbers(frequency=1.0e-5)
print(result["success"], world.love_number_k)
```

The world builder reads bundled or user TOML files, a file path, or a Python dictionary. See the [TOML schema](structures_x/config/toml_schema.md) and the `Demos (_x)` notebooks for more details.

### The standalone radial solver

The array-based `radial_solver` call is nearly identical between the two systems:

```python
# Classic
from TidalPy.RadialSolver import radial_solver

# New
from TidalPy.RadialSolver_x import radial_solver
```

The new solver adds dense radial evaluation at arbitrary radii (`get_radial_solution`), implicit CyRK integrators (`BDF`, `LSODA`, `Radau`), and a surface conditioning diagnostic (`surface_solve_amplification`). See [dense radial solutions](RadialSolver_x/dense_radial_solution.md).

The input builders keep the classic argument names but take `rheology_x` models (classic `TidalPy.rheology` models raise `TypeError`), and a single model can stand in for the per-layer tuple:

```python
# Classic
from TidalPy.RadialSolver import build_rs_input_homogeneous_layers
from TidalPy.rheology.models import Elastic, Maxwell
build_data = build_rs_input_homogeneous_layers(
    ..., shear_rheology_model_tuple=(Maxwell(), Maxwell(), Maxwell()),
    bulk_rheology_model_tuple=(Elastic(), Elastic(), Elastic()), ...)

# New (one model applies to every layer; per-layer tuples still work)
from TidalPy.RadialSolver_x import build_rs_input_homogeneous_layers
from TidalPy.rheology_x import Elastic, Maxwell
build_data = build_rs_input_homogeneous_layers(
    ..., shear_rheology_model_tuple=Maxwell(), bulk_rheology_model_tuple=Elastic(), ...)
```

See [input builders](RadialSolver_x/build_inputs.md).

### Rheology models

```python
# Classic
from TidalPy.rheology.models import Maxwell
rheology = Maxwell()
complex_modulus = rheology(frequency, shear_modulus, viscosity)

# New (note the argument order: modulus, viscosity, frequency)
from TidalPy.rheology_x import Maxwell, make_rheology
rheology = Maxwell()
complex_modulus = rheology.calc_complex_modulus(shear_modulus, viscosity, frequency)

# Or build by name (case-insensitive, alias-aware)
rheology = make_rheology("maxwell")
```

Every model also provides vectorized `calc_*` variants that accept numpy arrays and broadcast.

## Learning the new system

The best introduction is the `Demos (_x)` notebooks in the navigation: `Basics` notebooks (configuration, world building, save/load), `Physics` notebooks (orbits, tides, rheology, Love numbers, 3D heating, thermal/EOS), and `Systems` notebooks (multi-world systems, coupled thermal-orbital evolution). The `Benchmarks (_x)` pages validate the new system against published results and track its performance.
