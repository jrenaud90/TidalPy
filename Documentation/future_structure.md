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

Performance tests were run with the classic and new backend. Both are timed after warm up as the best of seven batches. The machine is an 8-core AMD desktop running Windows 11, Python 3.13, numpy 2.4, numba 0.67, scipy 1.18, and BurnMan 2.1. Ratios move with the machine and the problem size, so read them as rough magnitudes, and measure your own workload before relying on any of them.

The short version: the new backend is dramatically faster where the classic path called out to BurnMan or paid a numba compile, modestly faster on array work, about even on the radial solver, and **slower** on a few paths where the classic implementation collapses work that the new one does mode by mode. Those are listed too.

### Where it is faster

| Task | Classic | New | Change |
|---|---|---|---|
| Build a planet with its interior (Io, 3 layers) | 178 ms | 3.4 ms | **52x faster** |
| Radiogenic heating, one evaluation | 0.31 us | 0.033 us | **9.5x faster** |
| Radiogenic heating, 10k times | 0.085 ms | 0.019 ms | **4.4x faster** |
| Orbit-averaged 3D heating map (50 x 16 x 32) | 13.3 ms | 4.2 ms | **3.2x faster** |
| Rheology, 10k complex moduli | 0.175 ms | 0.054 ms | **3.2x faster** |
| Global tidal heating, e^2 truncation | 0.020 ms | 0.008 ms | **2.5x faster** |
| Build a world from config (2 layers, no interior solve) | 1.42 ms | 0.62 ms | **2.3x faster** |
| Homogeneous Love numbers (closed form) | 0.22 us | 0.10 us | **2.2x faster** |

Building a planet is the one that changes how the package feels to use. The classic path handed the interior to BurnMan, which does mineral-physics lookups and its own root finding; the new path integrates the equation of state in C++. A fresh Io went from 1.4 seconds to 4 milliseconds.

### Where it is about even

The standalone radial solver was already Cython calling CyRK, so there was little left to win, and the rewrite bought 10 to 35 percent on realistic problems while losing about 15 percent on a tiny one where call overhead dominates.

| Task | Classic | New | Change |
|---|---|---|---|
| `radial_solver`, 1 layer, 10 slices | 0.32 ms | 0.38 ms | 0.84x, slower |
| `radial_solver`, 1 layer, 200 slices | 0.47 ms | 0.41 ms | 1.14x faster |
| `radial_solver`, 3 layers, 300 slices | 0.85 ms | 0.72 ms | 1.17x faster |
| `radial_solver`, propagation matrix, 200 slices | 0.110 ms | 0.081 ms | 1.36x faster |
| Convective cooling, one evaluation | 0.20 us | 0.20 us | 1.0x, even |

### Where it is slower

| Task | Classic | New | Change |
|---|---|---|---|
| Instantaneous 3D heating map (50 x 16 x 32 x 8 times) | 13.6 ms | 28.5 ms | **0.48x, 2.1x slower** |
| Global tidal heating, degrees 2 to 4, e^6 | 0.038 ms | 0.122 ms | **0.32x, 3.2x slower** |
| Global tidal heating, e^6 truncation and above | 0.023 ms | 0.042 ms | **0.55x, 1.8x slower** |
| Convective cooling, 10k evaluations | 0.158 ms | 0.200 ms | 0.79x, 1.3x slower |
| Rheology, one complex modulus | 0.10 us | 0.14 us | 0.72x, 1.4x slower |

Three of these have a known cause. A single scalar rheology or cooling call is dominated by the Python-to-C++ boundary rather than by the arithmetic, and the numba path crosses a cheaper one; use the vectorized calls, where the new backend wins by 3x, whenever there is more than a handful of values. Global tidal heating scales with the mode count in the new backend, which evaluates 4, 11, and 35 modes at e^2, e^4, and e^6, while the classic implementation first collapses the modes onto their unique forcing frequencies and evaluates the rheology once per frequency, which is why its cost barely moves with the truncation level. The same accounts for the degree sweep. Collapsing by unique frequency is a real optimization the new backend has not adopted yet.

The instantaneous 3D map is the largest gap and does not have a tidy explanation; both paths evaluate the same grid, and the new one takes about twice as long. The orbit-averaged case reverses it, because the new backend has a secular form that never builds a time axis while the classic route has to average one.

### The first call

Steady-state timings hide something users feel immediately. The classic backend compiles its numba kernels on first use, once per session, and that cost is not small:

| First call, classic | First call, new |
|---|---|
| 3D heating map: 4.5 s | 29 ms |
| Build a planet with its interior: 1.4 s | 4.0 ms |
| Dual-body dissipation rates: 1.1 s | n/a |
| Tidal heating, degrees 2 to 4: 0.7 s | 0.7 ms |
| Build a world from config: 0.24 s | 1.2 ms |
| Convective cooling, 10k: 18 ms | 0.35 ms |

A script that computes one 3D map and exits spends four and a half seconds in the classic backend and under thirty milliseconds in the new one, whatever the steady-state ratio says.

### Two results that are not the same result

Timings only mean something when both sides compute the same thing, and in two places they do not.

Classic `quick_tidal_dissipation` returns 55 percent of the correct heating for a homogeneous degree-2 body. The closed form `(21/2) (-Im k2) G M^2 R^5 e^2 n / a^6` gives 2.649e10 W for the Io-like case used above; the new backend returns 2.649e10 W and the classic returns 1.464e10 W. The cause is in `TidalPy.tides.love1d.effective_rigidity_general`, whose coefficient is written `(2 l^2 + 4 l + 3 / l)` where it should be `(2 l^2 + 4 l + 3) / l`, giving 17.5 instead of 19/2 at degree 2. This is classic-only and is fixed by construction in the new backend.

The classic 3D grid output is not finite everywhere. On the grid timed above it carries NaN across the whole `r = 0` slice, which is the expected coordinate singularity, and also across about 8 percent of the interior points, along with infinities that make a plain volume integral diverge. The new output is finite everywhere except `r = 0`, and integrating it over the volume reproduces the 1D total to 0.05 percent.

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
