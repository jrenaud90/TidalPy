# The Future TidalPy Structure

TidalPy's internals are being rewritten in C++. The new implementation lives in modules that carry a `_x` suffix (`structures_x`, `Tides_x`, `RadialSolver_x`, `rheology_x`, and so on) and ships side by side with the classic modules today. In a future major release the `_x` modules become the only TidalPy: the classic modules are removed and the new ones drop their suffix. Nothing about the classic API changes until then, but all new development happens in the `_x` modules, and the classic modules receive bug fixes only. New projects should start with the `_x` modules, and existing ones should plan to switch.

The 0.8.X series is the last to include the classic modules. It will continue to receive bug fixes, but no new features, until the end of 2026, and support for 0.8.X after 2026 is not guaranteed. Plan to finish porting before then.

This page explains what is different, maps the classic modules to their replacements, and shows how to port common workflows. The <a href="code_map.html">interactive code map</a> draws the new backend's main classes and functions, the calls between them, and the purpose, inputs, and outputs of each call.

TidalPy announces this transition once per session when the package is imported. The notice can be silenced with:

```python
import warnings
from TidalPy.exceptions import TidalPyDeprecationWarning
warnings.filterwarnings("ignore", category=TidalPyDeprecationWarning)
```

## Reasons for the Rework

* **Performance.** All core physics runs in C++ (with the Eigen linear algebra library and CyRK integrators), wrapped by thin Cython layers. There is no numba JIT warmup, and hot paths avoid Python entirely.
* **Predictability.** The classic system stored planet state on Python objects and propagated changes through cascading updates. The new classes store configuration and return results from `calc_*` methods without mutating state.
* **Consistency.** Every physics module (rheology, cooling, radiogenics, viscosity, partial melting, equations of state, tides) follows one pattern: a C++ class hierarchy, a name-based factory (`make_<module>`), direct callable functions, vectorized variants, TOML configuration, and binary files that can be saved and loaded from disk for fast and accurate reproducibility.

## Performance

Performance tests were run with the classic and the new backend. Both are timed after warm up as the best of seven batches, and each figure is the lowest of three independent runs in fresh processes, with console logging limited to errors so terminal output is not timed. The machine is an 8-core AMD desktop running Windows 11, Python 3.13, numpy 2.4, numba 0.67, scipy 1.18, and BurnMan 2.1. Ratios move with the machine and the problem size, so read them as rough magnitudes and measure your own workload before relying on any of them.

The new backend is much faster where the classic path called out to BurnMan or paid a numba compile, 4 to 11 times faster on 3D heating maps, two to three times faster on array work, 1.6 to 2.4 times faster on global tidal heating, and slower on the standalone radial solver, which now runs through the world path, and on two scalar calls; all are listed with their causes.

### Where It Is Faster

| Task | Classic | New | Change |
|---|---|---|---|
| Build a planet with its interior (Io, 3 layers) | 189 ms | 1.8 ms | **103x faster** |
| Orbit-averaged 3D heating map (50 x 16 x 32) | 13.3 ms | 1.24 ms | **10.7x faster** |
| Radiogenic heating, one evaluation | 0.29 us | 0.031 us | **9.5x faster** |
| Radiogenic heating, 10k times | 0.088 ms | 0.020 ms | **4.5x faster** |
| Instantaneous 3D heating map (50 x 16 x 32 x 8 times) | 13.3 ms | 3.1 ms | **4.2x faster** |
| Rheology, 10k complex moduli | 0.176 ms | 0.056 ms | **3.2x faster** |
| Homogeneous Love numbers (closed form) | 0.15 us | 0.064 us | **2.4x faster** |
| Global tidal heating, e^2 truncation | 0.020 ms | 0.0084 ms | **2.4x faster** |
| Build a world from config (2 layers, no interior solve) | 1.48 ms | 0.64 ms | **2.3x faster** |
| Global tidal heating, degrees 2 to 4, e^10 | 0.048 ms | 0.023 ms | **2.1x faster** |
| Global tidal heating, e^4 truncation | 0.021 ms | 0.010 ms | **2.0x faster** |
| Global tidal heating, e^10 truncation | 0.025 ms | 0.016 ms | **1.6x faster** |

The planet-building row is the largest change. The classic path handed the interior to BurnMan, which does mineral-physics lookups and its own root finding: the new path integrates the equation of state in C++. A fresh Io went from 1.4 seconds to about 3 milliseconds.

The global tidal heating rows use the homogeneous Love method, which solves the same problem as the classic `quick_tidal_dissipation`, at eccentricity truncations both backends tabulate. At level n both keep the eccentricity terms through $e^n$: the classic tables hold the squared $G^2$, the new ones every term of the unsquared $G$, so the new side sums more modes at the same level. The new backend accepts e^1 through e^5, e^10, e^15, and e^20 and promotes any other requested level to the next tabulated one, so a request for e^6 or e^8 runs at e^10. Its cost follows the number of distinct forcing frequencies rather than the number of modes: Love numbers are solved once per frequency and degree, and the layer-averaged shear modulus the homogeneous methods need is formed once per frequency and shared by every degree, so adding degrees adds little. With the `radial_solver` Love method each frequency and degree is a full radial solve instead, and that solve dominates.

### Where It Is About Even

| Task | Classic | New | Change |
|---|---|---|---|
| Convective cooling, one evaluation | 0.16 us | 0.16 us | 1.0x, even |

### Where It Is Slower

| Task | Classic | New | Change |
|---|---|---|---|
| `radial_solver`, 1 layer, 10 slices | 0.30 ms | 0.51 ms | 0.59x, 1.7x slower |
| `radial_solver`, 1 layer, 200 slices | 0.45 ms | 0.59 ms | 0.77x, 1.3x slower |
| `radial_solver`, 3 layers (static liquid core), 300 slices | 0.44 ms | 0.65 ms | 0.68x, 1.5x slower |
| `radial_solver`, propagation matrix, 200 slices | 0.11 ms | 0.19 ms | 0.56x, 1.8x slower |
| Convective cooling, 10k evaluations | 0.160 ms | 0.208 ms | 0.77x, 1.3x slower |
| Rheology, one complex modulus | 0.057 us | 0.077 us | 0.75x, 1.3x slower |

The standalone radial solver was about even with the classic one until it became a wrapper over the world path (one code path for both entry points). It now builds a temporary world from the supplied arrays, solves that world's equation of state, and integrates the Love-number equations against the same dense structure the world path uses. The two solvers take identical integration steps on these problems (71, 68, and 90 for the three independent solutions of the one-layer body), so the whole gap is the cost of each right-hand-side read: the new path evaluates the equation-of-state interpolant for gravity, the interpolated material for density and both static moduli, and the two supplied complex-modulus arrays, where the classic path did four linear interpolations of its input arrays. The equation-of-state solve itself is a tenth of the time, the temporary world a twentieth, and the Python-side handling about the same as the world build. What the dense read buys is accuracy: against the closed-form homogeneous sphere the new solver's degree-2 k2 is 2.6 times closer (6.8e-5 against 1.8e-4), and the two solvers agree to 1e-15 on the one-layer rows and 6e-10 on the three-layer one at the tolerances timed here (`integration_rtol` 1e-8, `integration_atol` 1e-12, both solvers).

The knobs that move it, in order. `integration_rtol` and `integration_atol` set the step count and so the read count: at the `[radial_solver]` defaults (1e-6 and 1e-10, looser than the rows above) the new solver takes 0.31, 0.39, and 0.49 ms on the three shooting rows, even with or faster than the classic solver at its tighter setting, with k2 moving by 4e-9. The equation-of-state settings (`eos_rtol`, `eos_atol`, `eos_integration_method`) change the total by under 10 percent, and the slice count matters little once the searches are seeded. `RK45` for the Love integration is slower than `DOP853`, which reaches the tolerance in fewer steps. For repeated solves of one body, build a `LayeredWorld` instead: its equation of state is solved once, its Love solves are cached per degree and frequency, and `calc_tides` reuses them across modes, which is where the new backend's time went.

The two scalar rows have a known cause. A single scalar rheology call is dominated by the Python-to-C++ boundary rather than by the arithmetic, and the numba path crosses a cheaper boundary. Use the vectorized calls, where the new backend is about 3x faster, whenever there is more than a handful of values. The vectorized convective cooling gap has not been investigated.

### First Call

Steady-state timings leave out the startup cost. The classic backend compiles its numba kernels the first time they run and caches the machine code on disk, so the first session after installing or upgrading pays the full compile and every later session still pays to load and dispatch the cached kernels. The new backend has nothing to compile. Each figure below is the median first call in a fresh process:

| First call | Classic, first session after installing | Classic, later sessions | New |
|---|---|---|---|
| Tidal heating, degrees 2 to 4, e^10 | 7.1 s | 1.2 s | 0.12 ms |
| 3D heating map (50 x 16 x 32 x 8 times) | 4.7 s | 1.0 s | 3.6 ms |
| Build a planet with its interior (Io) | not measured | 1.4 s | 2.8 ms |
| Dual-body dissipation rates | not measured | 1.1 s | no single-call equivalent |
| Build a world from config | not measured | 0.24 s | 1.3 ms |

A script that computes one 3D map and exits spends about a second in the classic backend once its cache is warm, nearly five seconds the first time after installing, and about four milliseconds in the new one.

### Threads for 3D Grids

The 3D grid methods, `calc_3d_tides`, `calc_3d_stress_strain`, `calc_3d_displacements`, and `get_3d_tidal_heating_array`, take `num_threads`, which spreads the per-point evaluation over threads. The default of 1 leaves parallelism to the caller, such as a process pool, and every thread count returns identical values. The classic backend has no equivalent. The table times three grids of a homogeneous Io at degrees 2 to 3 with eccentricity, a non-synchronous spin, and obliquity, on 20 radii by 45 colatitudes by 90 longitudes, using the `tides_3d:*_1_thread` and `tides_3d:*_all_threads` tasks in `Benchmarks_x/Performance` on the same machine and its 16 hardware threads. Each figure is the lowest of three fresh processes, each taking the best of three batches.

| Grid | 1 thread | 16 threads | Change |
|---|---|---|---|
| Secular heating map | 166 ms | 40 ms | **4.1x faster** |
| Stress and strain, 4 times | 261 ms | 57 ms | **4.5x faster** |
| Displacements, 24 times | 291 ms | 69 ms | **4.2x faster** |

The gain stops well short of the thread count because the radial solves, about 20 ms of each call here, always run on one thread. The work after them grows with the grid while the solves do not, so larger grids gain more.

## Module Map

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

## Porting Examples

### Building a World and Getting Love Numbers

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

### Standalone Radial Solver

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

### Rheology Models

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

## Learning the New System

Start with the `Demos (_x)` notebooks in the navigation: `Basics` notebooks (configuration, world building, save/load), `Physics` notebooks (orbits, tides, rheology, Love numbers, 3D heating, thermal/EOS), and `Systems` notebooks (multi-world systems, coupled thermal-orbital evolution). The `Benchmarks (_x)` pages validate the new system against published results and track its performance.
