# The Future TidalPy Structure

_Updated: 2026-09-24_

TidalPy's internals were rewritten in C++ for v0.8.0. The new implementation lives in modules that carry a `_x` suffix (`structures_x`, `Tides_x`, `RadialSolver_x`, `rheology_x`, and so on) and lives side by side with the classic modules today. In a future major release the `_x` modules become the only TidalPy: the classic modules are removed and the new ones drop their suffix. Nothing about the classic API changes until then, but all new development to TidalPy will happen to `_x` modules, and the classic modules will only receive bug fixes. New projects should start with the `_x` modules, and existing ones should plan to switch. Support for the old modules will end no later than 2026-12-31.

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
* **Consistency.** Every physics module (rheology, cooling, radiogenics, viscosity, partial melting, equations of state, tides) follows a similar C++ class hierarchy, a name-based factory (`make_<module>`), direct callable functions, vectorized variants, TOML configuration, and binary files that can be saved and loaded from disk for fast and accurate reproducibility.

## Performance

Performance tests were run with the classic and the new backend. Ratios move with the machine and the problem size, so read them as rough magnitudes and measure your own workload before relying on any of them.

The new backend is much faster where the classic path called out to BurnMan or paid a numba compile, 2.5 to 7.3 times faster on 3D heating maps, about three times faster on vectorized rheology and world building, 3.2 to 3.5 times faster on global tidal heating, and slower on the standalone radial solver at tight tolerances, which now runs through the world path, on two vectorized sweeps, and on scalar rheology calls; all are listed with their causes.

### Where It Is Faster

| Task | Classic | New | Change |
|---|---|---|---|
| Build a planet with its interior (Io, 3 layers) | 180 ms | 0.72 ms | **250x faster** |
| Orbit-averaged 3D heating map (50 x 16 x 32) | 27.1 ms | 3.69 ms | **7.3x faster** |
| Radiogenic heating, one evaluation | 0.30 us | 0.045 us | **6.8x faster** |
| Global tidal heating, degrees 2 to 4, e^10 | 0.142 ms | 0.040 ms | **3.5x faster** |
| Global tidal heating, e^2 truncation | 0.061 ms | 0.018 ms | **3.4x faster** |
| Global tidal heating, e^4 truncation | 0.063 ms | 0.019 ms | **3.3x faster** |
| Global tidal heating, e^10 truncation | 0.074 ms | 0.023 ms | **3.2x faster** |
| Build a world from config (2 layers, no interior solve) | 1.38 ms | 0.43 ms | **3.2x faster** |
| Rheology, 10k complex moduli | 0.174 ms | 0.054 ms | **3.2x faster** |
| Instantaneous 3D heating map (50 x 16 x 32 x 8 times) | 26.8 ms | 10.6 ms | **2.5x faster** |
| Homogeneous Love numbers (closed form) | 0.16 us | 0.068 us | **2.4x faster** |

The planet-building row is the largest change. The classic path handed the interior to BurnMan, which does mineral-physics lookups and its own root finding: the new path integrates the equation of state in C++.

The 3D maps are timed at eccentricity truncation level 2 on both sides, where both backends keep the same heating terms (see below).

The global tidal heating and 3D map rows were timed again on 2026-09-24, after the eccentricity tables changed, on a busier machine that ran both backends about two to three times slower than for the other rows: compare their ratios, not their times, with the rest of the table.

The global tidal heating rows use the homogeneous Love method, which solves the same problem as the classic `quick_tidal_dissipation`, at eccentricity truncations both backends tabulate. At level n both keep every term of the heating through $e^n$: the classic tables hold the squared $G^2$ cut at $e^n$, and the new ones the unsquared $G$, whose products the new tide engines cut at $e^n$, so both sum the same heating terms. The new backend tabulates levels 2, 4, 6, 8, 10, 20, and 50 and promotes any other requested level to the next tabulated one, so a request for e^5 runs at e^6. Its cost follows the number of distinct forcing frequencies rather than the number of modes: Love numbers are solved once per frequency and degree, and the layer-averaged shear modulus the homogeneous methods need is formed once per frequency and shared by every degree, so adding degrees adds little. With the `radial_solver` Love method each frequency and degree is a full radial solve instead, and that solve dominates.

### Where It Is About Even

| Task | Classic | New | Change |
|---|---|---|---|
| Convective cooling, one evaluation | 0.16 us | 0.17 us | 0.97x, even |

### Where It Is Slower

| Task | Classic | New | Change |
|---|---|---|---|
| `radial_solver`, 1 layer, 10 slices | 0.30 ms | 0.45 ms | 0.66x, 1.5x slower |
| `radial_solver`, 1 layer, 200 slices | 0.45 ms | 0.50 ms | 0.89x, 1.1x slower |
| `radial_solver`, 3 layers (static liquid core), 300 slices | 0.45 ms | 0.54 ms | 0.82x, 1.2x slower |
| `radial_solver`, propagation matrix, 200 slices | 0.11 ms | 0.14 ms | 0.80x, 1.25x slower |
| Radiogenic heating, 10k times | 0.087 ms | 0.099 ms | 0.88x, 1.1x slower |
| Convective cooling, 10k evaluations | 0.157 ms | 0.204 ms | 0.77x, 1.3x slower |
| Rheology, one complex modulus | 0.056 us | 0.076 us | 0.73x, 1.4x slower |

The standalone radial solver was about even with the classic one until it became a wrapper over the world path (one code path for both entry points). It now builds a temporary world from the supplied arrays, solves that world's equation of state, and integrates the Love-number equations against the same dense structure the world path uses. The two solvers take identical integration steps on these problems (71, 68, and 90 for the three independent solutions of the one-layer body), so the whole gap is the cost of each right-hand-side read: the new path evaluates the equation-of-state interpolant for gravity, the interpolated material for density and both static moduli, and the two supplied complex-modulus arrays, where the classic path did four linear interpolations of its input arrays. The equation-of-state solve itself is a tenth of the time, the temporary world a twentieth, and the Python-side handling about the same as the world build. What the dense read buys is accuracy: against the closed-form homogeneous sphere the new solver's degree-2 k2 is 2.6 times closer (6.8e-5 against 1.8e-4), and the two solvers agree to 1e-15 on the one-layer rows and 6e-10 on the three-layer one at the tolerances timed here (`integration_rtol` 1e-8, `integration_atol` 1e-12, both solvers).

The knobs that move it, in order. `integration_rtol` and `integration_atol` set the step count and so the read count: at the `[radial_solver]` defaults (1e-6 and 1e-10, looser than the rows above) the new solver takes 0.25, 0.30, and 0.36 ms on the three shooting rows, faster than the classic solver at its tighter setting, with k2 moving by 4e-9. The equation-of-state settings (`eos_rtol`, `eos_atol`, `eos_integration_method`) change the total by under 10 percent, and the slice count matters little once the searches are seeded. `RK45` for the Love integration is slower than `DOP853`, which reaches the tolerance in fewer steps. For repeated solves of one body, build a `LayeredWorld` instead: its equation of state is solved once, its Love solves are cached per degree and frequency, and `calc_tides` reuses them across modes, which is where the new backend's time went.

The scalar rheology row has a known cause. A single scalar rheology call is dominated by the Python-to-C++ boundary rather than by the arithmetic, and the numba path crosses a cheaper boundary. Use the vectorized calls, where the new backend is about 3x faster, whenever there is more than a handful of values. The radiogenic sweep evaluates one exponential per isotope and time on both sides; numba vectorizes those exponentials and the C++ loop does not. The vectorized convective cooling gap has not been investigated.

### First Call

Steady-state timings leave out the startup cost. The classic backend compiles its numba kernels the first time they run and caches the machine code on disk, so the first session after installing or upgrading pays the full compile and every later session still pays to load and dispatch the cached kernels. The new backend has nothing to compile. Each figure below is the median first call in a fresh process:

| First call | Classic, first session after installing | Classic, later sessions | New |
|---|---|---|---|
| Tidal heating, degrees 2 to 4, e^10 | 6.5 s | 1.1 s | 0.13 ms |
| 3D heating map (50 x 16 x 32 x 8 times) | 4.4 s | 0.93 s | 5.5 ms |
| Build a planet with its interior (Io) | not measured | 1.4 s | 2.1 ms |
| Dual-body dissipation rates | not measured | 1.1 s | no single-call equivalent |
| Build a world from config | not measured | 0.21 s | 2.1 ms |

A script that computes one 3D map and exits spends about a second in the classic backend once its cache is warm, more than four seconds the first time after installing, and about six milliseconds in the new one.

### Threads for 3D Grids

The 3D grid methods, `calc_3d_tides`, `calc_3d_stress_strain`, `calc_3d_displacements`, and `get_3d_tidal_heating_array`, take `num_threads`, which spreads the per-point evaluation over threads. The default of 1 leaves parallelism to the caller, such as a process pool, and every thread count returns identical values. The classic backend has no equivalent. The table times three grids of a homogeneous Io at degrees 2 to 3 with eccentricity, a non-synchronous spin, and obliquity, on 20 radii by 45 colatitudes by 90 longitudes, using the `tides_3d:*_1_thread` and `tides_3d:*_all_threads` tasks in `Benchmarks_x/Performance` on the same machine and its 16 hardware threads. Each figure is the lowest of three fresh processes, each taking the best of three batches.

| Grid | 1 thread | 16 threads | Change |
|---|---|---|---|
| Secular heating map | 169 ms | 42 ms | **4.0x faster** |
| Stress and strain, 4 times | 265 ms | 57 ms | **4.7x faster** |
| Displacements, 24 times | 298 ms | 72 ms | **4.1x faster** |

The gain stops well short of the thread count because the radial solves, about 22 ms of each call here, always run on one thread. The work after them grows with the grid while the solves do not, so larger grids gain more.

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
