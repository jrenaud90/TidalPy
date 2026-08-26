# TidalPy Major Change Log

## Version 0.8.X

### Version 0.8.0 (2026-NNN)
_The `_x` in module and function names indicates experimental versions. This suffix will be dropped in future releases._

#### Fixes
* `RadialSolver`: Fixed bug in Takeuchi starting conditions where y6 was pulling the incorrect value. Kamata starting conditions were not affected.
* `TidalPy.Utilities_x.arrays.interp`: Fixed an out-of-bounds read in the C++ `c_binary_search_with_guess` (`interp_.hpp`) for short domains. The NumPy-derived cache-window guard `guess < length - LIKELY_IN_CACHE_SIZE - 1` underflowed under `std::size_t` (NumPy uses signed `npy_intp`), so for arrays of length <= 9 it read `array[guess + 8]` past the end and returned a stack-garbage-dependent wrong value. Rewritten as an addition (`guess + LIKELY_IN_CACHE_SIZE + 1 < length`) so it is correct for unsigned indices. This resolves the intermittently-failing `Test_Utilities_x/Test_Arrays/test_interp_01.py::test_scalar_matches_numpy` (a scalar interpolation on a 6-element domain that could return the wrong node) and removes a latent bug in any short-array `c_interp` use.
* `TidalPy.constants`: A user override of `debug.test_constant` supplied through `reinit()` was being lost. `update_constants_x` (run after the legacy `update_constants` during initialization) re-read `test_constant` from the `_x` config and clobbered the override. `update_constants_x` no longer touches `test_constant`; the debug knob is owned by the legacy config path (`debug.test_constant`), so `reinit({'debug': {'test_constant': ...}})` again propagates to `TidalPy.constants.test_constant`.
* `TidalPy.Material_x.eos`: Fixed the pressure-to-density inversion (`eos_invert_eta` in `material_eos_.hpp`) used by the Birch-Murnaghan and Vinet material equations of state. The root was bracketed over a fixed, over-wide compression range (`eta` in `[1e-3, 1e3]`) on the assumption that the pressure law is monotonic in `eta`. That holds for 3rd-order Birch-Murnaghan only when `K0' == 4`; for `K0' != 4` the finite-strain term `1 + (3/4)(K0'-4)(eta^(2/3)-1)` turns the pressure over and drives it negative at large `eta`, so `P(eta=1000)` was negative and the guard `if (target >= p_hi) return hi` fired for every positive target, returning `eta = 1000` (density = `rho0 * 1000`). `calc_density(P)` therefore returned garbage for any `P > 0` whenever `K0' != 4`, which in turn made `solve_eos` fail to converge and produced nonsensical `planet_mass_eos` for compressible worlds. The inversion now brackets the root by expanding outward from `eta = 1` (where `P = 0`) while the pressure moves monotonically toward the target, stopping at the turning point, so it stays inside the valid range. Validated against BurnMan 2.1.0: the `calc_density(P)` relation agrees to about one part in 10^8, and the full self-gravitating single-material planet (density, pressure, gravity, mass, central pressure) to a few parts in 10^5 (grid resolution), across a range of materials and `K0'` values (see `Benchmarks_x/EOS/`). `ConstantDensityEOS` and `InterpolatedEOS` were never affected.
* `TidalPy.Material_x.eos`: `InterpolatedEOS` never validated that its optional shear/bulk modulus and viscosity tables match the radius table's length; a longer table read past the end of the radius array during interpolation. Mismatched lengths now raise a `ValueError` at construction. The EOS solution's array-interpolation path also left the two viscosity output slots uninitialized (the dense path fills 11 outputs, the array path filled 9); they are now interpolated when stored and NaN otherwise. The EOS solver's integrator-failure message now includes the integrator's own message (previously unreachable).
* `TidalPy.RadialSolver_x`: The package-level `TidalPy.RadialSolver_x.radial_solver` re-export pointed at the original `TidalPy.RadialSolver.solver.radial_solver` rather than the `_x` solver in `TidalPy.RadialSolver_x.solver`. The two share the same call signature, but code importing from the package root silently exercised the original implementation. The re-export now points at the `_x` solver. (The `build_rs_input_*` helper re-exports remain shared with the original module by design.)
* `TidalPy.RadialSolver_x` / `TidalPy.RadialSolver`: `radial_solver` wrote one boundary-condition model per `solve_for` entry into a fixed 5-slot buffer without a length check, so a `solve_for` tuple with more than 5 entries overflowed the stack. Both solvers now reject more than 5 entries with a clear error.
* `TidalPy.RadialSolver_x` / `TidalPy.RadialSolver`: Neither solver validated that the density and complex modulus arrays match the radius array's length (or that the per-layer tuples/arrays agree with the layer count). A shorter density array was written past its end during non-dimensionalization, corrupting the heap (crash 0xc0000374 on Windows). Both solvers now validate every input length up front with a clear error.
* `TidalPy.RadialSolver_x`: A radial solution exported from a world (`release_radial_storage`) kept a non-owning pointer to the world's own EOS solution as its dense structure source; calling the exported solution's `eos_call`/`eos_call_si` after the world re-solved its EOS (or was destroyed) dereferenced freed memory. The pointer is now cleared on release and the exported solution answers EOS queries from its stored dimensional arrays. `eos_call` at exactly the planet's surface radius now resolves to the top layer instead of raising (interior interfaces still resolve to the layer above).
* `TidalPy.structures_x`: The world's cached frequency-independent radial-solver inputs only checked the layer count, grid size, degree, and non-dimensionalization flag, so mutating a layer's `is_solid`/`is_static`/`is_incompressible` flag between Love-number solves silently reused the previous assumptions. The cache now also verifies the per-layer flags and rebuilds when any changed.
* `TidalPy.structures_x`: `LayeredWorld.add_layer` accepted a non-owning layer view (obtained from another world's layer access) and moved its pointer into the receiving world, giving two worlds ownership of one C++ layer (a double free waiting to happen). Views are now rejected with a clear error. `solve_love_numbers_supplied` now requires contiguous input arrays (a strided view would previously be read as if dense) and rejects empty inputs.
* `TidalPy.Tides_x`: The tide entry points (`set_tide_config`, `collapse_global_tides`, `global_potential`, `tidal_potential_3d_modes`, `eccentricity_func`) now validate the eccentricity truncation up front with a clear message listing the tabulated levels (1..5, 10, 15, 20); an untabulated level previously surfaced as an opaque "global potential failed" error at solve time. `eccentricity_func` also no longer raises on its own default argument (the old default `'gen'` had no string mapping; the default is now the tabulated level 3), and the drivers' error messages name the real supported sets. `calc_love_numbers`/`calc_neg_imk` type their optional `solver_love` argument as `LoveNumbers` (an arbitrary object was previously cast unchecked, which could crash the process).
* `TidalPy.structures_x.configs`: The world `[tides]` table now gets the same unknown-key rejection as every other part of the schema (a typo like `eccentricity_trunc_levl` previously fell back to the defaults silently). The world builder promotes an untabulated eccentricity truncation from a stale configuration file once per session instead of warning on every build, its obliquity fallback now matches the shipped config default (`"off"`; the mirrors had diverged, so a stripped config silently enabled obliquity terms), an unknown layer `class` raises a named `ValueError` instead of a raw `KeyError`, and a `data_file` world missing `radius_m` gets a named error. The generated `TidalPy_Configs_x.toml` header wrote the literal text `# " + "=" * 117 + "` (a refactor artifact) instead of a divider line.
* `TidalPy.utilities.arrays`: The interpolation helpers in `utilities/arrays/interp_.hpp` are now a thin pointer-style front end over the shared implementation in `Utilities_x/arrays/interp_.hpp`, closing out the short-array hazards of the old standalone copy (unsigned-underflow clamps, cache-window reads past the end for arrays shorter than 9, no length 0..2 handling, non-`inline` definitions in a header, and an initial-index-guess formula that always produced 0). The caller-provided index now seeds the search rather than being trusted blindly, so a stale index can no longer select the wrong interval.
* `TidalPy.radiogenics_x`: The `llri_and_slri` isotope dataset used the present-epoch reference time (4600 Myr) even though its abundances (including the canonical formation ratio 26Al/27Al = 5e-5 and the elevated undecayed long-lived concentrations) are quoted at solar-system formation. Evaluating it in its documented use window (the first few Myr) therefore overflowed the short-lived isotopes' heating to infinity. Its reference time is now 0 (formation), and the isotope decay exponential is guarded so an overflowing back-extrapolation returns NaN instead of infinity.
* `TidalPy.stellar_x`: The mass-to-luminosity piecewise relation switched to its linear high-mass branch at 20 solar masses, producing a factor ~13 luminosity discontinuity there (inherited from the original implementation). The branch point is now ~55 solar masses, where the adjacent power-law and linear branches meet continuously.
* `TidalPy.viscosity_x`: The reference and Arrhenius viscosity laws now treat a non-positive temperature as the cold limit and return an explicit infinite viscosity (a rigid, purely elastic response) instead of relying on IEEE division; an overflowing Arrhenius exponent saturates to the same limit. The Spohn and Henning partial-melt power laws and exponentials are guarded against overflow.

* `TidalPy.rheology_x`: The Voigt, Burgers, and Sundberg docstrings described `voigt_modulus_frac` as a "compliance fraction" with a default of `0.2`; the parameter is a modulus fraction (the Voigt element modulus as a multiple of the layer modulus) defaulting to `5.0`, so a user following the docstring got a 25x weaker Voigt element than intended. The docstrings now describe the real parameter. `check_schema_compatibility`'s Cython declaration also typed its return as a Python `bool` instead of the C++ `cpp_bool` (a latent mis-generation trap; the method was unused from Cython).

#### Refactors

##### Codebase-wide (`_x` cleanup pass)
* Consolidated the shared numerical helpers into `TidalPy/Utilities_x/math_x/numerics_.hpp`: the NaN-safe `c_isclose` (previously duplicated between `utilities/math/numerics_.hpp` and a local copy in the EOS solution) plus new guarded `c_safe_pow` / `c_safe_exp` helpers that return NaN instead of overflowing to infinity. The old `utilities/math/numerics_.hpp` forwards to it. The guarded helpers are used by the growth-law hot spots (Birch-Murnaghan/Vinet pressure laws, Arrhenius stress/grain-size power laws, Spohn/Henning partial melt, isotope decay).
* Near-zero divide guards use `TidalPyConstants::d_EPS` instead of bare `== 0.0` / `<= 0.0` comparisons in the partial-melt fraction, the viscosity laws, the EOS density inversion, the convective-cooling edge cases, and the layered world's bulk-density estimate.
* All `_x` Cython now types booleans as `cpp_bool` (from `libcpp`) instead of `bint`.
* American-English spelling pass over all `_x` comments, docstrings, and messages, and a sweep replacing internal development-plan references with timeless descriptions of current behavior.

#### Demos, Benchmarks, and Integration Tests
* Added `Demos_x/` — a 12-notebook tutorial curriculum for the `_x` Python API, numbered and grouped so a newcomer walks them in order: `Basics/` (configuration, world building, save/load), `Physics/` (orbits and insolation, tidal-heating basics, Io rheology, gas-giant fixed-Q vs fixed-dt, 1D Love numbers, 3D tidal heating, thermal/EOS), and `Systems/` (multi-world systems, a coupled thermal-orbital evolution integrated over 1 Gyr with CyRK). Each notebook is executed and committed with inline figures. Wired into the docs (`conf.py` copies them; `index.md` lists them).
* Added `Benchmarks_x/` — validation and performance benchmarks on the `_x` API. Validation: Earth load Love numbers versus Guo et al. (2004) and Farrell (1972) (`RadialSolver/Earth_Love_Numbers.ipynb`); Enceladus interior tidal `y`-functions versus Tobie et al. (2005) and Roberts & Nimmo (2008) (`RadialSolver/Enceladus_Tobie_Roberts.ipynb`); and an equation-of-state comparison against BurnMan (`EOS/EOS_vs_BurnMan.ipynb`), which agrees to machine precision and runs several times faster. Performance: a dependency-light harness (`Performance/`) that times the common tasks the demos teach and records each run as JSON tagged with date, TidalPy version, git commit, OS, machine, and CPU, so performance can be tracked over time (`Perf_Trends.ipynb` plots the trends). `pandas` is now a core dependency (it drives the trend views).
* Added end-to-end integration tests (`Tests/Test_E2E_x/`) that build worlds and a system from schema-0.2.0 TOML fixtures and run the full pipeline: terrestrial with a fixed-Q tide, terrestrial with a Maxwell rheology through the radial solver, gas giant with a fixed-time-lag tide, and a star-planet system's one-step orbital/spin evolution and binary round-trip.
* Performance-harness robustness: a failing task records an error row instead of discarding the whole run's timings, results files carry a `schema_version` stamp and the run label in their name (with a collision counter instead of silent overwrite), unreadable results files are skipped with a message, and the trend notebook exits cleanly when no results exist yet. Several demo/benchmark notebook claims were corrected against their own outputs (the constant-density vs Birch-Murnaghan central-pressure comparison, the BurnMan agreement figures, synchronous-spin wording) and the multi-world demo now actually reassigns the system host in its role-changing example.

#### `TidalPy.RadialSolver_x`
* All integration entry points now accept CyRK's implicit (stiff) methods `BDF`, `LSODA`, and `Radau` alongside the explicit `RK23`/`RK45`/`DOP853` (CyRK >= 0.18). This covers the standalone `radial_solver` (both `integration_method` and `eos_integration_method`), the original `TidalPy.RadialSolver.radial_solver` (same two arguments), and the world-level `LayeredWorld.solve_eos` / `solve_love_numbers` / `solve_love_numbers_supplied`. Method names are case-insensitive; an unknown name lists the supported set. Tested against the DOP853 reference for both solvers and the world solves (`Tests/Test_RadialSolver_x/test_i_implicit_integrators.py`, `Tests/Test_Structures_x/Test_Worlds/test_world_implicit_integrators_01.py`). One caveat: the whole-planet EOS integration starts at the planet's singular center (r ~ 0), where LSODA's startup can fail to take its first step (it reports a clean failure); BDF and Radau handle the singular start and match the explicit reference (Radau to machine precision).
* Reworked the shooting method to use a **dense radial-solution calling system** instead of a fixed grid. Each layer's independent solutions (3 solid / 2 dynamic-liquid / 1 static-liquid) are now integrated with CyRK dense output and the whole interpolant is retained; the surface boundary condition and the downward interface propagation produce a set of collapse constants that are stored per layer. The radial solution at any radius is then obtained by evaluating the dense interpolants there and collapsing on the fly (`c_RadialSolutionStorage::get_radial_solution`), rather than collapsing onto a fixed grid and linearly interpolating between slices. This is materially more accurate between grid points (the previous approach could be inaccurate for coarse grids) and removes the radial solver's dependence on the gridded `full_solution_vec`. The propagation-matrix method is unchanged: it still constructs its solution on the provided grid, and `get_radial_solution` serves it via linear interpolation of that grid (the documented, less-accurate behavior intrinsic to that method).
  * New Python API on `RadialSolverSolution`: `get_radial_solution(radius)` and the vectorized `get_radial_solution_array(radius_array)` return the complex `y1..y6` (SI) at arbitrary radii, evaluated in fast C++. The legacy `.result` / `.love` arrays are unchanged (the standalone solver samples the dense interpolants back onto the EOS grid for backward compatibility).
  * The surface y-solution (and the Love numbers derived from it) is computed once per solve and cached for both the shooting and matrix methods; interior radii remain evaluated on demand.
  * The structural arrays (radius, gravity, density, ...) continue to be owned by the EOS class; the radial solver holds no solution grid of its own for the shooting method.
  * Added a dense, SI-aware EOS readout at arbitrary radius (`RadialSolverSolution.eos_call_si(radius)`, C++ `c_RadialSolutionStorage::get_eos_si`): it maps the SI radius into the solver's non-dimensional domain and evaluates the solution's own dense EOS interpolant, returning SI gravity, pressure, mass, moment of inertia, density, and the complex shear/bulk moduli (the dense analogue of `get_radial_solution` for the structural/material state, where `eos_call` only accepts a raw non-dim radius). To keep this valid after the standalone solve returns, the per-layer EOS interpolation inputs are now persisted in the solution storage in non-dimensional solve units; previously they were `c_radial_solver` stack locals that dangled after return (and were re-dimensionalized in place), so a post-solve dense readout of the extra outputs produced doubly-scaled garbage. The 3D-tides kernel uses this to read complex moduli on radius via the dense path.
  * Added a frozen benchmark (`Tests/Test_RadialSolver_x/test_dense_benchmark/`) that compares the dense system against results frozen from the original `RadialSolver` (stored as `.npz`, since that solver will eventually be removed). Configurations: 1-layer solid, 2-layer solid/solid, and solid / dynamic-liquid / solid at a short (well-conditioned) forcing period. The dense and original solutions agree very well at the surface and at layer boundaries and the Love numbers agree tightly; the dense evaluation is the more accurate one between grid slices.
* The shooting method now measures how strongly its surface boundary condition solve amplifies error: the collapse constants can grow enormous and cancel (deep starting radii, high harmonic degrees), amplifying roundoff and integration error into the surface solution and the Love numbers. Both `TidalPy.RadialSolver_x.radial_solver` and the original `TidalPy.RadialSolver.radial_solver` record the factor on the returned solution (`surface_solve_amplification`) and log a warning when the resulting roundoff floor (amplification times machine epsilon) exceeds the requested `integration_rtol`, or when the amplification is severe (> 1e8). The layered world's Love solves emit the same warning and expose the diagnostic as `LayeredWorld.love_surface_amplification`. The diagnostic is computed only when the solve runs with `warnings` enabled, so performance-critical loops that disable warnings pay nothing for it. A healthy automatic-starting-radius solve sits near 1e6 (floor ~1e-10); the degree-3 deep-starting-radius pathology reports ~1e14 (floor ~5e-2), quantitatively matching its observed Love-number error. Tested in `Tests/Test_RadialSolver_x/test_j_surface_conditioning.py`.
* Numerical note - manual starting radii deep in the planet. Starting the shooting integration essentially at the center (e.g. `starting_radius=0.1` m on a 6000 km planet) at degree 3 with a dynamic incompressible layer leaves the independent solutions spanning so much dynamic range that the surface boundary-condition constants reach ~1e11: the surface solve then amplifies integration error into the leading digits of the Love numbers (RK45 at rtol 1e-7 reports success with k = 0.045 while the converged value is 0.1374, and a 0.001% rtol change moves k by 50%), and the surface `y2`/`y4` values pinned to zero by the boundary conditions come out as pure cancellation roundoff. The original and `_x` solvers behave identically here (with converged tolerances, DOP853, or the automatic starting radius they agree to 8+ digits, so this is inherent conditioning, not a regression); the old-vs-new comparison test now tightens tolerances in that regime so it compares converged solutions. Prefer the automatic starting radius (`starting_radius=0`); if a deep manual start is required, tighten `integration_rtol` and take the "large number of steps" warning seriously.
* Numerical note - dynamic liquid layers and long forcing periods. A **dynamic** liquid layer sandwiched between solids (e.g. solid inner core / dynamic-liquid outer core / solid mantle) is well-conditioned only at **short forcing periods** (high frequency). At long periods (low frequency) the dynamic formulation's `1/omega^2` terms make the independent solutions grow exponentially through the liquid, the surface boundary-condition matrix becomes near-singular, and the solve is unstable - the original and `_x` solvers then diverge from each other and eventually fail (NaN). This is inherent to the dynamic-liquid assumption, not a regression: the `_x` forward integration is bit-identical to the original, and at short periods the two agree exactly (validated by the benchmark). Use the **static** liquid assumption for long-period forcing (it is stable and consistent across all periods).

#### Package
* Added helper function `TidalPy.get_include` to get paths to cpp/hpp source files so they can be included in the build process of dependent packages (similar to `numpy.get_include`).
* The TidalPy data/config directories are now scoped to the package's `<major>.<minor>.X` version (e.g. `.../TidalPy/0.8.X/`) instead of the full patch version. Every patch release of a given major.minor now shares one directory, so user configs and downloaded data are not duplicated (or lost) on each bugfix release. New helper `TidalPy.paths.get_data_version()` returns the scoped label; `get_config_dir`, `get_log_dir`, `get_worlds_dir`, and `get_worlds_x_dir` all use it.
* Drops support for Python 3.9

#### `TidalPy.Utilities_x`
* Added an associated-Legendre utility (`TidalPy.Utilities_x.legendre`) providing the unnormalized associated Legendre functions `P_lm(cos theta)` with their first and second colatitude derivatives (Condon-Shortley phase, matching `scipy.special.assoc_legendre_p` with `branch_cut=2`). Two paths: `legendre(l, m, colatitude)` uses fast precomputed closed-form tables for degrees `l = 2..10` (hard-coded polynomials in `cos`/`sin`, no `std::pow` and no `1/sin theta` interior singularities; generated by `legendre/codegen/gen_legendre.py` via the Ferrers construction), and `legendre_generic(l, m, colatitude)` handles any degree via the vendored `xsf` library (value + derivatives through `xsf` dual numbers). Both agree with `scipy` and each other to ~1e-11 across all supported `(l, m)`; the `l = 2` table reproduces the tidal kernel's historical hand-coded forms exactly. Header-only C++ (`legendre_driver_.hpp` / `legendre_generic_.hpp`, namespace `tidalpy`) plus a Cython wrapper. These are the angular building blocks for the forthcoming class-free, truncation-driven 3D tidal potential engine.

#### `TidalPy.Tides_x`
* Created a new module called "Tides_x" where cythonized C++ code related to tide functions will be stored.
  * In a future release we will remove the old `TidalPy.tides` module in favor of this one (refactoring it to `Tides` to follow the same capitalization scheme as `RadialSolver`).
* Created C++ obliquity functions in `TidalPy.Tides_x.obliquity` a helper function is available in Python to call these `TidalPy.Tides_x.obliquity.obliquity_func`. See documentation for more details.
* Created C++ eccentricity functions in `TidalPy.Tides_x.eccentricity` a helper function is available in Python to call these `TidalPy.Tides_x.eccentricity.eccentricity_func`. See documentation for more details.
* Added a global (1D) tidal dissipation model hierarchy (`TidalPy.Tides_x.classes`) following the rheology/cooling/viscosity class pattern. Each model maps a per-mode Love number to the dissipation multiplier `-Im[k_l]` that the mode collapse multiplies into the per-mode potential terms produced by `c_global_potential`. All quantities MKS; supported degrees `l = 2..10`. The full complex Love-number suite (`k`, `h`, `l`) is the transport type everywhere (the `LoveNumbers` / `c_LoveNumbers` container) so the displacement Love numbers from the radial solver are retained even though only `k` drives heating/dynamics; the analytic models return `h`, `l` as `NaN` (no radial solution). Each model exposes `calc_love_numbers(degree_l, frequency, solver_love=None) -> LoveNumbers` and `calc_neg_imk(...)`.
  * `RheologyTide` (alias `rheology`) — `k_l` supplied by the radial solver (frequency dependent; driven by the world's tide solve).
  * `FixedQTide` (aliases `cpl`/`fixed_q`) — constant phase lag, `k_l*(1 - i/Q_l)` → `-Im[k_l] = k_l/Q_l`.
  * `FixedLagTide` (aliases `ctl`/`fixed_dt`) — constant time lag, `k_l*(1 - i*omega*dt_l)`.
  * `CTLQTide` (aliases `ctl_q`/`fixed_dt_q`) — `k_l*(1 - i*omega*dt_l/Q_l)`.
  * Name factory `make_tide(name, config)`, enum + binary-dispatch factory (`c_tide_from_binary`), and binary serialization (class ids 901–904). Per-degree parameters (`fixed_k`, `fixed_q`, `fixed_dt`) are lists indexed from `l = 2`.
* Added the standalone global mode collapse `TidalPy.Tides_x.classes.collapse_global_tides(...)`: runs the global-potential engine for an orbital/spin state and collapses the per-mode terms with an analytic tide model into the total tidal heating [W] and the three orbital potential derivatives (`dUdM`, `dUdw`, `dUdO`). For a synchronously rotating, low-eccentricity body the `fixed_q` result reproduces the classic CPL heating rate `(21/2)(k2/Q) G M_host^2 R^5 n e^2 / a^6` exactly. The `rheology` model needs the radial solver and is driven by the world instead.
* Wired the global tidal solve into the world: `LayeredWorld.calc_tides(orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass)` attaches-and-runs a tide model against the world's stored `[tides]` config, populating `get_tidal_heating()` [W], `get_tidal_potential_derivatives()` (`dUdM`, `dUdw`, `dUdO`), `get_num_tidal_modes()`, and per-layer `get_layer_tidal_heating(index)` (= world heating × the layer's `tidal_scale`). The tide model is attached via `set_tide_model(make_tide(...))` and the truncation/degree via `set_tide_config(...)`; the TOML world builder wires the `[tides]` table automatically (per-family default models: star→`fixed_q`, gasgiant→`fixed_dt`, terrestrial→`rheology`). The heavy global-potential engine is compiled into the world extension only (the orchestration lives in `world_tides_.hpp`, out-of-line from the lightweight `layered_.hpp`).
* Wired the `rheology` global tide model to the world radial solver. When `calc_tides` runs with the `rheology` model it solves the world's complex Love numbers once per unique tidal frequency (reusing the frequency-independent radial-solver cache), feeds the per-mode `k_l` into the mode collapse, and retains the full per-mode `k`/`h`/`l` suite. The EOS must be solved first (`solve_eos`); `calc_tides` raises a clear `RuntimeError` if it has not been, or if a per-frequency radial solve fails. The analytic models (`cpl`/`ctl`/`ctl_q`) keep their fixed-parameter path. Added `LayeredWorld.get_tidal_love_k(l, m, p, q)` to read a mode's radial-solver `k_l` (NaN for the analytic models). The world also writes each layer's tidal heating onto the layer object (`c_BaseLayer::get_tidal_heating()`, Cython `BaseLayer.get_tidal_heating()`), still readable from a built world via `get_layer_tidal_heating(index)`.
* Lifted the analytic global tide pipeline up to `BaseWorld` so **every** world type can dissipate tidally, including a layerless star. The tide-model holder, `[tides]` config, results, and the analytic `calc_tides` path (cpl/ctl/ctl_q) + the `set_tide_model`/`set_tide_config`/`tide_model_set`/`tides_solved`/`get_tidal_heating`/`get_tidal_potential_derivatives`/`get_num_tidal_modes`/`get_tidal_love_k` accessors now live on `BaseWorld` (inherited by `LayeredWorld`/`GasGiantWorld`/`StarWorld`). `LayeredWorld` extends `calc_tides` with the rheology (radial-solver) path and the per-layer heating distribution. A star supports the analytic models only; `calc_tides` raises if the `rheology` model is selected on a non-layered world (it needs the radial solver). The TOML builder now wires the `[tides]` table onto stars as well (default model `fixed_q`).
* Added a per-layer `tidal_scale_method` controlling how a layer's share of the world's global tidal heating is set (used by `LayeredWorld.calc_tides` to distribute the heat). `user_provided` (default) uses the layer's `tidal_scale` field; `volume_fraction` uses the layer-to-planet volume ratio; `tidal_timescale` is a log-Gaussian bell curve in the layer's Maxwell time `τ = η/μ` (from its static shear modulus + viscosity) peaking where `τ` equals the orbital forcing period `2π/n`, with a configurable width (`set_tide_config(tidal_timescale_width_decades=...)`, default 1 decade; in `defaultc_x.py`'s `[tides]`). It is a constructor argument on every layer class (`BaseLayer`/`PhysicsLayer`/`SolidLiquidLayer`/`GasLayer`, default `"user_provided"`), a settable Cython property (`layer.tidal_scale_method`, alias-aware), and an accepted layer TOML key, so a layer can pick its rule. Methods may differ per layer; non-tidal layers always get 0. The field is serialized in every layer's binary record.
* Added global tidal-dissipation defaults to the `_x` config (`defaultc_x.py` → `TidalPy_Configs_x.toml` → `TidalPy.config_x['tides']`): the default dissipation model per world family (`[tides.default_model]`: star→`fixed_q`, gasgiant→`fixed_dt`, terrestrial→`rheology`), the default truncation/degree (`min_degree_l`, `max_degree_l`, `eccentricity_trunc_lvl`, `obliquity_trunc_lvl`), and per-degree analytic parameters (`fixed_k`/`fixed_q`/`fixed_dt`, lists indexed from `l = 2`). The world builder reads these from `config_x` (world `[tides]` overrides them) rather than hard-coding, so the config file is the single source of truth.
* Added an on-demand 3D tidal stress/strain/heating kernel (`TidalPy.Tides_x.multilayer`). Rather than materializing a full `(radius × colatitude × longitude × time)` grid (the legacy `collapse_multilayer_modes` does, which blows up memory and cache), it computes the complex strain/stress tensors and the volumetric heating at a single point on demand by combining the radial solver's dense `y1..y6` at that radius with a 2D tidal potential. The strain kernel uses the exact Kervazo+2021-corrected forms with complex (viscoelastic) moduli at the mode frequency, a layer-type-dependent `dy1/dr`, and is restricted to solid layers (liquids contribute no shear dissipation). The complex moduli at radius are read through the radial solver's dense EOS path (`eos_call_si`). The compiled kernel (`Tides_x.multilayer.stress_strain`) is validated to machine precision against the legacy `calculate_strain_stress` / `calculate_volumetric_heating`; it is consumed by the world 3D path below.
  * The tidal potential is built by a class-free dynamic **Kaula engine** (`TidalPy.Tides_x.potential.potential_3d`, C++ `c_tidal_potential_3d_modes`) driven purely by the user's truncation levels (max degree `l`, eccentricity truncation, obliquity truncation) — no per-scenario potential-model objects. Following Kaula (1964) / Efroimsky & Williams (2009) Eq. 18, it enumerates the active `(l, m, p, q)` modes from the same eccentricity (`G_lpq`) and inclination/obliquity (`F_lmp`) functions the global (1D) path uses, times the associated Legendre functions `P_lm` (new `TidalPy.Utilities_x.legendre`), returning each mode's degree, signed forcing frequency `omega_lmpq = (l-2p+q) n - m*spin`, and potential angular factor (`U` + its colatitude/longitude derivatives). It is linear in `F_lmp`, `G_lpq`, `P_lm` (the global path squares `F`, `G` because global heating goes as the potential squared). This replaces and removes the earlier per-scenario potential-model class hierarchy (`c_TidalPotentialBase` + `SyncLowEPotential`/`NSRModesPotential`/`NSRMedObliquityPotential`, `make_tidal_potential`, binary ids 1001-1003) so any eccentricity/obliquity truncation is handled by one code path.
  * **World-level, all-C++ secular 3D heating:** `LayeredWorld.get_3d_tidal_heating(orbital/spin state, radius, colatitude)` returns the **secular** (cycle/orbit-averaged) 3D tidal volumetric heating for a built world — the physically time-averaged dissipated-power density `h_bar(r, theta) = sum_modes (omega_mode/2) Im(sigma_c : conj(eps_c))`, using the mode's complex potential amplitudes (a single `omega/2`, no `abs`; inherently non-negative for a dissipative material). It is longitude- and time-independent. The orchestration (build the active modes from the world's `[tides]` truncation config, solve the world radial response once per mode, accumulate each mode's secular contribution) lives on the rheology tide model (`c_RheologyTide::calc_3d_tidal_heating`), which the world delegates to; it calls the world's radial-solver/EOS members directly (no Python, no callbacks, no potential-model object). Requires the rheology tide model (`set_tide_model`) and a solved EOS; the analytic tide models are rejected (no depth-resolved solution). The truncation levels are read from `set_tide_config` (the TOML builder wires them from `[tides]`). **The volume integral of the secular 3D heating equals the 1D global tidal heating (`get_tidal_heating`) — both are the same total dissipated power — validated to ~0.1% for a homogeneous Maxwell sphere** (`Tests/Test_Structures_x/Test_Worlds/test_world_1d_vs_3d_tides_01.py`). (The earlier instantaneous-snapshot 3D path, which used `(omega/2)^2` and an `abs()` fudge and could not be reconciled with the 1D total, was replaced by this secular formulation.)
  * **Vectorized batch form for building maps:** `LayeredWorld.get_3d_tidal_heating_array(orbital/spin state, radii, colatitudes)` returns the secular 3D heating at an array of paired `(radius, colatitude)` points. The Kaula engine is split into a position-independent mode list (degree, order, parity, signed frequency, scalar amplitude) plus an on-demand angular-factor evaluator, so the world radial (Love-number) solve — which depends only on `(degree l, frequency)`, not on the query position — is built and solved once per unique `(l, frequency)` and reused across every point, rather than re-solved per point as a scalar loop would. It reproduces the scalar `get_3d_tidal_heating` point-for-point (to machine precision) and returns `NaN` where a radius has no depth-resolved solution (liquid layer / center / below the solver start). This is the efficient way to build a heating map (`Tests/Test_Structures_x/Test_Worlds/test_world_3d_tides_batch_01.py`).
  * **Full 3D grid + collapsed flavors + instantaneous power:** `LayeredWorld.calc_3d_tides(orbital/spin state, radii=None, colatitudes=None, longitudes=None, times=None, orbit_averaged=True, latitude_summed=False, longitude_summed=False, radial_summed=False)` produces the 3D tidal heating as a full grid over `(radius, colatitude, longitude[, time])` or reduced along any spatial dimension. With `orbit_averaged=True` (default) the quantity is the secular density `h_bar` [W m-3] (longitude/time-independent, so the default grid is `(nr, ncolat, nlon)` and constant along longitude); with `orbit_averaged=False` it is the instantaneous mechanical power density `sigma_ij(t) eps_dot_ij(t)` [W m-3] at each supplied time (a 4th axis), which depends on longitude and time and orbit-averages back to `h_bar` (modes with signed frequency `< 0` use the conjugated phasor at `+|omega|`, reproducing the true real field including cross-mode terms). Any spatial dimension can be integrated out (`latitude_summed`/`longitude_summed`/`radial_summed`): if any axis is summed the surviving spatial axes carry their Jacobian (`r^2`, `sin theta`, `1`) so a plain integral over them recovers the total; if none is summed the output is the raw density (it matches the scalar path). The colatitude integral uses an internal Gauss-Legendre grid (`latitude_nodes`), the radial integral an internal per-layer trapezoid (`radial_slices`), and the longitude integral the analytic `2*pi` when averaged or a `longitude_nodes` trapezoid when instantaneous. The returned dict carries the surviving axes plus either `heating` (the grid over the surviving axes, ordered radius, colatitude, longitude, time) or, when all three spatial axes are summed, `total` [W] and `per_layer` [W] (per-layer totals, innermost first; each an array over time when instantaneous). The fully collapsed `total` equals the 1D global `get_tidal_heating` (per-layer totals sum to it, replacing the `tidal_scale` distribution for the depth-resolved rheology path); the profiles integrate to it. The `latitude_nodes`/`radial_slices` defaults (16 each) were tuned so the collapsed total is within ~1% of the 1D heating for a homogeneous degree-2 body (~4 colatitude nodes and ~8 radial slices per layer suffice; the defaults add margin for higher degree l and layered bodies). Validated in `Tests/Test_Structures_x/Test_Worlds/test_world_3d_tides_collapse_01.py` (including the instantaneous power orbit-averaging back to the secular density and total).
  * **Analytic colatitude collapse (pre-integrated angular tables):** the secular `latitude_summed` colatitude integral is now done analytically instead of by Gauss-Legendre quadrature. The six angular functions the strain/stress needs form a bounded basis (`f1=P_lm, f2=dP/dtheta, f3=d2P/dtheta2, f4=P/sin, f5=-m^2 P/sin^2 + cot dP, f6=(dP - cot P)/sin`) whose sphere integrals `G_ij(l,m) = int_0^pi f_i f_j sin(theta) dtheta` are precomputed once (`TidalPy/Tides_x/multilayer/angular_gram_.hpp`, generated by `codegen/gen_angular_gram.py`; `l = 2..10`, `m = 0..l`), so the colatitude integral of each mode's heating reduces to `sum_k w_k sum_ij Im(sigma-coeff_ki conj(eps-coeff_kj)) G_ij` — no colatitude grid or Legendre evaluation. It is exact (machine-precision agreement with the Gauss-Legendre path) and, for a large radius grid (a radial profile or map), a few times faster than the quadrature (the Love-number solves otherwise dominate, so the fully radially-collapsed `total` is near break-even). It is the default; pass `latitude_analytic=False` to `calc_3d_tides` for the quadrature fallback. The Gram table is exposed for testing as `TidalPy.Tides_x.multilayer.stress_strain.angular_gram(l, m)`. Validated in `Tests/Test_Tides_x/Test_Multilayer/test_angular_gram_01.py` (table vs independent high-precision integration) and `Tests/Test_Structures_x/Test_Worlds/test_world_3d_tides_analytic_01.py` (analytic == numerical to machine precision for the total, per-layer, and radial-profile outputs, degrees `l = 2..6`).
  * Current limitations: obliquity modes and a vectorized materialized-map path are not implemented.

#### `TidalPy.RadialSolver_x`

#### `TidalPy.Material_x`
* Created a Cython-wrapped C++ module that duplicates the functionality of `TidalPy.Material` (TidalPy's EOS solver).
* Added a material equation-of-state (EOS) model hierarchy (`TidalPy.Material_x.eos.material_eos`) following the rheology/cooling/radiogenics class pattern. Each model returns density [kg/m³] from local pressure (analytic models) or radius (interpolated model) and is attached to a layer as the per-layer density source for the whole-planet EOS solve:
  * `ConstantDensityEOS` (alias `constant`) — incompressible.
  * `BirchMurnaghanEOS` (alias `bm`) — 3rd-order Birch-Murnaghan; density inverted from pressure (safeguarded Newton).
  * `VinetEOS` (alias `vinet`) — Vinet/UBER EOS; density inverted from pressure.
  * `InterpolatedEOS` (alias `interp`) — linear `density(radius)` lookup table (e.g. PREM profiles). It now optionally also carries radius-varying static shear modulus, bulk modulus, shear viscosity, and bulk viscosity tables (`shear_modulus_pa`, `bulk_modulus_pa`, `shear_viscosity_pas`, `bulk_viscosity_pas`), interpolated via the new `c_MaterialEOSBase::calc_static_shear_modulus` / `calc_static_bulk_modulus` / `calc_shear_viscosity` / `calc_bulk_viscosity` virtuals (analytic models return NaN). These viscoelastic quantities are now produced as **extra outputs of the EOS ODE**: the per-layer EOS pre-eval (`c_preeval_material_eos`) asks the model for them at each radius, so they ride the CyRK solution (and its dense output) alongside density/gravity/pressure rather than being re-interpolated in a separate pass. The EOS ODE gained two extra outputs for the static shear and bulk viscosity (`C_EOS_EXTRA_VALUES` 5 → 7, `C_EOS_DY_VALUES` 9 → 11; the dense-output buffers in `eos_data_.hpp`, `RadialSolver_x/shooting_.hpp`, and `RadialSolver_x/rs_solution.pyx` were widened to match). The world's per-layer viscoelastic pass reads these solution arrays (`complex_shear_array_vec`/`complex_bulk_array_vec` + new `shear_viscosity_array_vec`/`bulk_viscosity_array_vec`), falling back to the layer constant / viscosity model where the EOS returns NaN. So an interpolated (e.g. PREM) layer's moduli and viscosities vary with radius.
  * Name factory `make_material_eos(name, config)`, free pressure laws `birch_murnaghan_pressure` / `vinet_pressure`, enum + binary-dispatch factory (`c_material_eos_from_binary`), and binary serialization (class ids 601–604).
* Wired the material EOS models into the layer/world classes for the whole-planet EOS solve:
  * `BaseLayer.set_eos(model)` attaches a material EOS model to a layer (transfers ownership, like `set_cooling`); `BaseLayer.eos_set` reports whether one is attached.
  * `LayeredWorld.solve_eos(...)` integrates the planet's radial structure (gravity, pressure, enclosed mass, moment of inertia) from center to surface using each layer's EOS model as the local density source, reusing the `Material_x.eos` `c_solve_eos` machinery via a new `c_preeval_material_eos` pre-eval that dispatches to `model.calc_density(pressure, temperature, radius)`. It populates every layer's EOS profile and returns the radial profile arrays plus the converged surface/planet scalars. The full orchestration lives in C++ (`c_LayeredWorld::solve_eos(const c_WorldEOSSolveConfig&)`) so it can be called directly from other C++; the Cython method is a thin wrapper.
  * After a solve, `LayeredWorld.get_density(r)` / `get_gravity(r)` / `get_pressure(r)` (and the individual layers' getters) return the interpolated profile at any radius [m]; `eos_solved`, `surface_gravity_eos`, `central_pressure`, `planet_mass_eos`, and `planet_moi_eos` expose the scalar results.
  * Removed the superseded functional `Material_x.eos.solver.solve_eos` Python wrapper (its `c_solve_eos` C++ core is retained and reused) in favor of this class-based world method.

#### New Features

##### `TidalPy.dynamics_x` Module (new)
* Added a `TidalPy.dynamics_x` module for the spin and orbital rate equations that turn a tidal solve into the instantaneous evolution of a body's rotation and orbit (rates only; the System class integrates them). All quantities MKS.
* `Spin` — the spin-dynamics model, **attached to a world** and driven with the world's own EOS-based moment of inertia. `world.set_spin_model(Spin(moment_of_inertia_factor=...))`; then `world.get_moment_of_inertia()` returns the EOS-integrated value (`planet_moi_eos`) once the EOS is solved, else the model's uniform-density fallback `factor (2/5) M (R_o^5 - R_i^5)/(R_o^3 - R_i^3)`; `world.calc_spin_derivative(host_mass) = M_host dU/dO / I` (Ferraz-Mello et al. 2008) using the world's stored `dU/dO` (from `calc_tides`); and `world.calc_synchronous_spin(orbital_frequency)` (a tidally locked body spins at the orbital mean motion). The `Spin` model can also be used standalone (`calc_moment_of_inertia`/`calc_dspin_dt`/`calc_synchronous_spin`). Validated in `Tests/Test_Dynamics_x/Test_Spin/test_spin_01.py` and `Tests/Test_Structures_x/Test_Worlds/test_world_spin_01.py`.
* `OrbitSolver` — the orbital rate engine (Boué & Efroimsky 2019, CMDA, Eqs. 116-117): `calc_da_dt = (2/(n a)) dR/dM`, `calc_de_dt = (sqrt(1-e^2)/(n a^2 e))(sqrt(1-e^2) dR/dM - dR/dw)` (0 for a circular orbit), and `calc_dn_dt = -(3/2)(n/a) da/dt`, with `dR/dX = -((M_target + M_host)/M_target) dU/dX` from the dissipating body's tidal-potential derivatives. It is the low-level engine that the System class attaches and drives from its worlds. Validated in `Tests/Test_Dynamics_x/Test_OrbitSolver/test_orbit_solver_01.py`, and the world spin-rate + orbital rates together conserve energy (`tidal_heating = -(dE_orbit/dt + dE_spin/dt)`, machine precision) in `test_world_spin_01.py`.

##### `TidalPy.stellar_x` Module (new)
* Added a `TidalPy.stellar_x` module for stellar physics models. All quantities MKS.
* Added a luminosity model hierarchy (`TidalPy.stellar_x.luminosity`) following the rheology/cooling/radiogenics class pattern. Each model returns a star's luminosity [W] from its mass and shares the Stefan-Boltzmann effective-temperature conversions (`calc_luminosity_from_temperature`, `calc_temperature_from_luminosity`, both a function of the stellar radius) plus `calc_effective_temperature(mass, radius)` (mass → L → T). `calc_luminosity(mass)` accepts a float or an ndarray. The solar anchors come from the TidalPy constants (`d_MASS_SOLAR`, `d_LUMINOSITY_SOLAR`).
  * `FixedLuminosity` (alias `constant`) — luminosity supplied directly, independent of mass.
  * `MassToLuminosity` (aliases `cuntz_wang`, `cw`) — piecewise main-sequence `L(M)`: `0.23 (M/Msun)^2.3` below 0.2 Msun, the Cuntz & Wang (2018) polynomial exponent up to 0.85 Msun, `(M/Msun)^4` up to 2 Msun, `1.4 (M/Msun)^3.5` up to 20 Msun, and `3.2e4 (M/Msun)` above. Mirrors the legacy `TidalPy/stellar/stellar.py::luminosity_from_mass`, with the low-mass polynomial-exponent branch corrected to use the mass ratio `M/Msun` rather than the raw kilogram mass.
  * `PowerLawLuminosity` (alias `power_law`) — single power law `L = Lsun * coeff * (M/Msun)^exponent` (defaults `coeff=1`, `exponent=3.5`).
  * Name factory `make_luminosity(name, config)`, enum + binary-dispatch factory (`c_luminosity_from_binary`), and binary serialization (class ids 1001–1003).
* `StarWorld` can hold a luminosity model to derive its luminosity and effective temperature from its own mass and radius: `star.set_luminosity_model(MassToLuminosity())` (transfers ownership), then `star.calc_luminosity_from_mass()` [W], `star.calc_effective_temperature_from_mass()` [K], and `star.update_luminosity_from_mass()` (writes both onto the star's scalar fields). The model is optional; the mass-derived calls raise `RuntimeError` if none is attached. Validated in `Tests/Test_Stellar_x/Test_Luminosity/test_luminosity_01.py` and `Tests/Test_Structures_x/Test_Worlds/test_world_luminosity_01.py`.

##### `TidalPy.viscosity_x` Module (new)
* Added a viscosity model hierarchy (`TidalPy.viscosity_x.viscosity`) following the rheology/cooling/radiogenics class pattern. Each model returns the dynamic viscosity [Pa·s] as a function of temperature [K] and pressure [Pa] (the pre-melt "solid" viscosity that the partial-melt step weakens). The molar gas constant comes from the shared TidalPy config.
  * `ArrheniusViscosity` (alias `arr`) — Arrhenius flow law `η = A·σ^(1−n)·d^m·exp((E_a+P·V_a)/(R·T))` (optional extra factor of `T`).
  * `ReferenceViscosity` (alias `ref`) — relative-activation law `η = η_ref·exp(((E_a+P·V_a)/R)·(1/T − 1/T_ref))`.
  * `ConstantViscosity` (alias `const`) — temperature/pressure independent.
  * Name factory `make_viscosity(name, config)`, enum + binary-dispatch factory (`c_viscosity_from_binary`), and binary serialization (class ids 801–803). Math mirrors the validated legacy `TidalPy/rheology/viscosity/viscosity_models.py`.

##### `TidalPy.partial_melt_x` Module (new)
* Added a partial-melt (melt-weakening) model hierarchy (`TidalPy.partial_melt_x.partial_melt`) following the rheology/cooling/radiogenics class pattern. Each model maps a material's pre-melt viscosity + shear modulus (and temperature) to the post-melt viscosity + shear modulus and reports the volumetric melt fraction `φ = clip((T − solidus)/(liquidus − solidus), 0, 1)`. These quantities are frequency-independent and feed the downstream rheology (complex modulus) step of the whole-planet love-number pipeline.
  * `OffPartialMelt` (alias `none`) — no weakening (returns pre-melt strengths).
  * `SpohnPartialMelt` (alias `fischer`) — Fischer & Spohn (1990) temperature-based viscosity/shear law.
  * `HenningPartialMelt` — Henning (2009/2010) three-regime weakening (sub-critical, breakdown band, liquid-like), gated by `crit_melt_frac`/`crit_melt_frac_width`.
  * Name factory `make_partial_melt(name, config)`, enum + binary-dispatch factory (`c_partial_melt_from_binary`), and binary serialization (class ids 701–703). Math mirrors the validated legacy `TidalPy/rheology/partial_melt/melting_models.py`.

##### `TidalPy.structures_x` Module (new)
* Created `TidalPy/structures_x/` as the new C++/Cython module for world, layer, and system classes.
* Added a TOML-driven world builder and configuration system (`TidalPy.structures_x.configs`, schema version `0.2.0`). A world (the world object, its inner-to-outer layer stack, and each layer's attached physics models) is built from a single TOML file, a bundled world name, or a configuration `dict`. As elsewhere in `structures_x`, C++ never touches TOML: files are read/written at the Python level with the `toml` package and validated before the layer/world constructors and the `make_*` physics-model factories are called.
  * `build_world(source, force=False)` resolves `source` (bundled name, `.toml` path, or dict), validates it, and returns the built Cython world directly (a `BaseWorld` subclass: `LayeredWorld` / `GasGiantWorld` / `StarWorld`, per the world `type`). It is a thin wrapper over the new `BaseWorld.build(source, force=False)` static factory, which holds the build logic on the world class itself (there is no separate Python `World` wrapper). `build_world` and `BaseWorld.build` are re-exported from / available via `TidalPy.structures_x`.
  * The returned world retains the normalized build configuration on `world.source_config` (also exposed as the `world.config` property). `world.save_to_toml(path)` (and the lower-level `save_world_to_toml(config, path)`) write it back to TOML, stamping the current `schema_version`; so a build → save → reload cycle reproduces the same world. (`save_to_toml` falls back to `get_config_dict()` for a world constructed directly rather than via `build_world`.)
  * Lower-level builders `construct_world(config)` and `construct_layer(name, cfg, layer_index)` turn a validated dict into the underlying Cython objects; `available_worlds()` lists the bundled examples.
  * A `LayeredWorld` now hands its layers back as wrappers: `world.get_layer(i)` (0 = innermost, negative indices allowed), `world.layers` (inner-to-outer list), attribute access `world.<layer_name>` (e.g. `world.mantle`), and the sequence protocol — `for layer in world:`, `len(world)`, `world[i]`, and slicing `world[a:b]`. Each returns a **non-owning view** dispatched to the matching subclass (`PhysicsLayer`/`SolidLiquidLayer`/`GasLayer`/`BaseLayer`), so the layer's full Cython API is reachable (e.g. `world.mantle.shear_modulus_static`, `world.core.calc_complex_shear_modulus(r, ω)`, `world.mantle.get_tidal_heating()`). The world still owns the C++ layer; the view keeps the world alive (so it never dangles) and never double-frees it. The views are **built once and cached** (invalidated when a layer is added), so repeated `world.<layer>` access returns the same wrapper rather than re-allocating it. Dispatch uses a new `c_BaseLayer::get_layer_class_id()` virtual; attribute-by-name lookup is consulted only after normal attribute resolution (defined members win) and ignores names starting with `_`.
  * Loader/validation helpers in `configs.toml_loader`: `load_toml`, `validate_schema_version`, `validate_world_config`, `validate_layer_config`, and `merge_with_defaults`. The schema-version check is graded against `SCHEMA_VERSION`: a patch (`0.0.X`) difference is allowed silently, a minor (`0.X.0`) difference is allowed with a warning that functionality may break, and a major (`X.0.0`) difference raises a `ValueError` (use `force=True` to bypass). Validation rejects unknown keys (typo protection) and physics models attached to an incompatible layer `class` (e.g. cooling/radiogenics on a non-`solidliquid` layer).
  * Layers are built inner-to-outer, so a layer's inner radius is never written by the user (it is the previous layer's outer radius; 0 for the innermost, and supplying `radius_inner_m` is an error). Each layer sets its outer radius with exactly one of `radius_outer_m` (absolute), `radius_fraction` (× world radius), or `volume_fraction` (the layer's spherical-shell volume as a fraction of the whole-world volume, from which the outer radius is solved).
  * Each layer carries a `class` (`base` / `physics` / `solidliquid` / `gas`, selecting the Cython layer class) and an optional material `type` (`gas` / `mantle_rock` / `ice` / `hp_ice` / `iron`, selecting the per-material default block). Parameter and physics-model defaults resolve in three tiers: (1) the user's world TOML/dict, (2) the matching `[layers.<type>]` block of the new `_x` config (`TidalPy_Configs_x.toml`), filtered to what the layer `class` can hold, then (3) the C++/Cython constructor or physics-model-factory default. So a world can be specified by just `class` + `type` + geometry (as the bundled worlds now are), with all physics pulled from the material defaults.
  * Added the `_x` configuration system: `TidalPy/defaultc_x.py` builds `TidalPy_Configs_x.toml` (written to the user's TidalPy `Config` directory on first use, next to the legacy config) holding the `[numerical]` settings for the `_x` C++ config singleton plus per-material `[layers.<type>]` default blocks for `gas`/`mantle_rock`/`ice`/`hp_ice`/`iron`. It is loaded into `TidalPy.config_x` during initialization, and `TidalPy.constants.update_constants_x()` populates the shared C++ config singleton from its `[numerical]` section so all `_x` modules use the `_x` config. New default `_x` configuration should be added here, not to the legacy `defaultc.py`.
  * Added support for building a world from a **PREM-like radial data file** via a top-level `data_file` key (e.g. `data_file = "PREM.csv"`). The delimited file (comma/tab/whitespace) has columns radius [km], density [kg/m³], Vp [m/s], Vs [m/s] (optionally + shear and bulk viscosity). On build it is loaded and split into layers automatically by shear modulus (`Vs = 0` → liquid; each solid↔liquid transition starts a new layer, named `layer_0`, `layer_1`, … from the center out), and each layer gets an interpolated EOS carrying its radius-varying density and static shear/bulk moduli (and viscosities). The fast load/detection is in a new Cython module `TidalPy.structures_x.configs.prem` (`load_prem_arrays`, `detect_layer_boundaries`). Optional `[layers.layer_N]` tables refine the auto-detected layers (one per detected layer, matched in order; a provided radius must match the detected boundary; a user constant modulus/viscosity overrides the interpolated array; other keys like rheology sub-tables are merged in). A bundled `earth_prem` world (with `PREM.csv`) demonstrates this.
  * Bundled v0.2.0 example worlds (`earth_simple`, `jupiter_simple`, `sol`, `earth_prem`) ship in `TidalPy/WorldPack_x/` (alongside the legacy `WorldPack`). Mirroring the legacy world-pack behavior, they are copied into a version-scoped, user-editable data directory (`.../TidalPy/<version>/Worlds_x` via the new `TidalPy.paths.get_worlds_x_dir`) on first use; a bare-name lookup (`build_world("earth_simple")`) prefers the data-directory copy over the packaged one, so users can edit the installed TOML instead of the repository/package file. `install_worldpack_x(force=False)` performs the copy-if-absent install (or re-copies with `force=True`); `available_worlds()` lists the data-directory and packaged worlds combined.
  * A world-level `[tides]` table is preserved across load/save and wires the world's tide model (model selection, degree range, truncation levels); multi-world systems are built with `build_system` (see the System entries below).
  * New documentation pages under `Documentation/structures_x/config/`: `toml_schema.md` (the full v0.2.0 world/layer schema, default-resolution chain, and Python API) and `worldpack.md` (the world TOML file at a glance plus how WorldPack_x installs/resolves worlds).
* Added the `TidalPy.structures_x.worlds` world-class hierarchy (`c_BaseWorld` → `c_LayeredWorld` → `c_GasGiantWorld`, and `c_BaseWorld` → `c_StarWorld`):
  * `BaseWorld` — world identity plus albedo/emissivity/obliquity/spin scalars; `calc_surface_gravity`, `calc_escape_velocity`, `calc_mean_density`, and `calc_equilibrium_temperature(flux)` (fast-rotator radiative equilibrium). Binary class id 200.
  * `LayeredWorld` — owns an ordered (inner-to-outer) stack of layers via `add_layer` (transfers layer ownership, validates boundary continuity), with `num_layers`, `calc_total_mass`, `calc_internal_heating(time)` (sums radiogenic heating over `SolidLiquidLayer`s), and `validate_layers`. Binary class id 201.
  * `LayeredWorld.solve_love_numbers(...)` computes the whole-planet viscoelastic-gravitational Love numbers (k, h, l) for a tidal forcing frequency from the world's already-solved EOS (it never re-solves the EOS). Results are read via `love_k2`/`love_h2`/`love_l2`, `get_love_number_k/h/l(ytype_idx)`, `get_love_surface_y(...)`, and the `love_solved`/`love_success`/`love_error_code`/`love_message` status properties.
    * Two radial-solve methods are available behind the one call, selected by `use_prop_matrix`: the **shooting method** (default; arbitrary multi-layer worlds) and the **propagation-matrix method** (`use_prop_matrix=True`, with `core_model` selecting the core starting condition), the latter restricted to a single solid/static/incompressible layer. An incompatible world fails the solve gracefully (`love_success=False`, non-zero `love_error_code`) instead of raising.
    * The solve is backed by a cached, reusable C++ helper (`c_WorldRadialSolver`, owned by `c_LayeredWorld::p_radial_solver`) that separates the **frequency-independent** setup from the **frequency-dependent** work. The non-dimensionalization is itself frequency-independent (the `c_NonDimensionalScales` time scale is `1/(π·G·ρ_bulk)`, not `1/ω`), so the per-layer metadata, slice partitioning, non-dim structure arrays, non-dim scalars, and the reused `c_RadialSolutionStorage` are all built once and reused; only the complex moduli and the radial integration are recomputed per call. This makes repeated/swept solves (the hot loop for tidal and orbital-evolution calculations) allocation-free. The cache is invalidated automatically whenever `solve_eos` re-runs.
    * The long positional argument lists of the shooting (`c_shooting_solver`, 24 args) and propagation-matrix (`c_matrix_propagate`) solvers are replaced, on the world path, by two frequency-independent input structs (`c_ShootingInputs`, `c_MatrixInputs`) that the cache fills once; thin `c_shooting_solve`/`c_matrix_solve` wrappers expand a struct + frequency into the underlying call.
    * `LayeredWorld.solve_love_numbers_supplied(complex_shear, complex_bulk, radius_array, ...)` solves the Love numbers from externally-supplied complex moduli arrays (defined at `radius_array` [m], interpolated onto the world's EOS grid) instead of the layer rheology. It is the world-level entry point that the standalone array-based `RadialSolver_x.radial_solver` API delegates to.
  * `PhysicsLayer` exposes settable `is_solid` / `is_static` / `is_incompressible` properties (previously read-only) so the radial-solver layer assumptions (e.g. selecting the incompressible case required by the propagation-matrix method) can be set from Python; backed by new `c_PhysicsLayer::set_is_solid/set_is_static/set_is_incompressible`.
  * `GasGiantWorld` — a `LayeredWorld` for gas giants (default type `"gasgiant"`; binary class id 202).
  * `StarWorld` — a star with no layers/EOS; effective temperature ↔ luminosity via Stefan-Boltzmann (`L = 4·π·R²·σ·T⁴`). Binary class id 203. It can also hold a `TidalPy.stellar_x` luminosity model to derive its luminosity and effective temperature from its own mass (`set_luminosity_model`, `calc_luminosity_from_mass`, `calc_effective_temperature_from_mass`, `update_luminosity_from_mass`).
* Added a `System` class (`TidalPy.structures_x.system`) that links two or more worlds into a gravitationally bound group. It tracks two roles independently: the **tidal host** (the body a world orbits/raises tides on) and the **star** (the insolation source), which need not be the same body. Each world carries two two-body orbits — one about the tidal host and one about the star — because, for example, in the Earth-Moon-Sun system the Moon's tidal host is the Earth but its insolation comes from the Sun, whereas for an exoplanet the star is also the tidal host and the two orbits coincide. Worlds interact only with the host / the star, not with one another. The system co-owns each world with its Python wrapper (the world wrappers now hold the C++ world through a `shared_ptr`), so an added world stays fully usable and is the same object the system hands back.
  * Membership + roles: `add_world(world, is_host=False, is_star=False, semi_major_axis=None, eccentricity=0.0)` (the a/e are the tidal-host orbit); host via `set_host`/`host`/`host_index`/`has_host`; star via `set_star`/`star`/`star_index`/`has_star`/`get_star_luminosity`. Worlds are identified by index, name, or the world object, and the system is a sequence over its members (`for world in system`, `len(system)`, `system[i]`, slicing, and `system.<world_name>`).
  * Orbit about the tidal host: `get/set_semi_major_axis`, `get/set_eccentricity`, and Kepler helpers `calc_gravitational_parameter` (`G(M_host+M_world)`), `calc_orbital_frequency` (`sqrt(mu/a^3)`), `calc_semi_major_axis_from_frequency`.
  * Orbit about the star: `get/set_stellar_semi_major_axis`, `get/set_stellar_eccentricity`, `calc_stellar_gravitational_parameter`, `calc_stellar_orbital_frequency`.
  * Insolation: `calc_insolation_flux` (orbit-averaged incident flux `L_star / (4 π a² √(1-e²))` [W/m²], using the world's orbit about the star) and `calc_equilibrium_temperature` (the world's gray-body radiative balance `((1-A) F / (4 ε σ))^(1/4)` applied to that flux). The star's luminosity is read from the star world (a `StarWorld`).
  * Orbital and spin evolution: `calc_world_evolution(world)` solves an orbiting world's tides in the current system state (mean motion from Kepler's third law, spin/obliquity from the world, a/e from its orbit about the host, host mass from the host) and returns the orbital rates (`da_dt`, `de_dt`, `dn_dt`) from the `dynamics_x` rate engine plus the world's spin rate (`dspin_dt`) and the energy-balance terms (`dE_orbit_dt`, `dE_spin_dt`, `energy_residual`), verifying `heating = -(dE_orbit/dt + dE_spin/dt)`. In this single-body form only the orbiting world dissipates (the host is a point mass). `calc_system_evolution()` returns one such result per world.
  * Dual-body evolution: `calc_pair_evolution(world)` evolves an orbiting world together with its host, with both bodies raising a tide on the shared orbit. Each body's tides are solved with the other as the tide raiser (masses swapped); their orbital-rate contributions add and each body evolves its own spin, so the combined balance is `heating_world + heating_host = -(dE_orbit/dt + dE_spin_world/dt + dE_spin_host/dt)`. The result carries each body's full single-body contribution plus the combined shared-orbit rates and energy. A body with no tide model is rigid and contributes nothing, so a rigid host reduces the pair to the single-body result.
  * TOML builder: `TidalPy.structures_x.configs.build_system(source)` (a thin wrapper over `System.build`, mirroring `build_world` / `BaseWorld.build`) builds a whole `System` from a bundled system name, a `.toml` path, or a dict. Each `[worlds.<name>]` table names a `world` (a bundled world name, a path, or an inline world config, built via `build_world`) plus its `is_host`/`is_star` roles and its orbital elements about the tidal host (`semi_major_axis_m`/`eccentricity`) and the star (`stellar_semi_major_axis_m`/`stellar_eccentricity`); the table key becomes the world's name within the system (so a world template can be reused under different names). A bundled `sol_system` example (Sun + Earth + Jupiter) ships in `WorldPack_x/`. `construct_system(dict)` and `validate_system_config(dict)` are also exposed. `BaseWorld` gained a `name` setter (`set_name`) used by the builder.
  * A built system retains its normalized configuration on `System.source_config` and can be written back out with `System.save_to_toml(path)` (a faithful round-trip of the description, keeping world references) or via `System.get_config_dict()` (the self-contained live-state expansion, inlining each world's full config with its roles + orbital elements, used as the `save_to_toml` fallback for a directly-built system). `save_system_to_toml(config, path)` is exposed alongside `save_world_to_toml`.
  * Binary I/O: the Cython `System` now inherits `TidalPyBaseClass`, so it gains `save_binary` / `load_binary` / `get_schema_version_str` / `save_config`. `save_binary` writes the container state (name, host/star roles, each world's orbital elements about the host and the star) followed by every world's own record; `load_binary` rebuilds the heterogeneous world list via a new world binary-dispatch factory (`c_world_from_binary`, `structures_x/worlds/factory_.hpp`) that reconstructs each world's concrete subclass from its `BinaryClassID`, then re-wraps them as their concrete Python types (each world class gained a `_wrap` view constructor; `c_world_kind` discriminates the loaded type). As for a directly-loaded world, physics sub-models a world does not serialize (layer EOS data, tide/spin/luminosity models) are reattached after load; a star's effective temperature is serialized, so insolation survives the round-trip.
  * Binary class id 210. Backed by `c_System : c_TidalPyBaseClass` (`structures_x/system/system_.hpp`). Documentation: `Documentation/structures_x/system/system.md`; tests `Tests/Test_Structures_x/Test_System/test_system_01.py` and `test_system_02.py` (orbital + spin evolution), `Tests/Test_Structures_x/Test_Config/test_system_builder_01.py` (TOML builder).
  * Worlds serialize recursively: `save_binary` / `load_binary` round-trips a layered world together with every layer and each layer's attached rheology/cooling/radiogenics models (via the new `c_layer_from_binary` layer binary-dispatch factory). EOS profile data is excluded and repopulated by the EOS solve.
* Added `TidalPy.structures_x.layers.GasLayer` — C++ ideal-gas/fluid layer class (`c_GasLayer`, inherits `c_PhysicsLayer`).
  * Constructor adds `mean_molecular_weight_kg_mol`, `adiabatic_index`, `reference_temperature_k`, `reference_density_kg_m3`.
  * `calc_adiabatic_lapse_rate(g)` — dry adiabatic lapse rate: g·(γ−1)·M/(γ·R) [K/m].
  * `calc_scale_height(T, g)` — barometric scale height: R·T/(g·M) [m].
  * `calc_pressure_ideal_gas(T, rho)` — ideal-gas pressure: ρ·R·T/M [Pa].
  * `calc_sound_speed(T)` — adiabatic sound speed: sqrt(γ·R·T/M) [m/s].
  * No phase changes, cooling, or radiogenics sub-models.
  * Binary serialization includes all parent fields plus the 4 gas-property doubles; binary class ID 103.
* Added `TidalPy.structures_x.layers.SolidLiquidLayer` — C++ thermo-mechanical layer with phase-change tracking (`c_SolidLiquidLayer`, inherits `c_PhysicsLayer`).
  * Constructor adds 11 thermal/melt parameters: `thermal_conductivity_ref_w_mk`, `thermal_expansion_ref_1_k`, `heat_capacity_ref_j_kgk`, `activation_energy_j_mol`, `activation_volume_m3_mol`, `solidus_temperature_k`, `liquidus_temperature_k`, `melt_fraction_exponent`, `reference_density_kg_m3`, `reference_temperature_k`, `melt_viscosity_reduction`.
  * `calc_melt_fraction(T, P)` — power-law interpolation between solidus and liquidus [0, 1].
  * `calc_viscosity(T, P)` — Arrhenius temperature/pressure dependence with partial-melt exponential reduction.
  * `calc_shear_modulus(T, P)` — melt-fraction-reduced shear modulus: G_static·(1−φ).
  * `calc_thermal_conductivity(T)`, `calc_thermal_diffusivity(T)` — thermal transport properties.
  * `calc_adiabatic_temperature_gradient(T, P)` — uses EOS surface gravity when available.
  * `calc_heat_flux_conductive(T_base, T_top)` — conductive heat flux through the layer.
  * `calc_radiogenic_heating(time_s, mass_kg)` — delegates to optional `c_RadiogenicsBase` sub-model.
  * Optional cooling (`c_CoolingBase`) and radiogenics (`c_RadiogenicsBase`) sub-models attached via C++-level `set_cooling` / `set_radiogenics`; the concrete models are provided by `cooling_x` and `radiogenics_x`.
  * Binary serialization includes all parent fields plus the 11 thermal doubles, and recursively serializes the attached rheology, viscosity, partial-melt, cooling, and radiogenics sub-models; EOS profile data is excluded (re-solve after load). Binary class ID 102.
* Added `TidalPy.structures_x.layers.PhysicsLayer` — C++ mechanical-properties layer class (`c_PhysicsLayer`, inherits `c_BaseLayer`).
  * Constructor adds `shear_modulus_static_pa`, `bulk_modulus_static_pa`, `viscosity_static_pas`, `love_number_re`, `love_number_im` to the `BaseLayer` constructor parameters.
  * `love_number` property returns a Python `complex` (was `love_number_real: float`); `get_config_dict` now emits `love_number_re` and `love_number_im` keys.
  * `calc_tidal_susceptibility()` — geometrical tidal susceptibility (3/2)·r⁵/(G·m²) [m³].
  * `calc_complex_shear_modulus(freq)` / `calc_complex_bulk_modulus(freq)` — returns static modulus as real complex when no rheology is attached.
  * Rheology objects (`c_RheologyBase` subclasses) attach via C++-level `set_shear/bulk_rheology`.
  * Binary serialization includes all `BaseLayer` fields plus the four mechanical-property doubles, the Love numbers, the radial-solver classification flags, and (recursively) the attached rheology, viscosity, and partial-melt models; EOS profile data is excluded (re-solve after load).
* Added `TidalPy.structures_x.layers.BaseLayer` — C++ geometry-only layer class (`c_BaseLayer`, inherits `c_StructureBase`).
  * Constructor: `BaseLayer(name, layer_index, radius_inner_m, radius_outer_m, mass_kg, material_name, is_tidal, tidal_scale)`.
  * Read-only geometry properties: `radius_inner`, `radius_outer`, `thickness`, `volume`, `surface_area_inner`, `surface_area_outer`.
  * EOS profile: `eos_data_populated`, `update_eos_data(r, rho, g, p)`, `get_density(r)`, `get_gravity(r)`, `get_pressure(r)`.
  * Binary serialization and TOML config save inherited from `TidalPyBaseClass`.
  * Added `c_LayerEOSData` (header-only, `eos_data_.hpp`) for per-layer EOS interpolation data.

##### Recursive Binary Serialization (`_x` layers + physics models)
* Layers now serialize their attached physics sub-models recursively in `save_binary` / `load_binary`. A saved layer fully round-trips with its rheology, viscosity, partial-melt, cooling, and radiogenics models attached; no Python/Cython reconstruction step is needed after `load_binary`.
  * `PhysicsLayer` and `GasLayer` serialize their shear and bulk rheology models, their shear and bulk viscosity models, and their partial-melt model; `SolidLiquidLayer` additionally serializes its cooling and radiogenics models.
  * Every layer record also carries the `tidal_scale_method` selector and the radial-solver classification flags (`is_solid` / `is_static` / `is_incompressible`), so those survive a round trip on every layer class.
  * The world-attached tide, spin, and luminosity models (and the EOS profile data) are intentionally not serialized; reattach them (and re-run `solve_eos`) after loading.
  * Optional owned sub-objects are encoded as a one-byte presence flag followed (when present) by the sub-model's own binary record. New `binary_x` helpers: `write_optional_binary`, `read_optional_binary`, `optional_binary_flag_bytes` (alongside the existing `write_binary_string` / `read_binary_string`).
  * Each physics module gained a binary-dispatch factory (`c_rheology_from_binary`, `c_cooling_from_binary`, `c_radiogenics_from_binary`, `c_viscosity_from_binary`, `c_partial_melt_from_binary`) that peeks the record's `BinaryClassID` and reconstructs the correct concrete model. Every concrete model already carries a unique `BinaryClassID` (rheology 301–307, cooling 401–403, radiogenics 501–503, partial melt 701–703, viscosity 801–803).
* Added Python-level methods to attach sub-models to layers (ownership of the C++ model is transferred into the layer):
  * `PhysicsLayer.set_shear_rheology(rheology)` / `set_bulk_rheology(rheology)` (inherited by `SolidLiquidLayer` and `GasLayer`).
  * `SolidLiquidLayer.set_cooling(cooling)` / `set_radiogenics(radiogenics)`.

##### `TidalPy.cooling_x` Module (new)
* Created `TidalPy/cooling_x/` as the new C++/Cython module for cooling (heat-transport) models.
* Added the abstract base `c_CoolingBase` (inherits `c_PhysicsBase`) with pure-virtual `calc_cooling(c_CoolingInputs)` returning a `c_CoolingResult` (heat flux [W/m²], boundary-layer thickness [m], Rayleigh and Nusselt numbers). The eight physical inputs are bundled in a `c_CoolingInputs` struct (per the >5-argument style rule).
* Added the three cooling models, each exposed as a Cython/Python class:
  * `OffCooling` (alias `none`) — cooling disabled (zero flux; boundary layer = half thickness).
  * `ConductiveCooling` (alias `conductive`) — conduction: `q = k · ΔT / thickness`.
  * `ConvectiveCooling` (alias `convective`) — parameterized boundary-layer convection via the Rayleigh number; params `convection_alpha`, `convection_beta`, `critical_rayleigh` (Nusselt floor of 2). Uses the shared `minimum_layer_thickness` config floor.
  * Physics matches TidalPy's legacy `cooling.cooling_models` formulas.
* Added a rich `CoolingResult` container (`cooling_flux`, `boundary_layer_thickness`, `rayleigh`, `nusselt`; `to_dict`, iteration) whose fields are floats for scalar input or `float64` ndarrays for vectorized input.
* Added an enum-based C++ factory: `c_CoolingModel`, `c_cooling_model_from_name(name)` (alias-aware), and `c_find_cooling(model, config)` returning a `unique_ptr<c_CoolingBase>`. A name overload is also provided.
* Added `make_cooling(model_name, config=None)` — case-insensitive, alias-aware Python factory; unknown names raise `ValueError`.
* Added vectorized cooling methods on `c_CoolingBase` (inherited by all models): `calc_cooling_vectorize_temperature`, `..._vectorize_viscosity`, and `..._vectorize_all` (the two "live" inputs are the temperature drop and viscosity). The Cython wrappers return a `CoolingResult` of `float64` ndarrays.
* Added lower-case direct convenience functions (`cooling_off`, `conductive`, `convective`) that build a stack-allocated model, solve, and return a `CoolingResult`; `delta_temp_k` (and, for convection, `viscosity_pas`) accept floats or NumPy arrays (broadcast together).
* Each model supports `get_config_dict`, `save_config` (TOML), and `save_binary`/`load_binary` (binary class IDs 401–403).

##### `TidalPy.rheology_x` Module (new)
* Created `TidalPy/rheology_x/` as the new C++/Cython module for rheology (complex-compliance) models.
* Added the abstract base `c_RheologyBase` (inherits `c_PhysicsBase`) with pure-virtual `calc_complex_modulus(modulus, viscosity, frequency)` returning the complex (shear/bulk) modulus `μ*` [Pa] directly. Simple models are analytic; series composites (Burgers, Andrade, Sundberg) invert the sum of their element compliances internally (compliance is never exposed to Python).
* Added the seven rheology models, each exposed as a Cython/Python class:
  * `Elastic` (alias `off`) — purely elastic, no dissipation.
  * `Viscous` (alias `newton`) — purely viscous (Newtonian fluid).
  * `Voigt` (alias `voigt-kelvin`) — Voigt-Kelvin element; params `voigt_modulus_frac`, `voigt_viscosity_frac`.
  * `Maxwell` — standard Maxwell body.
  * `Burgers` — Maxwell + Voigt in series; params `voigt_modulus_frac`, `voigt_viscosity_frac`.
  * `Andrade` — Maxwell + Andrade transient term (∝ ω^{−α}); params `alpha`, `zeta`.
  * `Sundberg` (alias `sundberg-cooper`) — Andrade + Voigt; params `alpha`, `zeta`, `voigt_modulus_frac`, `voigt_viscosity_frac`.
  * Physics matches TidalPy's legacy `rheology.complex_compliance` formulas.
* Added an enum-based C++ factory: `c_RheologyModel` (one value per model), `c_rheology_model_from_name(name)` (alias-aware name → enum), and `c_find_rheology(model, config)` returning a `unique_ptr<c_RheologyBase>` to a heap-allocated model. A name overload `c_find_rheology(name, config)` is also provided.
* Added `make_rheology(model_name, config=None)` — case-insensitive, alias-aware Python factory; unknown names raise `ValueError`. It wraps the C++ enum factory (`c_rheology_model_from_name` → `c_find_rheology`) and adopts the returned `unique_ptr` into the matching rich Python wrapper.
* Added vectorized complex-modulus methods on `c_RheologyBase` (inherited by all models): `calc_complex_modulus_vectorize_modulus`, `..._vectorize_frequency`, and `..._vectorize_all`. Each writes into a caller-supplied `std::vector<std::complex<double>>&`; the Cython wrappers accept array-likes and return `complex128` NumPy arrays.
* Added lower-case direct convenience functions (`elastic`, `viscous`, `voigt`, `maxwell`, `burgers`, `andrade`, `sundberg`) that build a stack-allocated model, solve, and return. `frequency`, `modulus`, and `viscosity` accept floats or NumPy arrays (broadcast together); a scalar returns a Python `complex`, arrays return a `complex128` `ndarray`.
* Each model supports `get_config_dict`, `save_config` (TOML), and `save_binary`/`load_binary` (binary class IDs 301–307).
* Refactored the `PhysicsBase` Cython wrapper to be subclassable: subclasses own their most-derived C++ object via a `unique_ptr`, and `model_name` reads through the inherited `_ptr`.

##### `TidalPy.radiogenics_x` Module (new)
* Created `TidalPy/radiogenics_x/` as the new C++/Cython module for radiogenic-heating models.
* Added the abstract base `c_RadiogenicsBase` (inherits `c_PhysicsBase`) with pure-virtual `calc_heating(time_s, mass_kg)` returning the radiogenic heating `Q` [W] directly. All inputs/outputs are MKS (time and half-lives in seconds).
* Added a lightweight `c_Isotope` value type (no base class) describing one isotope (`name`, `heat_production_w_kg`, `half_life_s`, `mass_frac`, `concentration`) with `decay_constant()` and `specific_heating(time, ref_time)` helpers. `c_RadiogenicsConfig`/`c_IsotopeRadiogenics` carry a `std::vector<c_Isotope>`.
* Added the three radiogenics models, each exposed as a Cython/Python class:
  * `OffRadiogenics` (alias `none`) — radiogenics disabled, heating == 0.
  * `IsotopeRadiogenics` — sum of decaying isotopes; constructed from parallel arrays (`heat_production_w_kg`, `half_lives_s`, `mass_fracs`, `concentrations`, optional `names`) plus `ref_time_s`, or via `IsotopeRadiogenics.from_dataset(name)`.
  * `FixedRadiogenics` (alias `constant`) — single lumped rate with optional exponential decay; params `fixed_heat_production_w_kg`, `average_half_life_s` (≤0 disables decay), `ref_time_s`.
  * Physics matches TidalPy's legacy `radiogenics.radiogenic_models` formulas. The decay factor ln(0.5) is a module-level `constexpr` (`d_LN_HALF`).
* Added built-in literature isotope datasets (C++ `c_get_isotope_dataset` / `c_isotope_dataset_names`, exposed as Python `available_isotope_datasets()` and `isotope_dataset(name)`): `modern_day_chondritic` (Hussmann & Spohn 2004; Turcotte & Schubert 2001), `llri_and_slri` (Castillo-Rogez et al. 2007, adds short-lived Mn53/Fe60/Al26), and `bulk_silicate_earth` (McDonough & Sun 1995). All defined in MKS (Myr converted to seconds).
* Added an enum-based C++ factory: `c_RadiogenicsModel` (one value per model), `c_radiogenics_model_from_name(name)` (alias-aware name → enum), and `c_find_radiogenics(model, config)` returning a `unique_ptr<c_RadiogenicsBase>` to a heap-allocated model. A name overload is also provided.
* Added `make_radiogenics(model_name, config=None)` — case-insensitive, alias-aware Python factory; unknown names raise `ValueError`. For the isotope model it accepts a built-in dataset name (`isotopes`), explicit MKS arrays, or a global-config/inline dataset (stored in Myr and converted to seconds).
* Added vectorized heating methods on `c_RadiogenicsBase` (inherited by all models): `calc_heating_vectorize_time`, `..._vectorize_mass`, and `..._vectorize_all`. Each writes into a caller-supplied `std::vector<double>&`; the Cython wrappers accept array-likes and return `float64` NumPy arrays.
* Added lower-case direct convenience functions (`off`, `isotope`, `fixed`) that build a stack-allocated model, solve, and return. `time` and `mass` accept floats or NumPy arrays (broadcast together); a scalar returns a Python `float`, arrays return a `float64` `ndarray`.
* Each model supports `get_config_dict`, `save_config` (TOML), and `save_binary`/`load_binary` (binary class IDs 501–503; the Isotope model serializes its variable-length isotope list, including per-isotope names).

##### `TidalPy.Utilities_x` Module (new)
* Created `TidalPy/Utilities_x/` as the new C++/Cython foundation module housing base classes, logging, and binary I/O.
* Added `TidalPy.Utilities_x.arrays` — header-only 1-D linear interpolation utilities (`interp_.hpp`: `c_interp`, `c_interp_complex`, `c_binary_search_with_guess`) ported from `TidalPy.utilities.arrays` into the new scheme, with a `numpy.interp`-style Python wrapper `interp(x, xp, fp)`. The interpolated material EOS (`c_InterpolatedEOS`) now uses `c_interp` instead of a hand-rolled lookup.
* Added `TidalPy.Utilities_x.lookups` — integer-keyed lookup containers (`IntMapN` / `IntMapNComplex`, `N = 1..4`; C++ `c_IntMap<c_KeyN, ...>` + `keys_.hpp` key packing) that store a value against a packed key of up to four 16-bit integers. These back the tidal mode/frequency maps used throughout `Tides_x`. Relocated here from `TidalPy.utilities.lookups` so the `_x` tide code depends only on `_x` utilities.
* Added `TidalPy.Utilities_x.logging_x` — C++ logging via [spdlog v1.15.3](https://github.com/gabime/spdlog) with a Cython/Python wrapper.
  * `init_logger(config)`, `set_log_level(level)`, `shutdown_logger()` are available from Python.
  * C++ code uses `TIDALPY_LOG_DEBUG/INFO/WARN/ERROR/CRITICAL(...)` macros from `logger_.hpp`.
  * Added `spdlog` as a git submodule at `Dependencies/spdlog` (header-only, cross-platform).
* Added `TidalPy.Utilities_x.binary_x` — custom TidalPy binary file format with a fixed 20-byte header.
  * `check_binary_file(path)` and `get_current_schema_version()` available from Python.
  * C++ utilities in `binary_.hpp`: `write_binary_header`, `read_binary_header`, `check_binary_schema_version`, and `BinaryClassID` enum.
  * Added shared length-prefixed string serialization helpers `write_binary_string`, `read_binary_string`, and `binary_string_bytes` (used for model names and any variable-length text — one encoding in one place).
  * Schema version `0.2.0` (separate from package version); same `major.minor` required for compatibility.
* Added `TidalPy.Utilities_x.classes_x` — C++ base class hierarchy with Cython/Python wrappers.
  * `TidalPyBaseClass`: abstract base; provides `save_binary`, `load_binary`, `get_schema_version_str`, `get_config_dict`, `save_config`.
  * `StructureBase(radius_m, mass_kg)`: spherical geometry base with `calc_surface_area`, `calc_volume_sphere`, `calc_volume_shell`, `calc_surface_gravity`, `calc_mean_density`, `calc_escape_velocity` (all MKS).
  * `PhysicsBase(model_name)`: physics model base with `model_name` property and binary serialization. Provides shared `write_physics_binary(out, class_id, params)` / `read_physics_binary(in, force, n_params)` helpers so every physics model serializes uniformly (header + model name + scalar params) and a model's `write_binary`/`read_binary` reduce to one call each.

#### Utilities
* Converted `math.numerics` to c++.
* Implemented a new constant/parameter backend that can be accessed in C++ but modified in Python/cython.
  * Refactored `constants.d_DBL_MANT_DIG` to `constants.d_DBL_MANT_DIGITS` for readability.
  * Refactored `constants.d_PI_DBL` to `constants.d_PI`
  * Refactored `constants.d_NAN_DBL` to `constants.d_NAN`

#### Tests
* Added tests for `TidalPy.get_include`.
* Added tests for the new `IntMap`.
* Added tests for the new `obliquity` functions.
* Added tests for the new `eccentricity` functions.

#### Documentation
* Added documentation for `IntMap` in the Utilities module section.
* Added a note about `TidalPy.get_include` in the readme.
* Added documentation for `obliquity` functions in the Tides module section.
* Added documentation for `eccentricity` functions in the Tides module section.

#### Repository
* Fixes incorrect license url in codemeta.json.

#### Dependencies
* TidalPy now works on Python 3.14.
* Bumped version of `ipympl` to `<=0.11.0`.
* Bumped version of `CyRK` to `>=0.17.1, <0.18.1`.
* Increased numpy's max pinnings to `<2.5`.

##### `XSF` Submodule
* Adds [xsf]() package as a submodule to TidalPy. We use the spherical bessel function headers in various RadialSolver
  calculations.
  * Adds cython wrappers for spherical bessel functions to `TidalPy.utilities.math.special`.

##### `Eigen` Submodule
* Adds [Eigen]](https://gitlab.com/libeigen/eigen) package as a submodule to TidalPy.
  This provides much of the functionality of LAPACK without us having to compile it or find its symbols. Specifically
  we use it to do a LU-Decomp in `RadialSolver_x` module.

### Version 0.7.0 (2025-12-02)

#### Fixes
* Fixed top-level directories and files being mistakenly installed.
* Updated `np.trapz` to `scipy.integrate.trapezoid` due to deprecation of the former.
* Fixed bug in the henning melting model where shear modulus may not have had the correct value when fully melted.
* Fixed issues in some of the demo notebooks.
* Fixed notebooks to work with changes to Matplotlib.

#### GitHub Actions
* Added action to automatically update version numbers, dates, and commit hashs in `citation.cff`, `codemeta.json`, and `meta.yaml`.

#### Utilities
* Added `fontsize` kwarg to projection plotter in utilities.graphics.
* Changed the default of `yplot` so plotting imaginary parts of radial solution is by default off.

#### Documentation
* Added JOSS paper draft to Papers/2025-JOSS.
* Added documentation interface with Sphinx and ReadTheDocs.
* Added and updated readme badges.
* Extensive changes to all documentation files and structure.
* Changed how CHANGES.md (main change log) is structured. Major versions will appear as top level headers.

#### Other
* Change TidalPy's license to Apache 2.0.
* Added `CODE_OF_CONDUCT.md` with TidalPy's code of conduct.
* Added `CONTRIBUTING.md` with information on how others can contribute to TidalPy.
* Updated `citation.cff` with latest information.
* Added `codemeta.json` with citation and metadata.
* Removed MacOS-13 tests.
* Removed unused code, tests, and documentations.
* TidalPy now uses `ruff` to check for code issues.
  * Added ruff lint check to ubuntu tests.
* Fixed many code issues found with ruff (these are not bug fixes, just syntax compliance and best practices).
* Moved more module/file level constants into TidalPy.constants and in the configuration file.
* Added TidalPy logo.

## Version 0.6.X

### Version 0.6.11 (2025-11-06)

#### Dependencies
* Added more restrictive pinning for numpy's max version to work with latest CyRK.

### Version 0.6.10 (2025-11-05)

#### New
* RadialSolver solution can now provide quality factors `rs_solution.Q` and phase lag angles (in radians) `rs_solution.lag`.

#### Fixes
* Fixed issue where RadialSolver would incorrectly say a result was successful when the application of surface boundary conditions failed.
* Fixed issue where RadialSolver solution's diagnostics would print the incorrect parameter for a planet's mass.

#### Changes
* RadialSolver solution returns nan's instead of None for Love numbers when a solution fails.
* RadialSolver solution Love number attributes now return lists of np.nan if the user requested multiple solution types (e.g., "tidal, loading") and the solution failed. This allows the user to still subscript `rs_solution.k[0], rs_solution.k[1]` even if the solution was not successful.

#### Tests
* Added more tests to check RadialSolver Solution attributes.

### Version 0.6.9 (2025-09-19)

#### Fixes
* Fixed issue where `TidalPy.RadialSolver.shooting` would pick the incorrect starting index. If the starting layer (set by the starting radius) was not the first layer it could cause a int overflow and lead to access violation crashes. 

#### Dependencies
* Updates some GitHub action dependencies.

#### Tests
* Added more tests to check `RadialSolver` starting radius conditions to try to catch bugs like this patch fixed!

### Version 0.6.8 (2025-08-19)

#### Changes
* `TidalPy.tides.modes.collapse_multilayer_modes` now checks for ill-formed `radius_arrays` and will raise an error if it detects a problem.
* Improved how `projection_map` displays colorbar numbers.

#### Fixes
* Fixed issue where `TidalPy.tides.modes.collapse_multilayer_modes` was not treating macro-layer boundaries correctly.

### Version 0.6.7 (2025-08-19)

#### Changes
* `TidalPy.tides.modes.collapse_multilayer_modes` now requires a tuple of rheology instances for both shear and bulk rheology. One for each macro layer. This allows different rheologies to be used for different layers.

#### Fixes
* Fixed some warnings that were showing up in the multilayer demo notebook.

### Version 0.6.6 (2025-08-15)

#### Changes
* Updated to work with CyRK v0.15.1
  * This change results in about a 18% performance improvement for RadialSolver depending on the problem.
* Converted RadialSolver and EOSSolver status message to C++ strings.
* Updated some documentation.

#### Fixes
* Fixed potential issue where RadialSolver could use a very small or even negative maximum step size.

### Version 0.6.5 (2025-08-12)

#### New
* Added test to check if structure arrays have been changed.
* Added debug flag to installation files to help with cython debugging. 

#### Fixes
* Fixed issue where TidalPy structures (layers, planets, etc.) would return editable arrays instead of copies of arrays. This could lead to subsequent functions (like planet paint) changing the arrays. This fixes GitHub Issue [#74](https://github.com/jrenaud90/TidalPy/issues/74).
* Fixed issue where shooting method was corrupting data when starting radius was too large (which happened for higher orders of l) This fixes GitHub Issue [#72](https://github.com/jrenaud90/TidalPy/issues/72).
* Fixed broken demo notebooks.

### Version 0.6.4 (2025-04-10)

#### Dependencies
* Removed max version limit for platformdirs package

### Version 0.6.3 (2025-04-09)

#### Changes
* Changes GitHub actions to avoid testing when not needed.
* Now supports Python 3.13

#### Fixes
* Removed various files that were being included in the manifest

#### Dependencies
* Updated to work with CyRK 0.13.3
* Updated to work with Numpy 2.x

### Version 0.6.2 (2025-03-28)

* Bumping version to integrate into conda-forge.
* Fixed some build processes that could cause both numpy 1.X and 2.X to be installed and inconsistent use (only affected MacOS)

### Version 0.6.1 (2025-01-11)

#### Benchmarks and Performance
* Added performance benchmark for `TidalPy.RadialSolver.helpers`.

#### Fixes
* Fixed memory leak that was occurring in the EOS Solver (therefore also RadialSolver) due to CyRK's CySolverSolution not being dereferenced (this is a problem with CyRK, but hack applied to TidalPy until a fix can be made to CyRK).

#### Performance
* Improved `TidalPy.RadialSolver.helpers.build_rs_input_from_data` performance by factor of 10x to 150x depending on layer structure.
* Improved `TidalPy.RadialSolver.helpers.build_rs_input_homogeneous_layers` performance by factor of 15% to 3.5x depending on layer structure.
* Improved `TidalPy.RadialSolver.radial_solver` performance by around 10% depending on layer structure.

### Version 0.6.0 (2025-01-07)
#### Removed
* Removed support for the older non-cythonized `TidalPy.radial_solver` module in favor of `TidalPy.RadialSolver`

#### RadialSolver Changes
* Moved RadialSolver's Boundary Condition finder to its own function in `TidalPy.RadialSolver.boundaries.surface_bc.pyx` to allow it to be used by both the shooting and propagation matrix techniques.
* Decoupled radial solver from shooting method.
  * Moved the shooting method (formerly just called `cf_radial_solver`) to a dedicated file to prep for a different dedicated file for the prop matrix solver. 
  * Now `TidalPy.RadialSolver.solver` only contains driver functions and output structures.
* Added Propagation Matrix technique to RadialSolver
  * This is simplified for now. Only planets with 1 solid, static, incompressible layer are allowed. 
    * Other assumptions can be approximated, e.g., liquid layers use a small shear modulus.
    * Multiple layers should also work if you have discontinuities in density, shear, etc. within your "one layer".
  * A cythonized solid fundamental matrix implementation can be found in `TidalPy.RadialSolver.PropMatrix.solid_matrix`.
* New radial slicing scheme required:
  * TidalPy now requires `radial_solver` input arrays to be defined in a precise manner:
    * `radius_array` must start at 0.
    * Each layer's upper and lower radius must be in the `radius_array`. That means if there is more than one layer there will be two identical radius values!
      * E.g., if a planet has a ICB at 1000km and a CMB at 3500km. Then `radius_array` must be setup with 2 values of 1000km and 2 values of 3500km. 
      * Other parameters should be defined on a "as layer" basis. So shear modulus at the 1st 1000km would be the shear of the inner core, at the 2nd 1000km it would be the shear modulus of the outer core. Likewise shear modulus at the 1st 3500km would be for the outer core and at the 2nd 3500km would be the shear modulus for the mantle. Same goes for density and bulk modulus.
* Added warning to check for instabilities (based on large number of steps taken; requires `warnings=True`).
* Changes to `radial_solver` arguments:
  * Many changes to the order as well as additions and removals of arguments to `radial_solver` highly suggest looking through the updated documentation.
  * `radial_solver` no longer requires `gravity_array`.
    * New with this update is a self-consistent equation of state solver (EOSS) that is called from `radial_solver`. This EOSS is used to find gravity(r).
  * Bulk modulus must now be provided as a complex-valued array
    * If a non-zero imaginary portion is provided (e.g., found via a rheology) then bulk dissipation is now tracked.
    * If imaginary portions are zero, then bulk dissipation is ignored as in TidalPy v0.5.x and earlier. (this is actually dependent on your bulk rheology; it could also cause infinities...)
  * Several arguments have had slight name refactors which will break calls that used keyword arguments. Review the RadialSolver documentation or "TidalPy/RadialSolver/solver.pyx" for the new argument names.
  * `upper_radius_bylayer_array` must now be provided as a numpy array (previously a tuple of floats was acceptable).
  * Added new optional argument `surface_pressure` (default=0.0) used with EOSS to find pressure convergence.
  * Added new optional argument `core_model` (default=0) used to set the lower boundary condition when the propagation matrix technique is used.
  * Added new optional argument `starting_radius` (default=0.0) to allow the user to set the initial radius for the radial solution.
    * Setting the solution radius higher in the planet can improve solution stability when using the shooting method. Particularly if looking at higher harmonic `degree_l`s. There is a trade off with accuracy so advise testing.
    * The starting radius must be >= 0.0, if 0.0 is provided (the default) then TidalPy will use the Martens technique to find a suitable starting radius (function of `degree_l`, planet radius, and the new optional argument `start_radius_tolerance` which defaults to 1.0e-5).
  * Removed `limit_solution_to_radius` argument.
  * Added new optional argument `use_prop_matrix` (default=False) to use the propagation matrix technique over the shooting method.
  * Equation of State Solver arguments:
    * `eos_method_bylayer` - EOSS method to use for each layer (currently only "interpolate" is supported).
    * `eos_integration_method` Runge-Kutta method to use for EOSS (default="RK45"). `eos_rtol` and `eos_atol` can also be provided to control integration error.
    * `eos_pressure_tol` (default=1.0e-3) and `eos_max_iters` (default=40) control the pressure convergence of the EOSS.
  * Added optional argument `perform_checks` (default=True) that performs many checks on the user input before running the solution (small performance penalty, but highly recommend leaving on until your inputs are tested).
  * Added optional argument `log_info` (default=False) that will log key physical and diagnostic information to TidalPy's log (which can be set to be consol print, log file, or both via TidalPy's configurations). 
    * Note there is a performance hit when using this, particularly if logging to file is enabled.

**New RadialSolver Helpers**
* To assist with the generation of valid inputs to `radial_solver`, TidalPy now provides two helper functions:
  * For planets with homogeneous layers: `from TidalPy.RadialSolver import build_rs_input_homogeneous_layers` takes in attributes for a planet made of layers with constant physical properties and then provides the arrays and other required `radial_solver` inputs that conform to the new `radius_array` scheme.
  * For planets with inhomogeneous layers: `from TidalPy.RadialSolver import build_rs_input_from_data` which parses data arrays (loaded from a file or built elsewhere like using a more robust EOS than TidalPy offers) and makes changes to ensure they will work with `radial_solver`.

**Expanded RadialSolverSolution Attributes and Methods**
* The output of `radial_solver`, an instance of the `RadialSolverSolution`, has been greatly expanded to provide much more data and functionality to the user. Full details can be found in the new RadialSolver documentation. Highlights include:
  * EOSS results like `<solution>.mass`, `<solution>.moi`, `<solution>.moi_factor`, `<solution>.central_pressure`, `<solution>.surface_gravity`.
  * Diagnostic data like number of integration steps required per layer to find a solutions `<solution>.steps_taken`, or EOSS diagnostics: `<solution>.eos_iterations`, `<solution>.eos_pressure_error`, `<solution>.eos_success`, `<solution>.eos_message`, `<solution>.eos_steps_taken`.
    * Much of the new diagnostic data as well as key results can now be quickly printed using `<solution>.print_diagnostics()`.
  * Method to quickly plot the viscoelastic-gravitational solution y's `<solution>.plot_ys()`.
  * Method to quickly plot the EOS results `<solution>.plot_interior()`.
  * In addition to the previously provided attributes like `<solution>.love`, `<solution>.k`, `<solution>.h`, `<solution>.l`, `<solution>.result`.

#### Cython / C Changes
* Shifted away from `PyMem_Free` to `CyRK.utils.mem_free` to allow for consistency in future development.
  * Avoiding using these manual heap allocations whenever possible. Many new uses of smart pointers and C++ vectors.

#### Other Changes
* Updated GitHub actions.
* Moved to more consistent and robust types (e.g., int for degree_l vs. prior unsigned-char; unsigned-int).
* Added inverse function `cinv` in `TidalPy.utilities.math.complex`.
* New "TidalPy/constants.pyx" isolates all TidalPy constants. Available to both Python and Cython. Refactored all files to use the constants in this file.
* New numerics module `TidalPy.math.numerics` for low-level floating point functions.
  * New cythonized `isclose` function that matches functionality of python's `math.isclose`
* Cythonized radial sensitivity to shear/bulk functions in `TidalPy.tides.multilayer.sensitivity` (based on Tobie+2005)
* Cythonized radial heating functions that use the sensitivity to shear/bulk functions in `TidalPy.tides.multilayer.heating` (based on Tobie+2005)
* Improved logging so it is less spammy.
* Logger now logs all exceptions raised.
* Moved TidalPy's default config and world config dir to user's "Documents" folder (from system appdata folder). 
  * If upgrading from previous version of TidalPy, you can safely delete the old config directory.
    * On Windows the old dir was: "'C:\\Users\\<username>\\AppData\\Local\\TidalPy'"; The new dir is "'C:\\Users\\<username>\\Documents\\TidalPy'"
    * On Mac the old dir was: "'/Users/<username>/Library/Application Support/TidalPy'"; The new dir is "'/Users/<username>/Documents/TidalPy'"
    * On Linux the old dir was: "'/Users/<username>/.local/share/TidalPy'"; The new dir is "'/home/<username>/Documents/TidalPy'"
* New switch `TidalPy.log_to_file()` to quickly turn on saving log to file (this can also be adjusted in the TidalPy configurations).
* TidalPy now looks for an environment variable "TIDALPY_TEST_MODE" to turn on test mode during first initialization (can later be changed using the `TidalPy.test_mode()` command or setting `TidalPy._test_mode = False; TidalPy.reinit()`).
* Made use of more TidalPy-specific exceptions.
* Tweaked the `TidalPy.utilities.graphics -> yplot`.
* User can now override TidalPy.config using `TidalPy.reinit(<new config toml file path; or dictionary of configs>)`.
* Refactored and made improvements to `TidalPy.utilities.graphics.planet_plot`.

#### Dependencies
* Added support for CyRK v0.12.x
* Added support for Burnman v0.2.x

#### Documentation
* Reworked TidalPy's documentation structure in prep for a shift to Sphinx in the future.
* Greatly expanded and improved RadialSolver documentation which can be found in "TidalPy/Documentation/RadialSolver"

#### Fixes
* Fixed issue where `radial_solver` arrays could dealloc while references still pointed to them (hanging pointers).
* Missing Cython compile arguments in TidalPy's utilities, `nondimensional.pyx`.
* Fixed incorrect type in dynamic liquid layers that may have been causing some errors to propagate.
* Fixed issue where `TidalPy._config_path` was not being updated.
* Fixed issue where log files could not be created.

## Version 0.5.X

### Version 0.5.5 (2024-11-11)

#### Fixes:
* Fixed dependency compatibility issues.
* Fixed incorrect function signature type for scipy's `spherical_jn`. SciPy v.1.14.X uses a new signature which is breaking on MacOS. Limiting to "SciPy<1.14" for now. See GitHub Issue #65

### Version 0.5.4 (2024-04-30)

#### Fixes:
* Fixed RadialSolver frequency warning message for dynamic liquid layers not displaying for correct layer types.

#### Additions:
* Added way to suppress warning messages in RadialSolver.
  * To turn this suppression on, pass `warnings=False` to `RadialSolver.radial_solver`.

### Version 0.5.3 (2024-04-26)

#### Fixes:
* RadialSolver: Fixed bug where solutions between liquid and solid layers were not propagating correctly.

#### Additions:
* New Love number benchmarks for Earth provided by [Nick Wagner](https://github.com/nlwagner) in `Benchmarks & Performance\RadialSolver\Earth Love Numbers.ipynb` (Jupyter Notebook).

#### Changes:
* Pre-allocated several cythonized arrays to nans to help with debugging.
* Provided more error messages to improve user experience.
* Cythonized non-dim function now takes in the planet's density and radius as variables to change.
* Improved the Tobie and Roberts benchmarks for radial solver.

#### Other:
* Updated to work with CyRK 0.8.7

### Version 0.5.2 (2024-02-22)

#### Documentation
* Improved RadialSolver documentation regarding higher degree-l.
* Added info about issues that can arise from using non C-continguous arrays in cythonized functions.

#### Fixes:
* Added error message to `RadialSolver.radial_solver` if length of provided assumption tuples is not the same.
* Fixed issue where non C-continguous arrays were allowed in cythonized functions that required them.

### Version 0.5.1 (2024-02-14)

* Removed Python 3.8 support due to issues with building SciPy.

### Version 0.5.0 (2024-02-14)
_This version is likely to break code based on TidalPy v0.4.X and earlier_

#### Cythonizing TidalPy
* A major change starting with v0.5.0 is the switch from numba.njited functions to cython precompiled functions and extension classes. The reasons for doing this are numerous. This transition will be completed in stages with minor versions (v0.X.0) each bringing a new set of cythonized updates until all njited functions are retired.
* For this version: 
  * Converted `TidalPy.radial_solver.radial_solver` to cythonized `TidalPy.RadialSolver.radial_solver`.
    * The old radial solver method will be removed in TidalPy version 0.6.0.
  * Added new cython-based `TidalPy.utilities.classes.base_x` base cython extension class that other classes are built off of.
  * Converted `TidalPy.rheology.complex_compliances.compliance_models` to cythonized `TidalPy.rheology.models`.
    * Improved the new rheology methods to better handle extreme values of frequency and modulus.
    * The old rheology solvers will be removed in a future release of python.
  * Added several new cython-based helper functions in the utilities package.

#### Other Major Changes
* Added support for Python 3.11 and 3.12. TidalPy now runs on Python 3.8--3.12.
  * Note that currently the Burnman package does not support 3.12 so that functionality is limited to python 3.8-3.11.
* Removed support for `solver_numba` in the `radial_solver` module.
* Removed some imports from main package and sub modules' `__init__` to avoid slow load times.
* Moved conversion tools from `TidalPy.toolbox.conversions` to `TidalPy.utilities.conversions`.
* Changed setup files so that cython code can be compiled.
  * `special` - for high-performance, general, scientific functions.
* Moved TidalPy configs to a standard user directory by default. The specific location will depend on the OS.
  * Default configs will be installed on the first `import TidalPy` call after installation.
    * These defaults are stored in the `TidalPy.defaultc.py` as a string which is copy and pasted to the new `TidalPy_Configs.toml`.
  * There is a new `TidalPy.clear_data()` function to delete all data stored in these locations. Data will be rebuilt the next time TidalPy is imported.
  * New `TidalPy.set_config(config_path)` to change the active configuration file used by TidalPy.
    * Note that `TidalPy.reinit()` should be called after changing the configurations.
  * New `TidalPy.set_world_dir(world_dir_path)` to change which directory to pull world configs from. 
  * Moved away from the system of `default.py` configurations for sub modules. All default configs are stored in the same `TidalPy_Config.toml`
* Shifted from `json` to `toml` files for world configs.
  * Store all world configs to a zip file for easier distribution.
* TidalPy now requires:
  * CyRK>=0.8.6
  * Cython>=3.0.0
* Moved `BurnMan` 3rd party dependence to a more dedicated `Extending` folder for future development.
* To make TidalPy lighter weight we are starting to remove a lot of 3rd party packages.

#### Minor Changes and New Features
* `complex_compliance` configurations are now stored in the top level `rheology` in all configs.
  * For example, in prior TidalPy versions you would need to change the complex compliance model by editing `config['layers']['mantle']['rheology']['complex_compliance']['model'] = 'andrade'`. Now this would be: `config['layers']['mantle']['rheology']['model'] = 'andrade'`.
* Added unique frequency finder functions to the `modes` module.
* Moved most of the type hints behind the `typing.TYPE_CHECKING` flag.
* Moved non-critical files out of repository.
* Created a new `tides.heating` module and moved the volumetric heating calculations there.
* Expanded the performance suite to better track the `radial_solver` module.
* Moved `cache.py` to top-level.
* Turned off numba cacheing on several functions that may be used in the radial solver.
  * rheology
    * complex compliance functions
  * radial_solver.numerical
    * initial guess functions
    * interface functions
* Converted radial_solver.numerical initial guess and interface functions output to np.ndarrays rather than numba lists.
* Removed `config_helper.py` and the functions defined within.
* New RadialSolver class now supports more than just boolean inputs.
  * Future proofing to allow for a greater variety of layer types.
* Added exoplanet archive download functionality in `TidalPy.utilities.exoplanets`.
  
#### Bug Fixes
* Fixed floating point comparison bug in `multilayer_modes` solver.
* Fixed obliquity not being used issue in quick tides calculator.
* Fixed issue in incorrect TidalPy version being loaded into the package.

#### Performance Improvements
* Improved the performance of the stress and strain calculator by ~20%.
* Cythonize Performance Increases:
  * New `RadialSolver.radial_solver` leads to a ~50x performance boost.
  * New cythonized rheology models are 500% faster for arrays; 40,000% faster for scalars (not a typo!)

## Version 0.4.X

### Version 0.4.1 Alpha (Spring 2023)
#### Major Changes
* Moved `radial_solver` to a top-level module of TidalPy.
  * Added `find_love` function to the `radial_solver` module for the calculation of Love and Shida numbers.
  * Added `sensitivity_to_shear` function to the `radial_solver` module based on Tobie et al. (2005).
  * Added `sensitivity_to_bulk` function to the `radial_solver` module based on Tobie et al. (2005).
* Added `newton` and `elastic` complex compliances to rheology module for ALMA comparisons.

#### Minor Changes
* Updated and added to `radial_solver` test suite.
* Changed `radial_solver.interfaces` functions to only require the top most value of the radial solutions rather than the whole array.

#### Bug Fixes
* Fixed bug in `radial_solver` where interfaces between a dynamic and static liquid were not correctly handled.

### Version 0.4.0 Alpha (Winter 2022/2023)

#### Major Changes 
* Added a multilayer tidal potential that allows for arbitrary obliquity.
* Added in load Love number calculations to the multilayer code.
* Removed a lot of 3rd-party dependencies to make TidalPy's install more lean.
* Switched over to using the integrators from the new `CyRK` package
  * Changed the signature of the numerical-int multilayer solver.
    * **Breaks Old Code**
* Issue with Numba 0.55 and dictionary updates. This restricts TidalPy to Python version 3.9 or lower.
* Started a lot of prep work for a move to sphinx or similar for documentation (this will be a 0.5.0 feature).
* Created a generalized multi-layer collapse function in `TidalPy.tides.multilayer.numerical_int.collapse`.
  * This will likely break old code. It was introduced in v0.4.0.dev10.
  * Arbitrary layer structures can now be fed into the y-solver. e.g., liquid-solid, solid-liquid, liquid-solid-liquid-solid, etc.
    * Note: Multiple dynamic liquid layers will likely lead to stability problems.
  * y-solver no longer requires a `model_name` variable. The function will automatically utilize the correct model based on the `is_layer_solid` and `is_layer_static` lists.
  * Removed the old code that handled individual cases.
* Interfaces between multiple liquid layers are now supported (but not well tested).
  * Two liquid layers of different types can be next to one another (dynamic-static; static-dynamic).
* y-solvers functional argument signature has changed (both regular and numba versions)
  * **Breaks Old Code**
  * Numerical integrator declaration has changed.
  * Planet_bulk_density is now a required argument.
  * Modified several other TidalPy functions to match this new call signature for the y-solver.
* Removed the y-derivative return from the propagation matrix (to match the output format of the numerical int)
  * **Breaks Old Code**
* Refactored `tidal_y_solver` to `radial_solver` since non-tidal calculations can be made with it.
  * **Breaks Old Code (pre v0.4.0.dev11)**
* Switched from `setup.py` to a streamlined `pyproject.toml` installation process.
* Changes to radial ODE's
  * Input arguments and output diffeqs are now passed as numpy arrays rather than tuples.
  * Input and outputs are now passed as floats not complex (doubling the number of terms)
* Added `numba-scipy` dependence to allow the use of scipy's special functions. 
  * Removed the pre-calculated factorial method. Using scipy's gamma now.
  * TODO: Note the numba-scipy package on github is not updated to the newest version of scipy. Packaging numba-scipy with TidalPy for now.

#### Performance Improvements
* Improved the performance of the pure-numba radial solve by ~10%

#### Minor Changes
* Added newer functions to the performance recording suite.
* Improved performance of eccentricity and inclination functions for non-multilayer tidal calculations.
* Added TidalPy to Zenodo. DOI added to readme.
* Removed Gitter account for now.
* Added files and functions to quickly install additional 3rd party applications
* Improvements to GitHub workflows
* Switched over to using the 3rd party `cmcrameri` colormap package rather than trying to maintain it within tidalpy
* Added a delta time to HH:MM:SS converter to utilities.string_helper
* Made improvements to multiprocessing user info
* Cleaned up & added some docstrings and type hints
* Tidal y solver for the prop matrix method no longer returns the y-derivatives.
  * dy1/dr is now calculated directly in the `decompose()` function.
* Created a config helper function `TidalPy.test_mode()` to quickly setup TidalPy configs for pytest'ing
* Cleaned up comments and reordered items in `multilayer.numerical_int.collapse`.
* Fixed a bug when using SciPy's Radau integrator method.
* The version number is now checked with importlib in TidalPy.__init__. Version number should only be changed in the pyproject.toml.
* Updated the pure-numba version of the radial solver's argument signature to better match the python implementation.
* Updated Github Actions

#### Bug Fixes
* Fixed bug in GridPlot related to number of subplots.
* Fixed bug in global variable for world config loader.
* Fixed type hint bug in the numba-based tidal y solver.
* Fixed bug that caused numba-based tidal y solver to not compile.
* Fixed bug that was causing full TidalPy log to print while in a Jupyter Notebook environment.
  * If you would like the log to print in a notebook then use `TidalPy.toggle_log_print_in_jupyter()` or set the
`print_log_in_jupyter` to `True` in the "configurations.py" file.
* Fixed an error in the multimode volumetric heating calculation where "_rr", "_thth", "_phiphi" were being double counted
* Fixed an error in the stress/strain calculations where the static, instead of complex, shear was being used.

## Version 0.3.X

### Version 0.3.5 Alpha (Spring 2022)

#### Minor Changes
* Removed the soon-to-be deprecated np.float, np.complex, np.int references.
* Added a check to see if cartopy is installed before functions that depend on it are imported. 

#### Bug Fixes
* Fixed coverage problem with github actions.
* Fixed issue importing nbTuple.
* Fixed issues with GitHub actions and coverage.

### Version 0.3.4 Alpha (Winter/Spring 2022)

*Multilayer scripts based on 0.3.3 or earlier will likely break with this version!*

#### Major Changes
* Added `TidalPy.modes.multilayer_modes.py` module to offer simplified calculation of multilayer tidal
  heating.
* Added `GridPlot` class to quickly make grid-like matplotlib figures.
  Checkout `TidalPy.utilities.graphics.grid_plot.py`
* Added `Cartopy` dependence.
    * Can now make cool projection maps! Added basic functionality to `TidalPy.utilities.graphics.global_map.py`.
    * New jupyter notebooks to showcase map projects and GridPlot functionality.
* Improved performance on both mode and non-mode tidal potential functions by at least a factor of 3. If used
  correctly these can be nearly 100x faster.
* Added a new obliquity version of the mode version tidal potential.
* Stress and strain relationship for multi-layer tides now allows for arbitrary rheology.
* Created a single multilayer solver to handle an arbitrary layer structure.
  See `TidalPy.tides.multilayer.numerical_int.solver.py`
* Stress & Strain relationship now accounts for arbitrary rheology.
* Created a single multi-mode solver for multilayer problems. See `TidalPy.tides.modes.multilayer_modes.py`
* TidalPy now defaults to using the frequency dependent zeta versions of Andrade and Sundberg rheologies.
    * This was done to avoid issues with real(complex_comp) at zero frequency which happens in multi mode
      calculations.
* Added numba-safe version of multilayer calc

#### Minor Changes
* Added a helper function to quickly calculate masses, volumes, and gravity for spherical shells provided a radius
  and density array: `TidalPy.utilities.spherical_helper.calculate_mass_gravity_arrays`.
* Added more colormaps, updated how reserved versions are constructed, updated old maps.
* Better support for post-multiprocessing function inputs.
* Refactored tidal mode calculation functions into a new TidalPy.tides.modes module.
* Added a voxel calculator to TidalPy.utilities.spherical_helper. Also, added related tests & benchmarking tools.
* Added a dictionary to track known color maps. It can be imported at `TidalPy.utilities.cmaps.KNOWN_CMAPS`
* Made some improvements to the unique_path function in io_helper.py
* Rearranged the tidal potential argument order.
* Added some sanity checks on the various kinds on both mode and non-mode tidal potentials to compare with one
  another.
* Updated stress, strain, and displacement calculations to account for new low-memory calculation method.
* Greatly increased performance of multilayer stress, strain, and potential calculations.
* Refactored much of the multilayer functions from TidalPy.toolbox to TidalPy.tides.multilayer.numerical_int and sub
  modules
* Reworked numba-safe RK integrator. Does not reproduce scipy exactly for chaotic functions when numba is on.

#### Bug Fixes
* Fixed issue that was causing some multilayer boundary calculations to be 10x slower.
* Fixed issue where 2D arrays would not work on njited `neg_array_for_log_plot` function.
* Fixed issue where there could be negative love numbers with frequency is zero in multi-mode calculator

### Version 0.3.3 Alpha (Fall 2021)

#### Major Changes
* Updated to work with BurnMan v1.0 release.
* Updated the NSR tidal potential for multilayer code to a multi-modal version
* Improvements to the multilayer module

#### Minor Changes
* Added more tests
* Removed the remaining `_array` functions

### Version 0.3.2 Alpha (Summer 2021)

#### Major Changes
* Changes to setup to use released version of BurnMan.

### Version 0.3.1 Alpha (Summer 2021)

#### Major Changes
* Added a simplistic (not modal version) NSR tidal potential for use in the multilayer code.
* Added multiprocessor toolkit to `TidalPy.utilities.multiprocessing`

#### Minor Changes
* Did large scale code reformatting (just style, no refactoring) on nearly all files.
* Cleaned up some doc strings.

### Version 0.3.0 Alpha (Spring 2021)

*Scripts based on 0.2.x will likely break with this version!*

#### Major Changes:
* Added the first iteration of a multilayer tidal calculator module in `TidalPy.tides.multilayer` this module
  provides basic functionality to calculate tidal dissipation in a semi-homogeneous, shell-based approach. This is
  more accurate than the pure homogeneous model used throughout the rest of TidalPy. The downside with the current
  version is that it does not allow for NSR or high eccentricity / obliquity. A future version will attempt to add
  in a more robust Tidal Potential equation which will allow for additional physics.
* Setup.py has been revamped as has the installation process. This is in prep to allow for TidalPy to become
  available on PyPI.
* Did away with most of the `_array` functions. Found a way for njit to compile a function to handle either arrays
  or floats.
    * Left the `self._func_array` (in addition to `self._func`) in the `model.py` classes just in case we ever **
      do** need to define array functions in the future: all the infrastructure is still in place.
* Added a numba-safe Explicit Runge-Kutta integrator. This is fully wrapped in njit'd functions.
    * On its own this can be 5--20 times faster than `scipy.solve_ivp`.
    * This also allows the integration function to be used from within another njit'd function(s).
    * This is still very experimental so please use caution and testing when using it.

#### Minor Changes
* Numerous bug fixes.
* Removed the array versions of the dynamic functions.
    * `use_array` is still tracked in OOP and some quick calculation functions. These may all be not necessary now.
* New cookbook to showcase the multilayer calculations.
* Added surface area slices to base physical class.
* Fixed some issues with how radius slices are tracked within layers and worlds.
    * `<world/layer>.radii[0]` is never the radius at the bottom of the object, it is always one dx up.
    * An interpolation had to be added to Burnman layers since Burnman radii starts at 0.
* Improved various docstrings.
* Refactored the `TidalPy.tools` to `TidalPy.toolbox`.
* Refactored `Cookbooks` to `Demos`.
* conversions.semi_a2orbital_motion and orbital_motion2semi_a now always return np.nan where they used to return
  complex numbers.

## Version 0.2.X

### Version 0.2.1 Alpha (Fall 2020)

*Will break studies based on previous versions*

Note: TidalPy version of "0.2.0" was never made publicly available. This is the next version after 0.1.0.

#### Major Changes
* Major refactoring all over.
* Replaced the custom logging system with one based on Python's logging package.
    * Logger will not print to console if using a Jupyter notebook.
    * User can decide if the log is saved to disk or not.
* Modified OOP Backend
    * Models are now free to access the state properties of layers and planets. This eliminates the need to pass
      arguments into model classes for calculations. All the user has to do is change the layer and/or planet's
      state properties and then those changes will automatically propagate.
* New Rheology Scheme
    * To better follow real physics, all strength-based models have been moved under the redesigned `Rheology`
      class.
    * This includes: liquid and solid viscosity calculations, partial melting, and complex compliance
    * Love numbers are now calculated by a new `Tides` class instead of the `Rheology` (see next bullet point).
* New Tides Module
    * New `Tides` class and child classes have been implemented which handle all Love number and tidal calculations
      for both a rheology-based approach or for the CPL/CTL model.
    * This module contains functionality to calculate tidal heating and tidal potential derivatives based on a
      complex compliance function and the various thermal and orbital parameters.
    * The functions to calculate the Love numbers are also stored here.
    * New eccentricity functions have been added including terms up to and including e^20 and tidal order l=7.
    * New inclination functions have been added (with arbitrary accuracy) up to tidal order l=7.
* Other
    * Various functions and classes have been reworked or removed to simplify the code base.

#### Minor Changes
* Removed ".pyname" from classes.
* Utilities package has been reorganized.
* Changed the setup pipeline to only require one command.

#### QOL Improvements
* Many new docstrings, type hints, and overall clean up of functions and classes. All docstrings should now follow
  the numpy format.
* Many new tests for both the functional and OOP versions of TidalPy.
* setup.py no longer requires a separate command line call to install the Burnman package.
* More comments and spelling/typo fixes everywhere.
* Added CVD-friendly color maps made by Crameri (2018; http://doi.org/10.5281/zenodo.1243862) to the
  utilities.graphics module
* More log.debug() calls all over. This should hopefully help bugfixes in the future (especially OOP bugs).

#### Bug Fixes:
* Fixed bug in world_builder.py : build_from_world where nested dicts were not being overwritten as expected.
* Fixed bug in melt fraction checks for the float version of the Henning model.
* Fixed bug in the Henning melting model where viscosity and shear was not being calculated correctly during the
  breakdown band (critical melt fraction + ~5%). This only affected a phase space that was rarely important for
  tidal calculations (planets would pass through it *very* quickly).

## Version 0.1.X

### Version 0.1.0 Alpha (July 2019)

*Initial Release - Changes not tracked*
