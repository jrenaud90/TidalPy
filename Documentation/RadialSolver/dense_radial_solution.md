# Dense Radial Solutions

_Updated: 2026-10-06_

`TidalPy.RadialSolver` computes the viscoelastic-gravitational radial functions `y1..y6` with either a shooting method or a propagation matrix, and from them the one-dimensional tidal or loading Love numbers. This page describes how the shooting method retains its solution which can be reevaluated at any radius.

## Shooting Method

For the shooting method the planet is integrated layer by layer from a starting radius up to the surface. Each layer contributes a fixed number of independent solutions depending on its assumptions:

| Layer type            | Independent solutions |
|-----------------------|-----------------------|
| Solid                 | 3                     |
| Dynamic liquid        | 2                     |
| Static liquid         | 1                     |

The physical solution in each layer is a linear combination of that layer's independent solutions. The combination coefficients ("constants of integration") are fixed by:

1. Applying the surface boundary condition (a small matrix solve at the planet surface), then
2. Propagating those constants downward through every interface, one layer at a time, so each layer below inherits a consistent set of constants. The upward pass records each interface's map from the lower layer's solutions to the upper layer's (`c_solve_upper_y_at_interface` in `interfaces_.hpp`), and the collapse applies it in reverse (`c_collapse_through_interface`), so the two directions share one set of interface conditions.

## Dense Interpolants

Each independent solution is integrated with dense (continuous) output and the whole interpolant is retained, together with the per-layer collapse constants. The solution at any radius is then produced on demand:

1. Locate the layer containing the requested radius,
2. Evaluate that layer's 1-3 dense interpolants at the radius,
3. Multiply by the stored constants and sum (the collapse),
4. Reconstruct `y3` for dynamic-liquid layers from the other components and the EOS gravity/density,
5. Re-dimensionalize from the solver's internal (non-dimensional) units to SI.

The dense interpolant is the integrated solution itself, so a value between grid slices comes from the integrator's own polynomial rather than from a linear interpolation of a sampled grid. The radial solver therefore does not own a solution grid for the shooting method; the structural arrays (radius, gravity, density, moduli) remain owned by the EOS class and are queried as needed.

The surface solution (and the Love numbers derived from it) is always required, so it is computed once per solve and cached for both methods.

## Propagation Matrix Method

The propagation-matrix method builds its solution on a grid by propagating the fundamental matrix of each slice's material from one slice radius to the next. It exposes the same `get_radial_solution(r)` entry point, and that call continues the propagation to the radius asked for: inside a slice the solution is the slice's fundamental matrix at `r` applied to the coefficients the propagation carried into that slice, which reproduces the grid values at both ends of the slice and is exact between them for a uniform layer. Below the first propagated slice the regular solution of the innermost material continues the solution to the center, so `get_radial_solution` and the `result` grid are defined down to `r = 0`. With a `core_model` other than `0` the seed is not the regular solution and the solution stays NaN below the starting radius.

## Python API

```python
rs = radial_solver(...)            # RadialSolverSolution

# Arbitrary-radius evaluation (SI complex y1..y6):
y  = rs.get_radial_solution(300.0)           # length-6 complex128 array at r = 300 m, ytype 0
ys = rs.get_radial_solution_array(r_array)   # (N, 6) complex128, evaluated in vectorized C++

# Arbitrary-radius EOS / material state (SI), via the same dense interpolant:
eos = rs.eos_call(300.0)             # dict of scalars at r = 300 m
eos["density"], eos["pressure"]      # gravity, pressure, mass, moi, density, shear_modulus, bulk_modulus,
eos["complex_shear_modulus"]         # shear_viscosity, bulk_viscosity, temperature, heat_flow, melt_fraction,
profile = rs.eos_call(r_array)       # complex_shear_modulus, complex_bulk_modulus; arrays in, arrays out

# Gridded array API (the standalone solver samples the dense interpolants back onto the EOS grid):
rs.result    # gridded y-solution
rs.love      # k, h, l
```

`eos_call(radius)` is the dense analogue of `get_radial_solution`, mapping the SI radius into the solver's non-dimensional domain and evaluates the solution's own dense EOS interpolant. So a radius query uses the same dense evaluation the solver uses internally rather than a separate re-interpolation of the gridded modulus arrays. It answers with a dict of named fields, scalars for a float radius and arrays of the input's shape for an array of radii, NaN outside the body. Its fields, in layout order, are listed by `EOS_CALL_FIELDS`. The per-layer EOS interpolation inputs are persisted in the solution storage in non-dimensional solve units, so the dense evaluation stays valid after the standalone solve returns, and `c_EOSSolution::call_nondim` re-dimensionalizes the result to SI.

Everything in that layout is frequency-independent, so its `shear_modulus` and `bulk_modulus` are the unrelaxed ones. A viscoelastic response is a property of a rheology at a forcing frequency rather than of the equation of state, so the `complex_shear_modulus` and `complex_bulk_modulus` fields come from a second read, `get_complex_shear_modulus(radius)` and `get_complex_bulk_modulus(radius)` reproduce the moduli the solve actually used by calling a world-attached rheology (shared with the solution, so it can outlive the world if applicable) at the frequency it recorded in `love_frequency`, and a solve handed its moduli as arrays interpolates those arrays. At the C++ level the solver reads the same values in one call through `c_EOSSolution::call_material`, which returns a `c_EOSMaterialState` (gravity, density, and both complex moduli) whatever the solution was built from.

At the C++ level the same is available on `c_RadialSolutionStorage` (`get_radial_solution`, `get_radial_solution_array`, `get_surface_y`, `get_eos_si`, with `get_radial_solution_nondim` for a radius in solve units) and, for a built world, on `c_BaseWorld` (`get_radial_solution_y`, `get_love_surface_y`).

## Validation

`Tests/Test_RadialSolver/test_dense_benchmark/` holds frozen `.npz` reference results for a 1-layer solid planet, a 2-layer solid/solid planet, and a solid / dynamic-liquid / solid planet at a short forcing period, and checks that the dense system reproduces them. The dense and reference solutions agree closely at the surface and at layer boundaries and the Love numbers agree tightly; between grid slices the dense evaluation is the more accurate one.

## Dynamic Liquid Layers at Long Forcing Periods

The dynamic liquid equations carry no density-gradient term, so a liquid layer's stratification comes from how its density, gravity, and bulk modulus vary with radius. A liquid whose density follows its bulk modulus, $d\rho/dr = -\rho^2 g / K$, is neutral ($N^2 = 0$). A constant-density liquid with a finite bulk modulus is not: $N^2 = -\rho g^2 / K < 0$, which is unstable stratification. Its solutions then grow as $e^E$ through the layer, with

$$E = \frac{\sqrt{\ell(\ell+1)}}{\omega} \int \frac{\sqrt{-N^2}}{r}\, dr,$$

which rises with the forcing period. Once $e^E$ amplifies the integration tolerance to order one, the independent solutions are nearly dependent by the surface and the solve is wrong or fails. Every bundled liquid except PREM's is a constant-density liquid: solved dynamically, Earth-Simple's core breaks down by half a day, Mercury's by 3.5 days, the Moon's and PREM's by 10 days, and Pluto's ocean by 30 days.

Before every radial Love solve (`solve_love_numbers`, `calc_tides`, the 3D calls, and `radial_solver`), TidalPy estimates the error each dynamic liquid layer brings, $\mathrm{rtol}\, e^E$, from the layer's profile, and logs one warning per world when it passes 1%. It names the layer, the period, and $E$.

Guidance:

* Use a static liquid layer for long-period forcing. It is stable and accurate at every period for any density profile.
* To keep a constant-density liquid dynamic, make it incompressible (`is_incompressible = True`). An incompressible liquid of constant density is neutral: in the bundled worlds its dynamic solve stays within 1e-3 of the static one out to 100 days, except Earth-Simple's large core, which drifts to 6% there.
* To keep it dynamic and compressible, give it a pressure-dependent law (Birch-Murnaghan or Vinet), which keeps it neutral. The `luna_dynamic`, `mercury_dynamic`, `pluto_dynamic`, and `europa_dynamic` bundled worlds are built this way. A neutral compressible liquid still loses some accuracy at long periods in a large core (`luna_dynamic`'s yearly $k_2$ is 3e-4 off at the default `rtol = atol = 3e-8`, 4e-5 at `3e-9`, and was 0.8% off at `rtol = 1e-6`), since the equations recover its tangential displacement through a division by $\omega^2$. Thin oceans are unaffected.
* The dynamic liquids demo (`Demos/Physics/P12_dynamic_liquids.ipynb`) works through each case.
