# Dense Radial Solutions

_Updated: 2026-09-21_

`TidalPy.RadialSolver_x` computes the viscoelastic-gravitational radial functions `y1..y6` with either a shooting method or a propagation matrix, and from them the one-dimensional tidal or loading Love numbers. This page describes how the shooting method retains its solution and evaluates it at any radius.

## Shooting Method

For the shooting method the planet is integrated layer by layer from a starting radius up to the surface. Each layer contributes a fixed number of independent solutions depending on its assumptions:

| Layer type            | Independent solutions |
|-----------------------|-----------------------|
| Solid                 | 3                     |
| Dynamic liquid        | 2                     |
| Static liquid         | 1                     |

The physical solution in each layer is a linear combination of that layer's independent solutions. The combination coefficients ("constants of integration") are fixed by:

1. Applying the surface boundary condition (a small matrix solve at the planet surface), then
2. Propagating those constants downward through every interface (`reversed_.hpp`), one layer at a time, so each layer below inherits a consistent set of constants.

## Dense Interpolants

Each independent solution is integrated with dense (continuous) output and the whole interpolant is retained, together with the per-layer collapse constants. The solution at any radius is then produced on demand:

1. Locate the layer containing the requested radius,
2. Evaluate that layer's 1-3 dense interpolants at the radius,
3. Multiply by the stored constants and sum (the collapse),
4. Reconstruct `y3` for dynamic-liquid layers from the other components and the EOS gravity/density,
5. Re-dimensionalize from the solver's internal (non-dimensional) units to SI.

The dense interpolant is the integrated solution itself, so a value between grid slices comes from the integrator's own polynomial rather than from a linear interpolation of a sampled grid. The radial solver therefore does not own a solution grid for the shooting method; the structural arrays (radius, gravity, density, moduli) remain owned by the EOS class and are queried as needed.

The surface solution (and the Love numbers derived from it) is always required, so it is computed once per solve and cached for both methods. Interior radii are evaluated on demand.

## Propagation Matrix Method

The propagation-matrix method does not use a shooting integration and has no dense interpolant: it builds its solution on the provided grid and can only represent it there. For consistency it exposes the same `get_radial_solution(r)` entry point, but for the matrix method that call performs a linear interpolation of its constructed grid. The two methods therefore differ: the shooting method serves a dense evaluation, the matrix method a linear interpolation of its grid.

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

`eos_call(radius)` is the dense analogue of `get_radial_solution` for the structural and material state: it maps the SI radius into the solver's non-dimensional domain and evaluates the solution's own dense EOS interpolant, so an on-radius query uses the same dense evaluation the solver uses internally rather than a separate re-interpolation of the gridded modulus arrays. It answers with a dict of named fields, scalars for a float radius and arrays of the input's shape for an array of radii, NaN outside the body. `eos_call_nondim(radius)` takes a raw radius in the solve's own units and returns the bare row (its slot order is `EOS_CALL_FIELDS`); it is kept for internal use. The per-layer EOS interpolation inputs are persisted in the solution storage in non-dimensional solve units, so the dense evaluation stays valid after the standalone solve returns, and `c_EOSSolution::call_nondim` re-dimensionalizes the result to SI.

Everything in that layout is frequency-independent, so its `shear_modulus` and `bulk_modulus` are the unrelaxed ones. A viscoelastic response is a property of a rheology at a forcing frequency rather than of the equation of state, so the `complex_shear_modulus` and `complex_bulk_modulus` fields come from a second read: `get_complex_shear_modulus(radius)` and `get_complex_bulk_modulus(radius)` reproduce the moduli the solve actually used by calling a world-attached rheology (shared with the solution, so it can outlive the world if applicable) at the frequency it recorded in `love_frequency`, and a solve handed its moduli as arrays interpolates those arrays. At the C++ level the solver reads the same values in one call through `c_EOSSolution::call_material`, which returns a `c_EOSMaterialState` (gravity, density, and both complex moduli) whatever the solution was built from.

At the C++ level the same is available on `c_RadialSolutionStorage` (`get_radial_solution`, `get_radial_solution_array`, `get_surface_y`, `get_eos_si`, with `get_radial_solution_nondim` for a radius in solve units) and, for a built world, on `c_LayeredWorld` (`get_radial_solution_y`, `get_love_surface_y`).

## Validation

`Tests/Test_RadialSolver_x/test_dense_benchmark/` holds frozen `.npz` reference results for a 1-layer solid planet, a 2-layer solid/solid planet, and a solid / dynamic-liquid / solid planet at a short forcing period, and checks that the dense system reproduces them. The dense and reference solutions agree closely at the surface and at layer boundaries and the Love numbers agree tightly; between grid slices the dense evaluation is the more accurate one.

## Dynamic Liquid Layers at Long Forcing Periods

A dynamic liquid layer carries inertial ($1/\omega^2$) terms that are only significant at short forcing periods. When a dynamic liquid layer is sandwiched between solid layers and forced at a long period (low frequency), those terms make the layer's independent solutions grow exponentially through the liquid, so by the surface they are nearly linearly dependent and the surface boundary-condition matrix becomes near-singular. The solve is then unstable: small numerical differences (integration tolerance, grid-vs-dense sampling, the linear-solver implementation) change the result, and at sufficiently long periods the solve fails outright. This is inherent to the dynamic-liquid assumption.

Guidance:

* Use a dynamic liquid layer only for short-period forcing, where it is well-conditioned.
* Use a static liquid layer for long-period forcing; it is stable and consistent across all periods.
