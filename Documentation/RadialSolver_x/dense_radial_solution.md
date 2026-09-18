# Dense Radial Solutions

_Updated: 2026-09-16_

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
eos = rs.eos_call_si(300.0)  # length-11 float64 array at r = 300 m
#   [0] gravity [1] pressure [2] mass [3] moment of inertia [4] density
#   [5,6] complex shear (re, im) [7,8] complex bulk (re, im) [9,10] shear/bulk viscosity
shear = complex(eos[5], eos[6])

# Gridded array API (the standalone solver samples the dense interpolants back onto the EOS grid):
rs.result    # gridded y-solution
rs.love      # k, h, l
```

`eos_call_si(radius)` is the dense analogue of `get_radial_solution` for the structural and material state: it maps the SI radius into the solver's non-dimensional domain and evaluates the solution's own dense EOS interpolant, so an on-radius query (for example the complex moduli the 3D tides kernel needs) uses the same dense evaluation the solver uses internally rather than a separate re-interpolation of the gridded modulus arrays. The older `eos_call(radius)` accepts only a raw non-dimensional radius and is kept for internal use. The per-layer EOS interpolation inputs are persisted in the solution storage in non-dimensional solve units, so the dense evaluation stays valid after the standalone solve returns, and `c_EOSSolution::call` re-dimensionalizes the result to SI.

At the C++ level the same is available on `c_RadialSolutionStorage` (`get_radial_solution`, `get_radial_solution_array`, `get_surface_y`, `get_eos_si`) and, for a built world, on `c_LayeredWorld` (`get_radial_solution_y`, `get_love_surface_y`).

## Validation

`Tests/Test_RadialSolver_x/test_dense_benchmark/` holds frozen `.npz` reference results for a 1-layer solid planet, a 2-layer solid/solid planet, and a solid / dynamic-liquid / solid planet at a short forcing period, and checks that the dense system reproduces them. The dense and reference solutions agree closely at the surface and at layer boundaries and the Love numbers agree tightly; between grid slices the dense evaluation is the more accurate one.

## Dynamic Liquid Layers at Long Forcing Periods

A dynamic liquid layer carries inertial (`1/omega^2`) terms that are only significant at short forcing periods. When a dynamic liquid layer is sandwiched between solid layers and forced at a long period (low frequency), those terms make the layer's independent solutions grow exponentially through the liquid, so by the surface they are nearly linearly dependent and the surface boundary-condition matrix becomes near-singular. The solve is then unstable: small numerical differences (integration tolerance, grid-vs-dense sampling, the linear-solver implementation) change the result, and at sufficiently long periods the solve fails outright. This is inherent to the dynamic-liquid assumption.

Guidance:

* Use a dynamic liquid layer only for short-period forcing, where it is well-conditioned.
* Use a static liquid layer for long-period forcing; it is stable and consistent across all periods.
