# RadialSolver_x: the dense radial-solution calling system

This document explains how `TidalPy.RadialSolver_x` computes the viscoelastic-gravitational radial
solution (`y1..y6`) at a given radius, and how that differs from previous versions of the radial
solver.

## Background: how the shooting method builds a solution

For the shooting method the planet is integrated layer by layer from a starting radius up to the
surface. Each layer contributes a fixed number of independent solutions depending on its assumptions:

| Layer type            | Independent solutions |
|-----------------------|-----------------------|
| Solid                 | 3                     |
| Dynamic liquid        | 2                     |
| Static liquid         | 1                     |

The physical solution in each layer is a linear combination of that layer's independent solutions.
The combination coefficients ("constants of integration") are fixed by:

1. applying the surface boundary condition (a small matrix solve at the planet surface), then
2. propagating those constants downward through every interface (`reversed_.hpp`), one layer at a
   time, so each layer below inherits a consistent set of constants.

## What changed: dense interpolants instead of a fixed grid

**Previous behavior.** Each independent solution was integrated and **sampled onto a fixed radial
grid** (the EOS radius array). The constants were then applied at every grid slice to build a single
gridded `y1..y6` solution. To obtain the solution at a radius that was *not* a grid point - for
example, a stress evaluation at `r = 300 m` when the grid only stored `r = 100 m` and `r = 500 m` -
the gridded solution was **linearly interpolated**. For coarse grids this interpolation could be
noticeably inaccurate.

**Current behavior.** Each independent solution is integrated with **dense (continuous) output** and
the whole interpolant is retained. The per-layer collapse constants are stored alongside. The
solution at *any* radius is then produced on demand:

1. locate the layer containing the requested radius,
2. evaluate that layer's 1-3 dense interpolants at the radius,
3. multiply by the stored constants and sum (the collapse),
4. reconstruct `y3` for dynamic-liquid layers from the other components and the EOS gravity/density,
5. re-dimensionalize from the solver's internal (non-dimensional) units to SI.

Improved interpolation *between* grid slices is involved, the dense interpolant is the integrated
solution with complex interpolation (not linear). The radial solver therefore no
longer needs to own a solution grid for the shooting method; the structural arrays (radius, gravity,
density, moduli) remain owned by the EOS class and are queried as needed.

The **surface** solution (and the Love numbers derived from it) is always required, so it is computed
once per solve and cached for both methods. Interior radii are evaluated on demand.

## The propagation-matrix method is different by construction

The propagation-matrix method does not use a shooting integration and has no dense interpolant: by
definition it builds its solution on the **provided grid** and can only represent it there. For
consistency it exposes the same `get_radial_solution(r)` entry point, but for the matrix method that
call performs a **linear interpolation of its constructed grid** (the same, less-accurate behavior
the shooting method used to have). So the two methods deliberately differ: the shooting method serves
an accurate dense evaluation, the matrix method serves a linear interpolation of its grid.

## Python API

```python
rs = radial_solver(...)            # RadialSolverSolution

# Arbitrary-radius evaluation (SI complex y1..y6):
y = rs.get_radial_solution(300.0)              # length-6 complex128 array at r = 300 m, ytype 0
ys = rs.get_radial_solution_array(r_array)     # (N, 6) complex128, evaluated in vectorized C++

# Arbitrary-radius EOS / material state (SI), via the same dense interpolant:
eos = rs.eos_call_si(300.0)                    # length-11 float64 array at r = 300 m
#   [0] gravity [1] pressure [2] mass [3] moment of inertia [4] density
#   [5,6] complex shear (re, im) [7,8] complex bulk (re, im) [9,10] shear/bulk viscosity
shear = complex(eos[5], eos[6])

# Legacy array API is unchanged (the standalone solver samples the dense interpolants back onto the
# EOS grid so .result / .love behave exactly as before):
rs.result    # gridded y-solution
rs.love      # k, h, l
```

`eos_call_si(radius)` is the dense analogue of `get_radial_solution` for the structural / material
state: it maps the SI radius into the solver's non-dimensional domain and evaluates the solution's own
dense EOS interpolant, so on-radius queries (e.g. the complex moduli the 3D-tides kernel needs) use the
same dense evaluation the solver uses internally rather than a separate re-interpolation of the gridded
modulus arrays. (The older `eos_call(radius)` accepts only a raw non-dimensional radius and is kept for
internal use.) For this to stay valid after the standalone solve returns, the per-layer EOS
interpolation inputs are persisted in the solution storage in non-dimensional solve units; the dense
re-invoke of the EOS extra outputs then reads live, correctly-scaled data and `c_EOSSolution::call`
re-dimensionalizes it back to SI.

At the C++ level the same is available on `c_RadialSolutionStorage`
(`get_radial_solution` / `get_radial_solution_array` / `get_surface_y` / `get_eos_si`) and, for a built
world, on `c_LayeredWorld` (`get_radial_solution_y`, `get_love_surface_y`).

## Validation

`Tests/Test_RadialSolver_x/test_dense_benchmark/` freezes results from the original `RadialSolver`
(to `.npz`, since that solver will eventually be removed) and checks that the dense system reproduces
them. The benchmark covers a 1-layer solid planet, a 2-layer solid/solid planet, and a solid /
dynamic-liquid / solid planet at a short forcing period. The dense and original solutions agree very
well at the surface and at layer boundaries and the Love numbers agree tightly; between grid slices
the dense evaluation is the more accurate one.

## Numerical limitation: dynamic liquid layers at long forcing periods

A **dynamic** liquid layer carries inertial (`1/omega^2`) terms that are only significant at short
forcing periods. When a dynamic liquid layer is sandwiched between solid layers and forced at a
**long period** (low frequency), those terms make the layer's independent solutions grow
exponentially through the liquid, so by the surface they are nearly linearly dependent and the
surface boundary-condition matrix becomes near-singular. The solve is then unstable: small numerical
differences (integration tolerance, grid-vs-dense sampling, the linear-solver implementation) change
the result, and at sufficiently long periods the solve fails outright. This is inherent to the
dynamic-liquid assumption.

Guidance:

* Use a **dynamic** liquid layer only for **short-period** forcing, where it is well-conditioned.
* Use a **static** liquid layer for **long-period** forcing; it is stable and consistent across all
  periods.

(A future robustness improvement would be to re-orthonormalize the independent solutions during the
integration, Godunov-style, so the boundary-condition solve stays well-conditioned even at long
periods. This would benefit both the dense and the original solvers.)
