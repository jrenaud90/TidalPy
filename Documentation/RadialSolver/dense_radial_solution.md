# Dense Radial Solutions

_Updated: 2026-10-06_

The shooting method keeps its whole integrated solution, so the radial functions `y1..y6` and the interior can be evaluated at any radius after the solve, not only on the input grid. This page covers how to do that, the propagation-matrix equivalent, and the limits of dynamic liquid layers at long forcing periods.

## Python API

```python
rs = radial_solver(...)            # RadialSolverSolution

# Radial functions at any radius (SI, complex y1..y6):
y  = rs.get_radial_solution(300.0)           # length-6 array at r = 300 m, first boundary condition
ys = rs.get_radial_solution_array(r_array)   # (N, 6), vectorized

# Interior and material state at any radius (SI):
eos = rs.eos_call(300.0)             # dict of scalars at r = 300 m
eos["density"], eos["pressure"], eos["complex_shear_modulus"]
profile = rs.eos_call(r_array)       # arrays in, arrays out

# Gridded results (the dense solution sampled on the input grid):
rs.result    # y-solution
rs.love      # k, h, l
```

`eos_call(radius)` is the interior counterpart of `get_radial_solution`. It returns a dict of named fields (listed by `EOS_CALL_FIELDS` and in [The Solution Class](solution_class.md#the-interior)): scalars for a float radius, arrays of the input's shape for an array, NaN outside the body. It evaluates the solver's own interpolant, so it gives the values the solve used rather than a re-interpolation of the gridded arrays, and it stays valid after a standalone solve returns.

Its `shear_modulus` and `bulk_modulus` are the frequency-independent, unrelaxed moduli. A viscoelastic response belongs to a rheology at a forcing frequency, not to the equation of state, so the `complex_shear_modulus` and `complex_bulk_modulus` fields, like `get_complex_shear_modulus(radius)` and `get_complex_bulk_modulus(radius)`, come from a second read. It reproduces the moduli the solve used: from the world's rheology (shared with the solution, so it can outlive the world) at the frequency recorded in `love_frequency`, or, for a solve handed its moduli as arrays, by interpolating those arrays.

## Shooting Method

The planet is integrated layer by layer from a starting radius up to the surface. Each layer contributes a fixed number of independent solutions:

| Layer type            | Independent solutions |
|-----------------------|-----------------------|
| Solid                 | 3                     |
| Dynamic liquid        | 2                     |
| Static liquid         | 1                     |

The physical solution in each layer is a linear combination of that layer's independent solutions. The coefficients ("constants of integration") are fixed by the surface boundary condition, a small matrix solve at the surface, then carried down through every interface one layer at a time, using the same interface conditions as the upward integration.

Each independent solution is integrated with dense (continuous) output, and the interpolants are kept with the per-layer constants. A value at any radius is the containing layer's 1 to 3 interpolants evaluated there, multiplied by the constants and summed, with `y3` of a dynamic liquid rebuilt from the other components and the EOS gravity and density, and returned in SI. A value between grid slices therefore comes from the integrator's own polynomial rather than from a linear interpolation of a sampled grid.

## Propagation Matrix Method

The propagation matrix builds its solution on a grid by propagating each slice's fundamental matrix from one slice radius to the next. `get_radial_solution(r)` continues that propagation to the radius asked for: the slice's fundamental matrix at `r` applied to the coefficients carried into the slice. This reproduces the grid values at both ends of the slice and is exact between them for a uniform layer. Below the first propagated slice, the regular solution of the innermost material continues the solution to the center, so `get_radial_solution` and `result` are defined down to `r = 0`. With a `core_model` other than `0` the seed is not the regular solution, and the solution stays NaN below the starting radius.

## Dynamic Liquid Layers at Long Forcing Periods

The dynamic liquid equations carry no density-gradient term, so a liquid layer's stratification comes from how its density, gravity, and bulk modulus vary with radius. A liquid whose density follows its bulk modulus, $d\rho/dr = -\rho^2 g / K$, is neutral ($N^2 = 0$). A constant-density liquid with a finite bulk modulus is not: $N^2 = -\rho g^2 / K < 0$, an unstable stratification. Its solutions then grow as $e^E$ through the layer, with

$$E = \frac{\sqrt{\ell(\ell+1)}}{\omega} \int \frac{\sqrt{-N^2}}{r}\, dr,$$

which rises with the forcing period. Once $e^E$ amplifies the integration tolerance to order one, the independent solutions are nearly dependent by the surface and the solve is wrong or fails. Every bundled liquid except PREM's has constant density: solved dynamically, Earth-Simple's core breaks down by half a day, Mercury's by 3.5 days, the Moon's and PREM's by 10 days, and Pluto's ocean by 15 days.

Before every radial Love solve (`solve_love_numbers`, `calc_tides`, the 3D calls, and `radial_solver`), TidalPy estimates each dynamic liquid layer's error, $\mathrm{rtol}\, e^E$, from its profile, and logs one warning per world when it passes 1%, naming the layer, the period, and $E$.

Guidance:

* Use a static liquid layer for long-period forcing. It is stable and accurate at every period for any density profile.
* To keep a constant-density liquid dynamic, make it incompressible (`is_incompressible = True`), which makes it neutral. In the bundled worlds its dynamic solve then stays within 1e-3 of the static one out to 100 days, except Earth-Simple's large core, which drifts to 6% there.
* To keep it dynamic and compressible, give it a pressure-dependent law (Birch-Murnaghan or Vinet), which keeps it neutral; the `luna_dynamic`, `mercury_dynamic`, `pluto_dynamic`, and `europa_dynamic` bundled worlds are built this way. A large neutral compressible core still loses some accuracy at long periods, since the equations recover its tangential displacement through a division by $\omega^2$: `luna_dynamic`'s yearly $k_2$ is 3e-4 off at the default `rtol = atol = 3e-8`, 4e-5 at `3e-9`, and 0.8% at `rtol = 1e-6`. Thin oceans are unaffected.
* The dynamic liquids demo (`Demos/Physics/P12_dynamic_liquids.ipynb`) works through each case.

## C++ API

`c_RadialSolutionStorage` provides `get_radial_solution`, `get_radial_solution_array`, `get_surface_y`, `get_eos_si`, and `get_radial_solution_nondim` (for a radius in solve units). A built `c_BaseWorld` provides `get_radial_solution_y` and `get_love_surface_y`. `c_EOSSolution::call_material` returns a `c_EOSMaterialState` (gravity, density, and both complex moduli) whatever the solution was built from.

## Validation

`Tests/Test_RadialSolver/test_dense_benchmark/` holds frozen reference results for a 1-layer solid planet, a 2-layer solid planet, and a solid / dynamic-liquid / solid planet at a short forcing period. The dense solution agrees closely with them at the surface and at layer boundaries, and the Love numbers agree tightly; between grid slices the dense evaluation is the more accurate.
