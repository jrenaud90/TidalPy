# The Solution Class

_Updated: 2026-09-12_

Every radial solve returns a `RadialSolverSolution`, the Cython class in `TidalPy.RadialSolver_x.rs_solution`. It holds three things: whether the solve worked, what the interior looks like (the equation-of-state result), and the viscoelastic-gravitational answer itself (the radial functions and the Love numbers). The same object comes back from `radial_solver`, from `homogeneous_love_numbers`, and from a world's released radial storage.

```python
solution = radial_solver(*build_data, degree_l=2, solve_for=("tidal", "loading"))
if not solution.success:
    raise RuntimeError(solution.message)
k2 = solution.k[0]      # tidal k2; index 1 would be the loading value
```

## Did it work

Check `success` before trusting anything else. A failed solve returns normally unless you passed `raise_on_fail=True`.

| Member | Type | Meaning |
|---|---|---|
| `success` | bool | The radial solve converged. |
| `error_code` | int | `0` when there was no error. |
| `message` | str | Status text, and the first thing to read after a failure. See the troubleshooting section of [Calculating Love Numbers](calculating_love_numbers.md). |
| `surface_solve_amplification` | float | Worst-case error amplification of the surface boundary-condition solve (shooting method). Near one is healthy. Large values mean the solution constants are cancelling, so roundoff and integration error are being amplified into the Love numbers, and the achievable accuracy is roughly this number times machine epsilon. |
| `steps_taken` | int array `(num_layers, 3)` | Integration steps per layer per independent solution. Solid layers use three solutions, dynamic liquid layers two, static liquid layers one; unused entries are zero. A few hundred per solution per layer is normal, a few thousand is tolerable, and ten thousand or more means the solve is unstable. |
| `print_diagnostics(print_diagnostics=True, log_diagnostics=False)` | method | Assemble a readable summary of the solve. Printing it is the fastest triage; logging it sends the same text to TidalPy's log. |

## The interior

The solver runs an equation of state before the deformation problem, and keeps the result. These are the profiles that produced the moduli the deformation solve used, so they are the right thing to inspect when a solve behaves oddly. See the [equation of state](../material_x/material_eos.md) documentation for the models themselves.

| Member | Meaning |
|---|---|
| `eos_success`, `eos_error_code`, `eos_message` | Outcome of the equation-of-state solve. |
| `eos_pressure_error`, `eos_iterations`, `eos_steps_taken` | Convergence of the interior pressure iteration. |
| `radius_array`, `gravity_array`, `pressure_array`, `mass_array`, `moi_array`, `density_array` | Real-valued profiles through the planet. |
| `shear_modulus_array`, `bulk_modulus_array` | Complex moduli through the planet. |
| `layer_upper_radius_array` | Upper radius of each layer [m]. |
| `radius`, `volume`, `mass`, `moi`, `density_bulk` | Whole-planet scalars. |
| `moi_factor` | The moment of inertia normalized by the uniform-sphere value `0.4 M R^2`, so exactly 1 for a uniform body and below 1 for a centrally condensed one. Note that this is not the conventional moment of inertia factor `C / (M R^2)`, which is 0.4 for a uniform sphere; multiply by 0.4 to get that. |
| `central_pressure`, `surface_pressure`, `surface_gravity` | Boundary values. |
| `eos_call(radius)` | Dense equation-of-state outputs at any radius, evaluated from the solver's own interpolant rather than re-interpolating the gridded arrays. |
| `eos_call_si(radius)` | The same, in SI units. |

## The radial functions

`result` is the raw block of radial functions, shaped `(num_ytypes * 6, num_slices)`: the six functions of the first boundary condition, then the six of the next, and so on. Index it by name instead when you have more than one. TidalPy follows the Takeuchi and Saito (1972) convention, so in a solid layer these are the familiar y1 through y6.

Liquid layers do not define all six. A dynamic liquid layer has no y4, and a static liquid layer additionally has no y2, y3, y5, or y6; the Saito (1974) variable "y7" takes the y6 slot there. Undefined entries are NaN, which is a feature rather than a fault: it keeps the array shape uniform and makes an accidental use obvious.

```python
solution.result                  # (num_ytypes * 6, num_slices)
solution["tidal"]                # the tidal block by name, always (6, num_slices)
solution["loading"]              # the loading block, if it was solved for
solution.get_radial_solution(1.5e6)              # complex y1..y6 at one radius [m]
solution.get_radial_solution_array(radii)        # (n, 6) at many radii, all in C++
```

The two dense getters evaluate the shooting method's per-layer interpolants, so they are accurate anywhere, including between grid slices. They are the recommended way to ask for values at a radius; see [Dense Radial Solutions](dense_radial_solution.md).

## The Love numbers

`love` is a complex array shaped `(num_solve_for, 3)`: the first axis follows the order you passed to `solve_for`, and the second is k, h, l. The shortcuts return scalars when you solved for one boundary condition and arrays when you solved for several.

| Member | Meaning |
|---|---|
| `love` | The full block, `(num_solve_for, 3)`. |
| `k`, `h`, `l` | Potential, radial displacement, and tangential displacement Love numbers. |
| `Q_k`, `Q_h`, `Q_l`, `Q` | Effective dissipation quality factor, defined as the magnitude of the real part over the negative imaginary part. `Q` follows k. |
| `lag_k`, `lag_h`, `lag_l`, `lag` | Phase lag [rad], defined as the arctangent of the negative imaginary part over the real part. `lag` follows k. |
| `degree_l` | The harmonic degree that was solved. |

## Plotting

Both methods wrap [`Utilities_x.graphics_x`](../utilities_x/graphics_x.md) and return the matplotlib figure and axes, so you can keep adjusting them. Extra keyword arguments pass straight through.

```python
solution.plot_ys()                              # the six radial functions against radius
solution.plot_ys(plot_imaginary=True, benchmarks="tobie2005")
solution.plot_interior(planet_name="Enceladus")  # gravity, density, pressure, moduli
```

`plot_ys` is the fastest instability check there is. A converged solve gives smooth curves; large spikes, ringing, or kinks that do not follow the layer structure mean the integration did not converge. `plot_interior` needs a successful equation-of-state solve, and both raise an informative error rather than returning nothing when the underlying solve failed.
