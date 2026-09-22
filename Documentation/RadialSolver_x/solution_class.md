# The Solution Class

_Updated: 2026-09-21_

Every radial solve returns a `RadialSolverSolution`, the Cython class in `TidalPy.RadialSolver_x.rs_solution`. It holds the solve status, the equation-of-state result, and the radial functions and Love numbers. The same object comes back from `radial_solver`, from `homogeneous_love_numbers`, and from a world's released radial storage.

```python
solution = radial_solver(*build_data, degree_l=2, solve_for=("tidal", "loading"))
if not solution.success:
    raise RuntimeError(solution.message)
k2 = solution.k[0]      # tidal k2; index 1 would be the loading value
```

## Success Check

Check `success` before trusting anything else. A failed solve returns normally unless you passed `raise_on_fail=True`.

| Member | Type | Meaning |
|---|---|---|
| `success` | bool | The radial solve converged. |
| `error_code` | int | `0` when there was no error. |
| `message` | str | Status text, and the first thing to read after a failure. See the troubleshooting section of [Calculating Love Numbers](calculating_love_numbers.md). |
| `surface_solve_amplification` | float | Worst-case error amplification of the surface boundary-condition solve (shooting method). Near one is healthy. Large values mean the solution constants are cancelling, so roundoff and integration error are being amplified into the Love numbers, and the achievable accuracy is roughly this number times machine epsilon. |
| `steps_taken` | int array `(num_layers, 3)` | Integration steps per layer per independent solution. Solid layers use three solutions, dynamic liquid layers two, static liquid layers one; unused entries are zero. A few hundred per solution per layer is normal, a few thousand is tolerable, and ten thousand or more means the solve is likely unstable. |
| `print_diagnostics(print_diagnostics=True, log_diagnostics=False)` | method | Assemble a readable summary of the solve. Printing it is a quick triage; logging sends the same text to TidalPy's log. |

## The Interior

The solver runs an equation of state before the deformation problem, and keeps the result. These are the profiles that produced the moduli the deformation solve used, so they should be inspected when a solve behaves oddly. See the [equation of state](../material_x/material_eos.md) documentation for the models themselves.

| Member | Meaning |
|---|---|
| `eos_success`, `eos_error_code`, `eos_message` | Outcome of the equation-of-state solve. |
| `eos_pressure_error`, `eos_iterations`, `eos_steps_taken` | Convergence of the interior pressure iteration. |
| `get_gravity(r)`, `get_pressure(r)`, `get_mass(r)`, `get_moi(r)`, `get_density(r)` | The interior at radius `r` [m], a float or an array of radii. The solution keeps the solved EOS it came from and evaluates it, so nothing is tabulated and any radius may be asked for. |
| `get_shear_modulus(r)`, `get_bulk_modulus(r)` | The material's **static** (unrelaxed) moduli [Pa] at `r`. These are frequency-independent, which is why they are what the dense equation-of-state readout carries. |
| `get_complex_shear_modulus(r)`, `get_complex_bulk_modulus(r)` | The **complex** moduli [Pa] at `r`, as this solve used them. |
| `love_frequency` | The forcing frequency [rad s-1] the complex moduli above are evaluated at; NaN when the solve carried no rheology. |
| `get_shear_viscosity(r)`, `get_bulk_viscosity(r)` | Viscosities [Pa s] at `r`; NaN when the material names none. |
| `sample_radii(num_points=0)` | A radius grid [m] spanning the body, for a caller that wants one (plotting, tabulating). Nothing in the solve uses it and the solution keeps no copy; it defaults to the slice count the solve was configured with. |
| `layer_upper_radius_array` | Upper radius of each layer [m]. |
| `radius`, `volume`, `mass`, `moi`, `density_bulk` | Whole-planet scalars. |
| `moi_factor` | The moment of inertia factor $C/(MR^2)$, with $C$ the moment of inertia `moi`: 0.4 for a uniform sphere, 0.3307 for Earth, and smaller the more mass sits near the center. |
| `moi_sphere_ratio` | The same moment of inertia measured against a uniform sphere of equal mass and radius, $C/(0.4\, MR^2)$: exactly 1 when uniform, below 1 when centrally condensed. It is 2.5 times `moi_factor`. |
| `central_pressure`, `surface_pressure`, `surface_gravity` | Boundary values. |
| `eos_call(radius)` | The dense equation-of-state and material state at an SI radius [m] (a float, or an array for arrays out) as a dict of named fields: `gravity`, `pressure`, `mass`, `moi`, `density`, `shear_modulus`, `bulk_modulus`, `shear_viscosity`, `bulk_viscosity`, `temperature`, `heat_flow`, `melt_fraction`, and the `complex_shear_modulus` and `complex_bulk_modulus` the solve used at `love_frequency`. Evaluated from the solver's own interpolant; NaN outside the body. |
| `eos_call_nondim(radius)` | The raw dense row at a radius in the solve's own units (internal; `EOS_CALL_FIELDS` gives the slot order). |

## Radial Functions

`result` is the raw block of radial functions, shaped `(num_ytypes * 6, num_slices)`: the six functions of the first boundary condition, then the six of the next, and so on. Index it by name instead when you have more than one. TidalPy follows the Takeuchi and Saito (1972) convention, so in a solid layer these are the familiar y1 through y6.

Liquid layers do not define all six. A dynamic liquid layer has no y4, and a static liquid layer has no y2, y3, y4, y5, or y6; the Saito (1974) variable "y7" takes the y6 slot there. Undefined entries are NaN, which keeps the array shape uniform and makes an accidental use obvious.

```python
solution.result                            # (num_ytypes * 6, num_slices)
solution["tidal"]                          # the tidal block by name, always (6, num_slices)
solution["loading"]                        # the loading block, if it was solved for
solution.get_radial_solution(1.5e6)        # complex y1..y6 at one radius [m]
solution.get_radial_solution_array(radii)  # (n, 6) at many radii, all in C++
```

The two dense getters evaluate the shooting method's per-layer interpolants, so they are accurate anywhere, including between grid slices. They are the recommended way to ask for values at a radius; see [Dense Radial Solutions](dense_radial_solution.md).

## Love Numbers

`love` is a complex array shaped `(num_solve_for, 3)`: the first axis follows the order you passed to `solve_for`, and the second is k, h, l. The shortcuts return scalars when you solved for one boundary condition and arrays when you solved for several.

| Member | Meaning |
|---|---|
| `love` | The full block, `(num_solve_for, 3)`. |
| `k`, `h`, `l` | Potential, radial displacement, and tangential displacement Love numbers. |
| `Q_k`, `Q_h`, `Q_l`, `Q` | Effective dissipation quality factor, defined as the magnitude of the real part over the negative imaginary part. `Q == Q_k`. |
| `lag_k`, `lag_h`, `lag_l`, `lag` | Phase lag \[rad\], defined as the arctangent of the negative imaginary part over the real part. `lag` follows k. |
| `degree_l` | The harmonic degree that was solved. |

## Plotting

Both methods wrap [`Utilities_x.graphics_x`](../utilities_x/graphics_x.md) and return the matplotlib figure and axes, so you can keep adjusting them. Extra keyword arguments pass straight through.

```python
solution.plot_ys()                                             # Plot of the six radial functions against radius
solution.plot_ys(plot_imaginary=True, benchmarks="tobie2005")  # Used to compare to Tobie et al. (2005)
solution.plot_interior(planet_name="Enceladus")                # Plot of EOS results: gravity, density, pressure, moduli
```

`plot_ys` is a fast instability check. A converged solve gives smooth curves; large spikes, ringing, or kinks that do not follow the layer structure mean the integration did not converge. Keep in mind that it is not unusual to get spikes near liquid-solid layer boundaries even in stable solutions. `plot_interior` needs a successful equation-of-state solve, and both raise an informative error rather than returning nothing when the underlying solve failed.
