# The Solution Class

_Updated: 2026-10-09_

Every radial solve returns a `RadialSolverSolution` (from `TidalPy.RadialSolver.rs_solution`), holding the solve status, the equation-of-state result, the radial functions, and the Love numbers. `radial_solver`, `homogeneous_love_numbers`, and a world's released radial storage all return it.

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
| `message` | str | Status text; read it first after a failure. See the troubleshooting section of [Calculating Love Numbers](calculating_love_numbers.md). |
| `surface_solve_amplification` | float | Worst-case error amplification of the surface boundary-condition solve (shooting method). Near one is healthy; large values mean the solution constants cancel, and the achievable accuracy is roughly this number times machine epsilon. It measures cancellation only, so it misses other semi-failures (like singular matrices). |
| `surface_solve_rcond` | float | Reciprocal condition number of the surface system (shooting method; NaN for the propagation matrix or a solve that stopped early), independent of units and of the starting normalization. At most one; healthy solves measure about `4e-5` to `1`, and a near-fluid layer stayed accurate at `3e-8`. Below the integration `rtol` is logged as poorly conditioned, below `[numerical] minimum_surface_rcond` (default `1e-12`) fails the solve (error code `-13`). Degree-1 loading measures the system with its frame row (y5 in place of y6). |
| `surface_frame_residual` | float | Degree-1 loading only (NaN otherwise): how far the surface condition replaced by the frame row is from met, relative to the y6 condition, worst over the boundary conditions. About roundoff when every layer is static, about $\omega^2 R / g$ with mixed static and dynamic layers; warned above `1e-2`. See [Degree-1 Load Love Numbers](calculating_love_numbers.md#degree-1-load-love-numbers). |
| `steps_taken` | int array `(num_layers, 3)` | Integration steps per layer, repeated for each of its independent solutions (three in a solid, two in a dynamic liquid, one in a static liquid; unused entries are zero), since a layer's solutions are integrated together. Tens to a few hundred per layer is normal; tens of thousands mean the solutions change on a scale far shorter than the layer, as in a strongly stratified dynamic liquid at a long period. |
| `orthonormalizations` | int array `(num_layers,)` | How many times the shooting method replaced each layer's solutions with an orthonormal basis of their span because they had become nearly dependent (see [Dense Radial Solutions](dense_radial_solution.md#shooting-method)). Zero is typical at short periods. Near-static solids take tens (a solid TRAPPIST-1b core took 81 at $10^{-16}$ rad/s), and a strongly stratified dynamic liquid takes more the longer the period (Earth-Simple's constant-density core made dynamic took 184, 1857, and 18,582 at 10, 100, and 1000 days). |
| `print_diagnostics(print_diagnostics=True, log_diagnostics=False)` | method | A readable summary of the solve, printed, logged, or both. |

## The Interior

The solution keeps the equation-of-state result behind the moduli the solve used; inspect it when a solve behaves oddly (models: [equation of state](../Material/material_eos.md)).

| Member | Meaning |
|---|---|
| `eos_success`, `eos_error_code`, `eos_message` | Outcome of the equation-of-state solve. |
| `eos_pressure_error`, `eos_iterations`, `eos_steps_taken` | Convergence of the interior pressure iteration. |
| `get_gravity(r)`, `get_pressure(r)`, `get_mass(r)`, `get_moi(r)`, `get_density(r)` | The interior at any radius `r` [m] (a float or an array). |
| `get_shear_modulus(r)`, `get_bulk_modulus(r)` | The **static** (unrelaxed, frequency-independent) moduli [Pa]. |
| `get_complex_shear_modulus(r)`, `get_complex_bulk_modulus(r)` | The **complex** moduli [Pa] as this solve used them. |
| `love_frequency` | The forcing frequency [rad s-1] of those complex moduli; NaN when the solve carried no rheology. |
| `get_shear_viscosity(r)`, `get_bulk_viscosity(r)` | Viscosities [Pa s]; NaN when the material names none. |
| `sample_radii(num_points=0)` | The radii [m] `result` was sampled on (the radius array passed to `radial_solver`, or the world's solve grid for a released solution). With `num_points`, that many evenly spaced radii across the body. |
| `layer_upper_radius_array` | Upper radius of each layer [m]. |
| `radius`, `volume`, `mass`, `moi`, `density_bulk` | Whole-planet scalars. |
| `moi_factor` | Moment of inertia factor $C/(MR^2)$, with $C$ = `moi`: 0.4 for a uniform sphere, 0.3307 for Earth, smaller the more mass sits near the center. |
| `moi_sphere_ratio` | $C/(0.4\, MR^2)$, the moment of inertia against a uniform sphere of equal mass and radius: 1 when uniform, below 1 when centrally condensed (2.5 times `moi_factor`). |
| `central_pressure`, `surface_pressure`, `surface_gravity` | Boundary values. |
| `eos_call(radius)` | Every field of the interior at an SI radius [m] (float or array) as a dict: `gravity`, `pressure`, `mass`, `moi`, `density`, `shear_modulus`, `bulk_modulus`, `shear_viscosity`, `bulk_viscosity`, `temperature`, `heat_flow`, `melt_fraction`, `buoyancy_frequency_squared` ($N^2 = -g\,(\rho'/\rho + \rho g / K)$ [s-2] along the solved structure, with $K$ the reported `bulk_modulus`, which the dynamic liquid equations read; positive where stably stratified and exactly 0 where the density follows the bulk modulus), `density_gradient` (d rho / dr [kg m-4] along the solved structure, which an incompressible dynamic liquid reads), and the `complex_shear_modulus` and `complex_bulk_modulus` used at `love_frequency`. NaN outside the body. See [Dense Radial Solutions](dense_radial_solution.md). |

## Radial Functions

`result` is the raw block of radial functions, shaped `(num_ytypes * 6, num_slices)`: the six functions (Takeuchi and Saito 1972 convention) of the first boundary condition, then the next six, and so on. With more than one, index by name.

Liquid layers do not define all six; undefined entries are NaN, which keeps the shape uniform. A dynamic liquid layer has no y4 (its y3 is rebuilt from its pressure variable, $y_3 = -P / (\rho \omega^2 r)$ with $P = y_2 - \rho g y_1 + \rho y_5$). A static liquid layer defines only y5; the Saito (1974) variable $y_7 = y_6 + (4 \pi G / g)\, y_2$ it integrates is not stored. At the free surface of a static liquid top layer, y2 is the surface boundary condition and $y_2 = \rho (g y_1 - y_5)$, so y1 and y2 are defined and h is finite (for a tidal solve h = 1 + k), while y3, and with it l, stays NaN.

```python
solution.result                            # (num_ytypes * 6, num_slices)
solution["tidal"]                          # the tidal block by name, always (6, num_slices)
solution["loading"]                        # the loading block, if it was solved for
solution.get_radial_solution(1.5e6)        # complex y1..y6 at one radius [m]
solution.get_radial_solution_array(radii)  # (n, 6) at many radii
```

The two dense getters are accurate at any radius, including between grid slices, and are the recommended way to ask for values at a radius (see [Dense Radial Solutions](dense_radial_solution.md)). A solve run with `love_only=True` (reported by `love_only`) keeps only the Love numbers and surface values: `result`, indexing by name, the dense getters, and `plot_ys` raise `ValueError`, while the interior members still answer.

## Love Numbers

`love` is a complex array shaped `(num_solve_for, 3)`: the first axis follows the order of `solve_for`, the second is k, h, l. The shortcuts return scalars for one boundary condition and arrays for several.

| Member | Meaning |
|---|---|
| `love` | The full block, `(num_solve_for, 3)`. |
| `k`, `h`, `l` | Potential, radial displacement, and tangential displacement Love numbers. |
| `Q_k`, `Q_h`, `Q_l`, `Q` | Effective quality factor, $-s \lvert k \rvert / \mathrm{Im}\, k$ with $s$ the sign of $\mathrm{Re}\, k$ ($\lvert k \rvert / (-\mathrm{Im}\, k)$ for the positive tidal k). Positive for a dissipative response of either sign (the loading k' and h' are negative, with a positive imaginary part when they lag); infinite when the imaginary part is zero. `Q == Q_k`. |
| `lag_k`, `lag_h`, `lag_l`, `lag` | Phase lag \[rad\], $\mathrm{atan2}(-s\, \mathrm{Im}\, k, \lvert \mathrm{Re}\, k \rvert)$ with $s$ as above, positive for a dissipative response of either sign. `lag` follows k. |
| `degree_l` | The harmonic degree solved. |

## Plotting

Both methods wrap [`Utilities.graphics`](../Utilities/graphics.md), pass extra keyword arguments through, and return the matplotlib figure and axes.

```python
solution.plot_ys()                                             # The six radial functions against radius
solution.plot_ys(plot_imaginary=True, benchmarks="tobie2005")  # Compare to Tobie et al. (2005)
solution.plot_interior(planet_name="Enceladus")                # EOS results: gravity, density, pressure, moduli
```

`plot_ys` is a fast instability check: large spikes, ringing, or kinks that do not follow the layer structure mean the integration did not converge, though spikes near liquid-solid boundaries can occur in stable solutions. `plot_interior` needs a successful equation-of-state solve; both raise an informative error when the solve failed.
