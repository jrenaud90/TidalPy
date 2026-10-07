# Calculating Love Numbers

_Updated: 2026-10-07_

`TidalPy.RadialSolver.radial_solver` is the array-based entry point to the viscoelastic-gravitational solve. You hand it a radial grid with density and complex moduli on it, a forcing frequency, and a description of the layers; it returns a [`RadialSolverSolution`](solution_class.md) carrying the radial functions and the Love numbers. If you already have a built world, prefer `BaseWorld.solve_love_numbers`, which fills these arrays from the layer rheologies for you.

## Example

```python
import numpy as np
from TidalPy.RadialSolver import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.Rheology import Elastic, Maxwell

# Use a helper function to build our radial solver inputs (see "build_inputs.md" for details on these helpers)
build_data = build_rs_input_homogeneous_layers(
    6000.0e3,                       # planet radius [m]
    2.0 * np.pi / (86400.0 * 7.5),  # forcing frequency [rad s-1]
    density_tuple                 = (8000.0, 5000.0, 3300.0),
    static_bulk_modulus_tuple     = (2.0e11, 1.5e11, 1.0e11),
    static_shear_modulus_tuple    = (1.0e11, 0.0, 5.0e10),
    bulk_viscosity_tuple          = (1.0e18, 1.0e18, 1.0e18),
    shear_viscosity_tuple         = (1.0e20, 1.0e3, 1.0e19),
    layer_type_tuple              = ("solid", "liquid", "solid"),
    layer_is_static_tuple         = (False, True, False),
    layer_is_incompressible_tuple = (False, False, False),
    shear_rheology_model_tuple    = (Maxwell(), Elastic(), Maxwell()),
    bulk_rheology_model_tuple     = Elastic(),
    radius_fraction_tuple         = (0.2, 0.55, 1.0),
    slices_tuple                  = (10, 12, 20))

solution = radial_solver(*build_data, degree_l=2, solve_for=("tidal",))
print(solution.k, solution.h, solution.l)
print(solution)                     # One line: success, degree, and k (the message when the solve failed)
```

The [input builders](build_inputs.md) exist because the solver's array requirements are strict. Build the arrays yourself only if you know they already satisfy the rules in the next section.

### Input Argument Rules

- Every array must be C-contiguous and of the stated dtype. The radius, density, and moduli arrays must all have the same size.
- The radius array starts at `r = 0` and strictly increases.
- Each layer needs at least 5 slices, so the arrays hold at least `5 * num_layers` points.
  - This entry point always builds an interpolated equation of state from the supplied arrays (`eos_method_bylayer` accepts only `"interpolate"`, and any other name raises), so the number of slices a layer needs depends on its profiles, as the last two rules explain.
- Every interface radius appears twice: once as the top of the lower layer and once as the base of the upper layer. The last radius is the planet radius and must equal the last entry of `upper_radius_bylayer_array`.
- A radius given twice inside a layer marks a discontinuity of its profile (PREM's 220, 400, and 670 km jumps, say). The solver accepts it, but the integrations then shorten their steps to cross the jump. Declaring a layer boundary at each such radius describes the same profile and is faster and more accurate: for PREM about half the time and about 20 times closer to the converged k2. The solver logs a warning naming the first such radius, once per session (`warn_if_internal_discontinuity` returns them all); `warnings=False` silences it. A world built from a data file is split at these radii already.
- The solver runs its own equation of state. Gravity, pressure, mass, and moment of inertia are integrated and read through the integrator's dense output at the exact integration radius. Density and the complex shear and bulk moduli are linearly interpolated from the supplied arrays at every integration radius, and the interpolated density also drives the structural integration, so the slices control these three.
- The error depends on how curved those profiles are within a layer, since linear interpolation reproduces a straight line exactly. For a single solid layer at degree 2, k2 from 5 slices is within 1e-9 of the converged value when density and the moduli are constant and within 1e-8 when they vary linearly, but a density falling exponentially by a factor of two across the layer is off by 1.3% at 5 slices, 0.2% at 10, and 0.03% at 25. Give a layer enough slices to follow the curvature of its profiles, and test the count against a refined run for your own problem.

## Choosing a Method

`love_method` selects how the Love numbers are obtained. Names are case-insensitive and the aliases in parentheses work everywhere.

| `love_method` | What it does |
|---|---|
| `radial_solver` (`shooting`, `rs`) | Default. Integrates the radial ODEs from the starting radius to the surface. Handles arbitrary multi-layer, solid or liquid, static or dynamic, compressible or incompressible interiors. |
| `propagation_matrix` (`prop_matrix`, `pm`, `prop`) | Quasi-analytic matrix propagation, valid only for a single solid, static, incompressible layer. It is more sensitive to the number of slices than the shooting method, since that sets the matrix dimension. Included mainly for comparison. |

The analytic methods (`homogeneous`, `cpl`, `ctl`) are not available here because this API takes moduli arrays rather than a layered world. Use `BaseWorld.solve_love_numbers(love_method=...)` for those, or the closed-form functions in [`TidalPy.Tides.love`](../Tides/love/love_numbers.md).

## Arguments

Every solver setting whose default is `None` takes its value from the TidalPy configuration: the `[radial_solver]` section for the shooting method and the `[eos_solver]` section for the equation of state (see [Configurations](../Overview/2_TidalPy_Configurations.md) for the packaged values and how they were chosen). These are the same defaults the world-attached `solve_love_numbers` and `solve_eos` use, so a configuration file plus a world file reproduces a result. An explicit argument always wins.

**Common arguments**

| Argument | Default | Meaning |
|---|---|---|
| `degree_l` | `2` | Spherical-harmonic degree. Stability degrades as the degree rises: expect trouble beyond `l = 10`, though stable solutions exist to `l = 40` with a higher starting radius and tighter tolerances. |
| `solve_for` | `None` | Tuple of surface boundary conditions: `'tidal'`, `'loading'`, `'free'`. `None` means `('tidal',)`. Note the trailing comma in a one-element tuple. Solving several at once is cheaper than separate calls, and sets the first dimension of `solution.love`. |
| `love_method` | `'radial_solver'` | The radial technique; see the table above. |
| `nondimensionalize` | `None` (config) | Non-dimensionalize the EOS and the shooting solve internally and restore SI before returning. Leave it on: the tolerances then mean the same thing for every planet. |

**Shooting method**

| Argument | Default | Meaning |
|---|---|---|
| `starting_radius` | `0.0` | Radius where integration begins [m]. `0.0` picks one automatically using the Martens (2016) criterion and `start_radius_tolerance`. Starting very deep at high degree makes the surface boundary solve ill-conditioned: the solution constants grow enormous and cancel, amplifying Love numbers error. The solver measures this on every solve and warns when the achievable accuracy drops below the requested tolerance; prefer the automatic radius when that warning appears. A manual radius above `[numerical] max_start_radius_fraction` of the planet radius (default 0.9) is refused, here with a `ValueError` and on the world path with a failed solve. |
| `start_radius_tolerance` | `None` (config) | Tolerance for that automatic choice: the start is at $R \cdot \mathrm{tol}^{1/l}$. |
| `starting_method` | `None` (config) | Starting conditions at the starting radius: `'takeuchi'` (Takeuchi and Saito 1972; alias `'ts'`), `'kamata'` (Kamata et al. 2015), `'power_series'` (Martens 2016; aliases `'powerseries'`, `'ps'`, `'martens'`), or `'unity'` (unit vectors). Every method covers every layer type. See [Starting Conditions](starting_conditions.md) for when to use each. |
| `degree1_frame` | `None` (config) | Reference frame of degree-1 load Love numbers and radial functions: `'CE'`, `'CM'`, `'CF'`, `'CL'`, or `'CH'`. Read only by a degree-1 solve for loading. See [Degree-1 Load Love Numbers](#degree-1-load-love-numbers). |
| `integration_method` | `None` (config) | `'RK23'`, `'RK45'`, `'DOP853'`, or the implicit methods `'BDF'`, `'LSODA'`, `'Radau'` for stiff problems. |
| `integration_rtol`, `integration_atol` | `None` (config) | Relative and absolute integration tolerances. |
| `scale_rtols_bylayer_type` | `None` (config) | Scale the relative tolerance by layer type; liquid layers generally want a tighter value. Experimental. |
| `max_num_steps` | `None` (config) | Step ceiling per integration. A healthy solve needs a few hundred steps per solution per layer. |
| `expected_size` | `None` (config) | Hint for the integrator's initial allocation. Overshooting costs little. |
| `max_ram_MB` | `None` (config) | Memory ceiling for the integrator. Real usage runs somewhat higher. |
| `max_step` | `0` | Largest allowed step [m]; `0` lets the integrator choose. |
| `love_only` | `False` | Keep only what the Love numbers need. The integration builds no dense output, so it is faster. The solution's radial functions (`result`, `get_radial_solution`, `plot_ys`) raise `ValueError`. |

**Propagation matrix**

| Argument | Default | Meaning |
|---|---|---|
| `core_model` | `0` | Inner-core starting condition. `0` Henning and Hurford (2014) seed matrix, `1` Roberts and Nimmo (2008) very small liquid core, `2` Henning and Hurford (2014) solid inner core, `3` Tobie et al. (2005) liquid inner core, `4` Sabadini and Vermeersen (2004) interface matrix. Only `0` is the regular solution of the modeled layer; the others change k2 of a uniform body by about 3 (r_start / R)^3, a few 1e-6 at the automatic starting radius, so they require the automatic starting radius and a manual `starting_radius` fails the solve (error code -22). |

**Equation of state**

The radial solver must have a EOS solution before it can solve the viscoelastic-gravitational problem. These arguments tell the solver what parameters to use for that additional solve.

| Argument | Default | Meaning |
|---|---|---|
| `eos_method_bylayer` | `None` | Per-layer EOS method. `"interpolate"` is the only accepted value and `None` selects it everywhere, so the only reason to pass this is explicitness; any other name raises. |
| `surface_pressure` | `0.0` | Pressure at the surface [Pa], the outer boundary condition for the interior pressure solve. |
| `eos_integration_method` | `None` (config) | As `integration_method`, for the EOS solve. LSODA's startup can fail cleanly at the singular center at tight tolerances; DOP853, BDF, and Radau handle it. |
| `eos_rtol`, `eos_atol` | `None` (config) | EOS integration tolerances. |
| `eos_pressure_tol` | `None` (config) | Convergence tolerance on the surface-pressure mismatch, relative to the central-pressure scale $(2/3) \pi G \rho^2 R^2$. Keep it above `eos_rtol`, the integrator's own noise on the surface pressure. |
| `eos_max_iters` | `None` (config) | Ceiling on the central-pressure iterations. The iteration is a secant method, so a compressible planet converges in a few steps; the cap is reported through `eos_iterations` and the message. |

**Reporting**

| Argument | Default | Meaning |
|---|---|---|
| `verbose` | `False` | Print solver status while running. |
| `warnings` | `True` | Emit solver warnings, including the surface-conditioning diagnostic. The diagnostics themselves (`surface_solve_amplification` and `surface_solve_rcond`) are recorded on every solve. |
| `raise_on_fail` | `False` | Raise instead of failing quietly. By default a failed solve returns with `success = False` and an explanatory `message`. |
| `perform_checks` | `True` | Accepted and ignored: the solver always validates its inputs. |
| `log_info` | `False` | Log the solution's key diagnostics. There is a cost, more so with file logging enabled. |

### Keyword Names on a World

`radial_solver` keeps the keyword names of the 0.7.X solver, while a world's `solve_love_numbers` and `solve_eos` take the shorter names of the `[radial_solver]` and `[eos_solver]` configuration keys. The settings are the same; only the names differ.

| `radial_solver` argument | World method argument | Configuration key |
|---|---|---|
| `start_radius_tolerance` | `solve_love_numbers(start_radius_tol=...)` | `[radial_solver] start_radius_tolerance` |
| `integration_method` | `solve_love_numbers(integration_method=...)` | `[radial_solver] integration_method` |
| `integration_rtol`, `integration_atol` | `solve_love_numbers(rtol=..., atol=...)` | `[radial_solver] rtol`, `atol` |
| `scale_rtols_bylayer_type` | `solve_love_numbers(scale_rtols=...)` | `[radial_solver] scale_rtols` |
| `starting_method`, `degree1_frame`, `max_num_steps`, `expected_size`, `nondimensionalize` | the same names | the same names |
| `max_ram_MB` | `solve_love_numbers(max_ram_MB=...)` | `[radial_solver] max_ram_mb` |
| `eos_integration_method` | `solve_eos(integration_method=...)` | `[eos_solver] integration_method` |
| `eos_rtol`, `eos_atol` | `solve_eos(rtol=..., atol=...)` | `[eos_solver] rtol`, `atol` |
| `eos_pressure_tol` | `solve_eos(pressure_tol=...)` | `[eos_solver] pressure_tol` |
| `eos_max_iters` | `solve_eos(max_iters=...)` | `[eos_solver] max_iters` |

The lower-level radial-solver functions that take Newton's constant (`find_starting_conditions` and the Kamata, Takeuchi, and power series starting conditions, `apply_surface_bc`, `solve_upper_y_at_interface`, `fundamental_matrix`), like the tidal-potential functions (`global_potential`, `tidal_potential_3d_modes`, `collapse_global_tides`) and the Kepler conversions, read a `G_to_use` of `None` as the TidalPy configuration's value (SciPy's G).

## Degree-1 Load Love Numbers

At degree 1 a rigid translation of the whole body (y1 = y3 = a constant, y5 = that constant times g, no stress) satisfies the equations and every surface condition of a load. The load Love numbers are therefore defined only once a reference frame is chosen (Farrell 1972; Blewitt 2003). A degree-1 tidal or free-surface response is only a translation, so only `solve_for=('loading',)` takes `degree_l=1`.

The solver finds the response in CE, the frame of the solid body's own center of mass, where k' = 0. It replaces the last surface condition (y6, or y7 for a static liquid surface layer) with y5 = 1 at the surface (Guo et al. 2004 Eq. 9; Martens 2016 Eq. 4.150). In a static body the replaced condition then holds by itself, because the body feels no net force (Saito 1974 Eq. 25). `degree1_frame` moves the result to another frame. The Love numbers and the radial functions move together: y1 and y3 shift by a constant c, y5 by c times g(r), and the stresses do not change.

| `degree1_frame` | Origin | Defined by | Shift alpha from CE |
|---|---|---|---|
| `'CE'` (default) | Center of mass of the solid body | k' = 0 | 0 |
| `'CM'` | Center of mass of the body and the load | 1 + k' = 0 | 1 |
| `'CF'` | Center of surface figure | h' + 2 l' = 0 | (h' + 2 l') / 3 |
| `'CL'` | Center of lateral figure | l' = 0 | l' |
| `'CH'` | Center of height figure | h' = 0 | h' |

In each frame h', l', and 1 + k' are their CE values minus alpha (Blewitt 2003 Eq. 17). h' - k' and l' - k' do not depend on the frame, and neither do the strain, the stress, or the dissipation. CE is the frame most published values use. CM and CF are the frames of satellite laser ranging and GNSS.

`Tests/Test_RadialSolver/test_degree_one_loading_01.py` checks the solver against two references:

- A homogeneous sphere with $V_P$ = 10 km/s, $V_S$ = 5 km/s, and $\rho$ = 5000 kg m$^{-3}$ gives h' = -0.206731 and l' = 0.161380. That is the exact value from the Takeuchi and Saito solutions; Martens (2016) Table C.4 gives -0.2069 and 0.1617. Every starting method and every tolerance from $10^{-6}$ to $10^{-12}$ agree, and a dynamic solve gives the same numbers from $\omega = 10^{-4}$ to $10^{-11}$ rad s$^{-1}$.
- The Guo et al. (2004) PREM gives h' = -0.28559 to -0.28564 and l' = 0.10372 to 0.10375 across static and dynamic layer choices. They give -0.285694 and 0.103633.

**Liquid surface layers.** A static liquid surface layer carries no horizontal displacement. Its l' is NaN, and the CF and CL frames, which read l', are refused with error code -16. Its surface is hydrostatic, which fixes h' = 1 - $\bar{\rho} / \rho_s$ in CE exactly, whatever lies beneath ($\bar{\rho}$ is the bulk density, $\rho_s$ the density at the surface). An incompressible solid lid of the same density tends to this value as its rigidity falls. A dynamic liquid surface layer solves like a solid one. Its l' is large, since the horizontal motion of an inviscid liquid grows as $1 / \omega^2$, and at long periods a dynamic compressible liquid fails as it does at any degree, or returns a wrong answer flagged by a large `surface_frame_residual` and a conditioning warning.

**Mixed static and dynamic layers.** A body with inertia in some layers but not others has no exact degree-1 frame: the replaced condition holds only to about $\omega^2 R / g$. `solution.surface_frame_residual` (`world.love_surface_frame_residual` on a world) reports how far it is from met, relative to the y6 condition, and the solver warns when it is above $10^{-2}$. Dynamic solid layers under a static ocean leave $5 \times 10^{-3}$ at a one-day period. Making every layer static, or every layer dynamic, removes it.

Other notes:

- k' is 0 in CE and -1 in CM, so its quality factor and lag mean nothing there.
- The propagation matrix refuses degree 1 (error code -13). The `homogeneous`, `cpl`, and `ctl` methods give tidal Love numbers only, at every degree.

## Troubleshooting

Start with `solution.message`, then `solution.steps_taken`, then plot. The messages below come from the shooting method.

**"Maximum number of steps (set by user) exceeded during integration."** The solve hit `max_num_steps`. Raising it is rarely the real fix: this normally means the solution is unstable, so work through the list below first.

**"Maximum number of steps (set by system architecture) exceeded during integration."** The integrator's arrays passed `max_ram_MB`. Same advice as above.

**"Error in step size calculation: Required step size is less than spacing between numbers."** The integrator could not make progress. Usually the tolerances are too tight, or counterintuitively too loose so that error compounds. A dynamic (non-static) compressible liquid layer is a common trigger; try the static assumption for it.

**Slow solves or huge step counts.** A healthy solve needs a few hundred steps per solution per layer; a few thousand happens in awkward cases; ten thousand or more means the solution is likely unstable. `solution.plot_ys()` is a quick check, since instability shows up as large spikes, ringing, or curves that do not vary smoothly with radius. Things to try, roughly in order: change the integration tolerances, change the integration method, switch the starting conditions with `starting_method` (see [Starting Conditions](starting_conditions.md)), lower `degree_l`, start higher in the planet with `starting_radius`, revisit the layer assumptions (dynamic compressible liquid layers especially), add a small solid core beneath a fully liquid one, or add more slices to the input arrays.

**NaN Love numbers from a successful solve.** The surface boundary condition solve was ill-conditioned. Check `solution.surface_solve_amplification`: values far above one mean the solution constants are cancelling catastrophically. Raise the starting radius, or let the solver choose it.

**"The surface boundary condition system is singular to working precision" (error code -13).** The reciprocal condition number of the surface system, `solution.surface_solve_rcond`, fell below `[numerical] minimum_surface_rcond` (default `1e-12`), so no set of solution constants is determined by the surface conditions. The independent solutions have become numerically dependent; start higher in the planet or use the automatic starting radius. A very weak solid starting layer does this for every start; see [Starting Conditions](starting_conditions.md).

**"The degree-1 frame asked for reads l'" (error code -16).** The `'CF'` and `'CL'` frames need l', which a static liquid surface layer does not define. Use `'CE'`, `'CM'`, or `'CH'`, or make the surface layer dynamic.

**"The surface condition its reference frame replaced is met only to ..." (warning).** A degree-1 loading solve on a body whose layers are partly static and partly dynamic. See [Degree-1 Load Love Numbers](#degree-1-load-love-numbers).

**"Layer ... is compressible but its bulk modulus is not positive" (error code -15).** A compressible solid or dynamic liquid layer needs a positive bulk modulus. An interpolated equation of state given no bulk-modulus table reports none (NaN), and a bulk modulus of zero with a shear modulus is a Poisson ratio of -1, which would give wrong Love numbers without any other sign of trouble. Give the material a bulk modulus, or mark the layer incompressible. A value between the threshold and the integration `rtol` is solved but logged as poorly conditioned, since the constants then carry the integration error divided by roughly `surface_solve_rcond`.

**"The starting radius ... is not inside the planet" or "No radial slice ... lies above the starting radius" (error code -5).** A manual starting radius must lie inside the planet with at least one radial slice above it in its layer. A starting radius exactly on an interface begins in the layer above it.

**A crash with no exception.** Rerun with the same inputs while watching memory. If it reproduces, record the inputs and open an issue on [GitHub](https://github.com/jrenaud90/TidalPy/issues).
