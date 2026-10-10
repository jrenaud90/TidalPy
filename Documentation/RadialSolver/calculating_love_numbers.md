# Calculating Love Numbers

_Updated: 2026-10-07_

`TidalPy.RadialSolver.radial_solver` is the array-based entry point to the viscoelastic-gravitational solve. Give it a radial grid with density and complex moduli, a forcing frequency, and the layer assumptions. It returns a [`RadialSolverSolution`](solution_class.md) with the radial functions and the Love numbers. With a built world, prefer `BaseWorld.solve_love_numbers`, which fills these arrays from the layer rheologies.

## Example

```python
import numpy as np
from TidalPy.RadialSolver import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.Rheology import Elastic, Maxwell

# A helper builds the solver inputs (see "build_inputs.md")
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

The [input builders](build_inputs.md) produce arrays that meet the rules below. Build the arrays yourself only if yours already do.

### Input Argument Rules

- Arrays are C-contiguous and of the stated dtype; the radius, density, and moduli arrays have one size.
- The radius starts at `r = 0`, strictly increases, and ends at the planet radius (the last entry of `upper_radius_bylayer_array`).
- Every interface radius appears twice: top of the lower layer, base of the upper one.
- Each layer has at least 5 slices.
- Density and the complex moduli are linearly interpolated between slices (`"interpolate"` is the only equation of state here), and that density drives the gravity and pressure integration. So the error depends on how curved each profile is within a layer. For one solid layer at degree 2, 5 slices give k2 within 1e-9 of converged for constant properties and 1e-8 for linear ones; a density falling exponentially by a factor of two across the layer gives 1.3% at 5 slices, 0.2% at 10, and 0.03% at 25. Test the count against a refined run.
- A radius repeated inside a layer marks a jump in its profile (PREM's 220, 400, and 670 km discontinuities). It is accepted, but a layer boundary there is faster and more accurate: for PREM, half the time and about 20 times closer to the converged k2. The solver warns once per session (`warn_if_internal_discontinuity` returns every such radius). Worlds built from data files are already split there.

## Choosing a Method

`love_method` selects how the Love numbers are obtained. Names are case-insensitive; the aliases in parentheses work everywhere.

| `love_method` | What it does |
|---|---|
| `radial_solver` (`shooting`, `rs`) | Default. Integrates the radial ODEs from the starting radius to the surface. Handles any multi-layer interior: solid or liquid, static or dynamic, compressible or incompressible. |
| `propagation_matrix` (`prop_matrix`, `pm`, `prop`) | Quasi-analytic, for a single solid, static, incompressible layer only. More sensitive to the slice count than the shooting method, since that sets the matrix dimension. Mainly for comparison. |

The analytic methods (`homogeneous`, `cpl`, `ctl`) need a layered world rather than moduli arrays: use `BaseWorld.solve_love_numbers(love_method=...)` or the closed forms in [`TidalPy.Tides.love`](../Tides/love/love_numbers.md).

## Degree-1 Load Love Numbers

At degree 1 a rigid translation of the body (y1 = y3 = a constant, y5 = that constant times g, no stress) satisfies the equations and every surface condition of a load, so load Love numbers are defined only in a chosen reference frame (Farrell 1972; Blewitt 2003). A degree-1 tidal or free-surface response is only a translation, so only `solve_for=('loading',)` takes `degree_l=1`.

The solver works in CE, the frame of the solid body's center of mass, where k' = 0. It replaces the last surface condition (y6, or y7 for a static liquid surface layer) with y5 = 1 at the surface (Guo et al. 2004 Eq. 9; Martens 2016 Eq. 4.150). In a static body the replaced condition then holds by itself, because the body feels no net force (Saito 1974 Eq. 25). `degree1_frame` moves the Love numbers and the radial functions together to another frame: y1 and y3 shift by a constant c, y5 by c times g(r), and the stresses do not change.

| `degree1_frame` | Origin | Defined by | Shift alpha from CE |
|---|---|---|---|
| `'CE'` (default) | Center of mass of the solid body | k' = 0 | 0 |
| `'CM'` | Center of mass of the body and the load | 1 + k' = 0 | 1 |
| `'CF'` | Center of surface figure | h' + 2 l' = 0 | (h' + 2 l') / 3 |
| `'CL'` | Center of lateral figure | l' = 0 | l' |
| `'CH'` | Center of height figure | h' = 0 | h' |

In each frame h', l', and 1 + k' are their CE values minus alpha (Blewitt 2003 Eq. 17). h' - k', l' - k', the strain, the stress, and the dissipation do not depend on the frame. Most published values are in CE; satellite laser ranging and GNSS use CM and CF. k' is 0 in CE and -1 in CM, so its quality factor and lag mean nothing there.

The test suite checks two references:

- A homogeneous sphere ($V_P$ = 10 km/s, $V_S$ = 5 km/s, $\rho$ = 5000 kg m$^{-3}$) gives h' = -0.206731 and l' = 0.161380, the exact value from the Takeuchi and Saito solutions (Martens 2016 Table C.4: -0.2069 and 0.1617). Every starting method, every tolerance from $10^{-6}$ to $10^{-12}$, and dynamic solves from $\omega = 10^{-4}$ to $10^{-11}$ rad s$^{-1}$ agree.
- The Guo et al. (2004) PREM gives h' = -0.28559 to -0.28564 and l' = 0.10372 to 0.10375 across static and dynamic layer choices, against their -0.285694 and 0.103633.

**Liquid surface layers.** A static liquid surface layer has no horizontal displacement: its l' is NaN, and the CF and CL frames, which read l', are refused with error code -16. Its hydrostatic surface fixes h' = 1 - $\bar{\rho} / \rho_s$ in CE exactly, whatever lies beneath ($\bar{\rho}$ is the bulk density, $\rho_s$ the surface density). An incompressible solid lid of that density tends to this value as its rigidity falls. A dynamic liquid surface layer solves like a solid one, with a large l', since an inviscid liquid's horizontal motion grows as $1 / \omega^2$. At long periods a dynamic compressible liquid fails as it does at any degree, or returns a wrong answer flagged by a large `surface_frame_residual` and a conditioning warning.

**Mixed static and dynamic layers.** A body with inertia in some layers but not others has no exact degree-1 frame: the replaced condition holds only to about $\omega^2 R / g$. `solution.surface_frame_residual` (`world.love_surface_frame_residual` on a world) reports the miss relative to the y6 condition, and the solver warns above $10^{-2}$. Dynamic solid layers under a static ocean leave $5 \times 10^{-3}$ at a one-day period. Making every layer static, or every layer dynamic, removes it.

The propagation matrix refuses degree 1 (error code -13). The `homogeneous`, `cpl`, and `ctl` methods give tidal Love numbers only, at every degree.

## Arguments

A setting whose default is `None` reads the `[radial_solver]` or `[eos_solver]` section of the [configuration](../Overview/2_TidalPy_Configurations.md), as a world's `solve_love_numbers` and `solve_eos` do. An explicit argument wins.

**Common arguments**

| Argument | Default | Meaning |
|---|---|---|
| `degree_l` | `2` | Spherical-harmonic degree. Expect trouble beyond `l = 10`; stable solutions exist to `l = 40` with a higher starting radius and tighter tolerances. |
| `solve_for` | `None` | Tuple of `'tidal'`, `'loading'`, `'free'`; `None` means `('tidal',)` (note the trailing comma). Several at once is cheaper than separate calls. Sets the first dimension of `solution.love`. |
| `love_method` | `'radial_solver'` | See [Choosing a Method](#choosing-a-method). |
| `nondimensionalize` | `None` (config) | Solve in non-dimensional units, return SI. Leave it on: the tolerances then mean the same for every planet. |

**Shooting method**

| Argument | Default | Meaning |
|---|---|---|
| `starting_radius` | `0.0` | Start of integration [m]; `0.0` chooses one (Martens 2016 criterion). A deep start at high degree makes the surface solve ill-conditioned, and the solver warns when accuracy drops below the requested tolerance. A manual radius above `[numerical] max_start_radius_fraction` of the planet radius (default 0.9) is refused: `ValueError` here, a failed solve on a world. |
| `start_radius_tolerance` | `None` (config) | The automatic start is at $R \cdot \mathrm{tol}^{1/l}$. |
| `starting_method` | `None` (config) | `'takeuchi'` (`'ts'`), `'kamata'`, `'power_series'` (`'powerseries'`, `'ps'`, `'martens'`), or `'unity'`. See [Starting Conditions](starting_conditions.md). |
| `degree1_frame` | `None` (config) | `'CE'`, `'CM'`, `'CF'`, `'CL'`, or `'CH'`; read only by a degree-1 loading solve. See [Degree-1 Load Love Numbers](#degree-1-load-love-numbers). |
| `integration_method` | `None` (config) | `'RK23'`, `'RK45'`, `'DOP853'`, or the implicit `'BDF'`, `'LSODA'`, `'Radau'` for stiff problems. |
| `integration_rtol`, `integration_atol` | `None` (config) | Integration tolerances. |
| `scale_rtols_bylayer_type` | `None` (config) | Scale the relative tolerance by layer type (liquids generally want it tighter). Experimental. |
| `max_num_steps` | `None` (config) | Step ceiling per integration. |
| `expected_size` | `None` (config) | Initial allocation hint for the integrator. Overshooting costs little. |
| `max_ram_MB` | `None` (config) | Integrator memory ceiling. Real usage runs somewhat higher. |
| `max_step` | `0` | Largest step [m]; `0` lets the integrator choose. |
| `love_only` | `False` | Keep only the Love numbers; faster. `result`, `get_radial_solution`, and `plot_ys` then raise `ValueError`. |

**Propagation matrix**

| Argument | Default | Meaning |
|---|---|---|
| `core_model` | `0` | Inner-core start: `0` Henning and Hurford (2014) seed matrix, `1` Roberts and Nimmo (2008) very small liquid core, `2` Henning and Hurford (2014) solid inner core, `3` Tobie et al. (2005) liquid inner core, `4` Sabadini and Vermeersen (2004) interface matrix. Only `0` is the regular solution; the others change k2 of a uniform body by about 3 (r_start / R)^3 (a few 1e-6 at the automatic start) and fail with a manual `starting_radius` (error code -22). |

**Equation of state** (solved before the deformation problem)

| Argument | Default | Meaning |
|---|---|---|
| `eos_method_bylayer` | `None` | Only `"interpolate"` (what `None` selects); other names raise. |
| `surface_pressure` | `0.0` | Surface pressure [Pa]. |
| `eos_integration_method` | `None` (config) | As `integration_method`. |
| `eos_rtol`, `eos_atol` | `None` (config) | EOS integration tolerances. |
| `eos_pressure_tol` | `None` (config) | Tolerance on the surface-pressure mismatch, relative to $(2/3) \pi G \rho^2 R^2$. Keep it above `eos_rtol`. |
| `eos_max_iters` | `None` (config) | Ceiling on the secant iterations for the central pressure; reported in `eos_iterations`. |

**Reporting**

| Argument | Default | Meaning |
|---|---|---|
| `verbose` | `False` | Print solver status. |
| `warnings` | `True` | Emit solver warnings. The diagnostics are recorded either way. |
| `raise_on_fail` | `False` | Raise instead of returning `success = False` with a `message`. |
| `perform_checks` | `True` | Ignored: inputs are always validated. |
| `log_info` | `False` | Log key diagnostics, at some cost. |

The radial-solver functions that take Newton's constant (the starting-condition functions, `apply_surface_bc`, `solve_upper_y_at_interface`, `fundamental_matrix`) read `G_to_use=None` as the configured value (SciPy's G).

### Keyword Names on a World

A world's `solve_love_numbers` and `solve_eos` take the shorter configuration-key names for the same settings.

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

## Troubleshooting

Read `solution.message`, then `solution.steps_taken`, then plot. These messages come from the shooting method.

**Slow solves, huge step counts, or "Maximum number of steps ... exceeded".** "(set by user)" means `max_num_steps` was hit, "(set by system architecture)" means `max_ram_MB`. Raising the limit is rarely the fix. A healthy solve takes a few hundred steps per solution per layer, a few thousand in awkward cases; ten thousand or more means it is likely unstable. `solution.plot_ys()` shows instability as spikes, ringing, or curves that are not smooth in radius. Try, roughly in order: other tolerances, another integration method, another `starting_method`, a lower `degree_l`, a higher `starting_radius`, other layer assumptions (dynamic compressible liquids especially), a small solid core under a fully liquid one, or more slices.

**"Error in step size calculation: Required step size is less than spacing between numbers."** The tolerances are usually too tight, or too loose so that error compounds. A dynamic compressible liquid layer is a common trigger; try it static.

**NaN Love numbers from a successful solve.** The surface solve was ill-conditioned; a `surface_solve_amplification` far above one means catastrophic cancellation. Raise the starting radius or use the automatic one.

**"The surface boundary condition system is singular to working precision" (error code -13).** `surface_solve_rcond` fell below `[numerical] minimum_surface_rcond` (default `1e-12`): the independent solutions became numerically dependent. Start higher or use the automatic radius. An rcond between that floor and the integration `rtol` is solved but logged as poorly conditioned. A very weak solid starting layer does this for every start; see [Starting Conditions](starting_conditions.md).

**"Layer ... is compressible but its bulk modulus is not positive" (error code -15).** A compressible solid or dynamic liquid needs a positive bulk modulus. An interpolated equation of state without a bulk-modulus table reports NaN, and a zero bulk modulus with a shear modulus (Poisson ratio -1) would silently give wrong Love numbers. Give the material a bulk modulus or mark the layer incompressible.

**"The degree-1 frame asked for reads l'" (error code -16).** `'CF'` and `'CL'` need l', which a static liquid surface layer lacks. Use `'CE'`, `'CM'`, or `'CH'`, or make the surface layer dynamic.

**"The surface condition its reference frame replaced is met only to ..." (warning).** A degree-1 loading solve with both static and dynamic layers; see [Degree-1 Load Love Numbers](#degree-1-load-love-numbers).

**"The starting radius ... is not inside the planet" or "No radial slice ... lies above the starting radius" (error code -5).** A manual starting radius needs at least one slice above it in its layer. A radius exactly on an interface starts in the layer above.

**A crash with no exception.** Rerun while watching memory. If it reproduces, open an issue on [GitHub](https://github.com/jrenaud90/TidalPy/issues) with the inputs.
