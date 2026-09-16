# Calculating Love Numbers

_Updated: 2026-09-12_

`TidalPy.RadialSolver_x.radial_solver` is the array-based entry point to the viscoelastic-gravitational solve. You hand it a radial grid with density and complex moduli on it, a forcing frequency, and a description of the layers; it returns a [`RadialSolverSolution`](solution_class.md) carrying the radial functions and the Love numbers. If you already have a built world, prefer `LayeredWorld.solve_love_numbers`, which fills these arrays from the layer rheologies for you.

## Example

```python
import numpy as np
from TidalPy.RadialSolver_x import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.rheology_x import Elastic, Maxwell

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
```

The [input builders](build_inputs.md) exist because the solver's array requirements are strict. Build the arrays yourself only if you know they already satisfy the rules in the next section.

### Input Argument Rules

- Every array must be C-contiguous and of the stated dtype. The radius, density, and moduli arrays must all have the same size.
- The radius array starts at `r = 0` and strictly increases.
- Each layer needs at least 5 slices, so the arrays hold at least `5 * num_layers` points.
  - This entry point always builds an interpolated equation of state out of the arrays you supply: `eos_method_bylayer` accepts only `"interpolate"`, and any other name raises. There is no analytic-EOS option here that would make the slice count irrelevant, so how many slices a layer needs depends on its profiles, as the last two rules below explain.
- Every interface radius appears **twice**: once as the top of the lower layer and once as the base of the upper layer. The last radius is the planet radius and must equal the last entry of `upper_radius_bylayer_array`.
- The solver runs its own equation of state, and the slice count does not affect everything it produces equally. Gravity, pressure, mass, and moment of inertia are integrated and read back through the integrator's dense output, so they are evaluated at the exact integration radius rather than interpolated between your slices. Density and the complex shear and bulk moduli are the ones the slices control: at every integration radius they are linearly interpolated from the arrays you supplied, and that interpolated density is also what the structural integration above is driven by.
- What that costs is set by how *curved* those profiles are within a layer, not by whether they vary at all, because linear interpolation reproduces a straight line exactly. For a single solid layer at degree 2, k2 from 5 slices lands within 1e-9 of the converged value when density and the moduli are constant and within 1e-8 when they vary linearly, but a strongly curved profile (density falling exponentially by a factor of two across the layer) is off by 1.3% at 5 slices, 0.2% at 10, and 0.03% at 25. Give a layer enough slices to follow the curvature of whichever property carries it, and test the count against a refined run for your own problem.

## Choosing a Method

`love_method` selects how the Love numbers are obtained. Names are case-insensitive and the aliases in parentheses work everywhere.

| `love_method` | What it does |
|---|---|
| `radial_solver` (`shooting`, `rs`) | Default. Integrates the radial ODEs from the starting radius to the surface. Handles arbitrary multi-layer, solid or liquid, static or dynamic, compressible or incompressible interiors. |
| `propagation_matrix` (`prop_matrix`, `pm`, `prop`) | Quasi-analytic matrix propagation, valid only for a single solid, static, incompressible layer. It is more sensitive to the number of slices than the shooting method, since that sets the matrix dimension. Included mainly for comparison. |

The analytic methods (`homogeneous`, `cpl`, `ctl`) are not available here because this API takes moduli arrays rather than a layered world. Use `LayeredWorld.solve_love_numbers(love_method=...)` for those, or the closed-form functions in [`TidalPy.Tides_x.love`](../Tides_x/love/love_numbers.md). Passing one of them raises `ValueError` with that pointer.

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
| `starting_radius` | `0.0` | Radius where integration begins [m]. `0.0` picks one automatically using the Martens (2016) criterion and `start_radius_tolerance`. Starting very deep at high degree makes the surface boundary solve ill-conditioned: the solution constants grow enormous and cancel, amplifying Love numbers error. The solver measures this on every solve and warns when the achievable accuracy drops below the requested tolerance; prefer the automatic radius when that warning appears. |
| `start_radius_tolerance` | `None` (config) | Tolerance for that automatic choice: the start is at $R \cdot \mathrm{tol}^{1/l}$. |
| `use_kamata` | `None` (config) | Use the Kamata et al. (2015) starting conditions instead of Takeuchi and Saito (1972). Kamata is the more stable choice for incompressible layers, and is required for an incompressible solid layer at the center, where the Takeuchi and Saito form is undefined. It does not cover a static incompressible solid layer. |
| `integration_method` | `None` (config) | `'RK23'`, `'RK45'`, `'DOP853'`, or the implicit methods `'BDF'`, `'LSODA'`, `'Radau'` for stiff problems. |
| `integration_rtol`, `integration_atol` | `None` (config) | Relative and absolute integration tolerances. |
| `scale_rtols_bylayer_type` | `None` (config) | Scale the relative tolerance by layer type; liquid layers generally want a tighter value. Experimental. |
| `max_num_steps` | `None` (config) | Step ceiling per integration. A healthy solve needs a few hundred steps per solution per layer. |
| `expected_size` | `None` (config) | Hint for the integrator's initial allocation. Overshooting costs little. |
| `max_ram_MB` | `None` (config) | Memory ceiling for the integrator. Real usage runs somewhat higher. |
| `max_step` | `0` | Largest allowed step [m]; `0` lets the integrator choose. |

**Propagation matrix**

| Argument | Default | Meaning |
|---|---|---|
| `core_model` | `0` | Inner-core starting condition. `0` Henning and Hurford (2014) seed matrix, `1` Roberts and Nimmo (2008) very small liquid core, `2` Henning and Hurford (2014) solid inner core, `3` Tobie et al. (2005) liquid inner core, `4` Sabadini and Vermeersen (2004) interface matrix. The choice matters more the higher you start. |

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
| `warnings` | `True` | Emit solver warnings, including the surface-conditioning diagnostic. The diagnostic itself (`surface_solve_amplification`) is recorded on every solve. |
| `raise_on_fail` | `False` | Raise instead of failing quietly. By default a failed solve returns with `success = False` and an explanatory `message`. |
| `perform_checks` | `True` | Accepted for compatibility with the classic solver; the new backend validates unconditionally. |
| `log_info` | `False` | Log the solution's key diagnostics. There is a cost, more so with file logging enabled. |

## Troubleshooting

Start with `solution.message`, then `solution.steps_taken`, then plot. The messages below come from the shooting method.

**"Maximum number of steps (set by user) exceeded during integration."** The solve hit `max_num_steps`. Raising it is rarely the real fix: this normally means the solution is unstable, so work through the list below first.

**"Maximum number of steps (set by system architecture) exceeded during integration."** The integrator's arrays passed `max_ram_MB`. Same advice as above.

**"Error in step size calculation: Required step size is less than spacing between numbers."** The integrator could not make progress. Usually the tolerances are too tight, or counterintuitively too loose so that error compounds. A dynamic (non-static) compressible liquid layer is a common trigger; try the static assumption for it.

**Slow solves or huge step counts.** A healthy solve needs a few hundred steps per solution per layer; a few thousand happens in awkward cases; ten thousand or more means the solution is likely unstable. `solution.plot_ys()` is the quickest check, since instability shows up as large spikes, ringing, or curves that do not vary smoothly with radius. Things to try, roughly in order: change the integration tolerances, change the integration method, switch the starting condition with `use_kamata`, lower `degree_l`, start higher in the planet with `starting_radius`, revisit the layer assumptions (dynamic compressible liquid layers especially), add a small solid core beneath a fully liquid one, or add more slices to the input arrays.

**NaN Love numbers from a successful solve.** The surface boundary condition solve was ill-conditioned. Check `solution.surface_solve_amplification`: values far above one mean the solution constants are cancelling catastrophically. Raise the starting radius, or let the solver choose it.

**A crash with no exception.** Rerun with the same inputs while watching memory. If it reproduces, record the inputs and open an issue on [GitHub](https://github.com/jrenaud90/TidalPy/issues).
