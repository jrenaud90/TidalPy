# Calculating Love Numbers

_Updated: 2026-09-12_

`TidalPy.RadialSolver_x.radial_solver` is the array-based entry point to the viscoelastic-gravitational solve. You hand it a radial grid with density and complex moduli on it, a forcing frequency, and a description of the layers; it returns a [`RadialSolverSolution`](solution_class.md) carrying the radial functions and the Love numbers. If you already have a built world, prefer `LayeredWorld.solve_love_numbers`, which fills these arrays from the layer rheologies for you.

## A first solve

```python
import numpy as np
from TidalPy.RadialSolver_x import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.rheology_x import Elastic, Maxwell

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

## What the input arrays must satisfy

- Every array must be C-contiguous and of the stated dtype. The radius, density, and moduli arrays share one length, the total number of slices.
- The radius array starts at `r = 0` and increases.
- Each layer needs at least 5 slices, so the arrays hold at least `5 * num_layers` points.
- Every interface radius appears **twice**: once as the top of the lower layer and once as the base of the upper layer. The last radius is the planet radius and must equal the last entry of `upper_radius_bylayer_array`.
- The solver runs its own equation of state to obtain gravity, pressure, mass, and moment of inertia. The current method interpolates the supplied profiles, so give a layer enough slices to resolve any property that varies within it. A layer with constant properties is fine with 5; a layer with a real density profile is not.

## Choosing a method

`love_method` selects how the Love numbers are obtained. Names are case-insensitive and the aliases in parentheses work everywhere.

| `love_method` | What it does |
|---|---|
| `radial_solver` (`shooting`, `rs`) | Default. Integrates the radial ODEs from the starting radius to the surface. Handles arbitrary multi-layer, solid or liquid, static or dynamic, compressible or incompressible interiors. |
| `propagation_matrix` (`prop_matrix`, `pm`, `prop`) | Quasi-analytic matrix propagation, valid only for a single solid, static, incompressible layer. It is more sensitive to the number of slices than the shooting method, since that sets the matrix dimension. Included mainly for comparison. |

The analytic methods (`homogeneous`, `cpl`, `ctl`) are not available here because this API takes moduli arrays rather than a layered world. Use `LayeredWorld.solve_love_numbers(love_method=...)` for those, or the closed-form functions in [`TidalPy.Tides_x.love`](../Tides_x/love/love_numbers.md). Passing one of them raises `ValueError` with that pointer.

## Arguments

The four arrays, the frequency, the bulk density, the three per-layer tuples, and the upper-radius array are positional and required. Everything below is optional, shown with its default.

**The problem**

| Argument | Default | Meaning |
|---|---|---|
| `degree_l` | `2` | Spherical-harmonic degree. Stability degrades as the degree rises: expect trouble beyond `l = 10`, though stable solutions exist to `l = 40` with a higher starting radius and tighter tolerances. |
| `solve_for` | `None` | Tuple of surface boundary conditions: `'tidal'`, `'loading'`, `'free'`. `None` means `('tidal',)`. Note the trailing comma in a one-element tuple. Solving several at once is cheaper than separate calls, and sets the first dimension of `solution.love`. |
| `love_method` | `'radial_solver'` | The radial technique; see the table above. |
| `nondimensionalize` | `True` | Non-dimensionalize internally and restore SI before returning. Generally more stable; leave it on. |

**Shooting method**

| Argument | Default | Meaning |
|---|---|---|
| `starting_radius` | `0.0` | Radius where integration begins [m]. `0.0` picks one automatically using the Martens (2016) criterion and `start_radius_tolerance`. Starting very deep at high degree makes the surface boundary solve ill-conditioned: the solution constants grow enormous and cancel, amplifying roundoff into the Love numbers. The solver measures this on every solve and warns when the achievable accuracy drops below the requested tolerance; prefer the automatic radius when that warning appears. |
| `start_radius_tolerance` | `1.0e-5` | Tolerance for that automatic choice. |
| `use_kamata` | `False` | Use the Kamata et al. (2015) starting conditions instead of Takeuchi and Saito (1972). Kamata is the more stable choice for incompressible layers, and is required for an incompressible solid layer at the center, where the Takeuchi and Saito form is undefined. It does not cover a static incompressible solid layer. |
| `integration_method` | `'DOP853'` | `'RK23'`, `'RK45'`, `'DOP853'`, or the implicit methods `'BDF'`, `'LSODA'`, `'Radau'` for stiff problems. |
| `integration_rtol`, `integration_atol` | `1.0e-5`, `1.0e-8` | Relative and absolute integration tolerances. |
| `scale_rtols_bylayer_type` | `False` | Scale the relative tolerance by layer type; liquid layers generally want a tighter value. Experimental. |
| `max_num_steps` | `500000` | Step ceiling per integration. A healthy solve needs a few hundred steps per solution per layer. |
| `expected_size` | `1000` | Hint for the integrator's initial allocation. Overshooting costs little. |
| `max_ram_MB` | `500` | Memory ceiling for the integrator. Real usage runs somewhat higher. |
| `max_step` | `0` | Largest allowed step [m]; `0` lets the integrator choose. |

**Propagation matrix**

| Argument | Default | Meaning |
|---|---|---|
| `core_model` | `0` | Inner-core starting condition. `0` Henning and Hurford (2014) seed matrix, `1` Roberts and Nimmo (2008) very small liquid core, `2` Henning and Hurford (2014) solid inner core, `3` Tobie et al. (2005) liquid inner core, `4` Sabadini and Vermeersen (2004) interface matrix. The choice matters more the higher you start. |

**Equation of state**

| Argument | Default | Meaning |
|---|---|---|
| `eos_method_bylayer` | `None` | Per-layer EOS method; `None` uses interpolation everywhere. |
| `surface_pressure` | `0.0` | Pressure at the surface [Pa], the outer boundary condition for the interior pressure solve. |
| `eos_integration_method` | `'DOP853'` | As `integration_method`, for the EOS solve. |
| `eos_rtol`, `eos_atol` | `1.0e-3`, `1.0e-5` | EOS integration tolerances. |
| `eos_pressure_tol` | `1.0e-3` | Convergence tolerance for the pressure iteration. |
| `eos_max_iters` | `40` | Iteration ceiling for that loop. |

**Reporting**

| Argument | Default | Meaning |
|---|---|---|
| `verbose` | `False` | Print solver status while running. |
| `warnings` | `True` | Emit solver warnings, including the surface-conditioning diagnostic. |
| `raise_on_fail` | `False` | Raise instead of failing quietly. By default a failed solve returns with `success = False` and an explanatory `message`. |
| `perform_checks` | `True` | Accepted for compatibility with the classic solver; the new backend validates unconditionally. |
| `log_info` | `False` | Log the solution's key diagnostics. There is a cost, more so with file logging enabled. |

## Troubleshooting

Start with `solution.message`, then `solution.steps_taken`, then plot. The messages below come from the shooting method.

**"Error in step size calculation: Required step size is less than spacing between numbers."** The integrator could not make progress. Usually the tolerances are too tight, or counterintuitively too loose so that error compounds. A dynamic (non-static) compressible liquid layer is a common trigger; try the static assumption for it.

**"Maximum number of steps (set by user) exceeded during integration."** The solve hit `max_num_steps`. Raising it is rarely the real fix: this normally means the solution is unstable, so work through the list below first.

**"Maximum number of steps (set by system architecture) exceeded during integration."** The integrator's arrays passed `max_ram_MB`. Same advice as above.

**Slow solves or huge step counts.** A healthy solve needs a few hundred steps per solution per layer; a few thousand happens in awkward cases; ten thousand or more means the solution is unstable. `solution.plot_ys()` is the quickest check, since instability shows up as large spikes, ringing, or curves that do not vary smoothly with radius. Things to try, roughly in order: change the integration tolerances, change the integration method, switch the starting condition with `use_kamata`, lower `degree_l`, start higher in the planet with `starting_radius`, revisit the layer assumptions (dynamic compressible liquid layers especially), add a small solid core beneath a fully liquid one, or add slices to the input arrays.

**NaN Love numbers from a successful solve.** The surface boundary condition solve was ill-conditioned. Check `solution.surface_solve_amplification`: values far above one mean the solution constants are cancelling catastrophically. Raise the starting radius, or let the solver choose it.

**A crash with no exception.** Rerun with the same inputs while watching memory. If it reproduces, record the inputs and open an issue on [GitHub](https://github.com/jrenaud90/TidalPy/issues).

## Migrating from the classic solver

The call signature is otherwise unchanged from `TidalPy.RadialSolver.radial_solver`, with three differences: the boolean `use_prop_matrix` is replaced by `love_method`, the input builders take `rheology_x` models (see [input builders](build_inputs.md)), and failures raise `ValueError` or `RuntimeError` rather than the classic `ArgumentException` and `UnknownModelError`. The solution object gains dense radial evaluation, the implicit integrators, and the conditioning diagnostic.
