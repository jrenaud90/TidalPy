# Dense Radial Solutions

_Updated: 2026-10-10_

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

A layer's independent solutions are integrated together as one system with dense (continuous) output, so each evaluation of the equations reads the layer's material once for all of them. Integrated outward, the solutions turn toward the fastest growing one, and the combinations the surface conditions need then cancel. Where their normalized Gram determinant (1 when they are orthogonal, 0 when they are dependent) has fallen by `[numerical] minimum_solution_independence` (default `1e-4`) from where the integration started, a terminal event stops the integration, the solutions are replaced by an orthonormal basis of their span (a QR factorization; Godunov 1961, Conte 1966), and a new integration continues from it. `RadialSolverSolution.orthonormalizations` counts these restarts per layer. The determinant and the basis are taken with each radial function divided by its characteristic size in the non-dimensional units (the displacement unit $1/(\pi G \bar{\rho} R)$, the stress unit $\bar{\rho}$, and 1 or $1/R$ for the potential terms), so a solve in SI units restarts where a non-dimensional one does, and the surface solve's `surface_solve_amplification` is measured in the same sizes. The equations are linear, so the solutions before a restart equal those after it times the triangular factor $R$ at every radius, and each segment's constants are $R^{-1}$ times those of the segment above it. The solutions a layer starts from are used as given, since near the center the starting conditions can be nearly dependent in their leading powers of $r$ and orthonormalizing them would only make the integrator resolve components the surface solve does not need. Solutions dependent to working precision, as from a start within meters of the center, are orthonormalized before the integration starts.

A value at any radius is the containing segment's interpolant evaluated there, multiplied by that segment's constants and summed, with `y3` of a dynamic liquid rebuilt from its pressure variable (below), and returned in SI. A value between grid slices therefore comes from the integrator's own polynomial rather than from a linear interpolation of a sampled grid.

## Propagation Matrix Method

The propagation matrix builds its solution on a grid by propagating each slice's fundamental matrix from one slice radius to the next. `get_radial_solution(r)` continues that propagation to the radius asked for: the slice's fundamental matrix at `r` applied to the coefficients carried into the slice. This reproduces the grid values at both ends of the slice and is exact between them for a uniform layer. Below the first propagated slice, the regular solution of the innermost material continues the solution to the center, so `get_radial_solution` and `result` are defined down to `r = 0`. With a `core_model` other than `0` the seed is not the regular solution, and the solution stays NaN below the starting radius.

## Dynamic Liquid Layers at Long Forcing Periods

In the variables of Takeuchi and Saito (1972), a dynamic liquid's $y_1$ and $y_2$ are coupled through $\ell(\ell+1)/(\omega^2 r)$. That coupling gives the equations a pair of rates near $\pm\sqrt{\ell(\ell+1) \rho g^2 / K}/(\omega r)$, which cancel in the physical solution but set an explicit integrator's step, so the step shrinks with the forcing period and the independent solutions lose their independence on the way to the surface. TidalPy therefore integrates a dynamic liquid in the pressure variable

$$P = y_2 - \rho g\, y_1 + \rho\, y_5,$$

which is minus the Eulerian pressure perturbation plus $\rho$ times the potential perturbation, in place of $y_2$. Then $y_3 = -P / (\rho \omega^2 r)$, and with $dg/dr = 4 \pi G \rho - 2 g / r$,

$$\frac{dP}{dr} = -\omega^2 \rho\, y_1 + \left(\rho' + \frac{\rho^2 g}{K}\right) \left(y_5 - g\, y_1\right) - \frac{\rho g}{K} P.$$

The large coupling now enters through $P$ alone, and the remaining rates are the buoyancy ones, $\sqrt{-\ell(\ell+1) N^2}/(\omega r)$, with $N^2 = -g\, (\rho'/\rho + \rho g / K)$ and $\rho' = d\rho/dr$ ($1/K = 0$ for an incompressible liquid). The change of variables is exact only when $\rho'$ is the derivative of the density the equations read. The solved interior reports it through `buoyancy_frequency_squared`, the $N^2$ of its own bulk modulus, which a compressible liquid reads, and through `density_gradient`, which an incompressible liquid reads (exactly 0 for a constant density). The chain rule gives both through the material's density with the structure's $dP/dr = -\rho g$ and $dT/dr$ (the table's slope for a density tabulated in radius, including at the table's own end radii). That $N^2$ is formed from the density's partial derivatives, so a material whose density follows its bulk modulus gives exactly 0. The difference $\rho' + \rho^2 g / K$ of two large terms would otherwise leave a roundoff $N^2$ of about $10^{-16} \rho g^2 / K$, a buoyancy period near $10^7$ days in Europa's ocean, which would set the solve at longer periods. The solver converts $y_2$ to $P$ where a liquid layer starts and back at its top, so the interface conditions, the surface solve, and every reported radial function use $y_2$.

A physical solution's $P = -\rho \omega^2 r\, y_3$ is of order $\omega^2$ in stress units, while the layer's independent solutions carry stress-sized $P$. Stored in stress units, the physical $P$ would be the cancellation of theirs, and $y_3 = -P / (\rho \omega^2 r)$, and through it $y_1$ and the Love numbers, would carry that cancellation's roundoff divided by $\omega^2$. The solver therefore stores $P$ in units of $\omega^2 / (\pi G \bar{\rho})$ of a stress below the frequency unit $\sqrt{\pi G \bar{\rho}}$, where a physical $P$ is of order $\rho r y_3$ and the integration tolerances apply to it. Where a dynamic liquid starts, its two solutions are recombined into one whose $P$ is exactly 0 and one that carries $P$, before the integration divides any roundoff in $P$ by $\omega^2$. In a stratified liquid the solutions carry $P$ about $N / \omega$ times $y_1$ in that unit, and the independence measures $P$ at that size.

The stratification then sets the cost of a long-period solve:

* A neutral liquid, whose density follows its bulk modulus ($\rho' = -\rho^2 g / K$), has $N^2 = 0$. The bundled `luna_dynamic`, `mercury_dynamic`, `pluto_dynamic`, and `europa_dynamic` worlds are built this way, on pressure-dependent laws. Their liquids take 4 to 9 steps at every period down to the $10^{-16}$ rad s$^{-1}$ of `[numerical] minimum_frequency`. At the default tolerances their $k_2$ agrees with `rtol = 1e-12` solves to 2e-9 or better over that whole range, and from $10^{-12}$ to $10^{-16}$ rad s$^{-1}$ it agrees with the static liquid's, the long-period limit, to 2e-9 in non-dimensional and SI solves alike. The radial functions inside the liquid hold the tolerance too, with Europa's ocean $y_1$ to $y_6$ within 1e-11 of an `rtol = 1e-12` solve at $10^{-10}$ rad s$^{-1}$.
* A stably stratified liquid ($N^2 > 0$, such as an incompressible liquid whose density rises with depth) carries internal gravity waves whose radial wavelength shrinks with the period, so its steps grow in proportion to the period. Europa's ocean made incompressible takes about 330 steps at 100 days and 33,000 at $10^4$ days. A stably stratified liquid at the surface itself, forced far below its buoyancy frequency, has a Shida $l$ of order $10^5$ whose error follows the integration `rtol` (about $10^{-4}$ at $10^{-8}$ rad/s at the default tolerance) while $k$ and $h$ stay accurate; tighten `rtol` for $l$ there.
* An unstably stratified liquid ($N^2 < 0$, such as a constant-density compressible liquid, where $N^2 = -\rho g^2 / K$) has solutions that grow through the layer as $e^E$, with $E = (\sqrt{\ell(\ell+1)}/\omega) \int \sqrt{-N^2}/r\, dr$. Re-orthonormalization keeps the solve accurate whatever $E$, but the steps grow in proportion to the period here too. Earth-Simple's constant-density core solved dynamically takes about 900 steps at 10 days, 9,300 at 100 days, and 93,000 at 1000 days, and reaches `max_num_steps` a few times further out. Every one of those steps and restarts is kept for the dense solution (about 200 MB at 1000 days and 600 MB at 3000 days), and `max_ram_MB` bounds that too. Its $k_2$ approaches the static liquid's at long periods.

A density tabulated in radius is linear between rows, so its $\rho'$ steps at each row, and the $P$ equation steps with it. PREM's outer core solved dynamically at a half-day period takes about 290 steps where TS72's $y_2$ form, which reads no $\rho'$, takes 124, and its $k_2$ is 2.2e-7 from an `rtol = 1e-12` solve where the $y_2$ form's is 1.6e-6. The table also leaves the core weakly stratified ($|N^2|$ up to about $6 \times 10^{-8}$ s$^{-2}$, of either sign), so its steps grow with the period beyond about $10^{-6}$ rad s$^{-1}$, and at $10^{-10}$ rad s$^{-1}$ the solve reaches `max_ram_MB`. TS72's $y_2$ form gives $k_2$ = 0.91 at $10^{-6}$ rad s$^{-1}$ (0.66 at `rtol = 1e-12`, against 0.2984) and fails from $10^{-7}$ rad s$^{-1}$.

Guidance:

* A static liquid layer is the long-period limit of a neutral dynamic one and is the cheapest choice at long periods for any density profile.
* To solve a liquid dynamically at long periods, give it a pressure-dependent law (Birch-Murnaghan or Vinet), which keeps it neutral.
* A start inside a compressible dynamic liquid at a long period fails with the Takeuchi and Saito, Kamata, and power series starts, which treat the liquid as a homogeneous constant-density sphere whose solutions overflow there (in `luna_dynamic`'s outer core at degree 5, from about 200 days). The automatic starting radius falls inside a liquid core at high degrees in a body with a small inner core (`luna_dynamic` from degree 5) or none. Start in the layer below with a manual `starting_radius`, or use `starting_method="unity"`, which solves there (see [Starting Conditions](starting_conditions.md#accuracy-and-cost)).

## C++ API

`c_RadialSolutionStorage` provides `get_radial_solution`, `get_radial_solution_array`, `get_surface_y`, `get_eos_si`, and `get_radial_solution_nondim` (for a radius in solve units). It keeps each layer's integration segments (`c_ShootingSegment`: the dense CyRK result, the radius span, and the basis change at its start) and their collapse constants. A built `c_BaseWorld` provides `get_radial_solution_y` and `get_love_surface_y`. `c_EOSSolution::call_material` returns a `c_EOSMaterialState` (gravity, density, the buoyancy frequency squared and the static bulk modulus it is measured against, and both complex moduli) whatever the solution was built from. The equations of each layer kind live in `RadialSolver/derivatives/odes_.hpp` and the re-orthonormalization in `RadialSolver/orthonormalize_.hpp`.

## Validation

`Tests/Test_RadialSolver/test_dense_benchmark/` holds frozen reference results for a 1-layer solid planet, a 2-layer solid planet, and a solid / dynamic-liquid / solid planet at a short forcing period. The dense solution agrees closely with them at the surface and at layer boundaries, and the Love numbers agree tightly; between grid slices the dense evaluation is the more accurate. `Tests/Test_RadialSolver/test_long_period_01.py` checks the pressure form against the $y_2$ form at short periods (to 1e-8), the long-period Love numbers against an independent model integrated with re-orthonormalization after every step (to 1e-8 at $10^3$ and $10^4$ days), and the dense radial functions across restarts.
