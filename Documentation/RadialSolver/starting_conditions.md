# Starting Conditions (`RadialSolver.starting`)

_Updated: 2026-10-07_

The shooting method integrates its independent solutions outward from a starting radius near the center (three in a solid, two in a dynamic liquid, one in a static liquid). The starting conditions are the values of those solutions at that radius. Every method below starts from solutions that are regular at the center of a homogeneous layer (constant density and moduli, gravity $g = \gamma r$ with $\gamma = 4 \pi G \rho / 3$), except unity, which makes no assumption. The `starting_method` argument chooses between them. It is taken by `radial_solver`, by `BaseWorld.solve_love_numbers`, and by the `[radial_solver]` configuration and world tables. The default is `"takeuchi"`.

| `starting_method` | Aliases | Form |
|---|---|---|
| `"takeuchi"` | `"ts"` | Takeuchi and Saito (1972) closed forms in spherical Bessel functions |
| `"kamata"` | | Kamata et al. (2015) closed forms |
| `"power_series"` | `"powerseries"`, `"ps"`, `"martens"` | Martens (2016) power series in $r^2$, summed to convergence |
| `"unity"` | | Unit vectors |

Every method covers every layer type: solid or liquid, static or dynamic, compressible or incompressible. Some of the forms are not in the papers and were derived for TidalPy (see each method below). Static liquids start from Saito's (1974) solution, $y_5 = r^l$ and $y_7 = 2 (l - 1) r^{l - 1}$, for every method but unity. That solution is exact for a homogeneous liquid, and the power series of the static liquid equations ends after its first term with the same values.

## Methods

### Takeuchi and Saito

`TidalPy.RadialSolver.starting.takeuchi` evaluates Eqs. 95-102 of Takeuchi and Saito (1972) for compressible solids and dynamic liquids. Two of a solid's solutions are spherical Bessel functions of the layer's two wavenumbers, and the third is a polynomial.

Takeuchi and Saito give no incompressible forms. TidalPy takes the limit of their solutions as the bulk modulus grows without bound ($1 / \alpha^2 \to 0$):

- One wavenumber tends to the shear wave's, $k^2 = \omega^2 / \beta^2$, and its solution keeps a finite limit.
- The other tends to zero, and its solution, rescaled, tends to an exact pressure-potential solution, $(y_1, \ldots, y_6) = (0, -\rho r^l, 0, 0, r^l, (2l + 1) r^{l - 1})$.
- The polynomial solution is unchanged, since its displacement has no divergence.

As $\omega \to 0$ the shear-wave solution converges onto the pressure-potential one, so it is replaced by its difference from it, written so that no division by $\omega^2$ remains. The static incompressible solid is that form at $\omega = 0$. The incompressible liquid's two solutions are the pressure-potential solution, $(y_1, y_2, y_5, y_6) = (0, -\rho r^l, r^l, (2l + 1) r^{l - 1})$, and the polynomial solution, both exact. The derivation is written out in `takeuchi_.hpp`.

At large $\omega r / \beta$ the shear-wave solution again nearly coincides with the pressure-potential one (at degree 8, a normalized smallest singular value of $4 \times 10^{-5}$ at $\omega r / \beta = 20$), as the published compressible forms do at large arguments. Kamata's normalization stays independent there.

### Kamata

`TidalPy.RadialSolver.starting.kamata` evaluates Kamata et al. (2015) Eqs. B1-B37, the same solutions as Takeuchi and Saito with incompressible forms added. For a dynamic incompressible solid, the first solution is replaced by its difference from the second times $\gamma / \omega^2$. The published pair converges as $\omega \to 0$, while the combination stays independent at tidal periods. That combination is finite at $\omega = 0$, and the static equations are the dynamic ones there, so the static incompressible solid, which Kamata et al. do not give, is that form at $\omega = 0$.

Kamata's forms carry a common factor of $r^{-l}$, so at the starting radius their vectors are orders of magnitude larger than the other methods'. This matters for cost and accuracy (see [Size of the Starting Vectors](#size-of-the-starting-vectors)).

### Power Series

`TidalPy.RadialSolver.starting.power_series` expands the solutions about the center (Crossley 1975; Smylie 2013; Martens 2016, Sec. 4.2.8; Martens et al. 2019). Martens defines $y_1$ to $y_5$ as TidalPy does (Takeuchi and Saito 1972), but her $y_6$ leaves out a term:

$$y_6^{\mathrm{TidalPy}} = y_6^{\mathrm{Martens}} + \frac{l + 1}{r} y_5$$

Her series is converted with this relation. With $y_i = r^{l - 2 + s_i + \nu_0} z_i$, every layer's equations become $r \, dz/dr = (B_0 + B_2 r^2) z$. A regular solution is then $z = \sum_k c_k r^{2k}$, with $(\nu_0 I - B_0) c_0 = 0$ and $((\nu_0 + 2k) I - B_0) c_k = B_2 c_{k - 1}$. The scalings are $s = (1, 0, 1, 0, 2, 1)$ for a solid's $y_1$ to $y_6$ and $s = (1, 2, 2, 1)$ for a dynamic liquid's $y_1$, $y_2$, $y_5$, $y_6$. The implementation differs from Martens' series in a few ways:

- Terms are added until the next one falls below machine precision. Martens stops at three terms ($r^0$, $r^2$, $r^4$). Summing to convergence keeps the series accurate at the larger starting radii of high degrees.
- Solids of every kind and dynamic liquids have a series, derived from TidalPy's own equations. Martens covers the compressible solid. Static liquids use Saito's solution, the first and only term of their series.
- The first solution is Takeuchi and Saito's polynomial solution, a combination of Martens' $A_{1,1}$ and $A_{6,1}$ solutions that is exact as its first term and written without a series. In an incompressible solid the $A_{6,1}$ solution is the exact pressure-potential solution and is also written without a series.
- In a compressible solid with gravity, the other two solutions are the series of Takeuchi and Saito's two wavenumber solutions. Each starts as $(\alpha^2 f - (l + 1) \beta^2)$ times Martens' $A_{6,1}$ vector, and the free part of its resonant $r^2$ step is set to Takeuchi and Saito's $z_4$. Each solution then carries one wavenumber. Elsewhere (incompressible solids, which have one wavenumber, and solids without gravity) the remaining solutions are Martens' $A_{6,1}$ and $A_{4,0}$ solutions, with the free part of the resonant step set to zero, as Martens sets $A_{4,2} = 0$.
- The recurrence is scaled by the layer's shear modulus and by $\gamma + \omega^2$. Its linear solves and its convergence test therefore hold in SI as well as in non-dimensional units.

The series refuses to start, with a message naming the closed-form starts, in four cases:

- It needs more than 100 terms.
- Cancellation would take more than half the digits of any component.
- A solution grows from its first terms by more than $1/\sqrt{\epsilon}$, which puts the start deep in a layer's exponential regime, such as a dynamic liquid at long periods (see [Dynamic Liquid Layers at Long Forcing Periods](dense_radial_solution.md#dynamic-liquid-layers-at-long-forcing-periods)).
- In a solid, the larger wavenumber times the starting radius exceeds $\tfrac{1}{2} \ln(1/\sqrt{\epsilon}) \approx 9$. The wavenumber solutions are summed alongside each other's roundoff, and the slower one's error grows as $\epsilon e^{2 |k| r}$. This is the case of a weak solid starting layer.

The closed forms carry that growth analytically and apply in all four cases. Where the series does start, its solutions match the Takeuchi and Saito solutions one by one to better than $10^{-6}$ in solids with $|\mu|$ down to $10^{-6} \rho g R$ (degrees 2 to 10). Below that, an accepted static start can differ by up to $10^{-5}$. At degree 1 a static solid has a rigid-translation solution, so single solutions are not defined there.

### Unity

The unit vectors are $y_1$, $y_4$, and $y_6$ in a solid, $y_1$ and $y_6$ in a dynamic liquid, and $y_7$ in a static liquid. These are the leading components of Martens' free-constant solutions. In a comparison of every set of unit vectors on homogeneous bodies, this set's errors were within a factor of two of the smallest. Sets that include $y_3$ can leave the regular solutions degenerate.

A unit vector holds singular content as well as regular. Integrating outward, the singular part decays relative to the regular part as $(r_0 / r)^{2l - 1}$ in a solid and $(r_0 / r)^{2l + 1}$ in a liquid, where $r_0$ is the starting radius. Its error at the surface therefore depends on the starting radius rather than on the integration tolerance. Unity is a check on the other methods, not a default. It is unreliable in a weak solid starting layer, with errors of $10^{-2}$ to $5$ in $k$ in the comparison below.

## Choosing a Method

- `"takeuchi"`, the default, is accurate wherever the starting layer is not a weak solid. Only the power series matched its cost in the comparison below.
- `"kamata"` had median errors equal to or below Takeuchi and Saito's in every comparison below: 3 times below over the bundled worlds and $10^{2}$ to $10^{5}$ times below in weak solid starting layers. It cost 1.0 to 3.1 times as much. Use it when the starting layer is a weak solid (a viscoelastic $|\mu^*|$ far below $\rho g R$, as in a low-viscosity Maxwell mantle at tidal periods), or when the extra digits matter.
- `"power_series"` costs about the same as Takeuchi and Saito and gives the same Love numbers to the integration's accuracy where it starts. It needs no special functions and, with the Takeuchi and Saito incompressible forms, it stays finite without gravity, where the compressible Takeuchi and Saito forms and every Kamata form divide by $\gamma$. It refuses rather than start in a weak solid or deep in an exponential regime.
- `"unity"` makes no assumption about the layer. Use it to test whether a result depends on the start.

> [!NOTE]
> Only the layer the integration starts in uses these conditions. Every layer above it starts from the interface conditions with the layer below.

> [!WARNING]
> In a very weak solid starting layer ($|\mu^*|$ near $10^{-4} \rho g R$ and below), the solve is ill-conditioned for every start. The Takeuchi and Saito and power series solves can then succeed with $k$ wrong by order one at the default tolerances: `surface_solve_rcond` falls to $10^{-13}$, just above the `[numerical] minimum_surface_rcond` refusal of $10^{-14}$. Check `surface_solve_rcond`, use `"kamata"`, and tighten the tolerances.

## Accuracy and Cost

These were one-time measurements, made on 2026-10-07 with TidalPy 0.8.0 on one machine. The scripts are not part of the repository.

### Bundled Worlds

The comparison covered every bundled world under five interior conditions:

- as bundled, which mostly means static layers;
- dynamic solids;
- every layer dynamic;
- every layer dynamic with incompressible liquids;
- incompressible solids.

It used degrees 2, 3, 5, and 10 and forcing periods of 0.5, 3, and 30 days, for 1367 tidal solves per method at the packaged `[radial_solver]` settings (DOP853, `rtol = atol = 3e-8`, `start_radius_tolerance = 1e-5`). The error in $k$ is relative to solves at `rtol = atol = 1e-12`. It is counted only where at least two of the Takeuchi, Kamata, and power series references succeeded and agree to $10^{-6}$ (1334 cases). Cost is the best of five solve times relative to the Takeuchi and Saito solve of the same case. The median Takeuchi and Saito solve took 0.36 ms.

| Method | Solves that succeed | Median error in $k$ | 90th percentile | Error above $10^{-5}$ | Median cost |
|---|---|---|---|---|---|
| `takeuchi` | 1360 | 5.5e-9 | 1.4e-6 | 4.1% | 1.00 |
| `kamata` | 1365 | 1.7e-9 | 7.7e-7 | 2.6% | 1.44 |
| `power_series` | 1355 | 9.6e-9 | 1.8e-6 | 3.5% | 1.04 |
| `unity` | 1365 | 1.6e-7 | 4.3e-5 | 16.2% | 1.38 |

Every method failed a few solves in dynamic liquids at long periods, where the solutions grow exponentially. The power series also refused 11 such starts.

### Synthetic Bodies

These were solved with `radial_solver` (R = 6371 km) at degrees 2, 3, 5, 10, and 20, with references at `rtol = 1e-12`, `atol = 1e-16`:

- homogeneous solids (static and dynamic, compressible and incompressible);
- a Maxwell solid ($\mu$ = 60 GPa, $\eta = 10^{9}$ to $10^{13}$ Pa s, at a 1.8 day period) under a 171 km elastic lid;
- a liquid core under a solid mantle (static, dynamic compressible, dynamic incompressible; 1 and 10 day periods).

The table gives the solves that succeed, the median error in $k$, and the median cost relative to Takeuchi and Saito.

| Body | `takeuchi` | `kamata` | `power_series` | `unity` |
|---|---|---|---|---|
| Homogeneous solid | 20/20, 2e-9, 1.00 | 20/20, 2e-9, 1.35 | 20/20, 2e-9, 1.04 | 20/20, 2e-7, 1.32 |
| Maxwell start, $\eta = 10^{12}$ to $10^{13}$ | 10/10, 2e-8, 1.00 | 10/10, 5e-11, 3.06 | 7/10, 2e-9, 1.08 | 10/10, 5e-7, 1.50 |
| Maxwell start, $\eta = 10^{11}$ | 5/5, 6e-7, 1.00 | 5/5, 8e-12, 2.64 | 3/5, 3e-9, 1.01 | 5/5, 3e-2, 1.44 |
| Maxwell start, $\eta = 10^{10}$ | 5/5, 6e-1, 1.00 | 5/5, 8e-5, 1.88 | 3/5, 3e-1, 1.03 | 5/5, 5e0, 1.38 |
| Liquid core, static | 10/10, 9e-9, 1.00 | 10/10, 9e-9, 1.00 | 10/10, 9e-9, 0.97 | 10/10, 2e-8, 1.14 |
| Liquid core, dynamic incompressible | 10/10, 7e-8, 1.00 | 10/10, 5e-9, 3.08 | 10/10, 7e-8, 0.99 | 10/10, 2e-8, 1.26 |
| Liquid core, dynamic compressible | 6/10, 7e-8, 1.00 | 6/10, 5e-8, 1.20 | 4/10, 7e-8, 1.12 | 6/10, 2e-8, 1.25 |

At $\eta = 10^{9}$ no two references agreed, so those solves are not scored. The dynamic compressible liquid core fails at a 10 day period for every method, since its solutions grow by many orders of magnitude. The power series' missing Maxwell solves are refusals past its wavenumber bound (degrees 10 and 20, and degree 5 at the weakest viscosity).

### Cost of the Start

Through `find_starting_conditions`, with a Python call overhead of about 0.3 µs included, a start took:

- 0.3 to 1.1 µs for the closed forms and unity;
- for the power series: 13 µs in a compressible solid, 2.5 to 3.4 µs in an incompressible solid, 7.5 µs in a compressible liquid, and 0.9 µs in an incompressible liquid.

Against a 0.2 to 0.5 ms solve, the start is a small fraction of the cost. The cost differences in the tables come from the integration.

### Size of the Starting Vectors

The methods start from the same solution spaces, but their vectors differ in size by orders of magnitude at the starting radius. Kamata's are the largest (by the factor $r^{-l}$). The integrator measures its absolute tolerance against those sizes, so Kamata's start is integrated to a tighter effective tolerance: that is its extra accuracy and its extra cost. Scaling every method's start to unit size equalized them: Takeuchi and Saito and the power series became 16 to 18 percent slower and more accurate, and Kamata 21 to 28 percent faster and less accurate. In the $\eta = 10^{10}$ Maxwell start, Takeuchi and Saito reached Kamata's value once `atol` was lowered to $10^{-14}$. Tightening `atol` is therefore an alternative to switching methods.

### Validation

- Martens' series as written (three terms, her variables), ported separately from TidalPy, matches LoadDef's `fundamental_solutions_powser` to $10^{-14}$. Integrated with her equations, it reproduces her homogeneous-Earth Love numbers (Martens 2016, Table C.4: $V_P$ = 10 km/s, $V_S$ = 5 km/s, $\rho$ = 5000 kg m$^{-3}$, at the M2 period with $G = 6.672 \times 10^{-11}$) to the table's rounding for $n$ = 2 to 100.
- In TidalPy's variables the series spans the Takeuchi and Saito and the Kamata solutions to better than $10^{-9}$ (about $10^{-11}$ in non-dimensional units) and reproduces Table C.4 through TidalPy's equations. Each of its compressible-solid wavenumber solutions matches the corresponding Takeuchi and Saito solution.
- The derived Takeuchi and Saito incompressible forms and Kamata's static incompressible solid satisfy TidalPy's equations, span the power series' solutions to about $10^{-11}$, and agree with Kamata's dynamic incompressible forms. Every regular start reproduces Love's (1911) closed form for a static incompressible solid, $k_2 = (3/2) / (1 + 19 \mu / (2 \rho g R))$.
- On `earth_prem`, all four methods give the same tidal and load Love numbers to four decimals for degrees 2 to 10.

`Tests/Test_RadialSolver/test_starting_methods_01.py` keeps the span, wavenumber, refusal, Table C.4, Love (1911), and weak-start checks. The LoadDef and `earth_prem` comparisons were one-time checks.

## Setting the Method

A call's argument wins over a world's pinned `[radial_solver]` table, which wins over the configuration.

```python
import TidalPy
from TidalPy.Structures import build_world

io = build_world("io")  # A bundled world
io.solve_eos()  # The Love solve reads the interior structure
io.solve_love_numbers(
    frequency=4.11e-5,
    starting_method="power_series",
)  # Start this solve from the power series
print(io.love_number_k)

io.set_solver_defaults(
    radial_solver={"starting_method": "kamata"},
)  # Every later solve of this world starts from Kamata's forms

TidalPy.reinit(
    provided_config={"radial_solver": {"starting_method": "ps"}},
)  # Every solve that sets nothing else starts from the power series
```

In a world or configuration file:

```toml
[radial_solver]
starting_method = "power_series"
```

A name TidalPy does not know is refused in a world file. In `TidalPy_Configs.toml` it falls back to the default with a warning.

The starting conditions can also be evaluated directly. Each function writes one row per independent solution into a complex array of shape (number of solutions, 2 × number of solutions), as the driver `find_starting_conditions(layer_type, is_static, is_incompressible, starting_method, ...)` does.

```python
import numpy as np
from TidalPy.RadialSolver.starting.power_series import power_series_solid_dynamic_compressible

starting_conditions = np.zeros((3, 6), dtype=np.complex128)  # Three solutions of y1 to y6
power_series_solid_dynamic_compressible(
    4.11e-5,
    1.0e5,
    3300.0,
    1.0e11 + 0.0j,
    6.0e10 + 1.0e8j,
    2,
    None,
    starting_conditions,
)  # Frequency [rad s-1], radius [m], density [kg m-3], bulk and shear moduli [Pa], degree, G (None: SciPy's)
```

## References

- Crossley, D. J. (1975). The free-oscillation equations at the centre of the Earth. *Geophysical Journal of the Royal Astronomical Society*, 41(2), 153-163.
- Kamata, S., Matsuyama, I., and Nimmo, F. (2015). Tidal resonance in icy satellites with subsurface oceans. *Journal of Geophysical Research: Planets*, 120(9), 1528-1542.
- Love, A. E. H. (1911). *Some Problems of Geodynamics*. Cambridge University Press.
- Martens, H. R. (2016). *Using Earth deformation caused by surface mass loading to constrain the elastic structure of the crust and mantle*. PhD thesis, California Institute of Technology.
- Martens, H. R., Rivera, L., and Simons, M. (2019). LoadDef: A Python-based toolkit to model elastic deformation caused by surface mass loading on spherically symmetric bodies. *Earth and Space Science*, 6(2), 311-323.
- Saito, M. (1974). Some problems of static deformation of the earth. *Journal of Physics of the Earth*, 22(1), 123-140.
- Smylie, D. E. (2013). *Earth Dynamics: Deformations and Oscillations of the Rotating Earth*. Cambridge University Press.
- Takeuchi, H., and Saito, M. (1972). Seismic Surface Waves. In *Methods in Computational Physics: Advances in Research and Applications*, 11, 217-295.
