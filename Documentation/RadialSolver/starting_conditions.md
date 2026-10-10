# Starting Conditions (`RadialSolver.starting`)

_Updated: 2026-10-10_

The shooting method integrates its independent solutions (three in a solid, two in a dynamic liquid, one in a static liquid) outward from a starting radius near the center. The starting conditions are their values there. `starting_method` chooses them; `radial_solver`, `BaseWorld.solve_love_numbers`, and the `[radial_solver]` configuration and world tables all take it. The default is `"takeuchi"`.

| `starting_method` | Aliases | Form |
|---|---|---|
| `"takeuchi"` | `"ts"` | Takeuchi and Saito (1972) closed forms in spherical Bessel functions |
| `"kamata"` | | Kamata et al. (2015) closed forms |
| `"power_series"` | `"powerseries"`, `"ps"`, `"martens"` | Martens (2016) power series in $r^2$, summed to convergence |
| `"unity"` | | Unit vectors |

Every method covers every layer type. All but unity start from solutions regular at the center of a homogeneous layer (constant density and moduli, gravity $g = \gamma r$ with $\gamma = 4 \pi G \rho / 3$). Some forms were derived for TidalPy (see [Methods](#methods)). Only the layer the integration starts in uses these conditions; every layer above starts from the interface conditions with the layer below.

## Choosing a Method

- `"takeuchi"`, the default, is accurate wherever it solves, weak solid starting layers included (a very weak one can make it fail; see the note below), and is the cheapest along with the power series.
- `"kamata"` was as accurate as Takeuchi and Saito at degrees 2 to 10 (median errors within a factor of two of each other) at about 1.1 times the cost (see [Size of the Starting Vectors](#size-of-the-starting-vectors)). At degrees 20 to 40 its errors were lower, up to 1000 times (`luna` at degree 40: 2.5e-9 against 1.9e-6).
- `"power_series"` costs about as much as Takeuchi and Saito and gives the same Love numbers where it starts. It needs no special functions and stays finite without gravity, where the compressible Takeuchi and Saito forms and every Kamata form divide by $\gamma$. It refuses to start in a weak solid or deep in an exponential regime.
- `"unity"` assumes nothing about the layer, so it also starts where the closed forms overflow, inside a dynamic liquid at long periods and in a very weak solid. Its start leaves an error that depends on the starting radius and grows large at periods of minutes (see [Unity](#unity)). Use it there, and to test whether a result depends on the start.

> [!NOTE]
> In a very weak solid starting layer the closed-form solutions grow steeply across the layer, and the integrator can fail to follow them. With an $\eta = 10^5$ Pa s Maxwell start ($|\mu^*| \approx 10^{-11} \rho g R$, the weak solid model of [Accuracy and Cost](#accuracy-and-cost)), Takeuchi and Saito failed from degree 5 and Kamata from degree 10 with error code -11 ("Required step size is less than spacing between numbers"), and the power series refused. Unity solved every degree from 2 to 20, and the starts that did solve agreed with it to 1e-10. Tighter tolerances fail sooner. Use `"unity"` there.

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
)  # This solve only
print(io.love_number_k)

io.set_solver_defaults(
    radial_solver={"starting_method": "kamata"},
)  # Every later solve of this world

TidalPy.reinit(
    provided_config={"radial_solver": {"starting_method": "ps"}},
)  # Every solve that sets nothing else
```

In a world or configuration file:

```toml
[radial_solver]
starting_method = "power_series"
```

An unknown name is refused in a world file; in `TidalPy_Configs.toml` it falls back to the default with a warning.

To evaluate the starting conditions directly, each function (and the driver `find_starting_conditions(layer_type, is_static, is_incompressible, starting_method, ...)`) fills a complex array with one row per independent solution, of shape (number of solutions, 2 × number of solutions).

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

## Methods

Static liquids start from Saito's (1974) solution, $y_5 = r^l$ and $y_7 = 2 (l - 1) r^{l - 1}$, for every method but unity. It is exact for a homogeneous liquid, and the power series of the static liquid equations ends after its first term with the same values.

### Takeuchi and Saito

`TidalPy.RadialSolver.starting.takeuchi` evaluates Eqs. 95-102 of Takeuchi and Saito (1972) for compressible solids and dynamic liquids. Two of a solid's solutions are spherical Bessel functions of the layer's two wavenumbers; the third is a polynomial.

Takeuchi and Saito give no incompressible forms. TidalPy takes the limit of their solutions as the bulk modulus grows without bound ($1 / \alpha^2 \to 0$):

- One wavenumber tends to the shear wave's, $k^2 = \omega^2 / \beta^2$, and its solution keeps a finite limit.
- The other tends to zero, and its solution, rescaled, tends to an exact pressure-potential solution, $(y_1, \ldots, y_6) = (0, -\rho r^l, 0, 0, r^l, (2l + 1) r^{l - 1})$.
- The polynomial solution is unchanged, since its displacement has no divergence.

As $\omega \to 0$ the shear-wave solution converges onto the pressure-potential one, so it is replaced by its difference from it, written with no division by $\omega^2$. The static incompressible solid is that form at $\omega = 0$. The incompressible liquid's two solutions are the pressure-potential solution, $(y_1, y_2, y_5, y_6) = (0, -\rho r^l, r^l, (2l + 1) r^{l - 1})$, and the polynomial solution, both exact.

At large $\omega r / \beta$ the shear-wave and pressure-potential solutions again nearly coincide (at degree 8, a normalized smallest singular value of $4 \times 10^{-5}$ at $\omega r / \beta = 20$), as the published compressible forms do. Kamata's normalization stays independent there.

### Kamata

`TidalPy.RadialSolver.starting.kamata` evaluates Kamata et al. (2015) Eqs. B1-B37: Takeuchi and Saito's solutions plus incompressible forms. For a dynamic incompressible solid, the first solution is replaced by its difference from the second times $\gamma / \omega^2$: the published pair converges as $\omega \to 0$, while the combination stays independent at tidal periods. It is finite at $\omega = 0$, where the static equations are the dynamic ones, so the static incompressible solid (not given by Kamata et al.) is that form at $\omega = 0$.

Kamata's forms carry a common factor of $r^{-l}$, which makes their starting vectors orders of magnitude larger than the other methods' and sets their cost and accuracy (see [Size of the Starting Vectors](#size-of-the-starting-vectors)).

### Power Series

`TidalPy.RadialSolver.starting.power_series` expands the solutions about the center (Crossley 1975; Smylie 2013; Martens 2016, Sec. 4.2.8; Martens et al. 2019). Martens defines $y_1$ to $y_5$ as TidalPy does (Takeuchi and Saito 1972), but her $y_6$ leaves out a term, and her series is converted with

$$y_6^{\mathrm{TidalPy}} = y_6^{\mathrm{Martens}} + \frac{l + 1}{r} y_5$$

With $y_i = r^{l - 2 + s_i + \nu_0} z_i$, every layer's equations become $r \, dz/dr = (B_0 + B_2 r^2) z$. A regular solution is then $z = \sum_k c_k r^{2k}$, with $(\nu_0 I - B_0) c_0 = 0$ and $((\nu_0 + 2k) I - B_0) c_k = B_2 c_{k - 1}$. The scalings are $s = (1, 0, 1, 0, 2, 1)$ for a solid's $y_1$ to $y_6$ and $s = (1, 2, 2, 1)$ for a dynamic liquid's $y_1$, $y_2$, $y_5$, $y_6$. TidalPy differs from Martens in these ways:

- Terms are added until the next falls below machine precision (Martens stops at $r^4$), which keeps the series accurate at the larger starting radii of high degrees.
- Every kind of solid and dynamic liquid has a series, derived from TidalPy's equations; Martens covers the compressible solid. Static liquids use Saito's solution.
- The first solution is Takeuchi and Saito's polynomial solution, a combination of Martens' $A_{1,1}$ and $A_{6,1}$ solutions that is exact as its first term. In an incompressible solid the $A_{6,1}$ solution is the exact pressure-potential solution. Both are written without a series.
- In a compressible solid with gravity, the other two solutions are the series of Takeuchi and Saito's two wavenumber solutions. Each starts as $(\alpha^2 f - (l + 1) \beta^2)$ times Martens' $A_{6,1}$ vector, with the free part of its resonant $r^2$ step set to Takeuchi and Saito's $z_4$, so each carries one wavenumber. Elsewhere (incompressible solids, which have one wavenumber, and solids without gravity) the remaining solutions are Martens' $A_{6,1}$ and $A_{4,0}$, with that free part set to zero, as Martens sets $A_{4,2} = 0$.
- The recurrence is scaled by the shear modulus and by $\gamma + \omega^2$, so it holds in SI as well as non-dimensional units.

The series refuses to start, naming the closed-form starts, when:

- it needs more than 100 terms;
- cancellation would take more than half the digits of any component;
- a solution grows from its first terms by more than $1/\sqrt{\epsilon}$: the start is deep in an exponential regime, such as a dynamic liquid at long periods (see [Dynamic Liquid Layers at Long Forcing Periods](dense_radial_solution.md#dynamic-liquid-layers-at-long-forcing-periods));
- in a solid, the larger wavenumber times the starting radius exceeds $\tfrac{1}{2} \ln(1/\sqrt{\epsilon}) \approx 9$, since the slower wavenumber solution's roundoff error grows as $\epsilon e^{2 |k| r}$. This is the weak solid starting layer.

The closed forms carry that growth analytically and apply in all four cases. Where the series starts, its solutions match Takeuchi and Saito's one by one to better than $10^{-6}$ in solids with $|\mu|$ down to $10^{-6} \rho g R$ (degrees 2 to 10); below that an accepted static start can differ by up to $10^{-5}$. At degree 1 a static solid has a rigid-translation solution, so single solutions are not defined there.

### Unity

The unit vectors are $y_1$, $y_4$, and $y_6$ in a solid, $y_1$ and $y_6$ in a dynamic liquid, and $y_7$ in a static liquid: the leading components of Martens' free-constant solutions. This set's errors were within a factor of two of the best set of unit vectors on homogeneous bodies; sets that include $y_3$ can leave the regular solutions degenerate.

A unit vector holds singular content as well as regular. Integrating outward, the singular part decays relative to the regular part as $(r_0 / r)^{2l - 1}$ in a solid and $(r_0 / r)^{2l + 1}$ in a liquid, where $r_0$ is the starting radius. Its surface error therefore depends on the starting radius, not the integration tolerance. Where the solutions oscillate at the start (a large $\omega r_0 / \beta$, at periods of minutes) the singular part need not decay: in a homogeneous solid at degrees 8, 10, and 20 with $\omega r_0 / \beta$ = 10, 20, and 40, unity was off by 1e-3 or more in 8 of the 9 cases (by 1 at degree 20 and 20) while the closed forms agreed to 1e-5. Unity is a check on the other methods and a start where they overflow, not a default.

## Accuracy and Cost

One-time measurements (2026-10-09, TidalPy 0.8.0), with errors in $k$ against references at `rtol = 1e-12` and costs relative to Takeuchi and Saito:

- **Bundled worlds** (every world as bundled, with dynamic solids, all dynamic, all dynamic with incompressible liquids, and with incompressible solids; degrees 2 to 10, periods 0.5 to 30 days, packaged settings; 1427 solves per method, median 0.24 ms): median errors of 7.4e-10 for Kamata, 1.1e-9 for Takeuchi and Saito and for the power series, and 2.3e-9 for unity, at median costs of 1.12, 1.00, 1.02, and 1.16. No solve was off by more than $10^{-5}$, and the largest errors were 8.0e-7, 9.6e-7, 1.9e-6, and 8.5e-6. Every failure was a start inside a dynamic liquid (below): 31 for Takeuchi and Saito and for Kamata, 78 for the power series, and none for unity.
- **Weak solid start** (a Maxwell solid, $\mu$ = 60 GPa, under a 171 km elastic lid, 1.8 day period, degrees 2 to 20): at $\eta = 10^{11}$, $10^{10}$, and $10^{9}$ Pa s every closed-form and unity start solved, with median errors between 3e-12 and 2e-10 (unity's the largest) and a `surface_solve_rcond` of at least 5e-4. The power series refused 2 or 3 of the 5 starts. Kamata cost 0.9 to 1.2 times Takeuchi and Saito.
- On homogeneous solids and a liquid core under a solid mantle, every median error was between 5e-10 and 2e-8, except unity's 3e-8 on homogeneous solids.
- Every closed-form and series start treats a dynamic liquid as a homogeneous liquid sphere. Where the liquid is compressible, that constant density is unstably stratified, so at long periods the regular solutions grow exponentially from the center and a start inside the liquid fails. An incompressible dynamic liquid of constant density is neutral, and starts inside it solve (`luna_dynamic`'s outer core made incompressible solves at 200, 1000, and $10^4$ days). In `luna_dynamic`'s outer core at degree 5, the Takeuchi and Saito and Kamata starts solve at 150 days and fail from 200 days, and the power series refuses from about 10 days. Unity's unit vectors carry no such growth, and there it solved out to 1000 days, within 4e-9 of a start in the solid inner core. In the bundled worlds it solved all 31 starts inside a dynamic liquid core that the closed forms failed, consistent to 2e-8 across starting radii where the deeper starts solved. A start in the layer below the liquid, with a manual `starting_radius`, also avoids the failure.
- The start itself takes 0.5 to 1.6 µs for the closed forms and up to 14 µs for the power series, against a 0.1 to 1 ms solve. The cost differences come from the integration.

### Size of the Starting Vectors

Kamata's starting vectors are larger than the others' by the factor $r^{-l}$. The integrator measures its absolute tolerance against those sizes, so Kamata's start is integrated to a tighter effective tolerance: that is its small edge in accuracy and its extra cost. A lower `atol` does the same for the other starts. On the bundled worlds as bundled, Takeuchi and Saito at `atol = 1e-10` had a median error of 3.9e-10 against Kamata's 1.0e-9 at the defaults, with 1.3 times the steps.

### Validation

- Martens' series as written matches LoadDef's `fundamental_solutions_powser` to $10^{-14}$ and reproduces her homogeneous-Earth Love numbers (Martens 2016, Table C.4: $V_P$ = 10 km/s, $V_S$ = 5 km/s, $\rho$ = 5000 kg m$^{-3}$, M2 period, $G = 6.672 \times 10^{-11}$) to the table's rounding for $n$ = 2 to 100. TidalPy's series reproduces Table C.4 too, and spans the Takeuchi and Saito and the Kamata solutions to better than $10^{-9}$.
- The derived incompressible forms satisfy TidalPy's equations, span the series' solutions to about $10^{-11}$, and agree with Kamata's dynamic incompressible forms. Every regular start reproduces Love's (1911) closed form for a static incompressible solid, $k_2 = (3/2) / (1 + 19 \mu / (2 \rho g R))$.
- On `earth_prem`, all four methods give the same tidal and load Love numbers to four decimals for degrees 2 to 10 (a one-time check, like the LoadDef one; the test suite keeps the rest).

## References

- Crossley, D. J. (1975). The free-oscillation equations at the centre of the Earth. *Geophysical Journal of the Royal Astronomical Society*, 41(2), 153-163.
- Kamata, S., Matsuyama, I., and Nimmo, F. (2015). Tidal resonance in icy satellites with subsurface oceans. *Journal of Geophysical Research: Planets*, 120(9), 1528-1542.
- Love, A. E. H. (1911). *Some Problems of Geodynamics*. Cambridge University Press.
- Martens, H. R. (2016). *Using Earth deformation caused by surface mass loading to constrain the elastic structure of the crust and mantle*. PhD thesis, California Institute of Technology.
- Martens, H. R., Rivera, L., and Simons, M. (2019). LoadDef: A Python-based toolkit to model elastic deformation caused by surface mass loading on spherically symmetric bodies. *Earth and Space Science*, 6(2), 311-323.
- Saito, M. (1974). Some problems of static deformation of the earth. *Journal of Physics of the Earth*, 22(1), 123-140.
- Smylie, D. E. (2013). *Earth Dynamics: Deformations and Oscillations of the Rotating Earth*. Cambridge University Press.
- Takeuchi, H., and Saito, M. (1972). Seismic Surface Waves. In *Methods in Computational Physics: Advances in Research and Applications*, 11, 217-295.
