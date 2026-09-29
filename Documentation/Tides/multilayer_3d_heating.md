# 3D Tidal Stress, Strain, and Heating (`Tides.multilayer`)

_Updated: 2026-09-29_

`TidalPy.Tides.multilayer` calculates the depth- and direction-resolved tidal response of a layered world: the complex strain and stress tensors and the volumetric heating. The response is evaluated point by point.

## Physics

The response at a point $(r, \theta, \phi, t)$ factorizes into a radial part and an angular and time part. The radial part is:

- The viscoelastic-gravitational functions $y_1 \dots y_6$ from the radial solver (Takeuchi and Saito 1972), evaluated at $r$ through the dense calling system (`RadialSolverSolution.get_radial_solution`).
- The complex shear and bulk moduli $\mu^{*}$ and $K^{*}$ at $r$, read at the solved frequency through the dense material path (`c_EOSSolution::call_material`, exposed in Python as `RadialSolverSolution.get_complex_shear_modulus`).

The angular and time part is the tidal potential $W(\theta, \phi, t)$ and its first and second $\theta$ and $\phi$ derivatives. The radial problem depends only on the degree $l$ and the frequency $|\omega|$, not on $m$, so it is solved once per $(l, |\omega|)$.

### Tidal Potential

The potential is the Kaula expansion of the [global tides](global_tides.md) page (Kaula 1964; Efroimsky and Williams 2009, Eq. 18), evaluated at the surface radius $R$. It is linear in $F_{lmp}$, $G_{lpq}$, and $P_{lm}$. The global path squares $F$ and $G$ because global heating scales with the potential squared. The engine (`c_tidal_potential_3d_modes`, wrapped as `TidalPy.Tides.potential.tidal_potential_3d_modes`) returns, for each active mode, the degree $l$, the signed forcing frequency, and the complex amplitude of $W$ and its derivatives, with $W(t) = \mathrm{Re}[W_c\,e^{i\omega t}]$:

$$W(\theta, \phi, t) = \sum_{l,m,p,q} A_{lmpq}\,P_{lm}(\cos\theta)\,\mathcal{T}_{lm}\!\left(\omega_{lmpq}t - m\phi\right), \qquad A_{lmpq} = \frac{G M_h}{a}\left(\frac{R}{a}\right)^{l}\frac{(l-m)!}{(l+m)!}\left(2-\delta_{0m}\right)F_{lmp}(I)\,G_{lpq}(e),$$

$$\omega_{lmpq} = (l - 2p + q)\,n - m\,\dot{\theta},$$

where $n$ is the mean motion, $\dot{\theta}$ the spin rate, and $\mathcal{T}_{lm}$ is $\cos$ for even $l - m$ and $\sin$ for odd $l - m$. Longitude $\phi$ is the body-fixed east longitude, increasing in the direction of rotation. $\phi = 0$ is the host's meridian at $t = 0$, when the host is at its ascending node and at periapse. As in Kaula, $P_{lm}$ carries no Condon-Shortley phase. TidalPy's Legendre tables (`TidalPy.Utilities.legendre`) include it, so each amplitude carries a factor $(-1)^{m}$ that cancels it. All radial dependence is in $y(r)$.

> [!WARNING]
> This assumes no periapse or node precession. It also assumes that the change in the mean anomaly can be approximated by the mean motion.

Three settings in the world's `[tides]` config set the truncation: `max_degree_l` (2..10), `eccentricity_trunc_lvl`, and `obliquity_trunc_lvl` (0 = off). Eccentricity level $N$ keeps the potential through $e^N$. As in the 1D path, it cuts every product of two eccentricity functions in the secular heating at $e^N$, so at zero obliquity the volume integral of the secular heating equals the 1D heating at any eccentricity and spin rate (see [Eccentricity Functions](eccentricity.md)). With both $e$ and $I$ nonzero the two differ by the cross terms described in [Coherent Waves](#coherent-waves). Obliquity level $N$ does the same in $I$: the potential through $I^N$, and every product of two obliquity functions cut at $I^N$ (see [Obliquity Functions](obliquity.md)). The instantaneous fields (displacements, stress, strain, and power) are linear in the potential and use the unsquared functions. A nonzero obliquity truncation turns on the odd-`m` harmonics (`P_21`, ...). A mode with $|\omega|$ at or below `[numerical] minimum_frequency` (1e-14 rad/s, `TidalPy.constants.min_frequency`) is later switched off, the same floor the 1D `calc_tides` path uses.

### Displacement, Strain, and Stress

For one mode, the displacement is (Tobie et al. 2005, Eq. 9)

$$\mathbf{u} = \left(y_{1}\,W,\;\; y_{3}\,\frac{\partial W}{\partial\theta},\;\; \frac{y_{3}}{\sin\theta}\,\frac{\partial W}{\partial\phi}\right),$$

and the strain is (Tobie et al. 2005, Eq. 10, with the correction of Kervazo et al. 2021, Appendix D, to $\varepsilon_{\phi\phi}$ and $\varepsilon_{\theta\phi}$)

$$\begin{aligned}
\varepsilon_{rr} &= \frac{dy_{1}}{dr}\,W, \\
\varepsilon_{\theta\theta} &= \frac{y_{1}}{r}\,W + \frac{y_{3}}{r}\,\frac{\partial^{2} W}{\partial\theta^{2}}, \\
\varepsilon_{\phi\phi} &= \frac{y_{1}}{r}\,W + \frac{y_{3}}{r}\left(\frac{1}{\sin^{2}\theta}\,\frac{\partial^{2} W}{\partial\phi^{2}} + \cot\theta\,\frac{\partial W}{\partial\theta}\right), \\
\varepsilon_{r\theta} &= \frac{y_{4}}{2\mu^{*}}\,\frac{\partial W}{\partial\theta}, \\
\varepsilon_{r\phi} &= \frac{y_{4}}{2\mu^{*}}\,\frac{1}{\sin\theta}\,\frac{\partial W}{\partial\phi}, \\
\varepsilon_{\theta\phi} &= \frac{y_{3}}{r}\,\frac{1}{\sin\theta}\left(\frac{\partial^{2} W}{\partial\theta\,\partial\phi} - \cot\theta\,\frac{\partial W}{\partial\phi}\right),
\end{aligned}$$

where $dy_{1}/dr$ comes from the radial equations of the layer's type. The stress follows the isotropic constitutive law (Takeuchi and Saito 1972)

$$\sigma_{ij} = 2\mu^{*}\varepsilon_{ij} + \lambda^{*}\varepsilon_{kk}\,\delta_{ij}, \qquad \lambda^{*} = K^{*} - \tfrac{2}{3}\mu^{*},$$

with both moduli taken at $r$ and at the mode's frequency. The isotropic part is taken from the radial stress $\sigma_{rr} = y_{2}\,W$:

$$\lambda^{*}\varepsilon_{kk} = \left(y_{2} - 2\mu^{*}\frac{dy_{1}}{dr}\right)W.$$

In a compressible layer this is the same quantity, since $y_{2} = \lambda^{*}\left(dy_{1}/dr + (2y_{1} - l(l+1)y_{3})/r\right) + 2\mu^{*}\,dy_{1}/dr$, and it stays accurate for a very large bulk modulus, where $\lambda^{*}$ times a nearly zero trace would amplify round-off. In a layer the radial solver treats as incompressible, the strain is traceless and this term is the pressure, which only $y_{2}$ carries.

The kernel applies to solid layers only. A liquid (a liquid layer, or a molten stretch of a solid layer that the radial solver treats as a static liquid) contributes no shear dissipation: its heating is 0 and its stress and strain are NaN. The poles ($\sin\theta$ within machine epsilon of 0, at $0$ and $\pi$) are singular points of the angular terms, so point-wise values there are NaN. Use colatitudes inside $(0, \pi)$, as the Gauss-Legendre nodes of the summed paths do.

### Coherent Waves

Before the heating is formed, each raw $(l, m, p, q)$ mode is mapped onto its non-negative frequency: a mode with $\omega < 0$ contributes the complex conjugate of its phasor at $|\omega|$, since $\mathrm{Re}[W_c e^{i\omega t}] = \mathrm{Re}[\overline{W_c}\,e^{-i\omega t}]$. It is then merged with every other mode that shares its real spatial function: the same degree $l$, order $m$, $|\omega|$, and azimuthal sign ($e^{-im\phi}$ or $e^{+im\phi}$, irrelevant for $m = 0$). The kernel works with the resulting coherent waves.

The $m = 0$ modes always come in pairs: $(l, 0, p, q)$ at $+\omega$ and $(l, 0, l-p, -q)$ at $-\omega$. The two carry equal amplitudes ($F_{l0p} = \pm F_{l0,l-p}$ with the parity sign, $G_{lpq} = G_{l,l-p,-q}$) and are the same function of time, since $\cos(-x) = \cos x$. Together they are one real sinusoid of twice the amplitude. Heating scales with the amplitude squared, so summing their cycle-averaged powers separately loses half of the zonal heating. For a homogeneous degree-2 body at zero obliquity the zonal terms are 9/84 of the total, so the loss is 4.5/84 = 5.36% of the heating of a synchronously rotating body, where only the eccentricity modes survive. The 1D formula counts the same pair through its $(2 - \delta_{0m})$ weighting, so it needs no merge.

At nonzero obliquity, modes of the same $(l, m)$ with different $(p, q)$ can also share a signed frequency. Their relative phase is set by the argument of periapsis $\omega$, which the engine takes as zero, so they also combine coherently. The 1D formula averages over $\omega$ and drops these cross terms. The 3D path keeps them. Its total is the heating of an orbit whose periapsis stays at the ascending node, while the 1D value is the secular heating of an orbit whose periapsis precesses. The difference needs both $e$ and $I$ nonzero, varies with $\cos 2\omega$ and $\cos 4\omega$, and changes sign at $\omega = \pi/2$. A direct calculation from the host's exact Kepler position confirms both values: the 1D heating equals its average over $\omega$ to 1e-14, and the 3D total equals its $\omega = 0$ value to the radial quadrature error. The difference is largest in synchronous rotation, where the semidiurnal tide is static, and grows with $e$ and $I$ at every spin rate. The table gives the 3D total over the 1D heating, minus one, at degree 2 (measured 2026-09-29).

| Body | $e$ | $I$ \[rad\] | Spin $n$ | Spin $1.5 n$ | Spin $2.5 n$ |
|---|---|---|---|---|---|
| Homogeneous Maxwell | 0.05 | 0.1 | -1.7e-3 | -6.9e-7 | 6.6e-8 |
| Homogeneous Maxwell | 0.05 | 0.3 | -3.8e-3 | -5.4e-6 | 5.4e-7 |
| Homogeneous Maxwell | 0.2 | 0.1 | -2.3e-3 | -1.7e-4 | -1.1e-6 |
| Homogeneous Maxwell | 0.2 | 0.3 | -1.6e-2 | -1.3e-3 | -8.0e-6 |
| Bundled Io | 0.05 | 0.1 | -1.6e-3 | 6.3e-6 | 1.4e-5 |
| Bundled Io | 0.05 | 0.3 | -3.6e-3 | 4.9e-5 | 1.1e-4 |
| Bundled Io | 0.2 | 0.1 | -2.2e-3 | -6.0e-5 | 1.9e-4 |
| Bundled Io | 0.2 | 0.3 | -1.5e-2 | -4.7e-4 | 1.6e-3 |

Use the 1D heating when the periapsis precesses over the time of interest.

### Secular Heating

The secular (cycle- or orbit-averaged) volumetric heating \[W m$^{-3}$\] is built from the complex amplitudes, so the cycle average is exact with no time grid (Tobie et al. 2005):

$$\bar{h}(r, \theta, \phi) = \sum_{\chi}\frac{\chi}{2}\sum_{k} w_{k}\,\mathrm{Im}\!\left(\sigma_{k}^{(\chi)}\,\overline{\varepsilon_{k}^{(\chi)}}\right), \qquad w_{k} = \begin{cases} 1, & k \in \{rr, \theta\theta, \phi\phi\}, \\ 2, & k \in \{r\theta, r\phi, \theta\phi\}, \end{cases}$$

where $\sigma_{k}^{(\chi)}$ and $\varepsilon_{k}^{(\chi)}$ are the total complex stress and strain amplitudes at the frequency $\chi = |\omega|$, summed over every wave at that frequency before the bilinear form. The form is non-negative for a dissipative material. Cross terms between waves at different frequencies average to zero over the orbit and are dropped, while those at the same frequency survive and are kept. When the waves at one frequency have different longitude structure, these cross terms make $\bar{h}$ depend on longitude. The zonal and sectoral waves of a synchronously rotating body, where every active mode sits at a multiple of $n$, produce the $\cos 2\phi$ and $\cos 4\phi$ patterns of a synchronous heating map, symmetric about the sub-host meridian. For a generic non-synchronous spin, every frequency carries one wave and $\bar{h}$ is a function of $r$ and $\theta$ only.

### Volume Integral and Collapse

The total power is the volume integral

$$\dot{E} = \int_{0}^{R}\int_{0}^{\pi}\int_{0}^{2\pi}\bar{h}\,r^{2}\sin\theta\;d\phi\,d\theta\,dr,$$

which `calc_3d_tides` evaluates with its summed arguments, one axis at a time:

- Longitude: waves with different signed azimuthal wavenumbers $\mu = \pm m$ integrate to zero over longitude, so the integral is $2\pi$ times the sum over $(\chi, \mu)$ groups evaluated at $\phi = 0$.
- Colatitude: every strain and stress component of a wave is a combination, with complex radial coefficients, of six real functions of $\theta$ for its $(l, m)$,

  $$f_{1} = P_{lm}, \quad f_{2} = \frac{dP_{lm}}{d\theta}, \quad f_{3} = \frac{d^{2}P_{lm}}{d\theta^{2}}, \quad f_{4} = \frac{P_{lm}}{\sin\theta}, \quad f_{5} = -\frac{m^{2}P_{lm}}{\sin^{2}\theta} + \cot\theta\,\frac{dP_{lm}}{d\theta}, \quad f_{6} = \frac{1}{\sin\theta}\left(\frac{dP_{lm}}{d\theta} - \cot\theta\,P_{lm}\right),$$

  so the colatitude integral of any pair of waves reduces to the Gram matrices

  $$\mathcal{G}_{ij}(l_{a}, l_{b}, m) = \int_{0}^{\pi} f_{i}^{(l_{a})}\,f_{j}^{(l_{b})}\,\sin\theta\;d\theta.$$

  They are tabulated for equal degrees and computed with 32-node Gauss-Legendre quadrature in $\cos\theta$ for different degrees, which is exact because the integrands are polynomials in $\cos\theta$ of degree at most $l_{a} + l_{b} + 2 \le 22$.
- Radius: Gauss-Legendre quadrature inside each layer with `radial_slices` nodes (default 16) and weight $r^{2}$; no node sits on a layer boundary. A node below the radial solver's starting radius, where no degree has a solution, is left out. The automatic starting radius keeps that region small, and the integral logs a warning when the nodes left out hold more of the body's volume than the radial solver's `rtol`.

At zero obliquity the volume integral equals the 1D global tidal heating (`get_tidal_heating`) to the radial quadrature error at every eccentricity and spin rate, including synchronous rotation. For the homogeneous Io of demo notebook 09 that error is below 1e-7 with the default 16 nodes per layer. With both $e$ and $I$ nonzero the two differ by the same-frequency cross terms of [Coherent Waves](#coherent-waves). The benchmark tests are `Tests/Test_Structures/Test_Worlds/test_world_1d_vs_3d_tides_01.py` and `test_world_3d_tides_coherent_01.py`.

## Python API

### World Method

A built world delegates these methods to its rheology tide model (`c_RheologyTide`), which runs in C++ and calls the world's radial-solver and EOS members directly.

```python
import numpy as np

from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.Rheology import Elastic, Maxwell
from TidalPy.Structures.layers.physics import PhysicsLayer
from TidalPy.Structures.worlds.layered import LayeredWorld
from TidalPy.Tides.classes import make_tide
from TidalPy.Viscosity import make_viscosity

radius  = 1.8e6    # [m]
density = 3500.0   # [kg m-3]
mass    = (4.0 / 3.0) * np.pi * radius**3 * density

layer = PhysicsLayer("mantle", 0, 0.0, radius, mass)
layer.set_eos(ConstantDensityEOS(reference_density=density, shear_modulus_static=6.0e10, bulk_modulus_static=1.0e11))
layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e19}))
layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e30}))
layer.set_shear_rheology(Maxwell())
layer.set_bulk_rheology(Elastic())

world = LayeredWorld("io_like", radius, mass)
world.add_layer(layer)
world.solve_eos()
world.set_tide_model(make_tide("rheology"))
# The tidal potential truncation comes from the tide config (or the [tides] TOML table):
world.set_tide_config(max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)

orbital_frequency = 4.1e-5      # [rad s-1]
spin_frequency = 5.6e-5         # non-synchronous
eccentricity = 0.004
obliquity = 0.0
semi_major_axis = 4.2e8         # [m]
host_mass = 1.9e27              # [kg]

h_bar = world.get_3d_tidal_heating(
    orbital_frequency, spin_frequency, eccentricity, obliquity,
    semi_major_axis, host_mass, 0.9 * world.radius, 0.8)   # [W m-3], secular
```

This requires the rheology tide model and a solved EOS. Analytic tide models such as fixed-Q have no depth-resolved solution and are rejected. Every 3D method checks the orbital state as `calc_tides` does: an eccentricity outside $[0, 1)$ or a semi-major axis that is not positive raises `ValueError`, and the truncation-range warnings apply. The radial solver's starting radius grows with degree, so with several degrees active the innermost region carries only the degrees that have a solution there. A radius is NaN only where no degree has a solution. This holds for every output with a radius axis, including a radial profile summed over colatitude and longitude.

#### Building the Map

`get_3d_tidal_heating` re-solves the world radial (Love-number) response at every point. `get_3d_tidal_heating_array` takes paired, equal-length `(radius, colatitude)` arrays. It builds the position-independent mode list once, solves the radial response once per unique `(l, frequency)`, and reuses it across all points. It returns an array of the same shape, NaN where a radius has no depth-resolved solution:

```python
radii = np.linspace(0.01 * world.radius, 0.999 * world.radius, 40)
colatitudes = np.linspace(1e-3, np.pi - 1e-3, 60)
radius_grid, colatitude_grid = np.meshgrid(radii, colatitudes, indexing="ij")
hbar = world.get_3d_tidal_heating_array(
    orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass,
    radius_grid.ravel(), colatitude_grid.ravel()).reshape(radius_grid.shape)   # [W m-3], secular
```

It matches the scalar `get_3d_tidal_heating` point for point to machine precision, since the scalar path is the batch path with one point. Both return the longitude mean of the secular density, so `get_3d_tidal_heating_array` is the efficient way to build a zonal-mean map. For a longitude-resolved map, use `calc_3d_tides` with a longitude grid.

#### Full Grid and Collapsed Forms: `calc_3d_tides`

`calc_3d_tides` returns the 3D heating on a full grid over `(radius, colatitude, longitude[, time])` or reduced along any spatial dimension. It computes one of two quantities:

* `orbit_averaged=True` (default): the secular heating density `h_bar` [W m-3], the time average of the instantaneous power at each point. Its longitude dependence is described in Secular Heating.
* `orbit_averaged=False`: the instantaneous mechanical power density `sigma_ij(t) * eps_dot_ij(t)` [W m-3] at each supplied time, on a 4th axis. It depends on longitude and time. Every wave is a real field at `+|omega|`, so every cross term is present and the time average equals `h_bar` at the same point through the truncation level's $e^N$. Past $e^N$, the instantaneous power, formed from the unsquared potential, keeps partial terms that `h_bar` cuts, so the two differ at large eccentricity and a low level.

Non-summed spatial axes take user arrays (`radii`, `colatitudes`, `longitudes`). The `times` array is required when `orbit_averaged=False`. If any spatial axis is summed, the surviving spatial axes carry their Jacobian (`r^2`, `sin theta`, `1`), so a plain integral over them recovers the total. If none is summed, the output is the raw density, whose longitude mean matches the scalar `get_3d_tidal_heating`.

The secular colatitude integral (`latitude_summed`) is analytic by default. It uses the Gram matrices of Volume Integral and Collapse, with the per-`(l, m)` table in `Tides.multilayer.stress_strain.angular_gram`, so no colatitude grid is needed. For a large radius grid (a radial profile or map), this is a few times faster than the numerical quadrature. Otherwise the Love-number solves dominate.

The analytic integral is of the longitude mean, so it is used only when longitude is summed too. A call that keeps its longitudes uses the Gauss-Legendre colatitude grid (`latitude_nodes`) whatever `latitude_analytic` says. Pass `latitude_analytic=False` to use that grid for a longitude-summed call as well. The two agree to machine precision. The longitude integral is the analytic `2*pi` times the longitude mean when averaged, or a `longitude_nodes` trapezoid when instantaneous.

With `latitude_summed`, `colatitude_min` and `colatitude_max` [rad] (defaults 0 and pi) restrict the colatitude integral to `colatitude_min <= theta <= colatitude_max`. Complementary bands add up to the full-sphere result, so zonal heating budgets (polar caps against an equatorial belt, for example) take a few banded calls. A band narrower than the full sphere always uses Gauss-Legendre quadrature, because the analytic Gram table covers the full sphere only. The band has no effect when colatitude is not summed.

The returned dict carries the surviving axes plus one of:

- `heating`: the grid over the surviving axes, in the order radius, colatitude, longitude, time.
- `total` [W] and `per_layer` [W], when all three spatial axes are summed.

The C++ code writes the heating directly into the returned arrays:

```python
longitudes = np.linspace(0.0, 2.0 * np.pi, 60)

# default: full secular density grid (nr, ncolat, nlon); longitude-dependent at synchronous rotation
grid = world.calc_3d_tides(
    orbital_frequency,
    spin_frequency,
    eccentricity,
    obliquity,
    semi_major_axis,
    host_mass,
    radii=radii,
    colatitudes=colatitudes,
    longitudes=longitudes)['heating']

# per-layer and total heating (all three spatial axes summed)
res = world.calc_3d_tides(
    orbital_frequency,
    spin_frequency,
    eccentricity,
    obliquity,
    semi_major_axis,
    host_mass,
    latitude_summed=True,
    longitude_summed=True,
    radial_summed=True)

res['total']      # [W], equals get_tidal_heating at zero obliquity
res['per_layer']  # [W] per layer (innermost first), sums to res['total']

# radial power profile dP/dr (integrates to the total)
prof = world.calc_3d_tides(
    orbital_frequency,
    spin_frequency,
    eccentricity,
    obliquity,
    semi_major_axis,
    host_mass,
    radii=np.linspace(0.01 * world.radius, world.radius, 400),   # A radius below the solver's start would be NaN
    latitude_summed=True,
    longitude_summed=True)
np.trapezoid(prof['heating'], prof['radii'])   # The total, less the small share inside 0.01 R

# time-resolved instantaneous power over one orbital period
period = 2.0 * np.pi / orbital_frequency
res = world.calc_3d_tides(
    orbital_frequency,
    spin_frequency,
    eccentricity,
    obliquity,
    semi_major_axis,
    host_mass,
    radii=radii,
    colatitudes=colatitudes,
    longitudes=longitudes,
    times=np.linspace(0.0, period, 400),
    orbit_averaged=False)
res['heating']   # shape (nr, ncolat, nlon, ntime)
```

The per-layer totals replace the `tidal_scale` distribution for the depth-resolved rheology path. The radial and colatitude profiles integrate to the total. The fallback `latitude_nodes` (16) adds margin over the ~4 colatitude nodes a homogeneous degree-2 body needs. The `radial_slices` nodes sit strictly inside each layer because a node on a layer boundary would read the modulus and radial solution of the layer below, so a layer sitting on a stiffer one would lose part of its own heating. With the nodes inside, the collapsed total at zero obliquity matches the 1D heating to better than 1e-4 from 4 nodes per layer, for a homogeneous body and for the bundled Io, whose thin asthenosphere has a sixtieth of the shear modulus of the mantle beneath it. The default adds margin for higher degree l.

#### Displacements: `calc_3d_displacements`

For each coherent wave, the radial functions $y_1$ (radial displacement) and $y_3$ (tangential displacement) set the complex displacement amplitude at a point, $\mathbf{u} = \left(y_1 U,\; y_3\, \partial U / \partial\theta,\; y_3\, (\partial U / \partial\phi) / \sin\theta\right)$ (Tobie et al. 2005, Eq. 9). The amplitude evolves in time as $\mathrm{Re}\left[\mathbf{u}\, e^{i|\omega| t}\right]$ and is summed over the waves. `calc_3d_displacements` returns the three instantaneous displacement components on the full `(radius, colatitude, longitude, time)` grid [m]:

```python
out = world.calc_3d_displacements(
    orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass,
    radii=[0.9 * world.radius, world.radius], colatitudes=np.linspace(0.1, np.pi - 0.1, 30),
    longitudes=np.linspace(0.0, 2 * np.pi, 60), times=np.linspace(0.0, period, 12))
out["radial"].shape        # (2, 30, 60, 12)   u_r [m]; also out["polar"], out["azimuthal"]
```

At the surface $y_1 = h / g$ and $y_3 = l / g$, so the surface radial displacement is $h U / g$ for a single mode. `calc_3d_displacements` requires the rheology tide model, a solved EOS, and a radial-solver Love-number method (the analytic `homogeneous`/`cpl`/`ctl` methods have no radial functions). A radius without a depth-resolved solution (the center, below the solver start) is NaN. The radial functions themselves are available at any radius through `world.get_love_radial_y(radius, ytype_idx, y_idx)` after a radial-solver Love solve.

#### Stress and Strain: `calc_3d_stress_strain`

At each point, every wave's complex stress and strain amplitudes are added into the total of its frequency. Each component at time $t$ is the sum over frequencies of $\mathrm{Re}\left[A\, e^{i|\omega| t}\right]$, with $A$ the complex amplitude. `calc_3d_stress_strain` returns both instantaneous tensors on the full `(radius, colatitude, longitude, time)` grid with the six components on the last axis, ordered `rr`, `theta_theta`, `phi_phi`, `r_theta`, `r_phi`, `theta_phi` (also returned as `components`). The strain is the symmetric gradient of the displacement field of `calc_3d_displacements`.

```python
tensors = world.calc_3d_stress_strain(
    orbital_frequency,
    spin_frequency,
    eccentricity,
    obliquity,
    semi_major_axis,
    host_mass,
    radii=[0.9 * world.radius, world.radius],
    colatitudes=np.linspace(0.1, np.pi - 0.1, 30),
    longitudes=np.linspace(0.0, 2.0 * np.pi, 60),
    times=np.linspace(0.0, period, 12))
tensors['stress'].shape                  # (2, 30, 60, 12, 6) [Pa]; tensors['strain'] has the same shape
radial_stress = tensors['stress'][..., 0]   # The rr component
peak_stress = np.abs(tensors['stress']).max(axis=3)   # Largest magnitude of each component over the times
```

Each tensor takes 48 bytes per grid point and time, and either can be skipped with `return_stress=False` or `return_strain=False`. The C++ code writes directly into the returned arrays. A point in a liquid layer, or at a radius without a depth-resolved solution, is NaN. Like the displacement grid, the stress and strain grids carry the time-varying tide only: modes at zero forcing frequency, the permanent tide, are not included. Over a common period of the modes, each component therefore averages to zero, even in a strongly dissipative body. Dissipation appears instead as a lag of the strain behind the stress and as the positive mean power that `calc_3d_tides` returns.

#### Threads

Every grid method (`get_3d_tidal_heating_array`, `calc_3d_tides`, `calc_3d_displacements`, and `calc_3d_stress_strain`) takes `num_threads`. The default, 0, uses the logical processors less 4, with a minimum of 1. The work is split as follows:

- The radial solves, one per degree and frequency, run first. They use up to `num_threads` threads once there are at least `[numerical]` `love_solve_min_parallel` of them (default 3).
- The per-point evaluation then runs over colatitude rows on up to `num_threads` threads. Rows that add into the same cells, as when colatitude is summed, are combined in row order.
- The analytic colatitude collapse of `calc_3d_tides`, the default when the secular heating is summed over colatitude, has no per-point grid. It spreads its radii over the threads instead, at least `[numerical]` `tides_3d_min_radii_per_thread` (default 8) to a thread.

The result is identical for any thread count.

Inside a process or thread pool whose workers already occupy the machine, pass `num_threads=1`. Notebook 13 times three grids on one thread and on every core.

```python
import os

grid = world.calc_3d_tides(
    orbital_frequency,
    spin_frequency,
    eccentricity,
    obliquity,
    semi_major_axis,
    host_mass,
    radii=radii,
    colatitudes=colatitudes,
    longitudes=longitudes,
    num_threads=os.cpu_count())['heating']   # The same values as num_threads=1, computed on every core
```

### Engine and Kernel Access

`tidal_potential_3d_modes` exposes the potential engine for custom pipelines. It returns one entry per raw `(l, m, p, q)` mode, before the coherent merge:

```python
from TidalPy.constants import G
from TidalPy.Tides.potential import tidal_potential_3d_modes

colatitude, longitude = 0.8, 0.3
degrees, freqs, pots = tidal_potential_3d_modes(
    world.radius, orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis,
    host_mass, G, colatitude, longitude,
    min_degree_l=2, max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)
# pots[i] = complex (U, dU/dtheta, dU/dphi, d2U/dtheta2, d2U/dphi2, d2U/dtheta_dphi) for mode i
```

`Tides.multilayer.stress_strain` exposes the compiled kernel of the world methods point by point:

- `strain_stress_heating_point` returns the six complex strain and six complex stress amplitudes at a point for one mode, and that mode's own heating [W m-3] at the frequency given.
- `displacement_point` returns the complex displacement amplitudes $u_r = y_1 U$, $u_\theta = y_3\, \partial U / \partial\theta$, and $u_\phi = y_3\, (\partial U / \partial\phi) / \sin\theta$ [m].
- `volumetric_heating(stress, strain, frequency)` returns the cycle-averaged heating `(|omega| / 2) |sum_k w_k Im(sigma_k conj(eps_k))|` [W m-3] of amplitudes at the frequency `omega`, with `w_k = 2` on the three off-diagonal components, the same factor the world path applies.

The first two take one row from `tidal_potential_3d_modes` together with the radial functions and complex moduli at the point, which the world provides after a radial-solver Love solve. A real row is treated as a phasor with zero phase. Assembling several modes follows the rules of the world methods:

- Solve the radial problem and evaluate the moduli at each mode's degree and `|omega|`, and conjugate the row of a mode with `omega < 0` so that its amplitudes sit at `+|omega|`.
- Sum the strain and stress amplitudes of every mode that shares one `|omega|` before forming the heating. `volumetric_heating` of those sums at `|omega|` is the secular heating at that frequency [W m-3], and the frequencies add.
- A field at time $t$ is $\mathrm{Re}\left[A\, e^{i|\omega| t}\right]$, with $A$ the complex amplitude, summed over the modes.

The pointwise secular density of `calc_3d_tides` can be rebuilt this way:

```python
from TidalPy.constants import min_frequency
from TidalPy.Tides.multilayer.stress_strain import strain_stress_heating_point, volumetric_heating

radius = 0.9 * world.radius
frequency_totals = dict()   # |omega| in units of the mean motion -> [|omega|, summed strain, summed stress]
for degree_l, frequency, row in zip(degrees, freqs, pots):
    magnitude = abs(frequency)
    if magnitude <= min_frequency:
        continue   # A static mode does not dissipate

    # Radial functions and moduli at this mode's degree and |omega|
    world.solve_love_numbers(
        frequency=magnitude,
        degree_l=int(degree_l))
    y = np.array([world.get_love_radial_y(radius, 0, y_index) for y_index in range(6)])
    strain, stress, _ = strain_stress_heating_point(
        y,
        world.calc_complex_shear_modulus(radius, magnitude),
        world.calc_complex_bulk_modulus(radius, magnitude),
        radius,
        float(degree_l),
        magnitude,
        True,    # The layer is solid
        False,   # The layer is compressible
        row if frequency > 0.0 else np.conj(row),   # Amplitudes at +|omega|
        colatitude)

    # Modes that share |omega| are summed before the heating is formed
    totals = frequency_totals.setdefault(round(magnitude / orbital_frequency, 9), [magnitude, 0.0, 0.0])
    totals[1] = totals[1] + strain
    totals[2] = totals[2] + stress

heating_from_kernel = sum(
    volumetric_heating(stress_sum, strain_sum, magnitude)
    for magnitude, strain_sum, stress_sum in frequency_totals.values())   # [W m-3]

# The same density from the world method
heating_from_world = world.calc_3d_tides(
    orbital_frequency,
    spin_frequency,
    eccentricity,
    obliquity,
    semi_major_axis,
    host_mass,
    radii=np.array([radius]),
    colatitudes=np.array([colatitude]),
    longitudes=np.array([longitude]))['heating'][0, 0, 0]
```

## References

- Kaula, W. M. (1964). Tidal dissipation by solid friction and the resulting orbital evolution. *Reviews of Geophysics*, 2(4), 661-685.
- Efroimsky, M., and Williams, J. G. (2009). Tidal torques: A critical review of some techniques. *Celestial Mechanics and Dynamical Astronomy*, 104, 257-289.
- Takeuchi, H., and Saito, M. (1972). Seismic Surface Waves. In *Methods in Computational Physics: Advances in Research and Applications*, 11, 217-295.
- Tobie, G., Mocquet, A., and Sotin, C. (2005). Tidal dissipation within large icy satellites: Applications to Europa and Titan. *Icarus*, 177(2), 534-549. The displacement, strain, and heating forms.
- Kervazo, M., Tobie, G., Choblet, G., Dumoulin, C., and Běhounková, M. (2021). Solid tides in Io's partially molten interior. *Astronomy and Astrophysics*, 650, A72. Appendix D, the corrected strain components.
