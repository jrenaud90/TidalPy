# 3D Tidal Stress, Strain, and Heating (`Tides.multilayer`)

_Updated: 2026-10-09_

`TidalPy.Tides.multilayer` calculates the depth- and direction-resolved tidal response of a layered world, point by point: strain, stress, displacement, and volumetric heating. Use it for heating maps, radial profiles, per-layer budgets, and stress or displacement fields. For the total heating alone the [global (1D) path](global_tides.md) is enough; at zero obliquity the two agree.

## Quick Example

```python
import numpy as np

from TidalPy.Material import Material, Phase
from TidalPy.Rheology import Elastic, Maxwell
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds import TerrestrialWorld
from TidalPy.Tides.classes import make_tide

radius  = 1.8e6    # [m]
density = 3500.0   # [kg m-3]
mass    = (4.0 / 3.0) * np.pi * radius**3 * density

rock = Phase(
    eos={"model": "constant", "reference_density_kg_m3": density, "bulk_modulus_pa": 1.0e11},
    shear_modulus={"model": "constant", "shear_modulus_pa": 6.0e10},
    shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e19},
    bulk_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e30})
layer = Layer(
    "mantle",
    0,
    0.0,
    radius,
    mass,
    material=Material(solid=rock),
    shear_rheology=Maxwell(),
    bulk_rheology=Elastic())

world = TerrestrialWorld("io_like", radius, mass)
world.add_layer(layer)
world.solve_eos()
world.set_tide_model(make_tide("rheology"))
# Potential truncation (also the [tides] TOML table)
world.set_tide_config(max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)

orbital_frequency = 4.1e-5      # [rad s-1]
spin_frequency = 5.6e-5         # non-synchronous
eccentricity = 0.004
obliquity = 0.0
semi_major_axis = 4.2e8         # [m]
host_mass = 1.9e27              # [kg]
orbit = (orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass)

h_bar = world.get_3d_tidal_heating(*orbit, 0.9 * world.radius, 0.8)   # radius, colatitude; [W m-3], secular
```

The 3D methods need the rheology tide model and a solved EOS; analytic tide models such as fixed-Q are rejected. Each takes the six orbital-state arguments of `calc_tides` first; any left as `None` comes from `get_tide_state()` when the world is in a `System` with a tidal host; then give the point or grid by keyword (`world.get_3d_tidal_heating(radius=0.9 * world.radius, colatitude=0.8)`). The orbital state is checked as in `calc_tides`: $e$ outside $[0, 1)$ or a non-positive semi-major axis raises `ValueError`, and the truncation-range warnings apply. The `[tides]` settings `max_degree_l` (2..10), `eccentricity_trunc_lvl`, and `obliquity_trunc_lvl` (0 = off) truncate the potential; see [Truncation](#truncation).

## Common Tasks

### Zonal-Mean Maps: `get_3d_tidal_heating_array`

`get_3d_tidal_heating` re-solves the radial (Love-number) response at every point. `get_3d_tidal_heating_array` takes paired, equal-length `(radius, colatitude)` arrays, or a scalar for one of them that pairs with every value of the other. It solves the radial response once per unique `(l, frequency)` and returns an array of the same shape, NaN where a radius has no depth-resolved solution:

```python
radii = np.linspace(0.01 * world.radius, 0.999 * world.radius, 40)
colatitudes = np.linspace(1e-3, np.pi - 1e-3, 60)
radius_grid, colatitude_grid = np.meshgrid(radii, colatitudes, indexing="ij")
hbar = world.get_3d_tidal_heating_array(
    *orbit,
    radius_grid.ravel(),
    colatitude_grid.ravel()
).reshape(radius_grid.shape)   # [W m-3], secular
```

It matches the scalar method to machine precision. Both return the longitude mean of the secular density. For a longitude-resolved map, use `calc_3d_tides`.

### Grids, Profiles, and Totals: `calc_3d_tides`

`calc_3d_tides` returns the heating on a full `(radius, colatitude, longitude[, time])` grid, or reduced along any spatial axis. It computes one of:

- `orbit_averaged=True` (default): the secular density $\bar{h}$ \[W m$^{-3}$\], the time average of the instantaneous power (see [Secular Heating](#secular-heating) for its longitude dependence).
- `orbit_averaged=False`: the instantaneous power density $\sigma_{ij}(t)\,\dot{\varepsilon}_{ij}(t)$ \[W m$^{-3}$\] at each of `times` (required), on a 4th axis. Every cross term is present, so its time average equals $\bar{h}$ through the truncation level's $e^N$. Past $e^N$ it keeps partial terms that $\bar{h}$ cuts, so the two differ at large eccentricity and a low truncation level.

Non-summed spatial axes take `radii`, `colatitudes`, and `longitudes`. When any axis is summed (`radial_summed`, `latitude_summed`, `longitude_summed`), the surviving axes carry their Jacobian ($r^2$, $\sin\theta$, $1$), so a plain integral over them gives the total. With none summed, the output is the raw density. The dict holds the surviving axes plus `heating` (radius, colatitude, longitude, time order), or, with all three summed, `total` and `per_layer` \[W\]. The per-layer totals replace the `tidal_scale` split of the 1D path.

```python
longitudes = np.linspace(0.0, 2.0 * np.pi, 60)

# Secular density grid (nr, ncolat, nlon); longitude-dependent at synchronous rotation
grid = world.calc_3d_tides(
    *orbit,
    radii=radii,
    colatitudes=colatitudes,
    longitudes=longitudes
)['heating']

# Per-layer and total heating
res = world.calc_3d_tides(
    *orbit,
    latitude_summed=True,
    longitude_summed=True,
    radial_summed=True
)

res['total']      # [W], equals get_tidal_heating at zero obliquity
res['per_layer']  # [W] per layer (innermost first), sums to res['total']

# Radial power profile dP/dr
prof = world.calc_3d_tides(
    *orbit,
    radii=np.linspace(0.01 * world.radius, world.radius, 400),   # Below the solver's start would be NaN
    latitude_summed=True,
    longitude_summed=True
)
np.trapezoid(prof['heating'], prof['radii'])   # The total, minus the share inside 0.01 R

# Instantaneous power over one orbital period
period = 2.0 * np.pi / orbital_frequency
res = world.calc_3d_tides(
    *orbit,
    radii=radii,
    colatitudes=colatitudes,
    longitudes=longitudes,
    times=np.linspace(0.0, period, 400),
    orbit_averaged=False
)
res['heating']   # shape (nr, ncolat, nlon, ntime)
```

The secular colatitude integral is analytic by default ([Volume Integral and Collapse](#volume-integral-and-collapse)), a few times faster than quadrature on a large radius grid. It needs longitude summed too. A call that keeps longitudes, sets `latitude_analytic=False`, or uses a colatitude band uses `latitude_nodes` Gauss-Legendre nodes (default 16; a homogeneous degree-2 body needs about 4). Both agree to machine precision. The longitude integral is analytic when averaged and a `longitude_nodes` trapezoid when instantaneous. The radial integral uses `radial_slices` nodes per layer (default 16).

With `latitude_summed`, `colatitude_min` and `colatitude_max` \[rad\] (defaults 0 and $\pi$) restrict the colatitude integral. Complementary bands add up to the full sphere, so a zonal budget (polar caps against an equatorial belt) takes a few calls.

### Displacements: `calc_3d_displacements`

Each coherent wave's displacement amplitude $\mathbf{u}$ ([Displacement, Strain, and Stress](#displacement-strain-and-stress)) evolves as $\mathrm{Re}\left[\mathbf{u}\, e^{i|\omega| t}\right]$, summed over the waves. `calc_3d_displacements` returns the three components \[m\] on the full `(radius, colatitude, longitude, time)` grid:

```python
out = world.calc_3d_displacements(
    *orbit,
    radii=[0.9 * world.radius, world.radius],
    colatitudes=np.linspace(0.1, np.pi - 0.1, 30),
    longitudes=np.linspace(0.0, 2 * np.pi, 60),
    times=np.linspace(0.0, period, 12)
)
out["radial"].shape        # (2, 30, 60, 12)   u_r [m]; also out["polar"], out["azimuthal"]
```

At the surface $y_1 = h / g$ and $y_3 = l / g$, so a single mode's surface radial displacement is $h U / g$. This needs a radial-solver Love-number method (the analytic `homogeneous`, `cpl`, and `ctl` methods have no radial functions). After such a solve, `world.get_love_radial_y(radius, ytype_idx, y_idx)` gives the radial functions at any radius.

### Stress and Strain: `calc_3d_stress_strain`

Each component at time $t$ is the sum over frequencies of $\mathrm{Re}\left[A\, e^{i|\omega| t}\right]$, with $A$ the summed complex amplitude of the waves at that frequency. Both tensors come on the full `(radius, colatitude, longitude, time)` grid, with the components `rr`, `theta_theta`, `phi_phi`, `r_theta`, `r_phi`, `theta_phi` (also returned as `components`) on the last axis. The strain is the symmetric gradient of the `calc_3d_displacements` field.

```python
tensors = world.calc_3d_stress_strain(
    *orbit,
    radii=[0.9 * world.radius, world.radius],
    colatitudes=np.linspace(0.1, np.pi - 0.1, 30),
    longitudes=np.linspace(0.0, 2.0 * np.pi, 60),
    times=np.linspace(0.0, period, 12)
)
tensors['stress'].shape                  # (2, 30, 60, 12, 6) [Pa]; tensors['strain'] has the same shape
radial_stress = tensors['stress'][..., 0]   # The rr component
peak_stress = np.abs(tensors['stress']).max(axis=3)   # Largest magnitude of each component over the times
```

Each tensor takes 48 bytes per point and time; skip one with `return_stress=False` or `return_strain=False`. Displacements, stress, and strain carry the time-varying tide only, without the permanent (zero-frequency) tide, so each component averages to zero over a common period of the modes, even in a dissipative body. Dissipation shows as the strain's lag behind the stress and as the mean power of `calc_3d_tides`.

### Threads

Every grid method takes `num_threads`. The default, 0, uses the logical processors minus 4, at least 1. The radial solves go parallel once there are at least `[numerical]` `love_solve_min_parallel` of them (default 3), then the points are split by colatitude row. The analytic colatitude collapse splits its radii instead, at least `[numerical]` `tides_3d_min_radii_per_thread` (default 8) per thread. The result is identical for any thread count. Inside a pool whose workers already occupy the machine, pass `num_threads=1`. Demo notebook P08 compares timings.

```python
import os

grid = world.calc_3d_tides(
    *orbit,
    radii=radii,
    colatitudes=colatitudes,
    longitudes=longitudes,
    num_threads=os.cpu_count()
)['heating']   # The same values as num_threads=1
```

## Physics

The response at $(r, \theta, \phi, t)$ factorizes. The radial part is the viscoelastic-gravitational functions $y_1 \dots y_6$ of the radial solver (Takeuchi and Saito 1972) and the complex moduli $\mu^{*}$ and $K^{*}$, at $r$ and the solved frequency (in Python, `RadialSolverSolution.get_radial_solution` and `get_complex_shear_modulus`). The angular and time part is the tidal potential $W(\theta, \phi, t)$ and its first and second $\theta$ and $\phi$ derivatives. The radial problem depends on $l$ and $|\omega|$ but not $m$, so it is solved once per $(l, |\omega|)$.

### Tidal Potential

The potential is the Kaula expansion of [Global Tides](global_tides.md) (Kaula 1964; Efroimsky and Williams 2009, Eq. 18) at the surface radius $R$. It is linear in $F_{lmp}$, $G_{lpq}$, and $P_{lm}$; the global path squares $F$ and $G$ because its heating scales with the potential squared. For each mode the engine returns $l$, the signed frequency, and the complex amplitude of $W$ and its derivatives, with $W(t) = \mathrm{Re}[W_c\,e^{i\omega t}]$:

$$W(\theta, \phi, t) = \sum_{l,m,p,q} A_{lmpq}\,P_{lm}(\cos\theta)\,\mathcal{T}_{lm}\!\left(\omega_{lmpq}t - m\phi\right), \qquad A_{lmpq} = \frac{G M_h}{a}\left(\frac{R}{a}\right)^{l}\frac{(l-m)!}{(l+m)!}\left(2-\delta_{0m}\right)F_{lmp}(I)\,G_{lpq}(e),$$

$$\omega_{lmpq} = (l - 2p + q)\,n - m\,\dot{\theta},$$

where $n$ is the mean motion, $\dot{\theta}$ the spin rate, and $\mathcal{T}_{lm}$ is $\cos$ for even $l - m$ and $\sin$ for odd $l - m$. $\phi$ is the body-fixed east longitude, increasing in the direction of rotation, and $\phi = 0$ is the host's meridian at $t = 0$, when the host is at its ascending node and at periapse. As in Kaula, $P_{lm}$ has no Condon-Shortley phase; TidalPy's Legendre tables (`TidalPy.Utilities.legendre`) include it, so each amplitude carries a canceling factor $(-1)^{m}$. All radial dependence is in $y(r)$.

> [!WARNING]
> This assumes no periapse or node precession. It also assumes that the change in the mean anomaly can be approximated by the mean motion.

### Truncation

Eccentricity level $N$ keeps the potential through $e^N$ and, as in the 1D path, cuts every product of two eccentricity functions in the secular heating at $e^N$ ([Eccentricity Functions](Eccentricity.md)). Obliquity level $N$ does the same in $I$ ([Obliquity Functions](Obliquity.md)). So at zero obliquity the volume integral equals the 1D heating at any eccentricity and spin rate; with $e$ and $I$ both nonzero they differ by the cross terms of [Coherent Waves](#coherent-waves). The instantaneous fields are linear in the potential and use the unsquared functions. A nonzero obliquity truncation turns on the odd-$m$ harmonics ($P_{21}$, ...). A mode of zero frequency is static and dropped. As in the 1D path, a mode of frequency $\omega$ near or below the world's continuation frequency $\omega_c$ (`world.calc_continuation_frequency()`, see [Numerical Settings](../Overview/2_TidalPy_Configurations.md#numerical-settings)) takes the radial solution at $\sqrt{\omega^2 + \omega_c^2}$ with its heating scaled by $|\omega|$ over that frequency, and its instantaneous power scales the strain rate the same way so that it averages to the secular heating.

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

where $dy_{1}/dr$ comes from the radial equations of the layer's type. The stress follows the isotropic law (Takeuchi and Saito 1972), with both moduli at $r$ and the mode's frequency,

$$\sigma_{ij} = 2\mu^{*}\varepsilon_{ij} + \lambda^{*}\varepsilon_{kk}\,\delta_{ij}, \qquad \lambda^{*} = K^{*} - \tfrac{2}{3}\mu^{*}.$$

The isotropic part is taken from the radial stress $\sigma_{rr} = y_{2}\,W$:

$$\lambda^{*}\varepsilon_{kk} = \left(y_{2} - 2\mu^{*}\frac{dy_{1}}{dr}\right)W.$$

In a compressible layer this is the same quantity, since $y_{2} = \lambda^{*}\left(dy_{1}/dr + (2y_{1} - l(l+1)y_{3})/r\right) + 2\mu^{*}\,dy_{1}/dr$, and it stays accurate for a very large bulk modulus, where $\lambda^{*}$ times a nearly zero trace would amplify round-off. In a layer the radial solver treats as incompressible, the strain is traceless and this term is the pressure, which only $y_{2}$ carries.

### Coherent Waves

Before the heating is formed, each raw $(l, m, p, q)$ mode is mapped onto its non-negative frequency: a mode with $\omega < 0$ contributes the conjugate of its phasor at $|\omega|$, since $\mathrm{Re}[W_c e^{i\omega t}] = \mathrm{Re}[\overline{W_c}\,e^{-i\omega t}]$. It is then merged with every mode that shares its real spatial function: the same $l$, $m$, $|\omega|$, and azimuthal sign ($e^{\mp im\phi}$, irrelevant for $m = 0$). The kernel works with these coherent waves.

The $m = 0$ modes come in pairs, $(l, 0, p, q)$ at $+\omega$ and $(l, 0, l-p, -q)$ at $-\omega$, with equal amplitudes ($F_{l0p} = \pm F_{l0,l-p}$ with the parity sign, $G_{lpq} = G_{l,l-p,-q}$) and the same time dependence, since $\cos(-x) = \cos x$. Together they are one sinusoid of twice the amplitude. Heating goes as amplitude squared, so summing their powers separately would lose half the zonal heating: for a homogeneous degree-2 body at zero obliquity the zonal terms are 9/84 of the total, so 4.5/84 = 5.36% of a synchronous body's heating (where only the eccentricity modes survive). The 1D formula counts the pair through its $(2 - \delta_{0m})$ weight instead.

At nonzero obliquity, modes of one $(l, m)$ with different $(p, q)$ can also share a signed frequency. Their relative phase is set by the argument of periapsis $\omega$, which the engine takes as zero, so they also add coherently. The 1D formula averages over $\omega$ and drops these cross terms. The 3D total is thus the heating of an orbit whose periapsis stays at the ascending node, and the 1D value that of an orbit whose periapsis precesses. The difference needs $e$ and $I$ both nonzero, varies with $\cos 2\omega$ and $\cos 4\omega$, and changes sign at $\omega = \pi/2$. Both values match a direct calculation from the host's exact Kepler position.

The difference is largest in synchronous rotation, where the semidiurnal tide is static. At degree 2, for a homogeneous Maxwell body and the bundled Io in synchronous rotation, the 3D total is 0.16 to 0.38% below the 1D heating at $e = 0.05$ ($I$ = 0.1 to 0.3 rad) and about 1.6% below at $e = 0.2$, $I = 0.3$ rad; at spins of $1.5n$ and $2.5n$ the difference is at most 0.16%, of either sign. Use the 1D heating when the periapsis precesses over the time of interest.

### Secular Heating

The secular (orbit-averaged) heating \[W m$^{-3}$\] is formed from the complex amplitudes, so the average is exact with no time grid (Tobie et al. 2005):

$$\bar{h}(r, \theta, \phi) = \sum_{\chi}\frac{\chi}{2}\sum_{k} w_{k}\,\mathrm{Im}\!\left(\sigma_{k}^{(\chi)}\,\overline{\varepsilon_{k}^{(\chi)}}\right), \qquad w_{k} = \begin{cases} 1, & k \in \{rr, \theta\theta, \phi\phi\}, \\ 2, & k \in \{r\theta, r\phi, \theta\phi\}, \end{cases}$$

where $\sigma_{k}^{(\chi)}$ and $\varepsilon_{k}^{(\chi)}$ are the stress and strain amplitudes summed over every wave at $\chi = |\omega|$. The form is non-negative for a dissipative material. Cross terms between different frequencies average to zero and are dropped; those at one frequency are kept, and make $\bar{h}$ depend on longitude when the waves differ in longitude structure. In synchronous rotation every mode sits at a multiple of $n$, and the zonal and sectoral waves give the $\cos 2\phi$ and $\cos 4\phi$ patterns of a synchronous heating map, symmetric about the sub-host meridian. For a generic non-synchronous spin each frequency carries one wave and $\bar{h}$ depends on $r$ and $\theta$ only.

### Volume Integral and Collapse

The total power is

$$\dot{E} = \int_{0}^{R}\int_{0}^{\pi}\int_{0}^{2\pi}\bar{h}\,r^{2}\sin\theta\;d\phi\,d\theta\,dr,$$

which `calc_3d_tides` evaluates one axis at a time:

- Longitude: waves with different signed azimuthal wavenumbers $\mu = \pm m$ integrate to zero, so the integral is $2\pi$ times the sum over $(\chi, \mu)$ groups at $\phi = 0$.
- Colatitude: every strain and stress component of a wave combines, with complex radial coefficients, six real functions of $\theta$ for its $(l, m)$,

  $$f_{1} = P_{lm}, \quad f_{2} = \frac{dP_{lm}}{d\theta}, \quad f_{3} = \frac{d^{2}P_{lm}}{d\theta^{2}}, \quad f_{4} = \frac{P_{lm}}{\sin\theta}, \quad f_{5} = -\frac{m^{2}P_{lm}}{\sin^{2}\theta} + \cot\theta\,\frac{dP_{lm}}{d\theta}, \quad f_{6} = \frac{1}{\sin\theta}\left(\frac{dP_{lm}}{d\theta} - \cot\theta\,P_{lm}\right),$$

  so the colatitude integral of any pair of waves reduces to the Gram matrices

  $$\mathcal{G}_{ij}(l_{a}, l_{b}, m) = \int_{0}^{\pi} f_{i}^{(l_{a})}\,f_{j}^{(l_{b})}\,\sin\theta\;d\theta.$$

  They are tabulated for equal degrees and computed with 32-node Gauss-Legendre quadrature in $\cos\theta$ for different degrees, which is exact because the integrands are polynomials in $\cos\theta$ of degree at most $l_{a} + l_{b} + 2 \le 22$.
- Radius: Gauss-Legendre quadrature with weight $r^{2}$ and `radial_slices` nodes strictly inside each layer. A node on a boundary would read the modulus and radial solution of the layer below, so a layer on a stiffer one would lose part of its heating. Nodes below the radial solver's starting radius, where no degree has a solution, are left out; a warning is logged when they hold more of the body's volume than the solver's `rtol`.

At zero obliquity the volume integral equals the 1D heating (`get_tidal_heating`) to the radial quadrature error, at every eccentricity and spin rate. That error is below 1e-4 from 4 nodes per layer, both for a homogeneous body and for the bundled Io (whose thin asthenosphere has a sixtieth of the shear modulus of the mantle beneath it), and below 1e-7 for the homogeneous Io of demo P06 at the default 16. The default adds margin for higher degrees.

## Engine and Kernel Access

`tidal_potential_3d_modes` exposes the potential engine. It returns one entry per raw $(l, m, p, q)$ mode, before the coherent merge:

```python
from TidalPy.constants import G
from TidalPy.Tides.potential import tidal_potential_3d_modes

colatitude, longitude = 0.8, 0.3
degrees, freqs, pots = tidal_potential_3d_modes(
    world.radius,
    orbital_frequency,
    spin_frequency,
    eccentricity,
    obliquity,
    semi_major_axis,
    host_mass,
    G,
    colatitude,
    longitude,
    min_degree_l=2,
    max_degree_l=2,
    eccentricity_truncation=6,
    obliquity_truncation=0
)
# pots[i] = complex (U, dU/dtheta, dU/dphi, d2U/dtheta2, d2U/dphi2, d2U/dtheta_dphi) for mode i
```

`Tides.multilayer.stress_strain` exposes the point kernel:

- `strain_stress_heating_point`: the six complex strain and six complex stress amplitudes of one mode at a point, and that mode's own heating \[W m$^{-3}$\].
- `displacement_point`: the complex displacements $u_r = y_1 U$, $u_\theta = y_3\, \partial U / \partial\theta$, and $u_\phi = y_3\, (\partial U / \partial\phi) / \sin\theta$ \[m\].
- `volumetric_heating(stress, strain, frequency)`: the cycle-averaged heating $(|\omega| / 2)\,\sum_k w_k\,\mathrm{Im}\left(\sigma_k\,\overline{\varepsilon_k}\right)$ \[W m$^{-3}$\], with $w_k = 2$ on the off-diagonal components. It is signed: a negative value flags a modulus with $\mathrm{Im}(\mu) < 0$ or a wrong-sign solution.

The first two take one row of `tidal_potential_3d_modes` (a real row is a phasor with zero phase) and the radial functions and complex moduli at the point. To combine modes as the world methods do, conjugate the row of a mode with $\omega < 0$, sum the amplitudes of all modes at one $|\omega|$ before calling `volumetric_heating`, and add the frequencies. This rebuilds a point of `calc_3d_tides`:

```python
from TidalPy.Tides.multilayer.stress_strain import strain_stress_heating_point, volumetric_heating

radius = 0.9 * world.radius
frequency_totals = dict()   # |omega| / n -> [|omega|, summed strain, summed stress]
for degree_l, frequency, row in zip(degrees, freqs, pots):
    magnitude = abs(frequency)
    if magnitude == 0.0:
        continue   # A static mode does not dissipate

    # Radial functions and moduli at this mode's degree and |omega|
    world.solve_love_numbers(
        frequency=magnitude,
        degree_l=int(degree_l)
    )
    y = np.array([world.get_love_radial_y(radius, 0, y_index) for y_index in range(6)])
    strain, stress, _ = strain_stress_heating_point(
        y,
        world.calc_complex_shear_modulus(radius, magnitude),
        world.calc_complex_bulk_modulus(radius, magnitude),
        radius,
        float(degree_l),
        magnitude,
        True,    # Solid layer
        False,   # Compressible layer
        row if frequency > 0.0 else np.conj(row),   # Amplitudes at +|omega|
        colatitude
    )

    # Sum modes that share |omega| before forming the heating
    totals = frequency_totals.setdefault(round(magnitude / orbital_frequency, 9), [magnitude, 0.0, 0.0])
    totals[1] = totals[1] + strain
    totals[2] = totals[2] + stress

heating_from_kernel = sum(
    volumetric_heating(stress_sum, strain_sum, magnitude)
    for magnitude, strain_sum, stress_sum in frequency_totals.values())   # [W m-3]

heating_from_world = world.calc_3d_tides(
    *orbit,
    radii=np.array([radius]),
    colatitudes=np.array([colatitude]),
    longitudes=np.array([longitude])
)['heating'][0, 0, 0]   # The same density
```

## Limits and Failure Modes

- Liquids: a liquid layer, or a liquid zone of a melting layer ([Pieces and Zones](../Structures/worlds/worlds.md#pieces-and-zones)), has no shear dissipation: heating 0, stress and strain NaN.
- Poles: values at $\sin\theta$ within machine epsilon of 0 are NaN. Use colatitudes inside $(0, \pi)$.
- Center: the solver's starting radius grows with degree, so the innermost region carries only the degrees solved there. A radius is NaN only where no degree has a solution, in every output with a radius axis.
- Colatitude bands have no effect unless colatitude is summed.

## References

- Kaula, W. M. (1964). Tidal dissipation by solid friction and the resulting orbital evolution. *Reviews of Geophysics*, 2(4), 661-685.
- Efroimsky, M., and Williams, J. G. (2009). Tidal torques: A critical review of some techniques. *Celestial Mechanics and Dynamical Astronomy*, 104, 257-289.
- Takeuchi, H., and Saito, M. (1972). Seismic Surface Waves. In *Methods in Computational Physics: Advances in Research and Applications*, 11, 217-295.
- Tobie, G., Mocquet, A., and Sotin, C. (2005). Tidal dissipation within large icy satellites: Applications to Europa and Titan. *Icarus*, 177(2), 534-549. The displacement, strain, and heating forms.
- Kervazo, M., Tobie, G., Choblet, G., Dumoulin, C., and Běhounková, M. (2021). Solid tides in Io's partially molten interior. *Astronomy and Astrophysics*, 650, A72. Appendix D, the corrected strain components.
