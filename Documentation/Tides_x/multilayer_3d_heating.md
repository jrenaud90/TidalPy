# 3D Tidal Stress, Strain, and Heating (`Tides_x.multilayer`)

_Updated: 2026-09-23_

This module computes the depth- and direction-resolved tidal response (the complex strain and stress tensors and the volumetric heating) of a layered world. The response is evaluated at a single point on demand, so a map is built only when the caller evaluates a set of points.

## Physics

The response at a point $(r, \theta, \phi, t)$ factorizes into a radial part as well as an angular and time part. The radial part is the set of viscoelastic-gravitational functions $y_1 \dots y_6$ from the radial solver (Takeuchi and Saito 1972), evaluated at $r$ through the dense calling system (`RadialSolverSolution.get_radial_solution`), together with the complex shear and bulk moduli $\mu^{*}$ and $K^{*}$ at $r$, read at the solved frequency through the dense material path (`c_EOSSolution::call_material`, exposed in Python as `RadialSolverSolution.get_complex_shear_modulus`). The angular and time part is the tidal potential $W(\theta, \phi, t)$ and its first and second $\theta$ and $\phi$ derivatives. The radial problem depends on the degree $l$ and the frequency $|\omega|$ only, not on $m$, so it is solved once per $(l, |\omega|)$.

### Tidal Potential

The potential is the Kaula expansion of the [global tides](global_tides.md) page (Kaula 1964; Efroimsky and Williams 2009, Eq. 18), evaluated at the surface radius $R$ and kept linear in $F_{lmp}$, $G_{lpq}$, and $P_{lm}$; the global path squares $F$ and $G$ because global heating goes as the potential squared. The engine (`c_tidal_potential_3d_modes`, wrapped as `TidalPy.Tides_x.potential.tidal_potential_3d_modes`) enumerates the active modes and returns, for each, the degree $l$, the signed forcing frequency, and the complex amplitude of $W$ and its derivatives, with $W(t) = \mathrm{Re}[W_c\,e^{i\omega t}]$:

$$W(\theta, \phi, t) = \sum_{l,m,p,q} A_{lmpq}\,P_{lm}(\cos\theta)\,\mathcal{T}_{lm}\!\left(\omega_{lmpq}t - m\phi\right), \qquad A_{lmpq} = \frac{G M_h}{a}\left(\frac{R}{a}\right)^{l}\frac{(l-m)!}{(l+m)!}\left(2-\delta_{0m}\right)F_{lmp}(I)\,G_{lpq}(e),$$

$$\omega_{lmpq} = (l - 2p + q)\,n - m\,\dot{\theta},$$

where $n$ is the mean motion, $\dot{\theta}$ the spin rate, and $\mathcal{T}_{lm}$ is $\cos$ for even $l - m$ and $\sin$ for odd $l - m$. Longitude $\phi$ is the body-fixed east longitude, increasing in the direction of rotation, and $\phi = 0$ is the host's meridian at $t = 0$, when the host is at its ascending node and at periapse. $P_{lm}$ carries no Condon-Shortley phase, as in Kaula; TidalPy's Legendre tables (`TidalPy.Utilities_x.legendre`) include it, so each amplitude carries $(-1)^{m}$ to cancel it. All radial dependence is carried by $y(r)$.

> [!WARNING]
> This assumes no periapse or node precession. It also assumes that the change in the mean anomaly can be approximated by the mean motion.

The user selects the truncation via three knobs (on the world's `[tides]` config): `max_degree_l` (2..10), `eccentricity_trunc_lvl`, and `obliquity_trunc_lvl` (0 = off). A nonzero obliquity truncation turns on the odd-`m` (`P_21`, ...) harmonics automatically. A mode whose $|\omega|$ does not exceed `[numerical] minimum_frequency` (1e-14 rad/s, `TidalPy.constants.min_frequency`) is switched off later, the same floor the 1D `calc_tides` path uses. (`[numerical] min_spin_orbit_diff` belongs to the classic backend; the new backend does not read it.)

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

with both moduli taken at $r$ and at the mode's frequency.

> [!NOTE]
> In a layer the radial solver treats as incompressible, the strain has no trace, so this law drops the pressure: the normal stresses ($\sigma_{rr}$, $\sigma_{\theta\theta}$, $\sigma_{\phi\phi}$) miss their isotropic part, and $\sigma_{rr}$ came out about 13 percent below that of a nearly incompressible compressible twin in a test. The deviatoric stress, the strain, and so the heating are unaffected (the 1D and 3D heating of such a layer agree to about 1e-4). Solve the layer as compressible when its normal stresses matter; carrying the pressure from $y_{2}$ is planned.

The kernel is a solid-layer computation: a liquid (a liquid layer, or a molten stretch of a solid layer that the radial solver treats as a static liquid) contributes no shear dissipation, so its heating is 0 while its stress and strain are NaN. The poles themselves ($\sin\theta$ within machine epsilon of 0, at both $0$ and $\pi$) are singular points of the angular terms, so point-wise values there are NaN; use colatitudes inside $(0, \pi)$, as the Gauss-Legendre nodes of the summed paths do.

### Coherent Waves

The heating paths do not consume the raw $(l, m, p, q)$ modes one at a time. Each mode is first mapped onto its non-negative frequency (a mode with $\omega < 0$ contributes the complex conjugate of its phasor at $|\omega|$, since $\mathrm{Re}[W_c e^{i\omega t}] = \mathrm{Re}[\overline{W_c}\,e^{-i\omega t}]$) and merged with every other mode that shares its real spatial function: the same degree $l$, order $m$, $|\omega|$, and azimuthal sign ($e^{-im\phi}$ against $e^{+im\phi}$; irrelevant for $m = 0$). The result is the list of coherent waves the kernel works with.

The merge matters because the $m = 0$ modes always come in pairs, $(l, 0, p, q)$ at $+\omega$ and $(l, 0, l-p, -q)$ at $-\omega$, that carry equal amplitudes ($F_{l0p} = \pm F_{l0,l-p}$ with the parity sign, $G_{lpq} = G_{l,l-p,-q}$) and are the same function of time, since $\cos(-x) = \cos x$. They are one real sinusoid of twice the amplitude, and the heating goes as the amplitude squared, so summing their cycle-averaged powers separately loses half of the zonal heating. For a homogeneous degree-2 body at zero obliquity the zonal terms are 9/84 of the total, so the loss is 4.5/84 = 5.36% of the heating of a synchronously rotating body, where only the eccentricity modes survive. The 1D formula counts the same pair through its $(2 - \delta_{0m})$ weighting, which is why it needs no merge.

At nonzero obliquity, modes of the same $(l, m)$ with different $(p, q)$ can also share a signed frequency. Their relative phase is set by the argument of periapse, which the engine takes as zero (no precession), so they too combine coherently. The 1D formula, being averaged over apsidal precession, does not carry that cross term, so the two paths agree only to the size of those terms at nonzero obliquity.

### Secular Heating

The secular (cycle- or orbit-averaged) volumetric heating \[W m$^{-3}$\] is built from the complex amplitudes, so the cycle average is exact with no time grid (Tobie et al. 2005):

$$\bar{h}(r, \theta, \phi) = \sum_{\chi}\frac{\chi}{2}\sum_{k} w_{k}\,\mathrm{Im}\!\left(\sigma_{k}^{(\chi)}\,\overline{\varepsilon_{k}^{(\chi)}}\right), \qquad w_{k} = \begin{cases} 1, & k \in \{rr, \theta\theta, \phi\phi\}, \\ 2, & k \in \{r\theta, r\phi, \theta\phi\}, \end{cases}$$

where $\sigma_{k}^{(\chi)}$ and $\varepsilon_{k}^{(\chi)}$ are the total complex stress and strain amplitudes at the frequency $\chi = |\omega|$, every wave at that frequency summed before the bilinear form. The form is non-negative for a dissipative material. Cross terms between waves at different frequencies average to zero over the orbit and are dropped; cross terms between waves at the same frequency survive the average and are kept. When the waves at one frequency have different longitude structure, as the zonal and sectoral waves of a synchronously rotating body do (every active mode then sits at a multiple of $n$), those cross terms make $\bar{h}$ depend on longitude: the familiar $\cos 2\phi$ and $\cos 4\phi$ patterns of a synchronous heating map, symmetric about the sub-host meridian. Away from such frequency coincidences (a generic non-synchronous spin) every frequency carries one wave and $\bar{h}$ is a function of $r$ and $\theta$ only.

At each point the machinery:
- Builds the active modes from the truncation config and merges them into coherent waves
- Solves the world radial response once per $(l, |\omega|)$
- Sums each frequency's waves into a total complex stress and strain, and accumulates $(\chi/2)\,\mathrm{Im}(\sigma : \overline{\varepsilon})$ per frequency.

### Volume Integral and Collapse

The total power is the volume integral

$$\dot{E} = \int_{0}^{R}\int_{0}^{\pi}\int_{0}^{2\pi}\bar{h}\,r^{2}\sin\theta\;d\phi\,d\theta\,dr,$$

which `calc_3d_tides` evaluates with its summed arguments, one axis at a time:

- Longitude: waves with different signed azimuthal wavenumbers $\mu = \pm m$ integrate to zero over longitude, so the integral is $2\pi$ times the sum over $(\chi, \mu)$ groups evaluated at $\phi = 0$. This is also the longitude mean that `get_3d_tidal_heating(radius, colatitude)` and its batch form return.
- Colatitude: every strain and stress component of a wave is a combination, with complex radial coefficients, of six real functions of $\theta$ for its $(l, m)$,

  $$f_{1} = P_{lm}, \quad f_{2} = \frac{dP_{lm}}{d\theta}, \quad f_{3} = \frac{d^{2}P_{lm}}{d\theta^{2}}, \quad f_{4} = \frac{P_{lm}}{\sin\theta}, \quad f_{5} = -\frac{m^{2}P_{lm}}{\sin^{2}\theta} + \cot\theta\,\frac{dP_{lm}}{d\theta}, \quad f_{6} = \frac{1}{\sin\theta}\left(\frac{dP_{lm}}{d\theta} - \cot\theta\,P_{lm}\right),$$

  so the colatitude integral of any pair of waves reduces to the Gram matrices

  $$\mathcal{G}_{ij}(l_{a}, l_{b}, m) = \int_{0}^{\pi} f_{i}^{(l_{a})}\,f_{j}^{(l_{b})}\,\sin\theta\;d\theta.$$

  They are tabulated for equal degrees and computed with 32-node Gauss-Legendre quadrature in $\cos\theta$ for different degrees, which is exact because the integrands are polynomials in $\cos\theta$ of degree at most $l_{a} + l_{b} + 2 \le 22$. With `latitude_analytic=False` the integral is taken numerically on `latitude_nodes` Gauss-Legendre nodes instead.
- Radius: Gauss-Legendre quadrature inside each layer with `radial_slices` nodes (default 16) and weight $r^{2}$; no node sits on a layer boundary.

The volume integral equals the 1D global tidal heating (`get_tidal_heating`): both describe the same total power. For a homogeneous body at zero obliquity the two agree to the radial quadrature error at every spin rate, synchronous rotation included; for the homogeneous Io of demo notebook 09 that is below 1e-7 with the default 16 nodes per layer. The benchmark tests are `Tests/Test_Structures_x/Test_Worlds/test_world_1d_vs_3d_tides_01.py` and `test_world_3d_tides_coherent_01.py`.

## Python API

### World Method

For a built world, the API lives on the rheology tide model (`c_RheologyTide`); the world delegates to it and everything runs in C++ (the tide model calls the world's radial-solver and EOS members directly).

```python
import numpy as np

from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
from TidalPy.rheology_x import Elastic, Maxwell
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.Tides_x.classes import make_tide
from TidalPy.viscosity_x import make_viscosity

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
world.set_tide_config(max_degree_l=2, eccentricity_truncation=3, obliquity_truncation=0)

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

Requires the rheology tide model and a solved EOS (the analytic tide models like fixed-Q are rejected as they have no depth-resolved solution). When a world is built from a TOML/dict config, the truncation levels flow automatically from the `[tides]` table (`max_degree_l`, `eccentricity_trunc_lvl`, `obliquity_trunc_lvl`). The radial solver's starting radius grows with degree, so with several degrees active the innermost region carries only the degrees that have a solution there; a radius is NaN only where no degree does.

#### Building the Map

Evaluating a heating map point by point through `get_3d_tidal_heating` re-solves the world radial response at every point, although the radial (Love-number) solve depends only on the tidal mode's `(degree l, frequency)` and not on the query position. `get_3d_tidal_heating_array` takes paired, equal-length `(radius, colatitude)` arrays, builds the position-independent tidal-mode list once, and solves the radial response once per unique `(l, frequency)`, reusing it across all points. It returns a same-shape array (NaN where a radius has no depth-resolved solution) and is the efficient way to build a map:

```python
radii = np.linspace(0.01 * world.radius, 0.999 * world.radius, 40)
colatitudes = np.linspace(1e-3, np.pi - 1e-3, 60)
radius_grid, colatitude_grid = np.meshgrid(radii, colatitudes, indexing="ij")
hbar = world.get_3d_tidal_heating_array(
    orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass,
    radius_grid.ravel(), colatitude_grid.ravel()).reshape(radius_grid.shape)   # [W m-3], secular
```

It reproduces the scalar `get_3d_tidal_heating` point-for-point (to machine precision): the scalar path is the batch path with one point. Both return the longitude mean of the secular density, so this is the efficient way to build a zonal-mean map; for a longitude-resolved map use `calc_3d_tides` with a longitude grid.

#### Full Grid and Collapsed Forms: `calc_3d_tides`

`calc_3d_tides` produces the 3D heating as a full grid over `(radius, colatitude, longitude[, time])` or reduced along any spatial dimension. Two quantities:

* `orbit_averaged=True` (default): the secular heating density `h_bar` [W m-3], the time average of the instantaneous power at each point. It depends on longitude wherever waves at one frequency have different longitude structure (synchronous rotation above all) and is constant along longitude otherwise.
* `orbit_averaged=False`: the instantaneous mechanical power density `sigma_ij(t) * eps_dot_ij(t)` [W m-3] at each supplied time (a 4th axis). It depends on longitude and time and time-averages to `h_bar` at the same point (every wave is a real field at `+|omega|`, so every cross term is present).

Non-summed spatial axes take user arrays (`radii`, `colatitudes`, `longitudes`); the `times` array is required when `orbit_averaged=False`. Reduction convention: if any spatial axis is summed, the surviving spatial axes carry their Jacobian (`r^2`, `sin theta`, `1`) so a plain integral over them recovers the total; if none is summed the output is the raw density (its longitude mean matches the scalar `get_3d_tidal_heating`). The colatitude integral (secular, `latitude_summed`) is done analytically by default: the six angular functions the strain/stress needs form a bounded basis whose sphere integrals are precomputed once into a per-`(l, m)` Gram table (`Tides_x.multilayer.stress_strain.angular_gram`), and the cross terms between coherent waves of different degree at one frequency use a cross-degree Gram matrix integrated by Gauss-Legendre quadrature (exact, the integrands being polynomials in `cos theta`), so no colatitude grid is needed on the collapse. This is exact and, for a large radius grid (a radial profile or map), a few times faster than the numerical quadrature (the Love-number solves otherwise dominate). The analytic integral is of the longitude mean, so it is used only when longitude is summed too; a call that keeps its longitudes uses the Gauss-Legendre colatitude grid (`latitude_nodes`) whatever `latitude_analytic` says. Pass `latitude_analytic=False` to use that grid for a longitude-summed call as well; the two agree to machine precision there. The radial integral uses `radial_slices` Gauss-Legendre nodes inside each layer, none on a layer boundary, and the longitude integral the analytic `2*pi` times the longitude mean when averaged or a `longitude_nodes` trapezoid when instantaneous.

The colatitude integral can also be restricted to a latitude band: with `latitude_summed`, pass `colatitude_min` / `colatitude_max` [rad] (defaults 0 and pi) to integrate only `colatitude_min <= theta <= colatitude_max`. Complementary bands add up to the full-sphere result, so zonal heating budgets (polar caps vs an equatorial belt, say) come from a few banded calls. A band narrower than the full sphere always uses the Gauss-Legendre quadrature (the analytic Gram table is full-sphere only); the band has no effect when colatitude is not summed.

The returned dict carries the surviving axes plus either `heating` (the grid over the surviving axes, in the order radius, colatitude, longitude, time) or, when all three spatial axes are summed, `total` [W] and `per_layer` [W]. The C++ code writes the heating straight into the returned arrays:

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

res['total']      # [W], the 1D get_tidal_heating
res['per_layer']  # [W] per layer (innermost first), sums to res['total']

# radial power profile dP/dr (integrates to the total)
prof = world.calc_3d_tides(
    orbital_frequency,
    spin_frequency,
    eccentricity,
    obliquity,
    semi_major_axis,
    host_mass,
    radii=np.linspace(1e3, world.radius, 400),
    latitude_summed=True,
    longitude_summed=True)
np.trapezoid(prof['heating'], prof['radii'])   # the total

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

The fully collapsed `total` equals the 1D global `get_tidal_heating`; per-layer totals sum to it (they replace the `tidal_scale` distribution for the depth-resolved rheology path); the radial and colatitude profiles integrate to it. The secular colatitude collapse is analytic (the angular Gram table), and the fallback `latitude_nodes` (16) adds margin over the ~4 colatitude nodes a homogeneous degree-2 body needs. The radial collapse keeps its `radial_slices` nodes (16 per layer by default) strictly inside each layer. A node on a layer boundary would read the modulus and radial solution of the layer below, so a layer sitting on a stiffer one would lose part of its own heating. With the nodes inside, the collapsed total matches the 1D heating to better than 1e-4 from 4 nodes per layer, for a homogeneous body and for the bundled Io, whose thin asthenosphere has a sixtieth of the shear modulus of the mantle beneath it. The default adds margin for higher degree l.

#### Displacements: `calc_3d_displacements`

The same machinery gives the instantaneous tidal displacements. For each coherent wave the radial functions $y_1$ (radial displacement) and $y_3$ (tangential displacement) of the wave's radial solution set the complex displacement amplitude at a point, $\mathbf{u} = \left(y_1 U,\; y_3\, \partial U / \partial\theta,\; y_3\, (\partial U / \partial\phi) / \sin\theta\right)$ (TB05 Eq. 9), which is evolved in time as $\mathrm{Re}\left[\mathbf{u}\, e^{i|\omega| t}\right]$ and summed over the waves (the phasor convention of the instantaneous heating). The world method returns the three components on the full `(radius, colatitude, longitude, time)` grid in metres:

```python
out = world.calc_3d_displacements(
    orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass,
    radii=[0.9 * world.radius, world.radius], colatitudes=np.linspace(0.1, np.pi - 0.1, 30),
    longitudes=np.linspace(0.0, 2 * np.pi, 60), times=np.linspace(0.0, period, 12))
out["radial"].shape        # (2, 30, 60, 12)   u_r [m]; also out["polar"], out["azimuthal"]
```

At the surface $y_1 = h / g$ and $y_3 = l / g$, so the surface radial displacement is $h U / g$ for a single mode. It requires the rheology tide model, a solved EOS, and a radial-solver Love-number method (the analytic `homogeneous`/`cpl`/`ctl` methods have no radial functions). A radius without a depth-resolved solution (the center, below the solver start) is NaN. The radial functions themselves are available at any radius through `world.get_love_radial_y(radius, ytype_idx, y_idx)` after a radial-solver Love solve.

#### Stress and Strain: `calc_3d_stress_strain`

The same waves give the instantaneous stress and strain tensors. At each point every wave's complex stress and strain amplitude is added into the total of its frequency, and each component at time $t$ is the sum over frequencies of $\mathrm{Re}\left[A\, e^{i|\omega| t}\right]$ for the complex amplitude $A$, the convention of the displacements and the instantaneous heating. The world method returns both tensors on the full `(radius, colatitude, longitude, time)` grid with the six components on the last axis, ordered `rr`, `theta_theta`, `phi_phi`, `r_theta`, `r_phi`, `theta_phi` (also returned as `components`). The strain is the symmetric gradient of the displacement field of `calc_3d_displacements`.

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

Each tensor takes 48 bytes per grid point and time, and either can be skipped with `return_stress=False` or `return_strain=False`. The C++ code writes directly into the returned arrays. The kernel applies to solid layers only, so a point in a liquid layer, or at a radius without a depth-resolved solution, is NaN. Like the displacement grid, the stress and strain grids carry the time-varying tide only: modes at zero forcing frequency, the permanent tide, are not included. Over a common period of the modes each component therefore averages to zero, even in a strongly dissipative body. Dissipation shows up instead as a lag of the strain behind the stress and as the positive mean power that `calc_3d_tides` returns.

#### Threads

Every grid method, `get_3d_tidal_heating_array`, `calc_3d_tides`, `calc_3d_displacements`, and `calc_3d_stress_strain`, takes `num_threads`, default 1. The radial solves run first on the calling thread. The per-point evaluation after them runs over colatitude rows on up to `num_threads` threads, and rows that add into the same cells, as when colatitude is summed, are combined in row order, so the result is identical for any thread count. The analytic colatitude collapse of `calc_3d_tides`, the default when the secular heating is summed over colatitude, has no per-point grid and always runs on one thread.

The default leaves parallelism to the caller: inside a process pool whose workers already occupy every core, keep `num_threads=1`. For one large grid in an interactive session, pass the number of cores. Notebook 13 times three grids both ways.

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

The potential engine is exposed directly for callers building their own pipelines. It returns each raw mode's degree, signed frequency, and the complex angular-factor amplitudes (one entry per `(l, m, p, q)`, before the coherent merge; a caller summing heating from these must combine the modes at each frequency first):

```python
from TidalPy.constants import G
from TidalPy.Tides_x.potential import tidal_potential_3d_modes

colatitude, longitude = 0.8, 0.3
degrees, freqs, pots = tidal_potential_3d_modes(
    world.radius, orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis,
    host_mass, G, colatitude, longitude,
    min_degree_l=2, max_degree_l=2, eccentricity_truncation=3, obliquity_truncation=0)
# pots[i] = complex (U, dU/dtheta, dU/dphi, d2U/dtheta2, d2U/dphi2, d2U/dtheta_dphi) for mode i
```

The compiled kernel the world methods use is also callable point by point from `Tides_x.multilayer.stress_strain`:

- `strain_stress_heating_point` returns the six complex strain and six complex stress amplitudes at a point for one mode, and that mode's own heating [W m-3] at the frequency given.
- `displacement_point` returns the complex displacement amplitudes $u_r = y_1 U$, $u_\theta = y_3\, \partial U / \partial\theta$, and $u_\phi = y_3\, (\partial U / \partial\phi) / \sin\theta$ [m].
- `volumetric_heating(stress, strain, frequency)` returns the cycle-averaged heating `(|omega| / 2) |sum_k w_k Im(sigma_k conj(eps_k))|` [W m-3] of amplitudes at the frequency `omega`, with `w_k = 2` on the three off-diagonal components, the same factor the world path applies.

The first two take one row from `tidal_potential_3d_modes` together with the radial functions and complex moduli at the point, which the world provides after a radial-solver Love solve. A real row is also accepted and is treated as a phasor with zero phase. Assembling several modes follows the rules the world methods use:

- Solve the radial problem and evaluate the moduli at each mode's degree and `|omega|`, and conjugate the row of a mode with `omega < 0` so that its amplitudes sit at `+|omega|`.
- Sum the strain and stress amplitudes of every mode that shares one `|omega|` before forming the heating. `volumetric_heating` of those sums at `|omega|` is the secular heating at that frequency [W m-3], and the frequencies add.
- A field at time $t$ is $\mathrm{Re}\left[A\, e^{i|\omega| t}\right]$, with $A$ the complex amplitude, summed over the modes.

The pointwise secular density of `calc_3d_tides` can be rebuilt this way:

```python
from TidalPy.constants import min_frequency
from TidalPy.Tides_x.multilayer.stress_strain import strain_stress_heating_point, volumetric_heating

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
