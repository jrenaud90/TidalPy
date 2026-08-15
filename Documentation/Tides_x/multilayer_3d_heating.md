# On-demand 3D tidal stress, strain, and heating (`Tides_x.multilayer`)

This module computes the depth- and direction-resolved tidal response (the complex strain and stress
tensors and the volumetric heating) of a layered world. The system here is **callable**: it returns
the response at a single point on demand, so a map is built only if the caller explicitly evaluates a set of points.

## What it combines

At a point `(r, colatitude θ, longitude φ, time t)` the response factorizes into

* the **radial** part — the viscoelastic-gravitational `y1..y6` from the radial solver, evaluated at `r`
  through the **dense calling system** (`RadialSolverSolution.get_radial_solution`), plus the complex
  (viscoelastic) shear/bulk moduli at `r` read through the dense EOS path
  (`RadialSolverSolution.eos_call_si`); and
* the **angular/time** part — a 2D tidal potential `U(θ, φ, t)` and its first/second θ,φ derivatives.

The strain kernel uses the exact Tobie+2005 forms with the Kervazo+2021 (A&A App. D) correction to the
θφ / φφ components, complex moduli, and a layer-type-dependent `dy1/dr`; it is a solid-layer computation
(liquids contribute no shear dissipation and return NaN). Stress follows the isotropic Takeuchi & Saito
constitutive law. The volumetric heating is
`h = | Σ_k Im(σ_k) Re(ε_k) − Re(σ_k) Im(ε_k) |` with a factor 2 on the three off-diagonal components
(Europa-book Eq. 42). The potential's `r²` coefficient is taken at the **surface** radius; all radial
dependence is carried by `y(r)` and the `1/r` factors in the kernel.

## Tidal potential — dynamic Kaula engine

The tidal potential is built by a **class-free dynamic engine** driven purely by the user's truncation
levels — there is no per-scenario potential-model object. Following Kaula (1964) / Efroimsky &
Williams (2009) Eq. 18, the engine (`c_tidal_potential_3d_modes`, wrapped as
`TidalPy.Tides_x.potential.tidal_potential_3d_modes`) enumerates the active `(l, m, p, q)` modes from
the same eccentricity (`G_lpq`) and inclination/obliquity (`F_lmp`) functions the global (1D) path
uses, times the associated Legendre functions `P_lm` (from `TidalPy.Utilities_x.legendre`). For each
active mode it returns the degree `l`, the signed forcing frequency

```
omega_lmpq = (l - 2p + q) n - m * spin
```

(`n` = orbital mean motion, `spin` = rotation rate; periapse/node precession dropped), and the
potential angular factor `U` with its first/second colatitude/longitude derivatives. The potential is
linear in `F_lmp`, `G_lpq`, and `P_lm` (the global 1D path squares `F`, `G` because global heating goes
as the potential squared). A mode whose `|frequency|` does not exceed `min_spin_orbit_diff` is switched
off downstream.

The user selects the truncation via three knobs (on the world's `[tides]` config): `max_degree_l`
(2..10), `eccentricity_trunc_lvl`, and `obliquity_trunc_lvl` (0 = off). A nonzero obliquity truncation
turns on the odd-`m` (`P_21`, ...) harmonics automatically.

## Secular (cycle/orbit-averaged) heating — the physical quantity

`get_3d_tidal_heating` returns the **secular** (cycle/orbit-averaged) tidal volumetric heating: the
physically time-averaged dissipated-power density. It is built from the mode's **complex** potential
amplitudes (the `e^{i omega t}` pulled out), so the cycle-average is exact with no time grid:

```
h_bar(r, theta) = sum_modes (omega_mode / 2) * Im( sigma_c : conj(eps_c) )
```

Per mode this uses a **single** `omega/2` (the cycle-average factor), the complex stress/strain
amplitudes, and **no** `abs()` — it is inherently non-negative for a dissipative material and modes sum
with sign (distinct-frequency cross terms average to zero over the orbit, so they are omitted). The
`e^{i m phi}` cancels per mode, so `h_bar` is **longitude- and time-independent** (a function of `r` and
`colatitude` only).

**Consistency with the 1D path:** the volume integral of `h_bar` equals the 1D global tidal heating
(`get_tidal_heating`) — both are the same total dissipated power. This is the authoritative correctness
check (validated to ~0.1% for a homogeneous Maxwell sphere; see
`Tests/Test_Structures_x/Test_Worlds/test_world_1d_vs_3d_tides_01.py`).

At each point the path: (1) builds the active modes from the truncation config; (2) solves the world
radial response once per mode `(l, frequency)` (the radial ODEs depend on `l` and `omega` only, not
`m`); (3) accumulates each mode's `(omega/2) Im(sigma_c : conj(eps_c))`.

## Python API

### World method

For a built world this is the canonical entry point. The orchestration lives on the rheology tide model
(`c_RheologyTide`); the world delegates to it and everything runs in C++ (the tide model calls the
world's radial-solver/EOS members directly).

```python
world = build_world(...)  # a LayeredWorld
world.solve_eos()
world.set_tide_model(make_tide("rheology"))
# The tidal potential truncation comes from the tide config (or the [tides] TOML table):
world.set_tide_config(max_degree_l=2, eccentricity_truncation=3, obliquity_truncation=0)
h_bar = world.get_3d_tidal_heating(
  orbital_frequency, spin_frequency, eccentricity, obliquity,
  semi_major_axis, host_mass, radius, colatitude)  # [W m-3], secular
```

Requires the rheology tide model and a solved EOS (the analytic tide models like fixed-Q are rejected
as they have no depth-resolved solution). When a world is built from a TOML/dict config, the
truncation levels flow automatically from the `[tides]` table (`max_degree_l`, `eccentricity_trunc_lvl`,
`obliquity_trunc_lvl`).

#### Building a map: the vectorized batch form

Evaluating a heating map point-by-point through `get_3d_tidal_heating` re-solves the world radial
response at every point, which is wasteful: the radial (Love-number) solve depends only on the tidal
mode's `(degree l, frequency)`, not on the query position. `get_3d_tidal_heating_array` takes paired,
equal-length `(radius, colatitude)` arrays, builds the position-independent tidal-mode list once, and
solves the radial response once per unique `(l, frequency)`, reusing it across all points. It returns a
same-shape array (NaN where a radius has no depth-resolved solution) and is the efficient way to build a
map:

```python
import numpy as np
rr = np.linspace(0.01 * world_radius, 0.999 * world_radius, 40)
th = np.linspace(1e-3, np.pi - 1e-3, 60)
rg, tg = np.meshgrid(rr, th, indexing="ij")
hbar = world.get_3d_tidal_heating_array(
  orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass,
  rg.ravel(), tg.ravel()).reshape(rg.shape)  # [W m-3], secular
```

It reproduces the scalar `get_3d_tidal_heating` point-for-point (to machine precision), so the two are
interchangeable in physics; only the performance differs.

#### Full grid + collapsed flavors: `calc_3d_tides`

`calc_3d_tides` produces the 3D heating as a full grid over `(radius, colatitude, longitude[, time])`
or reduced along any spatial dimension. Two quantities:

* `orbit_averaged=True` (default): the **secular** heating density `h_bar` [W m-3]. It is longitude- and
  time-independent (the per-mode `e^{i m phi}` cancels), so the longitude axis is constant and there is
  no time axis.
* `orbit_averaged=False`: the **instantaneous** mechanical power density `sigma_ij(t) * eps_dot_ij(t)`
  [W m-3] at each supplied time (a 4th axis). It depends on longitude and time and orbit-averages back
  to `h_bar` (modes with signed `omega < 0` use the conjugated phasor at `+|omega|`, reproducing the
  true real field including cross-mode terms).

Non-summed spatial axes take user arrays (`radii`, `colatitudes`, `longitudes`); the `times` array is
required when `orbit_averaged=False`. Reduction convention: if any spatial axis is summed, the surviving
spatial axes carry their Jacobian (`r^2`, `sin theta`, `1`) so a plain integral over them recovers the
total; if none is summed the output is the raw density (it matches the scalar `get_3d_tidal_heating`).
The colatitude integral (secular, `latitude_summed`) is done **analytically** by default: the six angular
functions the strain/stress needs form a bounded basis whose sphere integrals are precomputed once into a
per-`(l, m)` Gram table (`Tides_x.multilayer.stress_strain.angular_gram`), so no colatitude grid or
Legendre evaluation is needed. This is exact and, for a large radius grid (a radial profile or map), a few
times faster than the numerical quadrature (the Love-number solves otherwise dominate). Pass
`latitude_analytic=False` to fall back to a Gauss-Legendre colatitude grid (`latitude_nodes`); the two
agree to machine precision. The radial integral is an internal per-layer trapezoid (`radial_slices`), and
the longitude integral the analytic `2*pi` when averaged or a `longitude_nodes` trapezoid when instantaneous.

The returned dict carries the surviving axes plus either `heating` (the grid over the surviving axes,
in the order radius, colatitude, longitude, time) or, when all three spatial axes are summed, `total`
[W] and `per_layer` [W]:

```python
# default: full secular density grid (nr, ncolat, nlon) — constant along longitude
grid = world.calc_3d_tides(
    orbital_motion,
    spin,
    eccentricity,
    obliquity,
    semi_major_axis,
    host,
    radii=rr,
    colatitudes=cc,
    longitudes=ll)['heating']

# per-layer and total heating (all three spatial axes summed)
res = world.calc_3d_tides(
    orbital_motion,
    spin,
    eccentricity,
    obliquity,
    semi_major_axis,
    host,
    latitude_summed=True,
    longitude_summed=True,
    radial_summed=True)

res['total']      # [W], equals the 1D get_tidal_heating
res['per_layer']  # [W] per layer (innermost first), sums to res['total']

# radial power profile dP/dr (integrates to the total)
prof = world.calc_3d_tides(
    orbital_motion,
    spin,
    eccentricity,
    obliquity,
    semi_major_axis,
    host,
    radii=np.linspace(1e3, R, 400),
    latitude_summed=True,
    longitude_summed=True)
np.trapezoid(prof['heating'], prof['radii'])   # == total

# time-resolved instantaneous power over one orbital period
res = world.calc_3d_tides(
    orbital_motion,
    spin,
    eccentricity,
    obliquity,
    semi_major_axis,
    host,
    radii=rr,
    colatitudes=cc,
    longitudes=ll,
    times=np.linspace(0, 2*np.pi/n, 400),
    orbit_averaged=False)
res['heating']   # shape (nr, ncolat, nlon, ntime)
```

The fully collapsed `total` equals the 1D global `get_tidal_heating`; per-layer totals sum to it (they
replace the `tidal_scale` distribution for the depth-resolved rheology path); the radial and colatitude
profiles integrate to it. The secular colatitude collapse is analytic (the angular Gram table); the
`radial_slices` default (16 per layer) and the fallback `latitude_nodes` (16) were tuned so the collapsed
total is well within 1% of the 1D heating for a homogeneous degree-2 body (~8 radial slices per layer and
~4 colatitude nodes already suffice; the defaults add margin for higher degree l and layered bodies).

### Engine + kernel (raw) access

The dynamic potential engine is exposed directly for callers building their own pipelines. It returns
each mode's degree, signed frequency, and the **complex** angular-factor amplitudes:

```python
from TidalPy.Tides_x.potential import tidal_potential_3d_modes
degrees, freqs, pots = tidal_potential_3d_modes(
  orbital_frequency, spin_frequency, eccentricity, obliquity, host_mass, semi_major_axis,
  planet_radius, colatitude, longitude, G,
  max_degree_l=2, eccentricity_truncation=3, obliquity_truncation=0)
# pots[i] = complex (U, dU/dtheta, dU/dphi, d2U/dtheta2, d2U/dphi2, d2U/dtheta_dphi) for mode i
```

The compiled strain/stress/heating kernel is in `Tides_x.multilayer.stress_strain`
(`strain_stress_heating_point`, `volumetric_heating`) — low-level helpers that take a real potential
row (a snapshot at one time) and return the raw bilinear magnitude; the physical secular heating uses
the complex/signed form above.

## Feature summary

Implemented: the associated-Legendre tables (`Utilities_x.legendre`), the dynamic Kaula potential
engine (any eccentricity/obliquity truncation, `l = 2..10`), the point kernel, and the **world-level,
all-C++** secular `LayeredWorld.get_3d_tidal_heating` (orchestration on `c_RheologyTide`; the world
delegates), whose volume integral matches the 1D global heating to ~0.1%, and the **vectorized C++ batch
path** `get_3d_tidal_heating_array` (position-independent mode list built once, radial solve amortized
across points), and `calc_3d_tides` — the **full `(radius, colatitude, longitude[, time])` grid** with
its **collapsed flavors** (radial/colatitude profiles, per-layer and whole-planet totals; the total
equals the 1D global heating), the **analytic colatitude collapse** (precomputed angular Gram table, exact
and faster than the quadrature for large radius grids), and the **instantaneous `sigma:eps_dot` path**
(`orbit_averaged=False`), which orbit-averages back to the secular density.
