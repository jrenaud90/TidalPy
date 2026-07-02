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

## Multi-mode combination

Each active mode has its own frequency, so its radial response (complex moduli + `y`-solution) differs.
The on-demand path therefore, at each point:

1. builds the active modes from the truncation config;
2. solves the radial response for each mode (once per `(l, frequency)`; the radial ODEs depend on `l`
   and `omega` only, not `m`);
3. sums every mode's **stress and strain tensors**, each scaled by its own `freq_half = |omega| / 2`,
   then computes the heating **once** from the combined tensors.

Summing the tensors before heating (rather than heating per mode and summing) keeps the cross-mode terms.

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
h = world.get_3d_tidal_heating(
  orbital_frequency, spin_frequency, eccentricity, obliquity,
  semi_major_axis, host_mass, radius, colatitude, longitude, time)  # [W m-3]
```

Requires the rheology tide model and a solved EOS (the analytic tide models like fixed-Q are rejected
as they have no depth-resolved solution). When a world is built from a TOML/dict config, the
truncation levels flow automatically from the `[tides]` table (`max_degree_l`, `eccentricity_trunc_lvl`,
`obliquity_trunc_lvl`).

### Engine + kernel (raw) access

The dynamic potential engine is exposed directly for callers building their own pipelines:

```python
from TidalPy.Tides_x.potential import tidal_potential_3d_modes
degrees, freqs, pots = tidal_potential_3d_modes(
  orbital_frequency, spin_frequency, eccentricity, obliquity, host_mass, semi_major_axis,
  planet_radius, colatitude, longitude, time, G,
  max_degree_l=2, eccentricity_truncation=3, obliquity_truncation=0)
# pots[i] = (U, dU/dtheta, dU/dphi, d2U/dtheta2, d2U/dphi2, d2U/dtheta_dphi) for mode i
```

The compiled strain/stress/heating kernel is in `Tides_x.multilayer.stress_strain`
(`strain_stress_heating_point`, `volumetric_heating`).

## Status / follow-on work

Implemented: the associated-Legendre tables (`Utilities_x.legendre`), the dynamic Kaula potential
engine (any eccentricity/obliquity truncation, `l = 2..10`), the point kernel, and the **world-level,
all-C++** `LayeredWorld.get_3d_tidal_heating` (orchestration on `c_RheologyTide`; the world delegates).
**Under validation (in development):** the new engine's point heating agrees with the legacy
`collapse_multilayer_modes(use_modes=True)` to ~1%, versus <0.5% for the retired hand-coded provider
that shared the legacy convention; the faithful-Kaula normalization and the 1D-vs-3D total-heating
consistency are being validated against a direct/1D-consistent reference. Planned: pre-integrated
angular tables and the analytic dimension collapse (orbit-average and latitude/longitude/radial sums);
and a vectorized C++ batch path for building maps.
