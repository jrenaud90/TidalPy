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

## Tidal Potential Truncations

The tidal potential is a class hierarchy `c_TidalPotentialBase : c_PhysicsBase`
(`TidalPy.Tides_x.potential.tidal_potential`), built by `make_tidal_potential(name, config)` and
following the same model pattern as the rheology/tide hierarchies (enum factory, binary serialization,
config dict). Each model's `calc_modes(...)` returns every active mode's signed frequency and the
potential angular factor `U` with its colatitude/longitude derivatives.

| Model (alias) | Modes | Truncation |
|---------------|-------|------------|
| `SyncLowEPotential` (`sync_low_e`) | 1 (`n`) | synchronous rotation, low eccentricity (`e`), no obliquity |
| `NSRModesPotential` (`nsr_modes`) | up to 9 (`n, 2n, 3n, 2o±kn`) | moderate eccentricity (`e³`), non-synchronous rotation, no obliquity |
| `NSRMedObliquityPotential` (`nsr_modes_med_obliquity`) | up to 17 (the 9 no-obliquity modes plus the `m = 1` modes `o, o±kn` and `2o`) | moderate eccentricity (`e³`), moderate obliquity (`obliquity³`), non-synchronous rotation |

`o` is the rotation (spin) frequency and `n` the orbital frequency. A mode whose `|frequency|` does not
exceed `min_spin_orbit_diff` is switched off. The obliquity truncation adds the `m = 1` associated
Legendre harmonic (`P_21`); pass the planet `obliquity` (radians) to `calc_modes` when using it.

## Multi-mode combination

Each active mode has its own frequency, so its radial response (complex moduli + `y`-solution) differs.
The on-demand multi-mode path therefore:

1. solves the radial response **once per unique active mode frequency** (`prepare_nsr_modes`), caching it;
2. at each point, sums every mode's **stress and strain tensors**, each scaled by its own
   `freq_half = |ω_mode| / 2`, and then computes the heating **once** from the combined tensors
   (`tidal_heating_3d_point_nsr`).

Summing the tensors before heating (rather than heating per mode and summing) keeps the cross-mode terms,
matching the legacy collapse. This is validated to ~1e-12 (relative) against
`collapse_multilayer_modes(use_modes=True)` for a non-synchronously rotating homogeneous body.

## Python API

### World method

For a built world this is the canonical entry point. The orchestration lives on the rheology tide model
(`c_RheologyTide`); the world delegates to it and everything runs in C++ (the tide model calls the
world's radial-solver/EOS members directly).

```python
world = build_world(...)  # a LayeredWorld
world.solve_eos()
world.set_tide_model(make_tide("rheology"))
world.set_tidal_potential_model(make_tidal_potential("nsr_modes"))
h = world.get_3d_tidal_heating(
  orbital_frequency, spin_frequency, eccentricity, obliquity,
  semi_major_axis, host_mass, radius, colatitude, longitude, time)  # [W m-3]
```

Requires the rheology tide model, a tidal potential model, and a solved EOS (the analytic tide models like fixed-Q
are rejected as they have no depth-resolved solution).

### Standalone (raw-array) helpers

When you have radius/density/modulus arrays rather than a world, the standalone helpers can be
used to solve for the 3D tidal heating directly:

```python
from TidalPy.RadialSolver_x.solver import radial_solver
from TidalPy.Tides_x.multilayer.heating_3d import (
    tidal_heating_3d_point, prepare_nsr_modes, tidal_heating_3d_point_nsr)

# Single synchronous mode: solve the radial problem once, then query points.
rs = radial_solver(...)
h = tidal_heating_3d_point(
  rs, radius, longitude, colatitude, time,
  orbital_frequency, eccentricity, host_mass, semi_major_axis,
  layer_types, is_static_bylayer, is_incompressible_bylayer,
  upper_radius_bylayer, G)   # [W m-3]

# NSR multi-mode: prepare once (one radial solve per active mode frequency), then query points.
prepared = prepare_nsr_modes(
  orbital_frequency, spin_frequency, eccentricity, host_mass,
  semi_major_axis, radius_array, density_array, shear_array, bulk_array,
  shear_viscosity_array, bulk_viscosity_array,
  shear_rheology_bylayer, bulk_rheology_bylayer,
  upper_radius_bylayer, planet_bulk_density,
  layer_types, is_static_bylayer, is_incompressible_bylayer)
h = tidal_heating_3d_point_nsr(prepared, radius, longitude, colatitude, time)   # [W m-3]
```

The tidal-potential models live in the `Tides_x.potential.tidal_potential` extension
(`make_tidal_potential`, `SyncLowEPotential`, `NSRModesPotential`); the compiled strain/stress/heating
kernel is in `Tides_x.multilayer.stress_strain` (`strain_stress_heating_point`, `volumetric_heating`).

## Status / follow-on work

Implemented: the `c_TidalPotentialBase` model hierarchy (both truncations), the point kernel, the
on-demand single- and multi-mode heating, and the **world-level, all-C++** `LayeredWorld.get_3d_tidal_heating`
(orchestration on `c_RheologyTide`; the world delegates). Planned: obliquity-mode potentials;
pre-integrated angular tables and the analytic dimension collapse (orbit-average and
latitude/longitude/radial sums); and a vectorized C++ batch path for building maps.
