"""On-demand 3D tidal heating for the raw-array path, a convenience over the standalone
radial solver, when you have radius/density/modulus arrays rather than a built world.

For a built world the class-based, all-C++ path is ``LayeredWorld.get_3d_tidal_heating(...)`` (the
orchestration lives on ``c_RheologyTide``; the world delegates). This module is the array-driven
counterpart: given a solved ``RadialSolverSolution`` (the dense radial calling system) and the orbital
state, it returns the tidal volumetric heating [W m-3] at a single point (radius, longitude, colatitude,
time) - or, vectorized, at arrays of points. No multi-dimensional arrays are materialized; a full map is
built only if the caller passes the full set of points.

Two truncations are available: the synchronous / low-eccentricity / no-obliquity potential (one mode,
``tidal_heating_3d_point``) and the moderate-eccentricity, non-synchronous-rotation, no-obliquity
truncation (up to nine modes, ``prepare_nsr_modes`` + ``tidal_heating_3d_point_nsr``). The multi-mode
path solves the radial response once per active mode frequency, then sums each mode's stress/strain
tensors (each scaled by its own ``freq_half``) and heats once from the combined tensors (cross-mode
terms preserved). Obliquity modes and the analytic dimension-collapse paths build on this same kernel.

Normalization:
  * the tidal potential is evaluated at the surface radius (its r^2 factor uses the planet radius);
    the radial dependence is carried entirely by the radial-solver y-functions and the 1/r factors;
  * each mode's strain/stress is scaled by ``freq_half = mode_frequency / 2`` before the (bilinear)
    heating - for the single synchronous mode this is a ``(omega/2)**2`` factor on the heating.
"""
import numpy as np

from TidalPy.constants import min_spin_orbit_diff as _MIN_SPIN_ORBIT_DIFF
from TidalPy.RadialSolver_x.solver import radial_solver as _radial_solver
from TidalPy.Tides_x.multilayer.stress_strain import strain_stress_heating_point, volumetric_heating
from TidalPy.Tides_x.potential.tidal_potential import SyncLowEPotential, NSRModesPotential

# Stateless single-mode potential model, reused across point evaluations.
_SYNC_POTENTIAL = SyncLowEPotential()


def _layer_index_at(radius, upper_radius_bylayer):
    for i, upper_radius in enumerate(upper_radius_bylayer):
        if radius <= upper_radius:
            return i
    return len(upper_radius_bylayer) - 1


def _moduli_at(rs_solution, radius):
    """Complex (shear, bulk) [Pa] at a radius via the solution's own dense EOS interpolant.

    Routes through the radial solver's dense EOS output (``eos_call_si``) rather than re-interpolating
    the gridded modulus arrays, so an on-radius query gets the same dense evaluation the solver uses
    internally (accurate between grid slices, not just at them). The dense outputs are laid out
    ``[..., 5] shear_re, [6] shear_im, [7] bulk_re, [8] bulk_im`` (SI).
    """
    eos = rs_solution.eos_call_si(radius)
    shear = complex(eos[5], eos[6])
    bulk = complex(eos[7], eos[8])
    return shear, bulk


def tidal_heating_3d_point(
        rs_solution,
        radius,
        longitude,
        colatitude,
        time,
        orbital_frequency,
        eccentricity,
        host_mass,
        semi_major_axis,
        layer_types,
        is_static_bylayer,
        is_incompressible_bylayer,
        upper_radius_bylayer,
        G,
        degree_l=2,
        ytype_index=0):
    """Tidal volumetric heating [W m-3] at one point. Returns NaN in liquid layers / at the center."""
    y = rs_solution.get_radial_solution(radius, ytype_index)
    if not np.all(np.isfinite(np.asarray(y)[:4])):
        return np.nan

    shear, bulk   = _moduli_at(rs_solution, radius)
    layer_index   = _layer_index_at(radius, upper_radius_bylayer)
    is_solid      = (str(layer_types[layer_index]).lower() == 'solid')
    is_incomp     = bool(is_incompressible_bylayer[layer_index])
    planet_radius = float(rs_solution.radius_array[-1])

    _, potentials = _SYNC_POTENTIAL.calc_modes(
        orbital_frequency, 0.0, eccentricity, host_mass, semi_major_axis,
        planet_radius, colatitude, longitude, time)
    pot6 = tuple(potentials[0])
    _, _, h = strain_stress_heating_point(
        np.ascontiguousarray(y, dtype=np.complex128),
        shear,
        bulk,
        radius,
        float(degree_l),
        is_solid,
        is_incomp,
        pot6,
        colatitude)

    freq_half = orbital_frequency / 2.0
    return (freq_half * freq_half) * h


def _layer_slices(radius_array, upper_radius_bylayer):
    """Index ranges (start, end) of each layer in radius_array, matching the radial solver's split.

    A layer owns the slices with radius <= its upper radius that are not already claimed by a lower
    layer; a duplicated interface radius (the per-layer convention) is split between the two layers.
    """
    slices = []
    start = 0
    n = radius_array.size
    for upper in upper_radius_bylayer:
        end = start
        seen_interface = 0
        while end < n:
            r = radius_array[end]
            if abs(r - upper) <= 1.0e-9 * max(abs(upper), 1.0):
                seen_interface += 1
                if seen_interface > 1:
                    break
            elif r > upper:
                break
            end += 1
        slices.append((start, end))
        start = end
    return slices


def _complex_moduli_arrays(
        rheology_bylayer,
        modulus_array,
        viscosity_array,
        layer_slices,
        frequency):
    """Complex (viscoelastic) modulus array at `frequency`, built per layer from its rheology model."""
    out = np.empty(modulus_array.size, dtype=np.complex128)
    for layer_i, (start, end) in enumerate(layer_slices):
        segment = np.empty(end - start, dtype=np.complex128)
        rheology_bylayer[layer_i].vectorize_modulus_viscosity(
            frequency, modulus_array[start:end], viscosity_array[start:end], segment)
        out[start:end] = segment
    return out


def prepare_nsr_modes(
        orbital_frequency,
        spin_frequency,
        eccentricity,
        host_mass,
        semi_major_axis,
        radius_array,
        density_array,
        shear_array,
        bulk_array,
        shear_viscosity_array,
        bulk_viscosity_array,
        shear_rheology_bylayer,
        bulk_rheology_bylayer,
        upper_radius_bylayer,
        planet_bulk_density,
        layer_types,
        is_static_bylayer,
        is_incompressible_bylayer,
        degree_l=2,
        use_static=False,
        **radial_solver_kwargs):
    """Solve the radial response once per active NSR tidal mode, ready for on-demand 3D heating.

    The moderate-eccentricity, non-synchronous-rotation (no obliquity) truncation has up to nine modes
    (n, 2n, 3n, 2o+n, 2o-n, 2o-2n, 2o-3n, 2o-4n, 2o-5n). For each unique mode frequency above the
    spin-orbit threshold this computes the complex moduli at that frequency and runs the (dense) radial
    solver, so the per-mode y-solution and moduli are reused across every spatial/temporal query.

    Returns a dict consumed by `tidal_heating_3d_point_nsr` / `_grid_nsr`.
    """
    layer_slices = _layer_slices(radius_array, upper_radius_bylayer)
    potential = NSRModesPotential(use_static=use_static)

    # Mode frequencies are point-independent, so evaluate the model at a throwaway point to read them.
    signed_modes, _ = potential.calc_modes(
        orbital_frequency, spin_frequency, eccentricity, host_mass, semi_major_axis,
        radius_array[-1], 1.0, 0.0, 0.0)

    rs_kwargs = dict(
        degree_l=degree_l, solve_for=('tidal',), nondimensionalize=True,
        integration_method='DOP853', integration_rtol=1e-8, integration_atol=1e-10,
        max_num_steps=5_000_000, raise_on_fail=True)
    rs_kwargs.update(radial_solver_kwargs)

    # One radial solve per unique active frequency (modes that share |frequency| share the solution).
    rs_by_freq_key = {}
    active = []
    for mode_i in range(signed_modes.size):
        signed_mode = float(signed_modes[mode_i])
        abs_freq = abs(signed_mode)
        if abs_freq <= _MIN_SPIN_ORBIT_DIFF:
            continue   # switched-off mode (matches the potential's mode switch); no tidal response
        freq_key = round(abs_freq, 12)
        if freq_key not in rs_by_freq_key:
            cshear = _complex_moduli_arrays(
                shear_rheology_bylayer, shear_array, shear_viscosity_array, layer_slices, abs_freq)
            cbulk = _complex_moduli_arrays(
                bulk_rheology_bylayer, bulk_array, bulk_viscosity_array, layer_slices, abs_freq)
            rs = _radial_solver(
                radius_array.copy(), density_array.copy(), cbulk, cshear, abs_freq, planet_bulk_density,
                tuple(layer_types), tuple(is_static_bylayer), tuple(is_incompressible_bylayer),
                np.asarray(upper_radius_bylayer), **rs_kwargs)
            if not rs.success:
                raise RuntimeError(f"RadialSolver failed at mode frequency {abs_freq}: {rs.message}")
            rs_by_freq_key[freq_key] = rs
        active.append((mode_i, signed_mode, abs_freq, freq_key))

    return dict(
        active_modes=active, rs_by_freq_key=rs_by_freq_key, planet_radius=float(radius_array[-1]),
        potential=potential,
        orbital_frequency=orbital_frequency, spin_frequency=spin_frequency, eccentricity=eccentricity,
        host_mass=host_mass, semi_major_axis=semi_major_axis, use_static=use_static,
        layer_types=tuple(layer_types), is_incompressible_bylayer=tuple(is_incompressible_bylayer),
        upper_radius_bylayer=tuple(upper_radius_bylayer), degree_l=degree_l)


def tidal_heating_3d_point_nsr(prepared, radius, longitude, colatitude, time, ytype_index=0):
    """Tidal volumetric heating [W m-3] at one point for the NSR multi-mode truncation.

    Sums each active mode's stress and strain tensors (each scaled by its own ``freq_half = |omega|/2``)
    and computes the heating once from the combined tensors (keeping cross-mode terms, matching the
    legacy collapse). Returns NaN in liquid layers / at the center.
    """
    layer_index = _layer_index_at(radius, prepared['upper_radius_bylayer'])
    is_solid = (str(prepared['layer_types'][layer_index]).lower() == 'solid')
    is_incomp = bool(prepared['is_incompressible_bylayer'][layer_index])
    degree_l = prepared['degree_l']

    # Potential of every mode at this point (the model returns all nine; we use the active ones).
    # The potential's r^2 coefficient uses the SURFACE radius - the radial dependence is carried entirely
    # by the y-functions and the 1/r factors in the kernel (matches the single-mode path and the legacy).
    _, potentials = prepared['potential'].calc_modes(
        prepared['orbital_frequency'], prepared['spin_frequency'], prepared['eccentricity'],
        prepared['host_mass'], prepared['semi_major_axis'],
        prepared['planet_radius'], colatitude, longitude, time)

    strain_total = np.zeros(6, dtype=np.complex128)
    stress_total = np.zeros(6, dtype=np.complex128)
    any_mode = False
    for mode_i, _signed_mode, abs_freq, freq_key in prepared['active_modes']:
        rs = prepared['rs_by_freq_key'][freq_key]
        y = rs.get_radial_solution(radius, ytype_index)
        if not np.all(np.isfinite(np.asarray(y)[:4])):
            return np.nan
        shear, bulk = _moduli_at(rs, radius)
        pot6 = tuple(potentials[mode_i])
        strain, stress, _ = strain_stress_heating_point(
            np.ascontiguousarray(y, dtype=np.complex128), shear, bulk, radius, float(degree_l),
            is_solid, is_incomp, pot6, colatitude)
        freq_half = abs_freq / 2.0
        strain_total += freq_half * strain
        stress_total += freq_half * stress
        any_mode = True

    if not any_mode:
        return 0.0
    return volumetric_heating(
        np.ascontiguousarray(stress_total, dtype=np.complex128),
        np.ascontiguousarray(strain_total, dtype=np.complex128))


def tidal_heating_3d_grid(
        rs_solution,
        radius_array,
        longitude_array,
        colatitude_array,
        time_array,
        orbital_frequency,
        eccentricity,
        host_mass,
        semi_major_axis,
        layer_types,
        is_static_bylayer,
        is_incompressible_bylayer,
        upper_radius_bylayer,
        G,
        degree_l=2,
        ytype_index=0):
    """Build a (radius, longitude, colatitude, time) volumetric-heating map by point evaluation.

    This is the opt-in materialized map: the caller chooses the resolution and pays the memory. For a
    single point, prefer tidal_heating_3d_point. (Vectorized C++ + analytic collapse are follow-ons.)
    """
    out = np.empty((radius_array.size, longitude_array.size, colatitude_array.size, time_array.size),
                   dtype=np.float64)
    for ri, r in enumerate(radius_array):
        for lj, lon in enumerate(longitude_array):
            for ck, colat in enumerate(colatitude_array):
                for tm, t in enumerate(time_array):
                    out[ri, lj, ck, tm] = tidal_heating_3d_point(
                        rs_solution, r, lon, colat, t, orbital_frequency, eccentricity, host_mass,
                        semi_major_axis, layer_types, is_static_bylayer, is_incompressible_bylayer,
                        upper_radius_bylayer, G, degree_l, ytype_index)
    return out
