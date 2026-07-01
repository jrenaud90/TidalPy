"""End-to-end on-demand 3D heating vs the legacy collapse_multilayer_modes (homogeneous solid planet).

Validates the full grid-free pipeline (dense radial solver -> EOS moduli -> surface tidal potential ->
strain/stress/heating kernel, with the surface-potential + freq_half normalization) against the legacy
multilayer mode collapse, instantaneous (no orbit averaging). They should agree to the radial solver's
integration accuracy (the only difference is RadialSolver_x vs the original radial solver).
"""
import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.rheology import Maxwell, Elastic
from TidalPy.utilities.conversions import orbital_motion2semi_a
from TidalPy.utilities.spherical_helper.volume import calculate_voxel_volumes
from TidalPy.tides.modes.multilayer_modes import collapse_multilayer_modes
from TidalPy.RadialSolver_x.solver import radial_solver as rs_x
from TidalPy.Tides_x.multilayer.heating_3d import tidal_heating_3d_point


def test_on_demand_heating_matches_legacy_homogeneous_solid():
    Nr = 12
    planet_R = 1.0e6
    radius_array = np.linspace(0.0, planet_R, Nr)
    density_array = 5000.0 * np.ones(Nr)
    shear_array = 5.0e10 * np.ones(Nr)
    visc_array = 1.0e19 * np.ones(Nr)
    bulk_array = 1.0e11 * np.ones(Nr)
    n_orb = 2.0 * np.pi / 86400.0
    ecc = 0.05
    host = mass_trap1
    vol = (4.0 / 3.0) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
    planet_mass = float(np.sum(vol * density_array[1:]))
    bulk_density = planet_mass / np.sum(vol)
    a_sma = orbital_motion2semi_a(n_orb, host, planet_mass)

    cshear = np.empty(Nr, np.complex128); Maxwell().vectorize_modulus_viscosity(n_orb, shear_array, visc_array, cshear)
    cbulk = np.empty(Nr, np.complex128); Elastic().vectorize_modulus_viscosity(n_orb, bulk_array, visc_array, cbulk)

    colat = np.radians(np.linspace(0.1, 179.9, 4))
    lon = np.radians(np.linspace(0.0, 360.0, 5))
    tt = np.linspace(0.0, 2.0 * np.pi / n_orb, 4)
    vox = calculate_voxel_volumes(radius_array, lon, colat)
    LON, COLAT, TIME = np.meshgrid(lon, colat, tt, indexing='ij')

    leg = collapse_multilayer_modes(
        orbital_frequency=n_orb, spin_frequency=n_orb, semi_major_axis=a_sma, eccentricity=ecc,
        host_mass=host, planet_bulk_density=bulk_density, radius_array=radius_array,
        density_array=density_array, bulk_array=bulk_array, shear_array=shear_array,
        bulk_viscosity_array=visc_array, shear_viscosity_array=visc_array,
        shear_rheology_inst_bylayer=(Maxwell(),), bulk_rheology_inst_bylayer=(Elastic(),),
        upper_radius_bylayer_array=np.asarray((planet_R,)), longitude_matrix=LON,
        colatitude_matrix=COLAT, time_matrix=TIME, voxel_volume=vox, layer_types=('solid',),
        is_static_bylayer=(False,), is_incompressible_bylayer=(False,), degree_l=2,
        use_modes=False, use_simple_potential=True, orbit_average_results=False, nondimensionalize=True,
        integration_method='DOP853', integration_rtol=1e-8, integration_atol=1e-10)
    leg_volheat = leg[1]

    rs = rs_x(radius_array.copy(), density_array.copy(), cbulk.copy(), cshear.copy(), n_orb, bulk_density,
              ('solid',), (False,), (False,), np.asarray((planet_R,)), degree_l=2, solve_for=('tidal',),
              nondimensionalize=True, integration_method='DOP853', integration_rtol=1e-8,
              integration_atol=1e-10, max_num_steps=5_000_000, raise_on_fail=True)
    assert rs.success

    worst = 0.0; ncmp = 0
    for ri in range(1, Nr):
        for lj in range(lon.size):
            for ck in range(colat.size):
                for tm in range(tt.size):
                    mine = tidal_heating_3d_point(
                        rs, radius_array[ri], lon[lj], colat[ck], tt[tm], n_orb, ecc, host, a_sma,
                        ('solid',), (False,), (False,), (planet_R,), G, degree_l=2)
                    ref = abs(leg_volheat[ri, lj, ck, tm])
                    if ref > 1e-30 and np.isfinite(mine):
                        worst = max(worst, abs(mine - ref) / ref); ncmp += 1
    assert ncmp > 100
    assert worst < 1e-3, f"on-demand heating differs from legacy: {worst:.3e} over {ncmp} points"
