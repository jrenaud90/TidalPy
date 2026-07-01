"""On-demand NSR multi-mode 3D heating vs the legacy collapse_multilayer_modes (homogeneous solid planet).

Validates the multi-mode path (moderate eccentricity, non-synchronous rotation, no obliquity) of the
grid-free pipeline against the legacy mode collapse, instantaneous (no orbit averaging). Each active
tidal mode contributes its own radial solve + complex moduli; the stress/strain tensors are summed
(each scaled by its mode's freq_half) and the heating is taken once from the combined tensors. They
should agree to the radial solver's integration accuracy.

A non-synchronous spin (spin != orbital) is used so several modes have distinct, non-zero frequencies.
"""
import numpy as np
import pytest

from TidalPy.constants import mass_trap1
from TidalPy.rheology import Maxwell, Elastic
from TidalPy.utilities.conversions import orbital_motion2semi_a
from TidalPy.utilities.spherical_helper.volume import calculate_voxel_volumes
from TidalPy.tides.modes.multilayer_modes import collapse_multilayer_modes
from TidalPy.Tides_x.potential.tidal_potential import NSRModesPotential
from TidalPy.Tides_x.multilayer.heating_3d import prepare_nsr_modes, tidal_heating_3d_point_nsr


def _homogeneous_setup(spin_factor):
    n_slices = 12
    planet_radius = 1.0e6
    radius_array = np.linspace(0.0, planet_radius, n_slices)
    density_array = 5000.0 * np.ones(n_slices)
    shear_array = 5.0e10 * np.ones(n_slices)
    viscosity_array = 1.0e19 * np.ones(n_slices)
    bulk_array = 1.0e11 * np.ones(n_slices)
    orbital_frequency = 2.0 * np.pi / 86400.0
    spin_frequency = spin_factor * orbital_frequency
    eccentricity = 0.05
    host = mass_trap1
    shell_volume = (4.0 / 3.0) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
    planet_mass = float(np.sum(shell_volume * density_array[1:]))
    bulk_density = planet_mass / np.sum(shell_volume)
    semi_major_axis = orbital_motion2semi_a(orbital_frequency, host, planet_mass)
    return dict(
        radius_array=radius_array, density_array=density_array, shear_array=shear_array,
        viscosity_array=viscosity_array, bulk_array=bulk_array, orbital_frequency=orbital_frequency,
        spin_frequency=spin_frequency, eccentricity=eccentricity, host=host, planet_radius=planet_radius,
        bulk_density=bulk_density, semi_major_axis=semi_major_axis)


def test_potential_nsr_modes_frequencies():
    """The model returns the nine canonical mode frequencies (signed) n, 2n, 3n, 2o+-kn."""
    n = 2.0 * np.pi / 86400.0
    o = 1.5 * n
    nsr_potential = NSRModesPotential()
    assert nsr_potential.num_modes == 9
    modes, potentials = nsr_potential.calc_modes(n, o, 0.05, 1.0e27, 1.0e9, 1.0e6, 1.2, 0.3, 100.0)
    expected = np.array([n, 2 * n, 3 * n, 2 * o + n, 2 * o - n, 2 * o - 2 * n,
                         2 * o - 3 * n, 2 * o - 4 * n, 2 * o - 5 * n])
    assert np.allclose(modes, expected, rtol=1e-12)
    assert potentials.shape == (9, 6)


@pytest.mark.parametrize('spin_factor', (1.5, 2.0))
def test_nsr_multimode_heating_matches_legacy(spin_factor):
    cfg = _homogeneous_setup(spin_factor)
    Maxwell_rh = (Maxwell(),)
    Elastic_rh = (Elastic(),)
    planet_radius = cfg['planet_radius']

    colat = np.radians(np.linspace(0.1, 179.9, 4))
    lon = np.radians(np.linspace(0.0, 360.0, 5))
    tt = np.linspace(0.0, 2.0 * np.pi / cfg['orbital_frequency'], 4)
    vox = calculate_voxel_volumes(cfg['radius_array'], lon, colat)
    LON, COLAT, TIME = np.meshgrid(lon, colat, tt, indexing='ij')

    leg = collapse_multilayer_modes(
        orbital_frequency=cfg['orbital_frequency'], spin_frequency=cfg['spin_frequency'],
        semi_major_axis=cfg['semi_major_axis'], eccentricity=cfg['eccentricity'], host_mass=cfg['host'],
        planet_bulk_density=cfg['bulk_density'], radius_array=cfg['radius_array'],
        density_array=cfg['density_array'], bulk_array=cfg['bulk_array'], shear_array=cfg['shear_array'],
        bulk_viscosity_array=cfg['viscosity_array'], shear_viscosity_array=cfg['viscosity_array'],
        shear_rheology_inst_bylayer=Maxwell_rh, bulk_rheology_inst_bylayer=Elastic_rh,
        upper_radius_bylayer_array=np.asarray((planet_radius,)), longitude_matrix=LON,
        colatitude_matrix=COLAT, time_matrix=TIME, voxel_volume=vox, layer_types=('solid',),
        is_static_bylayer=(False,), is_incompressible_bylayer=(False,), degree_l=2,
        use_modes=True, use_simple_potential=False, orbit_average_results=False, nondimensionalize=True,
        integration_method='DOP853', integration_rtol=1e-8, integration_atol=1e-10)
    leg_volheat = leg[1]

    prepared = prepare_nsr_modes(
        cfg['orbital_frequency'], cfg['spin_frequency'], cfg['eccentricity'], cfg['host'],
        cfg['semi_major_axis'], cfg['radius_array'], cfg['density_array'], cfg['shear_array'],
        cfg['bulk_array'], cfg['viscosity_array'], cfg['viscosity_array'], Maxwell_rh, Elastic_rh,
        np.asarray((planet_radius,)), cfg['bulk_density'], ('solid',), (False,), (False,),
        degree_l=2, integration_rtol=1e-8, integration_atol=1e-10)

    # Several distinct, active mode frequencies (not a single-mode degenerate case).
    assert len(prepared['rs_by_freq_key']) >= 3

    worst = 0.0
    ncmp = 0
    for ri in range(1, cfg['radius_array'].size):
        for lj in range(lon.size):
            for ck in range(colat.size):
                for tm in range(tt.size):
                    mine = tidal_heating_3d_point_nsr(
                        prepared, cfg['radius_array'][ri], lon[lj], colat[ck], tt[tm])
                    ref = abs(leg_volheat[ri, lj, ck, tm])
                    if ref > 1e-20 and np.isfinite(mine):
                        worst = max(worst, abs(mine - ref) / ref)
                        ncmp += 1
    assert ncmp > 100
    assert worst < 1e-6, f"NSR multi-mode heating differs from legacy: {worst:.3e} over {ncmp} points"
