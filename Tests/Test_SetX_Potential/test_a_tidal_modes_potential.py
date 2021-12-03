""" Tests for extracting useful information out of a multilayer tidal propagation
"""

import numpy as np

import TidalPy
from TidalPy.constants import G
from TidalPy.tides.potential import tidal_potential_nsr_modes
from TidalPy.toolbox.conversions import orbital_motion2semi_a

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()

# Model planet - 2layers
density_array = 5000. * np.ones(7)
radius_array = np.linspace(0., 1.e6, 8)
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(7)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.complex)
host_mass = 10. * planet_mass
orbital_freq = (2. * np.pi / (86400. * 6.))
semi_major_axis = orbital_motion2semi_a(orbital_freq, host_mass, planet_mass)
eccentricity = 0.01

rad_mtx, time_mtx, long_mtx, lat_mtx = \
    np.meshgrid(radius_array, np.linspace(1000., 2000., 3),
                np.linspace(0., 2., 3), np.linspace(0., 1., 3), indexing='ij')


def test_tidal_potential_nsr_modes():
    """ Test the basic tidal potential assuming moderate eccentricity, no obliquity, and synchronous rotation """
    # Test arrays
    modes, potential_dict, potential_dtheta_dict, potential_dphi_dict, potential_d2theta_dict, \
        potential_d2phi_dict, potential_dtheta_dphi_dict = \
        tidal_potential_nsr_modes(
            rad_mtx, long_mtx, lat_mtx, orbital_freq, eccentricity, time_mtx,
            orbital_freq * 2., world_radius=radius_array[-1], host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False)

    # Check mode frequency types
    assert len(modes) == 9
    for mode_label, mode_freq in modes.items():
        assert type(mode_label) == str
        assert type(mode_freq) in [float, np.float, np.float64]

        assert type(potential_dict[mode_label]) == np.ndarray
        assert type(potential_dtheta_dict[mode_label]) == np.ndarray
        assert type(potential_dphi_dict[mode_label]) == np.ndarray
        assert type(potential_d2theta_dict[mode_label]) == np.ndarray
        assert type(potential_d2phi_dict[mode_label]) == np.ndarray
        assert type(potential_dtheta_dphi_dict[mode_label]) == np.ndarray

        assert potential_dict[mode_label].shape == rad_mtx.shape
        assert potential_dtheta_dict[mode_label].shape == rad_mtx.shape
        assert potential_dphi_dict[mode_label].shape == rad_mtx.shape
        assert potential_d2theta_dict[mode_label].shape == rad_mtx.shape
        assert potential_d2phi_dict[mode_label].shape == rad_mtx.shape
        assert potential_dtheta_dphi_dict[mode_label].shape == rad_mtx.shape

        assert type(potential_dict[mode_label][0, 0, 0, 0]) in [float, np.float, np.float64]
        assert type(potential_dtheta_dict[mode_label][0, 0, 0, 0]) in [float, np.float, np.float64]
        assert type(potential_dphi_dict[mode_label][0, 0, 0, 0]) in [float, np.float, np.float64]
        assert type(potential_d2theta_dict[mode_label][0, 0, 0, 0]) in [float, np.float, np.float64]
        assert type(potential_d2phi_dict[mode_label][0, 0, 0, 0]) in [float, np.float, np.float64]
        assert type(potential_dtheta_dphi_dict[mode_label][0, 0, 0, 0]) in [float, np.float, np.float64]
