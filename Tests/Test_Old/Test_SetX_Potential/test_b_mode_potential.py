""" Tests for extracting useful information out of a multilayer tidal propagation
"""
import numpy as np

import TidalPy


from TidalPy.constants import G
from TidalPy.tides.potential import (tidal_potential_nsr, tidal_potential_nsr_modes, tidal_potential_obliquity_nsr,
                                     tidal_potential_gen_obliquity_low_e_nsr_modes,
                                     tidal_potential_obliquity_nsr_modes, tidal_potential_gen_obliquity_nsr_modes)
from TidalPy.utilities.conversions import orbital_motion2semi_a

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array = np.linspace(0., 1.e6, 11)
longitude_array = np.radians(np.linspace(0., 360., 12))
colat_array = np.radians(np.linspace(0.5, 179.5, 13))

volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.complex128)
host_mass = 50000. * planet_mass
orbital_freq = (2. * np.pi / (86400. * 6.))
semi_major_axis = orbital_motion2semi_a(orbital_freq, host_mass, planet_mass)
eccentricity = 0.3
obliquity = np.radians(15.)

time_array = np.linspace(0., 2. * np.pi / orbital_freq, 5)

long_mtx, colat_mtx, time_mtx = np.meshgrid(longitude_array, colat_array, time_array)


def test_tidal_potential_nsr_modes():
    """ Test the modal tidal potential assuming moderate eccentricity, no obliquity, and synchronous rotation """
    # Test arrays - Static=False
    spin_freq = 1.5 * orbital_freq
    freqs, modes, potential_tuple_by_mode = \
        tidal_potential_nsr_modes(
            radius_array[-1], long_mtx, colat_mtx, time_mtx,
            orbital_freq, spin_freq,
            eccentricity,
            host_mass, semi_major_axis,
            use_static=False
            )

    # Check mode frequency types
    assert len(modes) == 9
    assert len(freqs) == 9
    for mode_label, mode_freq in modes.items():
        assert type(mode_label) == str
        assert type(mode_freq) in [float, np.float64]

        # Unpack potential tuple
        potential, potential_dtheta, potential_dphi, potential_d2theta, potential_d2phi, potential_dtheta_dphi = \
            potential_tuple_by_mode[mode_label]

        assert type(potential) == np.ndarray
        assert type(potential_dtheta) == np.ndarray
        assert type(potential_dphi) == np.ndarray
        assert type(potential_d2theta) == np.ndarray
        assert type(potential_d2phi) == np.ndarray
        assert type(potential_dtheta_dphi) == np.ndarray

        assert potential.shape == long_mtx.shape
        assert potential_dtheta.shape == long_mtx.shape
        assert potential_dphi.shape == long_mtx.shape
        assert potential_d2theta.shape == long_mtx.shape
        assert potential_d2phi.shape == long_mtx.shape
        assert potential_dtheta_dphi.shape == long_mtx.shape

        assert type(potential[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta[0, 0, 0]) in [float, np.float64]
        assert type(potential_dphi[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2theta[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2phi[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta_dphi[0, 0, 0]) in [float, np.float64]

    # Test arrays - Static=True
    spin_freq = 1.5 * orbital_freq
    freqs, modes, potential_tuple_by_mode = \
        tidal_potential_nsr_modes(
            radius_array[-1], long_mtx, colat_mtx, time_mtx,
            orbital_freq, spin_freq,
            eccentricity,
            host_mass, semi_major_axis,
            use_static=True
            )

    # Check mode frequency types
    assert len(modes) == 9
    assert len(freqs) == 9
    for mode_label, mode_freq in modes.items():
        assert type(mode_label) == str
        assert type(mode_freq) in [float, np.float64]

        # Unpack potential tuple
        potential, potential_dtheta, potential_dphi, potential_d2theta, potential_d2phi, potential_dtheta_dphi = \
            potential_tuple_by_mode[mode_label]

        assert type(potential) == np.ndarray
        assert type(potential_dtheta) == np.ndarray
        assert type(potential_dphi) == np.ndarray
        assert type(potential_d2theta) == np.ndarray
        assert type(potential_d2phi) == np.ndarray
        assert type(potential_dtheta_dphi) == np.ndarray

        assert potential.shape == long_mtx.shape
        assert potential_dtheta.shape == long_mtx.shape
        assert potential_dphi.shape == long_mtx.shape
        assert potential_d2theta.shape == long_mtx.shape
        assert potential_d2phi.shape == long_mtx.shape
        assert potential_dtheta_dphi.shape == long_mtx.shape

        assert type(potential[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta[0, 0, 0]) in [float, np.float64]
        assert type(potential_dphi[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2theta[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2phi[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta_dphi[0, 0, 0]) in [float, np.float64]


def test_tidal_potential_obliquity_nsr_modes():
    """ Test the modal tidal potential assuming moderate eccentricity, moderate obliquity, and synchronous rotation """
    # Test arrays - Static=False
    spin_freq = 1.5 * orbital_freq
    freqs, modes, potential_tuple_by_mode = \
        tidal_potential_obliquity_nsr_modes(
            radius_array[-1], long_mtx, colat_mtx, time_mtx,
            orbital_freq, spin_freq,
            eccentricity, obliquity,
            host_mass, semi_major_axis,
            use_static=False
            )

    # Check mode frequency types
    assert len(modes) == 17
    assert len(freqs) == 17
    for mode_label, mode_freq in modes.items():
        assert type(mode_label) == str
        assert type(mode_freq) in [float, np.float64]

        # Unpack potential tuple
        potential, potential_dtheta, potential_dphi, potential_d2theta, potential_d2phi, potential_dtheta_dphi = \
            potential_tuple_by_mode[mode_label]

        assert type(potential) == np.ndarray
        assert type(potential_dtheta) == np.ndarray
        assert type(potential_dphi) == np.ndarray
        assert type(potential_d2theta) == np.ndarray
        assert type(potential_d2phi) == np.ndarray
        assert type(potential_dtheta_dphi) == np.ndarray

        assert potential.shape == long_mtx.shape
        assert potential_dtheta.shape == long_mtx.shape
        assert potential_dphi.shape == long_mtx.shape
        assert potential_d2theta.shape == long_mtx.shape
        assert potential_d2phi.shape == long_mtx.shape
        assert potential_dtheta_dphi.shape == long_mtx.shape

        assert type(potential[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta[0, 0, 0]) in [float, np.float64]
        assert type(potential_dphi[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2theta[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2phi[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta_dphi[0, 0, 0]) in [float, np.float64]

    # Test arrays - Static=True
    spin_freq = 1.5 * orbital_freq
    freqs, modes, potential_tuple_by_mode = \
        tidal_potential_obliquity_nsr_modes(
            radius_array[-1], long_mtx, colat_mtx, time_mtx,
            orbital_freq, spin_freq,
            eccentricity, obliquity,
            host_mass, semi_major_axis,
            use_static=True
            )

    # Check mode frequency types
    assert len(modes) == 17
    assert len(freqs) == 17
    for mode_label, mode_freq in modes.items():
        assert type(mode_label) == str
        assert type(mode_freq) in [float, np.float64]

        # Unpack potential tuple
        potential, potential_dtheta, potential_dphi, potential_d2theta, potential_d2phi, potential_dtheta_dphi = \
            potential_tuple_by_mode[mode_label]

        assert type(potential) == np.ndarray
        assert type(potential_dtheta) == np.ndarray
        assert type(potential_dphi) == np.ndarray
        assert type(potential_d2theta) == np.ndarray
        assert type(potential_d2phi) == np.ndarray
        assert type(potential_dtheta_dphi) == np.ndarray

        assert potential.shape == long_mtx.shape
        assert potential_dtheta.shape == long_mtx.shape
        assert potential_dphi.shape == long_mtx.shape
        assert potential_d2theta.shape == long_mtx.shape
        assert potential_d2phi.shape == long_mtx.shape
        assert potential_dtheta_dphi.shape == long_mtx.shape

        assert type(potential[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta[0, 0, 0]) in [float, np.float64]
        assert type(potential_dphi[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2theta[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2phi[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta_dphi[0, 0, 0]) in [float, np.float64]

def test_tidal_potential_general_obliquity_low_e_nsr_modes():
    """ Test the modal tidal potential assuming low eccentricity, general obliquity, and synchronous rotation """
    # Test arrays - Static=False
    spin_freq = 1.5 * orbital_freq
    freqs, modes, potential_tuple_by_mode = \
        tidal_potential_gen_obliquity_low_e_nsr_modes(
            radius_array[-1], long_mtx, colat_mtx, time_mtx,
            orbital_freq, spin_freq,
            eccentricity, obliquity,
            host_mass, semi_major_axis,
            use_static=False
            )

    # Check mode frequency types
    assert len(modes) == 17
    assert len(freqs) == 17
    for mode_label, mode_freq in modes.items():
        assert type(mode_label) == str
        assert type(mode_freq) in [float, np.float64]

        # Unpack potential tuple
        potential, potential_dtheta, potential_dphi, potential_d2theta, potential_d2phi, potential_dtheta_dphi = \
            potential_tuple_by_mode[mode_label]

        assert type(potential) == np.ndarray
        assert type(potential_dtheta) == np.ndarray
        assert type(potential_dphi) == np.ndarray
        assert type(potential_d2theta) == np.ndarray
        assert type(potential_d2phi) == np.ndarray
        assert type(potential_dtheta_dphi) == np.ndarray

        assert potential.shape == long_mtx.shape
        assert potential_dtheta.shape == long_mtx.shape
        assert potential_dphi.shape == long_mtx.shape
        assert potential_d2theta.shape == long_mtx.shape
        assert potential_d2phi.shape == long_mtx.shape
        assert potential_dtheta_dphi.shape == long_mtx.shape

        assert type(potential[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta[0, 0, 0]) in [float, np.float64]
        assert type(potential_dphi[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2theta[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2phi[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta_dphi[0, 0, 0]) in [float, np.float64]

    # Test arrays - Static=True
    spin_freq = 1.5 * orbital_freq
    freqs, modes, potential_tuple_by_mode = \
        tidal_potential_gen_obliquity_low_e_nsr_modes(
            radius_array[-1], long_mtx, colat_mtx, time_mtx,
            orbital_freq, spin_freq,
            eccentricity, obliquity,
            host_mass, semi_major_axis,
            use_static=True
            )

    # Check mode frequency types
    assert len(modes) == 17
    assert len(freqs) == 17
    for mode_label, mode_freq in modes.items():
        assert type(mode_label) == str
        assert type(mode_freq) in [float, np.float64]

        # Unpack potential tuple
        potential, potential_dtheta, potential_dphi, potential_d2theta, potential_d2phi, potential_dtheta_dphi = \
            potential_tuple_by_mode[mode_label]

        assert type(potential) == np.ndarray
        assert type(potential_dtheta) == np.ndarray
        assert type(potential_dphi) == np.ndarray
        assert type(potential_d2theta) == np.ndarray
        assert type(potential_d2phi) == np.ndarray
        assert type(potential_dtheta_dphi) == np.ndarray

        assert potential.shape == long_mtx.shape
        assert potential_dtheta.shape == long_mtx.shape
        assert potential_dphi.shape == long_mtx.shape
        assert potential_d2theta.shape == long_mtx.shape
        assert potential_d2phi.shape == long_mtx.shape
        assert potential_dtheta_dphi.shape == long_mtx.shape

        assert type(potential[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta[0, 0, 0]) in [float, np.float64]
        assert type(potential_dphi[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2theta[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2phi[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta_dphi[0, 0, 0]) in [float, np.float64]

def test_tidal_potential_general_obliquity_nsr_modes():
    """ Test the modal tidal potential assuming moderate eccentricity, general obliquity, and synchronous rotation """
    # Test arrays - Static=False
    spin_freq = 1.5 * orbital_freq
    freqs, modes, potential_tuple_by_mode = \
        tidal_potential_gen_obliquity_nsr_modes(
            radius_array[-1], long_mtx, colat_mtx, time_mtx,
            orbital_freq, spin_freq,
            eccentricity, obliquity,
            host_mass, semi_major_axis,
            use_static=False
            )

    # Check mode frequency types
    assert len(modes) == 27
    assert len(freqs) == 27
    for mode_label, mode_freq in modes.items():
        assert type(mode_label) == str
        assert type(mode_freq) in [float, np.float64]

        # Unpack potential tuple
        potential, potential_dtheta, potential_dphi, potential_d2theta, potential_d2phi, potential_dtheta_dphi = \
            potential_tuple_by_mode[mode_label]

        assert type(potential) == np.ndarray
        assert type(potential_dtheta) == np.ndarray
        assert type(potential_dphi) == np.ndarray
        assert type(potential_d2theta) == np.ndarray
        assert type(potential_d2phi) == np.ndarray
        assert type(potential_dtheta_dphi) == np.ndarray

        assert potential.shape == long_mtx.shape
        assert potential_dtheta.shape == long_mtx.shape
        assert potential_dphi.shape == long_mtx.shape
        assert potential_d2theta.shape == long_mtx.shape
        assert potential_d2phi.shape == long_mtx.shape
        assert potential_dtheta_dphi.shape == long_mtx.shape

        assert type(potential[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta[0, 0, 0]) in [float, np.float64]
        assert type(potential_dphi[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2theta[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2phi[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta_dphi[0, 0, 0]) in [float, np.float64]

    # Test arrays - Static=True
    spin_freq = 1.5 * orbital_freq
    freqs, modes, potential_tuple_by_mode = \
        tidal_potential_gen_obliquity_nsr_modes(
            radius_array[-1], long_mtx, colat_mtx, time_mtx,
            orbital_freq, spin_freq,
            eccentricity, obliquity,
            host_mass, semi_major_axis,
            use_static=True
            )

    # Check mode frequency types
    assert len(modes) == 27
    assert len(freqs) == 27
    for mode_label, mode_freq in modes.items():
        assert type(mode_label) == str
        assert type(mode_freq) in [float, np.float64]

        # Unpack potential tuple
        potential, potential_dtheta, potential_dphi, potential_d2theta, potential_d2phi, potential_dtheta_dphi = \
            potential_tuple_by_mode[mode_label]

        assert type(potential) == np.ndarray
        assert type(potential_dtheta) == np.ndarray
        assert type(potential_dphi) == np.ndarray
        assert type(potential_d2theta) == np.ndarray
        assert type(potential_d2phi) == np.ndarray
        assert type(potential_dtheta_dphi) == np.ndarray

        assert potential.shape == long_mtx.shape
        assert potential_dtheta.shape == long_mtx.shape
        assert potential_dphi.shape == long_mtx.shape
        assert potential_d2theta.shape == long_mtx.shape
        assert potential_d2phi.shape == long_mtx.shape
        assert potential_dtheta_dphi.shape == long_mtx.shape

        assert type(potential[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta[0, 0, 0]) in [float, np.float64]
        assert type(potential_dphi[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2theta[0, 0, 0]) in [float, np.float64]
        assert type(potential_d2phi[0, 0, 0]) in [float, np.float64]
        assert type(potential_dtheta_dphi[0, 0, 0]) in [float, np.float64]

def test_tidal_potential_nsr_modes_vs_non_mode():
    """ Test the modal tidal potential assuming moderate eccentricity, no obliquity, and synchronous rotation vs.
        the non-modal version"""
    # First find the mode version
    # Test arrays - Static=False
    spin_freq = 1.5 * orbital_freq
    freqs, modes, potential_tuple_by_mode = \
        tidal_potential_nsr_modes(
            radius_array[-1], long_mtx, colat_mtx, time_mtx,
            orbital_freq, spin_freq,
            eccentricity,
            host_mass, semi_major_axis,
            use_static=False
            )

    # Then calculate the non-mode version.
    freqs_reg, modes_reg, potential_tuple_by_mode_reg = \
        tidal_potential_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=spin_freq,
            eccentricity=eccentricity,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode_reg['n']

    for mode_pot, reg_pot in [(0, potential),
                              (1, potential_partial_theta),
                              (2, potential_partial_phi),
                              (3, potential_partial2_theta2),
                              (4, potential_partial2_phi2),
                              (5, potential_partial2_theta_phi)]:

        # We need to collapse the mode version so that it matches the non-mode version.
        # They will only match if we make the CPL assumption (linearly add the modes together).
        collapsed_mode_pot = sum(
            [potential_[mode_pot] for potential_ in list(potential_tuple_by_mode.values())]
            )

        assert collapsed_mode_pot.shape == reg_pot.shape
        assert np.allclose(collapsed_mode_pot, reg_pot)


def test_tidal_potential_obliquity_nsr_modes_vs_non_mode():
    """ Test the modal tidal potential assuming moderate eccentricity, moderate obliquity, and synchronous rotation vs.
        the non-modal version"""
    # First find the mode version
    # Test arrays - Static=False
    spin_freq = 1.5 * orbital_freq
    freqs, modes, potential_tuple_by_mode = \
        tidal_potential_obliquity_nsr_modes(
            radius_array[-1], long_mtx, colat_mtx, time_mtx,
            orbital_freq, spin_freq,
            eccentricity, obliquity,
            host_mass, semi_major_axis,
            use_static=False
            )

    # Then calculate the non-mode version.
    freqs_reg, modes_reg, potential_tuple_by_mode_reg = \
        tidal_potential_obliquity_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=spin_freq,
            eccentricity=eccentricity, obliquity=obliquity,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode_reg['n']

    for mode_pot, reg_pot in [(0, potential),
                              (1, potential_partial_theta),
                              (2, potential_partial_phi),
                              (3, potential_partial2_theta2),
                              (4, potential_partial2_phi2),
                              (5, potential_partial2_theta_phi)]:

        # We need to collapse the mode version so that it matches the non-mode version.
        # They will only match if we make the CPL assumption (linearly add the modes together).
        collapsed_mode_pot = sum(
            [potential_[mode_pot] for potential_ in list(potential_tuple_by_mode.values())]
            )

        assert collapsed_mode_pot.shape == reg_pot.shape
        assert np.allclose(collapsed_mode_pot, reg_pot)
