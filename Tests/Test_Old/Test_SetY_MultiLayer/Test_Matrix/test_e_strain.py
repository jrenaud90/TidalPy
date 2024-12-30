""" Tests for calculating tidal strain tensor from the tidal potential and multi-layer solution.
"""

import numpy as np

import TidalPy



from TidalPy.constants import G
from TidalPy.tides.multilayer import calculate_displacements, calculate_strain_stress
from TidalPy.RadialSolver import radial_solver
from TidalPy.tides.potential import tidal_potential_nsr, tidal_potential_simple
from TidalPy.utilities.conversions import orbital_motion2semi_a

# Model planet - 2layers
density_array = 5000. * np.ones(10, dtype=np.float64, order='C')
radius_array = np.linspace(0., 1.e6, 10, dtype=np.float64)
longitude_array = np.radians(np.linspace(0., 360., 12, dtype=np.float64))
colat_array = np.radians(np.linspace(0.5, 179.5, 13, dtype=np.float64))

volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
volume_array = np.insert(volume_array, 0, 0.0)
mass_array = volume_array * density_array
planet_bulk_density = np.average(density_array)
planet_mass = sum(mass_array)
layer_types = ('solid',)
is_static_bylayer = (True,)
is_incompressible_bylayer = (True,)
upper_radius_bylayer_array = np.asarray((radius_array[-1],))
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(10)])
complex_shear_modulus_array = 5.e10 * np.ones(10, dtype=np.complex128, order='C')
complex_bulk_modulus_array = 10.e10 * np.ones(10, dtype=np.complex128, order='C')
host_mass = 50000. * planet_mass
orbital_freq = (2. * np.pi / (86400. * 6.))
semi_major_axis = orbital_motion2semi_a(orbital_freq, host_mass, planet_mass)
eccentricity = 0.05
obliquity = np.radians(15.)

time_array = np.linspace(0., 2. * np.pi / orbital_freq, 5, dtype=np.float64)

long_mtx, colat_mtx, time_mtx = np.meshgrid(longitude_array, colat_array, time_array)


def test_calc_displacements():
    # Calculate the fundamental matrix and its inverse
    radial_solution = radial_solver(
        radius_array,
        density_array,
        complex_bulk_modulus_array,
        complex_shear_modulus_array,
        orbital_freq,
        planet_bulk_density,
        layer_types,
        is_static_bylayer,
        is_incompressible_bylayer,
        upper_radius_bylayer_array,
        surface_pressure = 0.0,
        degree_l = 2,
        core_model = 1,
        use_prop_matrix = True,
        perform_checks = True
        )

    # Calculate tidal potential and its partial derivatives
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_simple(
                radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
                orbital_frequency=orbital_freq,
                eccentricity=eccentricity,
                host_mass=host_mass, semi_major_axis=semi_major_axis
                )

    potential, potential_partial_theta, potential_partial_phi, \
        potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    # Calculate displacements
    radial_displacement, polar_displacement, azimuthal_displacement = \
        calculate_displacements(
                potential, potential_partial_theta, potential_partial_phi, radial_solution.result, colatitude=colat_mtx
                )

    # Check shapes
    shape = (*radius_array.shape, *colat_mtx.shape)
    assert radial_displacement.shape == shape
    assert polar_displacement.shape == shape
    assert azimuthal_displacement.shape == shape

    # Check types
    assert radial_displacement.dtype in [np.complex128, complex]
    assert polar_displacement.dtype in [np.complex128, complex]
    assert azimuthal_displacement.dtype in [np.complex128, complex]


def test_calc_strains_simple():
    # Calculate the fundamental matrix and its inverse
    radial_solution = radial_solver(
        radius_array,
        density_array,
        complex_bulk_modulus_array,
        complex_shear_modulus_array,
        orbital_freq,
        planet_bulk_density,
        layer_types,
        is_static_bylayer,
        is_incompressible_bylayer,
        upper_radius_bylayer_array,
        surface_pressure = 0.0,
        degree_l = 2,
        core_model = 1,
        use_prop_matrix = True,
        perform_checks = True
        )

    # Calculate tidal potential and its partial derivatives
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_simple(
                radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
                orbital_frequency=orbital_freq,
                eccentricity=eccentricity,
                host_mass=host_mass, semi_major_axis=semi_major_axis
                )

    potential, potential_partial_theta, potential_partial_phi, \
        potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    # Calculate strain tensor
    strain_components, stress_components = \
        calculate_strain_stress(
                potential, potential_partial_theta, potential_partial_phi,
                potential_partial2_theta2, potential_partial2_phi2,
                potential_partial2_theta_phi,
                radial_solution.result,
                longitude_array=longitude_array, colatitude_array=colat_array, time_array=time_array,
                radius_array=radius_array, shear_moduli=radial_solution.shear_modulus_array,
                bulk_moduli=radial_solution.bulk_modulus_array, frequency=orbital_freq, degree_l=2
                )

    shape = (6, *radius_array.shape, longitude_array.size, colat_array.size, time_array.size)
    # Check shape
    assert strain_components.shape == shape
    assert stress_components.shape == shape

    # Check type
    assert strain_components.dtype in [np.complex128, complex]
    assert stress_components.dtype in [np.complex128, complex]


def test_calc_strains_nsr():
    # Calculate the fundamental matrix and its inverse
    radial_solution = radial_solver(
        radius_array,
        density_array,
        complex_bulk_modulus_array,
        complex_shear_modulus_array,
        orbital_freq,
        planet_bulk_density,
        layer_types,
        is_static_bylayer,
        is_incompressible_bylayer,
        upper_radius_bylayer_array,
        surface_pressure = 0.0,
        degree_l = 2,
        core_model = 1,
        use_prop_matrix = True,
        perform_checks = True
        )

    # Calculate tidal potential and its partial derivatives
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_nsr(
                radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
                orbital_frequency=orbital_freq, rotation_frequency=5. * orbital_freq,
                eccentricity=eccentricity,
                host_mass=host_mass, semi_major_axis=semi_major_axis, use_static=False
                )

    potential, potential_partial_theta, potential_partial_phi, \
        potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    # Calculate strain tensor
    strain_components, stress_components = \
        calculate_strain_stress(
                potential, potential_partial_theta, potential_partial_phi,
                potential_partial2_theta2, potential_partial2_phi2,
                potential_partial2_theta_phi,
                radial_solution.result,
                longitude_array, colat_array, time_array,
                radius_array, radial_solution.shear_modulus_array,
                radial_solution.bulk_modulus_array, frequency=orbital_freq, degree_l=2
                )

    shape = (6, *radius_array.shape, longitude_array.size, colat_array.size, time_array.size)
    # Check shape
    assert strain_components.shape == shape
    assert stress_components.shape == shape

    # Check type
    assert strain_components.dtype in [np.complex128, complex]
    assert stress_components.dtype in [np.complex128, complex]
