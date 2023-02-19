""" Tests for calculating tidal strain tensor from the tidal potential and multi-layer solution.
"""

import numpy as np

import TidalPy

TidalPy.test_mode()

from TidalPy.constants import G
from TidalPy.tides.multilayer import calculate_displacements, calculate_strain_stress
from TidalPy.radial_solver.sensitivity import sensitivity_to_shear as find_sensitivity_to_shear
from TidalPy.radial_solver.matrix import fundamental_matrix_orderl2
from TidalPy.radial_solver.matrix import matrix_propagate
from TidalPy.tides.potential import tidal_potential_nsr, tidal_potential_simple
from TidalPy.toolbox.conversions import orbital_motion2semi_a

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
bulk_array = 10.e10 * np.ones(10, dtype=np.complex128)
host_mass = 50000. * planet_mass
orbital_freq = (2. * np.pi / (86400. * 6.))
semi_major_axis = orbital_motion2semi_a(orbital_freq, host_mass, planet_mass)
eccentricity = 0.05
obliquity = np.radians(15.)

time_array = np.linspace(0., 2. * np.pi / orbital_freq, 5)

long_mtx, colat_mtx, time_mtx = np.meshgrid(longitude_array, colat_array, time_array)


def test_calc_displacements():
    # Calculate the fundamental matrix and its inverse
    F, F_inv, deriv_mtx = fundamental_matrix_orderl2(radius_array[1:], shear_array, density_array, gravity_array)

    # Central boundary condition
    ## From IcyDwarf: "They are inconsequential on the rest of the solution, so false assumptions are OK."
    core_condition = np.zeros((6, 3), dtype=np.complex128)
    # Roberts & Nimmo (2000): Liquid innermost zone.
    core_condition[2, 0] = 1.0
    core_condition[3, 1] = 1.0
    core_condition[5, 2] = 1.0

    # Find tidal solution
    tidal_y = matrix_propagate(F, F_inv, deriv_mtx, core_condition, world_radius=radius_array[-1], order_l=2)

    # Decompose the results
    sensitivity_to_shear = find_sensitivity_to_shear(tidal_y, radius_array[1:], shear_array,
                                                     bulk_modulus_array=200.e9 * np.ones_like(radius_array[1:]))

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
                potential, potential_partial_theta, potential_partial_phi, tidal_y, colatitude=colat_mtx
                )

    # Check shapes
    shape = (*radius_array[1:].shape, *colat_mtx.shape)
    assert radial_displacement.shape == shape
    assert polar_displacement.shape == shape
    assert azimuthal_displacement.shape == shape

    # Check types
    assert radial_displacement.dtype in [np.complex128, complex]
    assert polar_displacement.dtype in [np.complex128, complex]
    assert azimuthal_displacement.dtype in [np.complex128, complex]


def test_calc_strains_simple():
    # Calculate the fundamental matrix and its inverse
    F, F_inv, deriv_mtx = fundamental_matrix_orderl2(radius_array[1:], shear_array, density_array, gravity_array)

    # Central boundary condition
    ## From IcyDwarf: "They are inconsequential on the rest of the solution, so false assumptions are OK."
    core_condition = np.zeros((6, 3), dtype=np.complex128)
    # Roberts & Nimmo (2000): Liquid innermost zone.
    core_condition[2, 0] = 1.0
    core_condition[3, 1] = 1.0
    core_condition[5, 2] = 1.0

    # Find tidal solution
    tidal_y = matrix_propagate(F, F_inv, deriv_mtx, core_condition, world_radius=radius_array[-1], order_l=2)

    # Decompose the results
    sensitivity_to_shear = find_sensitivity_to_shear(
        tidal_y, radius_array[1:], shear_array,
        bulk_modulus_array=200.e9 * np.ones_like(radius_array[1:]))

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
                potential_partial2_theta_phi, tidal_y,
                colatitude=colat_mtx, radius=radius_array[1:], shear_moduli=shear_array,
                bulk_moduli=bulk_array, frequency=orbital_freq, order_l=2
                )

    shape = (6, *radius_array[1:].shape, *colat_mtx.shape)
    # Check shape
    assert strain_components.shape == shape
    assert stress_components.shape == shape

    # Check type
    assert strain_components.dtype in [np.complex128, complex]
    assert stress_components.dtype in [np.complex128, complex]


def test_calc_strains_nsr():
    # Calculate the fundamental matrix and its inverse
    F, F_inv, deriv_mtx = fundamental_matrix_orderl2(radius_array[1:], shear_array, density_array, gravity_array)

    # Central boundary condition
    ## From IcyDwarf: "They are inconsequential on the rest of the solution, so false assumptions are OK."
    core_condition = np.zeros((6, 3), dtype=np.complex128)
    # Roberts & Nimmo (2000): Liquid innermost zone.
    core_condition[2, 0] = 1.0
    core_condition[3, 1] = 1.0
    core_condition[5, 2] = 1.0

    # Find tidal solution
    tidal_y = matrix_propagate(F, F_inv, deriv_mtx, core_condition, world_radius=radius_array[-1], order_l=2)

    # Decompose the results
    sensitivity_to_shear = find_sensitivity_to_shear(
        tidal_y, radius_array[1:], shear_array,
        bulk_modulus_array=200.e9 * np.ones_like(radius_array[1:]))

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
                potential_partial2_theta_phi, tidal_y,
                colatitude=colat_mtx, radius=radius_array[1:], shear_moduli=shear_array,
                bulk_moduli=bulk_array, frequency=orbital_freq, order_l=2
                )

    shape = (6, *radius_array[1:].shape, *colat_mtx.shape)
    # Check shape
    assert strain_components.shape == shape
    assert stress_components.shape == shape

    # Check type
    assert strain_components.dtype in [np.complex128, complex]
    assert stress_components.dtype in [np.complex128, complex]
