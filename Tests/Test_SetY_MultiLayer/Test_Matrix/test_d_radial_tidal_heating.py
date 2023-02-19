""" Tests for extracting useful information out of a multilayer tidal propagation
"""

import numpy as np

import TidalPy
TidalPy.test_mode()

from TidalPy.constants import G
from TidalPy.tides.multilayer.heating import calc_radial_tidal_heating
from TidalPy.radial_solver.sensitivity import sensitivity_to_shear as find_sensitivity_to_shear
from TidalPy.radial_solver.matrix import fundamental_matrix_generic, fundamental_matrix_orderl2
from TidalPy.radial_solver.matrix import matrix_propagate
from TidalPy.toolbox.conversions import orbital_motion2semi_a

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array = np.linspace(0., 1.e6, 11)
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.complex128)
host_mass = 10. * planet_mass
orbital_freq = (2. * np.pi / (86400. * 6.))
semi_major_axis = orbital_motion2semi_a(orbital_freq, host_mass, planet_mass)
eccentricity = 0.01


def test_calc_fundamental_order2():
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

    # Calculate tidal heating as a function of radius
    radial_tidal_heating = calc_radial_tidal_heating(
        eccentricity, orbital_freq, semi_major_axis, host_mass,
        radius_array[1:], sensitivity_to_shear, shear_array,
        order_l=2
        )

    # Check shapes
    assert radial_tidal_heating.shape == (10,)

    # Check types
    assert type(radial_tidal_heating[0]) in [np.float64, float]


def test_calc_fundamental_order3():
    # Calculate the fundamental matrix and its inverse
    F, F_inv, deriv_mtx = fundamental_matrix_generic(
        radius_array[1:], shear_array,
        density_array, gravity_array, order_l=3
        )

    # Central boundary condition
    ## From IcyDwarf: "They are inconsequential on the rest of the solution, so false assumptions are OK."
    core_condition = np.zeros((6, 3), dtype=np.complex128)
    # Roberts & Nimmo (2000): Liquid innermost zone.
    core_condition[2, 0] = 1.0
    core_condition[3, 1] = 1.0
    core_condition[5, 2] = 1.0

    # Find tidal solution
    tidal_y = matrix_propagate(F, F_inv, deriv_mtx, core_condition, world_radius=radius_array[-1], order_l=3)

    # Decompose the results
    sensitivity_to_shear = find_sensitivity_to_shear(
        tidal_y, radius_array[1:], shear_array,
        bulk_modulus_array=200.e9 * np.ones_like(radius_array[1:]))

    # Calculate tidal heating as a function of radius
    radial_tidal_heating = calc_radial_tidal_heating(
        eccentricity, orbital_freq, semi_major_axis, host_mass,
        radius_array[1:], sensitivity_to_shear, shear_array,
        order_l=3
        )

    # Check shapes
    assert radial_tidal_heating.shape == (10,)

    # Check types
    assert type(radial_tidal_heating[0]) in [np.float64, float]
