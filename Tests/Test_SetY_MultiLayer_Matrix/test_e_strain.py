""" Tests for calculating tidal strain tensor from the tidal potential and multi-layer solution.
"""

import numpy as np

import TidalPy
from TidalPy.tides.multilayer.matrix.fundamental_solid import fundamental_matrix_orderl2
from TidalPy.tides.multilayer.matrix.propagate import propagate
from TidalPy.tides.multilayer.decompose import decompose
from TidalPy.tides.potential import tidal_potential_simple
from TidalPy.tides.multilayer.stress_strain import calculate_strain, calculate_displacements
from TidalPy.tools.conversions import orbital_motion2semi_a
from TidalPy.constants import G

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array =  np.linspace(0., 1.e6, 11)
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i+1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.complex)
host_mass = 10. * planet_mass
orbital_freq = (2. * np.pi / (86400. * 6.))
semi_major_axis = orbital_motion2semi_a(orbital_freq, host_mass, planet_mass)
eccentricity = 0.01

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
    tidal_y, tidal_y_deriv = propagate(F, F_inv, deriv_mtx, core_condition, world_radius=radius_array[-1], order_l=2)

    # Decompose the results
    sensitivity_to_shear, (k, h, l) = decompose(tidal_y, tidal_y_deriv, radius_array[1:], gravity_array,
                                                shear_array, bulk_modulus=200.0e9, order_l=2)

    # Calculate tidal potential and its partial derivatives
    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = \
        tidal_potential_simple(radius_array[-1], longitude=0.1, colatitude=0.1, orbital_frequency=orbital_freq,
                               eccentricity=eccentricity, time=1000.)

    # Calculate displacements
    radial_displacement, polar_displacement, azimuthal_displacement = \
        calculate_displacements(potential, potential_partial_theta, potential_partial_phi, tidal_y, colatitude=0.1)

    # Check shapes
    assert radial_displacement.shape == (10,)
    assert polar_displacement.shape == (10,)
    assert azimuthal_displacement.shape == (10,)

    # Check types
    assert type(radial_displacement[0]) in [np.complex, np.complex128, complex]
    assert type(polar_displacement[0]) in [np.complex, np.complex128, complex]
    assert type(azimuthal_displacement[0]) in [np.complex, np.complex128, complex]

def test_calc_strains():
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
    tidal_y, tidal_y_deriv = propagate(F, F_inv, deriv_mtx, core_condition, world_radius=radius_array[-1], order_l=2)

    # Decompose the results
    sensitivity_to_shear, (k, h, l) = decompose(tidal_y, tidal_y_deriv, radius_array[1:], gravity_array,
                                                shear_array, bulk_modulus=200.0e9, order_l=2)

    # Calculate tidal potential and its partial derivatives
    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = \
        tidal_potential_simple(radius_array[1:], longitude=0.1, colatitude=0.1, orbital_frequency=orbital_freq,
                               eccentricity=eccentricity, time=1000.)

    # Calculate strain tensor
    e_rr, e_thth, e_phph, e_rth, e_rph, e_thph = \
        calculate_strain(potential, potential_partial_theta, potential_partial_phi,
                         potential_partial2_theta2, potential_partial2_phi2,
                         potential_partial2_theta_phi, tidal_y, tidal_y_deriv,
                         colatitude=0.1, radius=radius_array[1:], shear_moduli=shear_array)

    for strain_component in [e_rr, e_thth, e_phph, e_rth, e_rph, e_thph]:
        # Check shape
        assert strain_component.shape == (10,)

        # Check type
        assert type(strain_component[0]) in [np.complex128, np.complex, complex]