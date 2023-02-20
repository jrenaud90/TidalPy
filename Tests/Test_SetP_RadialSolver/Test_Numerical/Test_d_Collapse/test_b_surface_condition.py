""" Test the `surface_condition` functions from the `TidalPy.radial_solver.numerical.collapse` module. """

import numpy as np
import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.constants import G
from TidalPy.radial_solver.numerical.collapse import solid_surface, static_liquid_surface, dynamic_liquid_surface

radius = 1.0e6
planet_bulk_density = 4000.
surface_gravity = 9.81

second2_conversion = 1. / (np.pi * G * 3000.)
mass_conversion = 3000. * radius**3
length_conversion = radius

fake_y_values = np.asarray([
    (1.38521163e+01-6.43810616e-01j) / (second2_conversion / length_conversion),
    (8.67690181e-14+4.33845091e-14j) / (mass_conversion / length_conversion**3),
    (3.17548877e+00-3.10732612e-01j) / (second2_conversion / length_conversion),
    (-1.22018932e-14-4.33845091e-14j) / (mass_conversion / length_conversion**3),
    (1.55681773e+00-2.61835144e-02j),
    (1.98412698e-05+1.65211760e-22j) / (1. / length_conversion)],
        np.complex128)


@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_solid_surface(order_l):
    """ Test the `surface_condition` function for a solid layer surface. """
    y_solutions_at_surface = (
        np.asarray(
            [
            fake_y_values[0],
            fake_y_values[1],
            fake_y_values[2],
            fake_y_values[3],
            fake_y_values[4],
            fake_y_values[5]
            ], dtype=np.complex128
        ),
        np.asarray(
            [
            3*fake_y_values[0],
            3*fake_y_values[1],
            3*fake_y_values[2],
            3*fake_y_values[3],
            3*fake_y_values[4],
            3*fake_y_values[5]
            ], dtype=np.complex128
        ),
        np.asarray(
            [
            7*fake_y_values[0],
            7*fake_y_values[1],
            7*fake_y_values[2],
            7*fake_y_values[3],
            7*fake_y_values[4],
            7*fake_y_values[5]
            ], dtype=np.complex128
        )
    )

    # Test tidal surface boundary condition.
    tidal_boundary_condition = np.zeros(3, dtype=np.complex128)
    tidal_boundary_condition[2] = (2. * order_l + 1.) / radius

    try:
        tidal_result = solid_surface(y_solutions_at_surface, tidal_boundary_condition)
    except np.linalg.LinAlgError:
        # TODO
        pytest.skip("Singular Matrix was found. This can happen on some machines since the y's here are fake. ")

    # Check type
    assert tidal_result.dtype == np.complex128

    # Check shape
    assert tidal_result.shape == (3,)

    # Test loading surface boundary condition.
    load_boundary_condition = np.zeros(3, dtype=np.complex128)
    load_boundary_condition[0] = -(2. * order_l + 1.) * planet_bulk_density / 3.
    load_boundary_condition[2] = (2. * order_l + 1.) / radius

    load_result = solid_surface(y_solutions_at_surface, load_boundary_condition)

    # Check type
    assert load_result.dtype == np.complex128

    # Check shape
    assert load_result.shape == (3,)


@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_dynamic_liquid_surface(order_l):
    """ Test the `surface_condition` function for a dynamic liquid layer surface. """
    y_solutions_at_surface = (
        np.asarray(
            [
            fake_y_values[0],
            fake_y_values[1],
            fake_y_values[4],
            fake_y_values[5]
            ], dtype=np.complex128
        ),
        np.asarray(
            [
            7*fake_y_values[0],
            7*fake_y_values[1],
            7*fake_y_values[4],
            7*fake_y_values[5]
            ], dtype=np.complex128
        )
    )

    # Test tidal surface boundary condition.
    tidal_boundary_condition = np.zeros(3, dtype=np.complex128)
    tidal_boundary_condition[2] = (2. * order_l + 1.) / radius

    try:
        tidal_result = dynamic_liquid_surface(y_solutions_at_surface, tidal_boundary_condition)
    except np.linalg.LinAlgError:
        # TODO
        pytest.skip("Singular Matrix was found. This can happen on some machines since the y's here are fake. ")

    # Check type
    assert tidal_result.dtype == np.complex128

    # Check shape
    assert tidal_result.shape == (2,)

    # Test loading surface boundary condition.
    load_boundary_condition = np.zeros(3, dtype=np.complex128)
    load_boundary_condition[0] = -(2. * order_l + 1.) * planet_bulk_density / 3.
    load_boundary_condition[2] = (2. * order_l + 1.) / radius

    load_result = dynamic_liquid_surface(y_solutions_at_surface, load_boundary_condition)

    # Check type
    assert load_result.dtype == np.complex128

    # Check shape
    assert load_result.shape == (2,)


@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_static_liquid_surface(order_l):
    """ Test the `surface_condition` function for a static liquid layer surface. """
    y_solutions_at_surface = (
        np.asarray(
            [
            fake_y_values[4],
            fake_y_values[5]
            ], dtype=np.complex128
        ),
    )

    # Test tidal surface boundary condition.
    tidal_boundary_condition = np.zeros(3, dtype=np.complex128)
    tidal_boundary_condition[2] = (2. * order_l + 1.) / radius

    try:
        tidal_result = static_liquid_surface(y_solutions_at_surface, tidal_boundary_condition, surface_gravity,
                                             G_to_use=G)
    except np.linalg.LinAlgError:
        # TODO
        pytest.skip("Singular Matrix was found. This can happen on some machines since the y's here are fake. ")

    # Check type
    assert tidal_result.dtype == np.complex128

    # Check shape
    assert tidal_result.shape == (1,)

    # Test loading surface boundary condition.
    load_boundary_condition = np.zeros(3, dtype=np.complex128)
    load_boundary_condition[0] = -(2. * order_l + 1.) * planet_bulk_density / 3.
    load_boundary_condition[2] = (2. * order_l + 1.) / radius

    load_result = static_liquid_surface(y_solutions_at_surface, tidal_boundary_condition, surface_gravity, G_to_use=G)

    # Check type
    assert load_result.dtype == np.complex128

    # Check shape
    assert load_result.shape == (1,)
