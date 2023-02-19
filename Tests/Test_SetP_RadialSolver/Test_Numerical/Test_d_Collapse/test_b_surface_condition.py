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


@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_solid_surface(order_l):
    """ Test the `surface_condition` function for a solid layer surface. """
    y_solutions_at_surface = (
        np.asarray(
            [
            (100. + 20j),
            (3000. + 50j),
            (1012. + 3j),
            (30. + 1.51j),
            (1.55 + 2j),
            (11. + 0.3420j)
            ], dtype=np.complex128
        ),
        np.asarray(
            [
            2*(100. + 20j),
            2*(3000. + 50j),
            2*(1012. + 3j),
            2*(30. + 1.51j),
            2*(1.55 + 2j),
            2*(11. + 0.3420j)
            ], dtype=np.complex128
        ),
        np.asarray(
            [
            3*(100. + 20j),
            3*(3000. + 50j),
            3*(1012. + 3j),
            3*(30. + 1.51j),
            3*(1.55 + 2j),
            3*(11. + 0.3420j)
            ], dtype=np.complex128
        )
    )

    # Test tidal surface boundary condition.
    tidal_boundary_condition = np.zeros(3, dtype=np.complex128)
    tidal_boundary_condition[2] = (2. * order_l + 1.) / radius

    tidal_result = solid_surface(y_solutions_at_surface, tidal_boundary_condition)

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
            (100. + 20j),
            (3000. + 50j),
            (1.55 + 2j),
            (11. + 0.3420j)
            ], dtype=np.complex128
        ),
        np.asarray(
            [
            2*(100. + 20j),
            2*(3000. + 50j),
            2*(1.55 + 2j),
            2*(11. + 0.3420j)
            ], dtype=np.complex128
        )
    )

    # Test tidal surface boundary condition.
    tidal_boundary_condition = np.zeros(3, dtype=np.complex128)
    tidal_boundary_condition[2] = (2. * order_l + 1.) / radius

    tidal_result = dynamic_liquid_surface(y_solutions_at_surface, tidal_boundary_condition)

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
            (1.55 + 2j),
            (11. + 0.3420j)
            ], dtype=np.complex128
        ),
    )

    # Test tidal surface boundary condition.
    tidal_boundary_condition = np.zeros(3, dtype=np.complex128)
    tidal_boundary_condition[2] = (2. * order_l + 1.) / radius

    tidal_result = static_liquid_surface(y_solutions_at_surface, tidal_boundary_condition, surface_gravity, G_to_use=G)

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
