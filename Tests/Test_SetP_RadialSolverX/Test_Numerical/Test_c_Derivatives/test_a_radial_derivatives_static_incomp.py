""" Test the `TidalPy.radial_solver.radial_derivatives_static_incomp` (cython extension) functionality. """

import numpy as np
import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.radial_solver.numerical.derivatives.radial_derivatives_static_incomp_wrap import (
    dy_solid_static_incompressible_wrap, dy_liquid_static_incompressible_wrap
    )

radius = 100.
gravity = 1.5
density = 3000.
shear_modulus = 1.0e10 + 200.j
bulk_modulus = 1.0e10
G = 6.67430e-11


@pytest.mark.parametrize('degree_l', (2, 3))
def test_dy_solid_static_incompressible(degree_l):

    # This function takes 6 complex-valued y variables, so we need to give it 12 floats.
    y = np.asarray(
        (100., 2., 100., 2., 300., 3., 400., 4., 600., 9., -33., 2.), dtype=np.float64
        )
    dy = np.empty(12, dtype=np.float64)

    # Check dimensions make sense
    dy_solid_static_incompressible_wrap(
        radius, y, dy, shear_modulus, density, gravity, degree_l, G_to_use=G
        )

    assert dy.shape == (12,)
    assert dy.dtype == np.float64
    # Check that the results match expectations
    if degree_l == 2:
        expected_dy = np.asarray(
            [1.600e+01, 1.400e-01, -9.600e+09, -8.400e+07, 2.000e+00, 1.000e-02,
                6.000e+09, 5.400e+07, -5.100e+01, 1.730e+00, -3.300e-01, 2.000e-02], dtype=np.float64
            )
    elif degree_l == 3:
        expected_dy = np.asarray(
            [3.400e+01, 3.200e-01, -2.040e+10, -1.920e+08, 2.000e+00, 1.000e-02,
                1.320e+10, 1.260e+08, -5.700e+01, 1.640e+00, -6.601e-01, 4.000e-02], dtype=np.float64
            )
    assert np.allclose(dy, expected_dy, rtol=0.001)

    # Check that results don't change if we use the default G
    dy2 = np.empty(12, dtype=np.float64)
    dy_solid_static_incompressible_wrap(
        radius, y, dy2, shear_modulus, density, gravity, degree_l
        )

    assert np.allclose(dy, dy2)

    # Check that they do change if we don't
    dy3 = np.empty(12, dtype=np.float64)
    dy_solid_static_incompressible_wrap(
        radius, y, dy3, shear_modulus, density, gravity, degree_l, G_to_use=1.
        )

    assert not np.allclose(dy, dy3)


@pytest.mark.parametrize('degree_l', (2, 3))
def test_dy_liquid_static_incompressible(degree_l):

    # This function takes 2 complex-valued y variables, so we need to give it 4 floats.
    y = np.asarray(
        (100., 2., 100., 2.), dtype=np.float64
        )
    dy = np.empty(4, dtype=np.float64)

    # Check dimensions make sense
    dy_liquid_static_incompressible_wrap(
        radius, y, dy, density, gravity, degree_l, G_to_use=G
        )

    assert dy.shape == (4,)
    assert dy.dtype == np.float64
    # Check that the results match expectations
    if degree_l == 2:
        expected_dy = np.asarray(
            [9.700e+01, 1.940e+00, 9.998e-01, 2.000e-02], dtype=np.float64
            )
    elif degree_l == 3:
        expected_dy = np.asarray(
            [9.600e+01, 1.920e+00, 2.000e+00, 4.000e-02], dtype=np.float64
            )
    assert np.allclose(dy, expected_dy, rtol=0.001)

    # Check that results don't change if we use the default G
    dy2 = np.empty(4, dtype=np.float64)
    dy_liquid_static_incompressible_wrap(
        radius, y, dy2, density, gravity, degree_l
        )

    assert np.allclose(dy, dy2)

    # Check that they do change if we don't
    dy3 = np.empty(4, dtype=np.float64)
    dy_liquid_static_incompressible_wrap(
        radius, y, dy3, density, gravity, degree_l, G_to_use=1.
        )

    assert not np.allclose(dy, dy3)