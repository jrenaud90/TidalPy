""" Test the `TidalPy.radial_solver.ode` (cython extension) functionality. """

import numpy as np
import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.radial_solver.numerical.derivatives.odes_x import (
    SolidStaticCompressible, SolidDynamicCompressible, SolidDynamicIncompressible, SolidStaticIncompressible,
    LiquidStaticCompressible, LiquidDynamicCompressible, LiquidDynamicIncompressible, LiquidStaticIncompressible
    )

N = 25
radius_array = np.linspace(0.1, 10., N)
shear_modulus_array = 1.0e10 * np.ones(radius_array.shape, dtype=np.complex128)
density_array = 3000. * np.ones(radius_array.shape, dtype=np.float64)
gravity_array = np.linspace(0.1, 2., N, dtype=np.float64)
bulk_modulus_array = 1.0e11 * np.ones(radius_array.shape, dtype=np.float64)

frequency = 1.0e-5
degree_l = 2
G_to_use = 1.0e-12

t_span = (0.1, 10.)

y0 = np.asarray((
    1.0, -0.2,
    2.0, -0.3,
    1.0, -0.2,
    2.0, -0.3,
    1.0, -0.2,
    2.0, -0.3), dtype=np.float64)

@pytest.mark.parametrize('rk_method', (1, 2))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_solid_static_compressible_solver(degree_l, rk_method):

    solver = SolidStaticCompressible(
        radius_array, shear_modulus_array, bulk_modulus_array, density_array, gravity_array, frequency, degree_l,
        G_to_use, t_span, y0, rk_method=rk_method
        )
    solver.solve()

    assert solver.solution_t.shape == (N,)
    assert solver.solution_y.shape == (12, N)
    assert solver.solution_t.dtype is np.float64
    assert solver.solution_t[0] == radius_array[0]
    assert solver.solution_t[-1] == radius_array[-1]
    assert np.all(np.isfinite(solver.solution_t))
    assert np.all(np.isfinite(solver.solution_y))


@pytest.mark.parametrize('rk_method', (1, 2))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_solid_dynamic_compressible_solver(degree_l, rk_method):
    solver = SolidDynamicCompressible(
        radius_array, shear_modulus_array, bulk_modulus_array, density_array, gravity_array, frequency, degree_l,
        G_to_use, t_span, y0, rk_method=rk_method
        )
    solver.solve()

    assert solver.solution_t.shape == (N,)
    assert solver.solution_y.shape == (12, N)
    assert solver.solution_t.dtype is np.float64
    assert solver.solution_t[0] == radius_array[0]
    assert solver.solution_t[-1] == radius_array[-1]
    assert np.all(np.isfinite(solver.solution_t))
    assert np.all(np.isfinite(solver.solution_y))

@pytest.mark.parametrize('rk_method', (1, 2))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_solid_dynamic_incompressible_solver(degree_l, rk_method):
    solver = SolidDynamicIncompressible(
        radius_array, shear_modulus_array, bulk_modulus_array, density_array, gravity_array, frequency, degree_l,
        G_to_use, t_span, y0, rk_method=rk_method
        )
    solver.solve()

    assert solver.solution_t.shape == (N,)
    assert solver.solution_y.shape == (12, N)
    assert solver.solution_t.dtype is np.float64
    assert solver.solution_t[0] == radius_array[0]
    assert solver.solution_t[-1] == radius_array[-1]
    assert np.all(np.isfinite(solver.solution_t))
    assert np.all(np.isfinite(solver.solution_y))

@pytest.mark.parametrize('rk_method', (1, 2))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_solid_static_incompressible_solver(degree_l, rk_method):
    solver = SolidStaticIncompressible(
        radius_array, shear_modulus_array, bulk_modulus_array, density_array, gravity_array, frequency, degree_l,
        G_to_use, t_span, y0, rk_method=rk_method
        )
    solver.solve()

    assert solver.solution_t.shape == (N,)
    assert solver.solution_y.shape == (12, N)
    assert solver.solution_t.dtype is np.float64
    assert solver.solution_t[0] == radius_array[0]
    assert solver.solution_t[-1] == radius_array[-1]
    assert np.all(np.isfinite(solver.solution_t))
    assert np.all(np.isfinite(solver.solution_y))

@pytest.mark.parametrize('rk_method', (1, 2))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_liquid_dynamic_compressible_solver(degree_l, rk_method):

    y0_ = np.copy(y0[:8])
    solver = LiquidDynamicCompressible(
        radius_array, shear_modulus_array, bulk_modulus_array, density_array, gravity_array, frequency, degree_l,
        G_to_use, t_span, y0_, rk_method=rk_method
        )
    solver.solve()

    assert solver.solution_t.shape == (N,)
    assert solver.solution_y.shape == (8, N)
    assert solver.solution_t.dtype is np.float64
    assert solver.solution_t[0] == radius_array[0]
    assert solver.solution_t[-1] == radius_array[-1]
    assert np.all(np.isfinite(solver.solution_t))
    assert np.all(np.isfinite(solver.solution_y))

@pytest.mark.parametrize('rk_method', (1, 2))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_liquid_dynamic_incompressible_solver(degree_l, rk_method):

    y0_ = np.copy(y0[:8])
    solver = LiquidDynamicIncompressible(
        radius_array, shear_modulus_array, bulk_modulus_array, density_array, gravity_array, frequency, degree_l,
        G_to_use, t_span, y0_, rk_method=rk_method
        )
    solver.solve()

    assert solver.solution_t.shape == (N,)
    assert solver.solution_y.shape == (8, N)
    assert solver.solution_t.dtype is np.float64
    assert solver.solution_t[0] == radius_array[0]
    assert solver.solution_t[-1] == radius_array[-1]
    assert np.all(np.isfinite(solver.solution_t))
    assert np.all(np.isfinite(solver.solution_y))

@pytest.mark.parametrize('rk_method', (1, 2))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_liquid_static_incompressible_solver(degree_l, rk_method):

    y0_ = np.copy(y0[:4])
    solver = LiquidStaticIncompressible(
        radius_array, shear_modulus_array, bulk_modulus_array, density_array, gravity_array, frequency, degree_l,
        G_to_use, t_span, y0_, rk_method=rk_method
        )
    solver.solve()

    assert solver.solution_t.shape == (N,)
    assert solver.solution_y.shape == (4, N)
    assert solver.solution_t.dtype is np.float64
    assert solver.solution_t[0] == radius_array[0]
    assert solver.solution_t[-1] == radius_array[-1]
    assert np.all(np.isfinite(solver.solution_t))
    assert np.all(np.isfinite(solver.solution_y))

@pytest.mark.parametrize('rk_method', (1, 2))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_liquid_static_compressible_solver(degree_l, rk_method):

    y0_ = np.copy(y0[:4])
    solver = LiquidStaticCompressible(
        radius_array, shear_modulus_array, bulk_modulus_array, density_array, gravity_array, frequency, degree_l,
        G_to_use, t_span, y0_, rk_method=rk_method
        )
    solver.solve()

    assert solver.solution_t.shape == (N,)
    assert solver.solution_y.shape == (4, N)
    assert solver.solution_t.dtype is np.float64
    assert solver.solution_t[0] == radius_array[0]
    assert solver.solution_t[-1] == radius_array[-1]
    assert np.all(np.isfinite(solver.solution_t))
    assert np.all(np.isfinite(solver.solution_y))
