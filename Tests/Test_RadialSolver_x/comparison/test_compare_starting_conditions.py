"""Compare RadialSolver_x starting conditions against original RadialSolver."""
import pytest
import numpy as np

from TidalPy.rheology.models import Maxwell

from TidalPy.RadialSolver.shooting import find_num_shooting_solutions as find_num_old
from TidalPy.RadialSolver.starting.driver import find_starting_conditions as find_starting_old

from TidalPy.RadialSolver_x.derivatives.odes import find_num_shooting_solutions as find_num_new
from TidalPy.RadialSolver_x.starting.driver import find_starting_conditions as find_starting_new

frequency    = 0.1
radius       = 0.1
density      = 7000.
bulk_modulus  = 100.0e9
shear         = 50.0e9
viscosity     = 1.0e20
rheo_inst     = Maxwell()
complex_shear = rheo_inst(frequency, shear, viscosity)
G_to_use      = 6.67430e-11


@pytest.mark.parametrize('layer_type', (0, 1))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
def test_compare_num_shooting_solutions(layer_type, is_static, is_incompressible):
    """Number of shooting solutions should match."""
    old_n = find_num_old(layer_type, is_static, is_incompressible)
    new_n = find_num_new(layer_type, is_static, is_incompressible)
    assert old_n == new_n


@pytest.mark.parametrize('layer_type', (0, 1))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('use_kamata', (True, False))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_compare_starting_conditions(layer_type, is_static, is_incompressible, use_kamata, degree_l):
    """Starting conditions should match between old and new implementations."""

    is_solid = layer_type == 0
    # Skip known-unimplemented combos.
    if ((not use_kamata) and is_incompressible and (is_solid or ((not is_solid) and (not is_static)))) \
            or (use_kamata and is_static and is_incompressible and is_solid):
        pytest.skip('Combination not yet implemented.')

    num_sols = find_num_old(layer_type, is_static, is_incompressible)
    num_ys = num_sols * 2

    old_arr = np.nan * np.ones((num_sols, num_ys), dtype=np.complex128, order='C')
    new_arr = np.nan * np.ones((num_sols, num_ys), dtype=np.complex128, order='C')

    find_starting_old(
        layer_type, is_static, is_incompressible, use_kamata,
        frequency, radius, density, bulk_modulus, complex_shear,
        degree_l, G_to_use, old_arr
    )
    find_starting_new(
        layer_type, is_static, is_incompressible, use_kamata,
        frequency, radius, density, bulk_modulus, complex_shear,
        degree_l, G_to_use, new_arr
    )

    np.testing.assert_allclose(new_arr, old_arr, rtol=1e-10)
