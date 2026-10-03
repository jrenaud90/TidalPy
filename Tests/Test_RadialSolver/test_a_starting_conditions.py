"""Starting-condition driver: which combinations are implemented. Their values are checked against the frozen classic
solver's in test_comparison/test_compare_starting_conditions."""
import numpy as np
import pytest

from TidalPy.Rheology import Maxwell
from TidalPy.RadialSolver.derivatives.odes import find_num_shooting_solutions
from TidalPy.RadialSolver.starting.driver import find_starting_conditions

FREQUENCY = 0.1
RADIUS = 0.1
DENSITY = 7000.
BULK_MODULUS = 100.0e9
COMPLEX_SHEAR = Maxwell().calc_complex_modulus(50.0e9, 1.0e20, FREQUENCY)
G_TO_USE = 6.67430e-11


def _run_starting(
        layer_type,
        is_static,
        is_incompressible,
        use_kamata,
        degree_l):
    """Starting conditions as a (num_solutions, 2 * num_solutions) array, NaN-initialized."""
    num_sols = find_num_shooting_solutions(layer_type, is_static, is_incompressible)
    initial_condition_array = np.nan * np.ones((num_sols, num_sols * 2), dtype=np.complex128, order='C')
    find_starting_conditions(
        layer_type,
        is_static,
        is_incompressible,
        use_kamata,
        FREQUENCY,
        RADIUS,
        DENSITY,
        BULK_MODULUS,
        COMPLEX_SHEAR,
        degree_l,
        G_TO_USE,
        initial_condition_array
    )
    return initial_condition_array


@pytest.mark.parametrize('layer_type', (0, 1))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('use_kamata', (True, False))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_initial_condition_driver(
        layer_type,
        is_static,
        is_incompressible,
        use_kamata,
        degree_l):
    """Implemented combinations fill every value; the rest raise NotImplementedError."""
    is_solid = layer_type == 0
    if ((not use_kamata) and is_incompressible and (is_solid or ((not is_solid) and (not is_static)))) \
            or (use_kamata and is_static and is_incompressible and is_solid):
        with pytest.raises(NotImplementedError):
            _run_starting(
                layer_type,
                is_static,
                is_incompressible,
                use_kamata,
                degree_l)
    else:
        initial_condition_array = _run_starting(
            layer_type,
            is_static,
            is_incompressible,
            use_kamata,
            degree_l)
        assert not np.any(np.isnan(initial_condition_array))

