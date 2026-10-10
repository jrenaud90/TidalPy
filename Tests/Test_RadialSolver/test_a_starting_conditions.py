"""Starting-condition driver: which combinations are implemented. Their values are checked against the frozen classic
solver's in test_comparison/test_compare_starting_conditions."""
import numpy as np
import pytest

from TidalPy.constants import STARTING_METHOD_NAMES, starting_method_from_name
from TidalPy.Rheology import Maxwell
from TidalPy.RadialSolver.derivatives.odes import find_num_shooting_solutions
from TidalPy.RadialSolver.starting.driver import find_starting_conditions

from starting_methods import STARTING_METHODS

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
        starting_method,
        degree_l):
    """Starting conditions as a (num_solutions, 2 * num_solutions) array, NaN-initialized."""
    num_sols = find_num_shooting_solutions(layer_type, is_static, is_incompressible)
    initial_condition_array = np.nan * np.ones((num_sols, num_sols * 2), dtype=np.complex128, order='C')
    find_starting_conditions(
        layer_type,
        is_static,
        is_incompressible,
        starting_method,
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
@pytest.mark.parametrize('starting_method', STARTING_METHODS)
@pytest.mark.parametrize('degree_l', (2, 3))
def test_initial_condition_driver(
        layer_type,
        is_static,
        is_incompressible,
        starting_method,
        degree_l):
    """Every method fills every value for every layer type."""
    initial_condition_array = _run_starting(
        layer_type,
        is_static,
        is_incompressible,
        starting_method,
        degree_l)
    assert not np.any(np.isnan(initial_condition_array))


@pytest.mark.parametrize('name, canonical', (
    ('takeuchi', 'takeuchi'), ('TS', 'takeuchi'), ('Kamata', 'kamata'), ('power_series', 'power_series'),
    ('PowerSeries', 'power_series'), ('ps', 'power_series'), ('martens', 'power_series'), ('UNITY', 'unity')))
def test_starting_method_names_and_aliases(name, canonical):
    """Every alias resolves, in any case, to its method's canonical name."""
    assert STARTING_METHOD_NAMES[starting_method_from_name(name)] == canonical


def test_unknown_starting_method_is_refused():
    """An unknown name raises ValueError listing the accepted names."""
    with pytest.raises(ValueError, match="Unsupported starting method"):
        _run_starting(0, False, False, "bessel", 2)
