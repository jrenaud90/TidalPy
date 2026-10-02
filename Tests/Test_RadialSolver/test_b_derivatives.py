"""Shooting-solution counts for each layer configuration."""
import pytest

from TidalPy.RadialSolver.derivatives.odes import find_num_shooting_solutions


@pytest.mark.parametrize('layer_type', (0, 1))
@pytest.mark.parametrize('is_static', (0, 1))
@pytest.mark.parametrize('is_incompressible', (0, 1))
def test_find_num_shooting_solutions(layer_type, is_static, is_incompressible):
    """Solid layers have 3 solutions, static liquids 1, and dynamic liquids 2."""
    num_solutions = find_num_shooting_solutions(layer_type, is_static, is_incompressible)
    if layer_type == 0:
        expected = 3
    else:
        expected = 1 if is_static else 2
    assert num_solutions == expected

