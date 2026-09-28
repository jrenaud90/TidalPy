"""Shooting-solution counts and ODE names for each layer configuration."""
import pytest

from TidalPy.RadialSolver.derivatives.odes import find_num_shooting_solutions, find_layer_diffeq_name


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


@pytest.mark.parametrize('layer_type', (0, 1))
@pytest.mark.parametrize('is_static', (0, 1))
@pytest.mark.parametrize('is_incompressible', (0, 1))
def test_find_layer_diffeq_name(layer_type, is_static, is_incompressible):
    """The ODE name is '<solid|liquid>_<static|dynamic>_<incompressible|compressible>'."""
    name = find_layer_diffeq_name(layer_type, is_static, is_incompressible)
    assert isinstance(name, str)

    layer_str = 'solid' if layer_type == 0 else 'liquid'
    static_str = 'static' if is_static == 1 else 'dynamic'
    incomp_str = 'incompressible' if is_incompressible == 1 else 'compressible'
    assert name == f'{layer_str}_{static_str}_{incomp_str}'
