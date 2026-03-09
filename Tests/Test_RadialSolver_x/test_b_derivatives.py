import pytest

from TidalPy.RadialSolver_x.derivatives.odes import find_num_shooting_solutions, find_layer_diffeq_name


@pytest.mark.parametrize('layer_type', (0, 1))
@pytest.mark.parametrize('is_static', (0, 1))
@pytest.mark.parametrize('is_incompressible', (0, 1))
def test_find_num_shooting_solutions(layer_type, is_static, is_incompressible):
    """Test the number of shooting solutions returned for each layer configuration."""
    n = find_num_shooting_solutions(layer_type, is_static, is_incompressible)

    if layer_type == 0:
        # Solid layers always have 3 independent solutions.
        assert n == 3
    else:
        # Liquid layers
        if is_static:
            # Static liquid: 1 solution.
            assert n == 1
        else:
            # Dynamic liquid: 2 solutions.
            assert n == 2


@pytest.mark.parametrize('layer_type', (0, 1))
@pytest.mark.parametrize('is_static', (0, 1))
@pytest.mark.parametrize('is_incompressible', (0, 1))
def test_find_layer_diffeq_name(layer_type, is_static, is_incompressible):
    """Test that find_layer_diffeq_name returns a valid name string."""
    name = find_layer_diffeq_name(layer_type, is_static, is_incompressible)
    assert isinstance(name, str)

    layer_str = 'solid' if layer_type == 0 else 'liquid'
    static_str = 'static' if is_static == 1 else 'dynamic'
    incomp_str = 'incompressible' if is_incompressible == 1 else 'compressible'
    expected = f'{layer_str}_{static_str}_{incomp_str}'
    assert name == expected
