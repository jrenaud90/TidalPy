"""radial_solver input validation errors."""
import numpy as np
import pytest

from TidalPy.exceptions import UnknownModelError, ArgumentException
from TidalPy.RadialSolver.solver import radial_solver as _radial_solver


def radial_solver(*args, **kwargs):
    """Force Kamata starts and map validation errors onto TidalPy's exception types."""
    kwargs.setdefault('starting_method', 'kamata')
    try:
        return _radial_solver(*args, **kwargs)
    except NotImplementedError as exc:
        # An unimplemented start is an input error here, except for the matrix method's own limitations.
        if kwargs.get('love_method', 'radial_solver') == 'propagation_matrix':
            raise
        raise ArgumentException(str(exc)) from exc
    except ValueError as exc:
        msg = str(exc)
        msg_lower = msg.lower()
        if 'frequency' in msg_lower:
            raise
        if 'layer type' in msg_lower:
            raise UnknownModelError(msg) from exc
        raise ArgumentException(msg) from exc


def _radii(*segments):
    """Concatenated linspace segments, each (start, stop, count)."""
    return np.concatenate([np.linspace(start, stop, count, dtype=np.float64) for start, stop, count in segments])


def _swapped_radii():
    radius_array = _radii((0, 1000, 10))
    radius_array[[4, 5]] = radius_array[[5, 4]]
    return radius_array


def _negative_radii():
    radius_array = _radii((0, 1000, 10))
    radius_array[4] = -radius_array[4]
    return radius_array


ONE_LAYER = dict(layer_types=("solid",), is_static=(True,), is_incompressible=(True,), upper_radii=(1000.0,))
TWO_LAYERS_LIQUID_BASE = dict(
    layer_types=("liquid", "solid"), is_static=(False, False), is_incompressible=(False, False),
    upper_radii=(500.0, 1000.0))
PROPAGATION_MATRIX = dict(love_method='propagation_matrix')

# id: (radius array, layer description, extra kwargs, expected error, match)
INVALID_INPUT_CASES = {
    "density_size": (
        _radii((0, 1000, 10)), dict(ONE_LAYER, density_size=9), {}, ArgumentException,
        "match the radius array length"),
    "upper_radius_order": (
        _radii((0, 500, 5), (500, 1000, 5)),
        dict(layer_types=("solid", "liquid"), is_static=(True, False), is_incompressible=(True, False),
             upper_radii=(1000.0, 500.0)),
        {}, ArgumentException, None),
    "frequency_too_low": (
        _radii((0, 1000, 10)), dict(ONE_LAYER, is_static=(False,), is_incompressible=(False,)),
        dict(frequency=1e-20), ValueError, None),
    "frequency_too_high": (
        _radii((0, 1000, 10)), dict(ONE_LAYER, is_static=(False,), is_incompressible=(False,)),
        dict(frequency=1e10), ValueError, None),
    "radius_not_starting_at_zero": (_radii((10, 1000, 10)), ONE_LAYER, {}, ArgumentException, None),
    "radius_not_ascending": (_swapped_radii(), ONE_LAYER, {}, ArgumentException, None),
    "radius_negative": (_negative_radii(), ONE_LAYER, {}, ArgumentException, None),
    "unknown_layer_type": (
        _radii((0, 500, 5), (500, 1000, 5)),
        dict(TWO_LAYERS_LIQUID_BASE, layer_types=("solid", "unknown")),
        {}, UnknownModelError, None),
    "interface_missing_upper_base": (
        _radii((0, 500, 5), (510, 1000, 5)), TWO_LAYERS_LIQUID_BASE, {}, ArgumentException, None),
    "interface_missing_lower_top": (
        _radii((0, 490, 5), (500, 1000, 5)), TWO_LAYERS_LIQUID_BASE, {}, ArgumentException, None),
    "interface_missing_both": (
        _radii((0, 490, 5), (510, 1000, 5)), TWO_LAYERS_LIQUID_BASE, {}, ArgumentException, None),
    "too_few_slices": (
        _radii((0, 500, 2), (510, 1000, 3)), TWO_LAYERS_LIQUID_BASE, {}, ArgumentException, None),
    "matrix_too_many_layers": (
        _radii((0, 500, 5), (510, 1000, 5)),
        dict(layer_types=("solid", "solid"), is_static=(True, True), is_incompressible=(True, True),
             upper_radii=(500.0, 1000.0)),
        PROPAGATION_MATRIX, ArgumentException, None),
    "matrix_liquid_layer": (
        _radii((0, 1000, 10)), dict(ONE_LAYER, layer_types=("liquid",)), PROPAGATION_MATRIX, ArgumentException,
        None),
    "matrix_dynamic_layer": (
        _radii((0, 1000, 10)), dict(ONE_LAYER, is_static=(False,)), PROPAGATION_MATRIX, ArgumentException, None),
    "matrix_compressible_layer": (
        _radii((0, 1000, 10)), dict(ONE_LAYER, is_incompressible=(False,)), PROPAGATION_MATRIX, ArgumentException,
        None),
    # The starting radius must be below 90% of the planet radius.
    "starting_radius_too_high": (
        _radii((0, 1000, 10)), ONE_LAYER, dict(starting_radius=0.91 * 1000), ArgumentException, None),
}


@pytest.mark.parametrize(
    "radius_array, layers, extra_kwargs, expected_error, match",
    INVALID_INPUT_CASES.values(),
    ids=INVALID_INPUT_CASES.keys())
def test_invalid_input_raises(
        radius_array,
        layers,
        extra_kwargs,
        expected_error,
        match):
    """Malformed radius, density, layer, frequency, or method inputs raise the expected error."""
    num_slices = radius_array.size
    extra_kwargs = dict(extra_kwargs)
    frequency = extra_kwargs.pop('frequency', 1.0)
    with pytest.raises(expected_error, match=match):
        radial_solver(
            radius_array.copy(),
            np.linspace(1000, 2000, layers.get('density_size', num_slices), dtype=np.float64),
            np.zeros(num_slices, dtype=np.complex128),
            np.zeros(num_slices, dtype=np.complex128),
            frequency,
            3000.0,
            layers['layer_types'],
            layers['is_static'],
            layers['is_incompressible'],
            np.array(layers['upper_radii'], dtype=np.float64),
            raise_on_fail=True,
            **extra_kwargs)


@pytest.mark.parametrize("layer_arguments, extra_kwargs, match", (
    pytest.param(((), (), (), np.array([])), dict(love_method='propagation_matrix'), "At least one layer",
                 id="no_layers"),
    pytest.param((("solid", "solid"), (False, False), (False, False), np.array([3.0e6, 6.0e6])),
                 dict(eos_method_bylayer=("interpolate",)), "one method per layer", id="short_eos_method_list"),
    pytest.param((("solid",), (False,), (False,), np.array([5.0e6])), {}, "planet radius",
                 id="top_layer_below_profile_top"),
))
def test_layer_description_errors_raise_value_error(layer_arguments, extra_kwargs, match):
    """An empty layer list, a short EOS method list, or a top layer below the last radius raise ValueError."""
    with pytest.raises(ValueError, match=match):
        _radial_solver(
            np.linspace(0.0, 6.0e6, 50),
            np.full(50, 5000.0),
            np.full(50, 1.0e11 + 0j),
            np.full(50, 5.0e10 + 1.0e8j),
            1.0e-5,
            5000.0,
            *layer_arguments,
            **extra_kwargs)
