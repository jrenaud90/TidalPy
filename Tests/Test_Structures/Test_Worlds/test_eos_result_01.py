"""BaseWorld.solve_eos returns an EOSResult: a dict in every respect whose repr is a short summary."""
import copy
import pickle

import numpy as np
import pytest

from TidalPy.Structures import build_world
from TidalPy.Structures.worlds import EOSResult


@pytest.fixture(scope="module")
def result():
    return build_world("io").solve_eos()


def test_the_result_is_a_dict(result):
    assert type(result) is EOSResult
    assert isinstance(result, dict)
    assert result["success"]
    assert result["planet_mass"] > 0.0


def test_the_result_equals_a_plain_dict_of_its_entries(result):
    plain = dict(result)
    assert type(plain) is dict
    assert plain.keys() == result.keys()
    assert all(plain[key] is result[key] for key in result)
    # Equality compares entries, so it holds against a dict without array entries, whose truth would be ambiguous.
    scalars = {key: value for key, value in result.items() if not isinstance(value, (np.ndarray, list))}
    assert EOSResult(scalars) == scalars


@pytest.mark.parametrize("make_copy", [copy.copy, copy.deepcopy, lambda value: pickle.loads(pickle.dumps(value))],
                         ids=["copy", "deepcopy", "pickle"])
def test_copies_keep_the_type_and_entries(result, make_copy):
    duplicate = make_copy(result)
    assert type(duplicate) is EOSResult
    assert duplicate.keys() == result.keys()
    np.testing.assert_array_equal(duplicate["radius"], result["radius"])
    assert duplicate["planet_mass"] == result["planet_mass"]


def test_the_repr_is_a_short_summary(result):
    text = repr(result)
    assert text.startswith("EOSResult(success=True, iterations=")
    assert f"planet_mass = {result['planet_mass']:.6g} kg" in text
    assert f"radius = {result['radius'][-1]:.6g} m" in text
    assert f"central_pressure = {result['central_pressure']:.6g} Pa" in text
    assert "layer_in_thermal_network" in text
    # A few lines, not the profile arrays.
    assert len(text.splitlines()) < 15
    assert "array(" not in text


def test_the_repr_of_an_emptied_result():
    assert "radius = nan m" in repr(EOSResult())
