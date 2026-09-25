"""Vectorized physics-model calls accept read-only arrays, and a model's name cannot be changed under it."""
import numpy as np
import pytest

from TidalPy.rheology_x.rheology import Maxwell, maxwell


def _read_only(values):
    array = np.array(values, dtype=np.float64)
    array.setflags(write=False)
    return array


def test_read_only_inputs_are_accepted():
    modulus = _read_only([5.0e10, 6.0e10])
    viscosity = _read_only([1.0e19, 1.0e20])
    frequency = _read_only([1.0e-5, 2.0e-5])
    direct = maxwell(modulus, viscosity, frequency)
    method = Maxwell().calc_complex_modulus_vectorize_all(modulus, viscosity, frequency)
    assert np.allclose(direct, method)


def test_model_name_is_read_only():
    model = Maxwell()
    with pytest.raises(AttributeError):
        model.model_name = "andrade"
    assert model.model_name == "maxwell"
