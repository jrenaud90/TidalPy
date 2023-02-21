""" Test the `TidalPy.utilities.performance.array` module. """

import numpy as np
import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.utilities.performance.array import interp

evenly_spaced_array = np.linspace(0., 100., 10000, dtype=np.float64)
unevenly_spaced_array = np.concatenate(
        (np.linspace(0., 10., 1000, dtype=np.float64),
        np.linspace(11., 40., 1000, dtype=np.float64),
        np.linspace(41., 50., 1000, dtype=np.float64),
        np.linspace(51., 100., 1000, dtype=np.float64))
        )


@pytest.mark.parametrize('value_to_check', (24., 56.2, 87.))
@pytest.mark.parametrize('array_to_use', (evenly_spaced_array, unevenly_spaced_array))
def test_interp(array_to_use, value_to_check):
    """ Test `utilities.performance.array` interp functionality. """

    # Find y-values
    y = np.sqrt(array_to_use)

    # Get Numpy version
    value_np = np.interp(value_to_check, array_to_use, y)

    # Get TidalPy version
    value_tpy = interp(value_to_check, array_to_use, y)

    assert np.isclose(value_np, value_tpy)
