""" Tests for TidalPy.dynamics.single_dissipation """
import numpy as np

import TidalPy


from TidalPy.dynamics.single_dissipation import (spin_rate_derivative, eccentricity_derivative,
                                                 semi_major_axis_derivative,
                                                 semia_eccen_derivatives)

semi_major_axis = 1.0e11
orbital_freq = 1.0e-5
eccentricity = 0.1

host_mass = 1.e27
host_radius = 1.e6
host_gravity = 10.
host_density = 5000.
host_spin_period = 1.2
host_moi = 1.e5
host_obliquity = 0.1

target_radius = 1.e6
target_mass = 1.e24
target_gravity = 10.
target_density = 5000.
target_moi = 1.e5
target_spin_period = 1.2
target_obliquity = 0.1

target_dRdM = 1.0e-5
target_dRdw = 5.0e-5
target_dRdO = 10.0e-5


def test_single_disp_floats():
    """ Test the three dynamic functions for float inputs """

    result = spin_rate_derivative(target_dRdO, target_moi, host_mass)
    assert type(result) in [float, np.float64]

    result = eccentricity_derivative(
        semi_major_axis, orbital_freq, eccentricity,
        target_mass, target_dRdM, target_dRdw, host_mass
        )
    assert type(result) in [float, np.float64]

    result = semi_major_axis_derivative(
        semi_major_axis, orbital_freq,
        target_mass, target_dRdM, host_mass
        )
    assert type(result) in [float, np.float64]

    result = semia_eccen_derivatives(
        semi_major_axis, orbital_freq, eccentricity,
        target_mass, target_dRdM, target_dRdw, host_mass
        )
    assert type(result) == tuple
    assert type(result[0]) in [float, np.float64]
    assert type(result[1]) in [float, np.float64]


def test_single_disp_arrays():
    """ Test the three dynamic functions for float inputs """

    sa = semi_major_axis * np.ones(10, dtype=np.float64)
    freq = orbital_freq * np.ones(10, dtype=np.float64)
    tdRdM = target_dRdM * np.ones(10, dtype=np.float64)
    tdRdw = target_dRdw * np.ones(10, dtype=np.float64)
    tdRdO = target_dRdO * np.ones(10, dtype=np.float64)

    result = spin_rate_derivative(tdRdO, target_moi, host_mass)
    assert type(result) == np.ndarray
    assert result.dtype in [float, np.float64]
    assert result.shape == (10,)

    result = eccentricity_derivative(
        sa, freq, eccentricity,
        target_mass, tdRdM, tdRdw, host_mass
        )
    assert type(result) == np.ndarray
    assert result.dtype in [float, np.float64]
    assert result.shape == (10,)

    result = semi_major_axis_derivative(
        sa, freq,
        target_mass, tdRdM, host_mass
        )
    assert type(result) == np.ndarray
    assert result.dtype in [float, np.float64]
    assert result.shape == (10,)

    result = semia_eccen_derivatives(
        sa, freq, eccentricity,
        target_mass, tdRdM, tdRdw, host_mass
        )
    assert type(result) == tuple
    assert type(result[0]) == np.ndarray
    assert result[0].dtype in [float, np.float64]
    assert result[0].shape == (10,)
    assert type(result[1]) == np.ndarray
    assert result[1].dtype in [float, np.float64]
    assert result[1].shape == (10,)
