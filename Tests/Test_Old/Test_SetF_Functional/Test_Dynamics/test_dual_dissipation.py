""" Tests for TidalPy.dynamics.dual_dissipation """
import numpy as np

import TidalPy


from TidalPy.dynamics.dual_dissipation import (eccentricity_derivative, semi_major_axis_derivative,
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

host_dRdM = 1.0e-8
host_dRdw = 5.0e-8

target_radius = 1.e6
target_mass = 1.e24
target_gravity = 10.
target_density = 5000.
target_moi = 1.e5
target_spin_period = 1.2
target_obliquity = 0.1

target_dRdM = 1.0e-5
target_dRdw = 5.0e-5


def test_dual_disp_floats():
    """ Test the three dynamic functions for float inputs """

    result = eccentricity_derivative(
        semi_major_axis, orbital_freq, eccentricity,
        host_mass, host_dRdM, host_dRdw,
        target_mass, target_dRdM, target_dRdw
        )
    assert type(result) in [float, np.float64]

    result = semi_major_axis_derivative(
        semi_major_axis, orbital_freq,
        host_mass, host_dRdM,
        target_mass, target_dRdM
        )
    assert type(result) in [float, np.float64]

    result = semia_eccen_derivatives(
        semi_major_axis, orbital_freq, eccentricity,
        host_mass, host_dRdM, host_dRdw,
        target_mass, target_dRdM, target_dRdw
        )
    assert type(result) == tuple
    assert type(result[0]) in [float, np.float64]
    assert type(result[1]) in [float, np.float64]


def test_dual_disp_arrays():
    """ Test the three dynamic functions for float inputs """

    sa = semi_major_axis * np.ones(10, dtype=np.float64)
    freq = orbital_freq * np.ones(10, dtype=np.float64)
    hdRdM = host_dRdM * np.ones(10, dtype=np.float64)
    hdRdw = host_dRdw * np.ones(10, dtype=np.float64)
    tdRdM = target_dRdM * np.ones(10, dtype=np.float64)
    tdRdw = target_dRdw * np.ones(10, dtype=np.float64)

    result = eccentricity_derivative(
        sa, freq, eccentricity,
        host_mass, hdRdM, hdRdw,
        target_mass, tdRdM, tdRdw
        )
    assert type(result) == np.ndarray
    assert result.dtype in [float, np.float64]
    assert result.shape == (10,)

    result = semi_major_axis_derivative(
        sa, freq,
        host_mass, hdRdM,
        target_mass, tdRdM
        )
    assert type(result) == np.ndarray
    assert result.dtype in [float, np.float64]
    assert result.shape == (10,)

    result = semia_eccen_derivatives(
        sa, freq, eccentricity,
        host_mass, hdRdM, hdRdw,
        target_mass, tdRdM, tdRdw
        )
    assert type(result) == tuple
    assert type(result[0]) == np.ndarray
    assert result[0].dtype in [float, np.float64]
    assert result[0].shape == (10,)
    assert type(result[1]) == np.ndarray
    assert result[1].dtype in [float, np.float64]
    assert result[1].shape == (10,)
