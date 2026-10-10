"""Every function that takes ``G_to_use`` reads ``None`` as the configured gravitational constant."""
import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.RadialSolver.boundaries.boundaries import apply_surface_bc
from TidalPy.RadialSolver.derivatives.odes import find_num_shooting_solutions
from TidalPy.RadialSolver.interfaces.interfaces import solve_upper_y_at_interface
from TidalPy.RadialSolver.matrix_types.solid_matrix import fundamental_matrix
from TidalPy.RadialSolver.starting.driver import find_starting_conditions
from TidalPy.RadialSolver.starting.kamata import (
    kamata_liquid_dynamic_compressible,
    kamata_liquid_dynamic_incompressible,
    kamata_solid_dynamic_compressible,
    kamata_solid_dynamic_incompressible,
    kamata_solid_static_compressible,
)
from TidalPy.RadialSolver.starting.power_series import (
    power_series_liquid_dynamic_compressible,
    power_series_solid_dynamic_compressible,
    power_series_solid_static_incompressible,
)
from TidalPy.RadialSolver.starting.takeuchi import (
    takeuchi_liquid_dynamic_compressible,
    takeuchi_solid_dynamic_compressible,
    takeuchi_solid_static_compressible,
)
from TidalPy.Tides.classes.collapse import collapse_global_tides
from TidalPy.Tides.potential import global_potential, tidal_potential_3d_modes
from TidalPy.Utilities.conversions import orbital_motion2semi_a, semi_a2orbital_motion

FREQUENCY = 2.0 * np.pi / 86400.0
RADIUS = 1.0e5
DENSITY = 7000.0
BULK = 200.0e9 + 0.0j
SHEAR = 100.0e9 + 1.0e8j
DEGREE_L = 2

# An Io-like orbit: radius [m], then the calc_tides orbital state.
ORBIT = (1.8216e6, 4.11e-5, 4.11e-5, 0.0041, 0.0, 4.217e8, 1.898e27)


def _empty(layer_type, is_static, is_incompressible):
    num_sols = find_num_shooting_solutions(layer_type, is_static, is_incompressible)
    return np.zeros((num_sols, 2 * num_sols), dtype=np.complex128, order="C")


STARTING = (
    (kamata_solid_dynamic_compressible, (FREQUENCY, RADIUS, DENSITY, BULK, SHEAR, DEGREE_L), (0, False, False)),
    (kamata_solid_static_compressible, (RADIUS, DENSITY, BULK, SHEAR, DEGREE_L), (0, True, False)),
    (kamata_solid_dynamic_incompressible, (FREQUENCY, RADIUS, DENSITY, SHEAR, DEGREE_L), (0, False, True)),
    (kamata_liquid_dynamic_compressible, (FREQUENCY, RADIUS, DENSITY, BULK, DEGREE_L), (1, False, False)),
    (kamata_liquid_dynamic_incompressible, (FREQUENCY, RADIUS, DENSITY, DEGREE_L), (1, False, True)),
    (takeuchi_solid_dynamic_compressible, (FREQUENCY, RADIUS, DENSITY, BULK, SHEAR, DEGREE_L), (0, False, False)),
    (takeuchi_solid_static_compressible, (RADIUS, DENSITY, BULK, SHEAR, DEGREE_L), (0, True, False)),
    (takeuchi_liquid_dynamic_compressible, (FREQUENCY, RADIUS, DENSITY, BULK, DEGREE_L), (1, False, False)),
    (power_series_solid_dynamic_compressible, (FREQUENCY, RADIUS, DENSITY, BULK, SHEAR, DEGREE_L), (0, False, False)),
    (power_series_solid_static_incompressible, (RADIUS, DENSITY, SHEAR, DEGREE_L), (0, True, True)),
    (power_series_liquid_dynamic_compressible, (FREQUENCY, RADIUS, DENSITY, BULK, DEGREE_L), (1, False, False)),
)


@pytest.mark.parametrize("wrapper, arguments, layer", STARTING, ids=[case[0].__name__ for case in STARTING])
def test_starting_conditions_take_none_for_the_configured_g(wrapper, arguments, layer):
    with_none = _empty(*layer)
    wrapper(*arguments, None, with_none)
    with_g = _empty(*layer)
    wrapper(*arguments, G, with_g)
    np.testing.assert_array_equal(with_none, with_g)


def test_starting_driver_takes_none_for_the_configured_g():
    results = []
    for G_to_use in (None, G):
        out = _empty(0, False, False)
        find_starting_conditions(
            0, False, False, "kamata", FREQUENCY, RADIUS, DENSITY, BULK, SHEAR, DEGREE_L, G_to_use, out)
        results.append(out)
    np.testing.assert_array_equal(*results)


def test_surface_conditions_take_none_for_the_configured_g():
    boundary = np.asarray([0.0, 0.0, (2.0 * DEGREE_L + 1.0) / RADIUS], dtype=np.float64)
    uppermost = np.zeros((3, 6), dtype=np.complex128, order="C")
    for index in range(3):
        uppermost[index, index] = 1.0 + 0.0j
        uppermost[index, 5] = 1.0e-6 * (index + 1)
    results = []
    for G_to_use in (None, G):
        constants = np.zeros(3, dtype=np.complex128)
        apply_surface_bc(constants, boundary, uppermost, 1.62, G_to_use, 0, 0, False, False)
        results.append(constants)
    np.testing.assert_array_equal(*results)


def test_interface_conditions_default_to_the_configured_g():
    lower_y = np.random.default_rng(2024).normal(size=(3, 6)) + 0.0j
    results = []
    for G_to_use in ("default", None, G):
        upper_y = np.full((3, 6), np.nan, dtype=np.complex128)
        if G_to_use == "default":
            solve_upper_y_at_interface(lower_y, upper_y, 0, False, 0, False, 1.62, 3000.0)
        else:
            solve_upper_y_at_interface(lower_y, upper_y, 0, False, 0, False, 1.62, 3000.0, G_to_use)
        results.append(upper_y)
    for result in results[1:]:
        np.testing.assert_array_equal(result, results[0])


def test_fundamental_matrix_default_is_the_configured_g():
    radius = np.linspace(1.0e5, 1.0e6, 5)
    density = np.full(5, 3000.0)
    gravity = np.linspace(0.5, 1.5, 5)
    shear = np.full(5, 5.0e10 + 1.0e8j)
    default = fundamental_matrix(radius, density, gravity, shear)
    explicit = fundamental_matrix(radius, density, gravity, shear, 2, G)
    for default_part, explicit_part in zip(default, explicit):
        np.testing.assert_array_equal(default_part, explicit_part)


def test_tidal_potentials_take_none_for_the_configured_g():
    default = global_potential(*ORBIT)[3]
    explicit = global_potential(*ORBIT, G)[3]
    assert default == explicit
    modes_none = tidal_potential_3d_modes(*ORBIT, None, 0.5, 0.3)
    modes_g = tidal_potential_3d_modes(*ORBIT, G, 0.5, 0.3)
    for none_part, g_part in zip(modes_none, modes_g):
        np.testing.assert_array_equal(none_part, g_part)
    tide_config = {"fixed_k": [0.3], "fixed_q": [100.0]}
    collapsed = [collapse_global_tides(*ORBIT, G_to_use, "fixed_q", tide_config) for G_to_use in (None, G)]
    assert collapsed[0] == collapsed[1]


def test_kepler_conversions_take_none_for_the_configured_g():
    assert semi_a2orbital_motion(4.217e8, 1.898e27) == semi_a2orbital_motion(4.217e8, 1.898e27, 0.0, G)
    assert orbital_motion2semi_a(4.11e-5, 1.898e27, G_to_use=None) == orbital_motion2semi_a(4.11e-5, 1.898e27, 0.0, G)
