"""The per-method starting-condition, interface, and surface wrappers.

``find_starting_conditions`` dispatches to one of these by (layer_type, is_static, is_incompressible,
use_kamata), and the rest of the suite only ever exercises the driver. These tests call each wrapper
directly and pin it against the driver's answer for the combination that selects it, so a wrapper cannot
drift from the path that uses it, and the surface and interface helpers get a shape-and-finiteness check.
"""
import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.RadialSolver_x.boundaries.boundaries import apply_surface_bc
from TidalPy.RadialSolver_x.interfaces.reversed import top_to_bottom_interface_bc
from TidalPy.RadialSolver_x.starting.common import takeuchi_phi_psi, z_calc
from TidalPy.RadialSolver_x.derivatives.odes import find_num_shooting_solutions
from TidalPy.RadialSolver_x.starting.driver import find_starting_conditions
from TidalPy.RadialSolver_x.starting.kamata import (
    kamata_liquid_dynamic_compressible,
    kamata_liquid_dynamic_incompressible,
    kamata_solid_dynamic_compressible,
    kamata_solid_dynamic_incompressible,
    kamata_solid_static_compressible,
)
from TidalPy.RadialSolver_x.starting.saito import saito_liquid_static_incompressible
from TidalPy.RadialSolver_x.starting.takeuchi import (
    takeuchi_liquid_dynamic_compressible,
    takeuchi_solid_dynamic_compressible,
    takeuchi_solid_static_compressible,
)
from TidalPy.Utilities_x.arrays.interp import partition_radius_by_layer

FREQUENCY = 2.0 * np.pi / 86400.0
RADIUS = 1.0e5
DENSITY = 7000.0
BULK = 200.0e9 + 0.0j
SHEAR = 100.0e9 + 1.0e8j
DEGREE_L = 2

# (wrapper, argument names it takes, layer_type, is_static, is_incompressible, use_kamata, num_solutions)
CASES = (
    (kamata_solid_dynamic_compressible, "full", 0, False, False, True, 3),
    (kamata_solid_static_compressible, "static_solid", 0, True, False, True, 3),
    (kamata_solid_dynamic_incompressible, "solid_incomp", 0, False, True, True, 3),
    (kamata_liquid_dynamic_compressible, "liquid", 1, False, False, True, 2),
    (kamata_liquid_dynamic_incompressible, "liquid_incomp", 1, False, True, True, 2),
    (takeuchi_solid_dynamic_compressible, "full", 0, False, False, False, 3),
    (takeuchi_solid_static_compressible, "static_solid", 0, True, False, False, 3),
    (takeuchi_liquid_dynamic_compressible, "liquid", 1, False, False, False, 2),
)


def _empty(layer_type, is_static, is_incompressible):
    """The driver sizes the view as (num_solutions, 2 * num_solutions); anything else it rejects."""
    num_sols = find_num_shooting_solutions(layer_type, is_static, is_incompressible)
    return np.zeros((num_sols, 2 * num_sols), dtype=np.complex128, order="C")


def _call(wrapper, kind, view):
    if kind == "full":
        wrapper(FREQUENCY, RADIUS, DENSITY, BULK, SHEAR, DEGREE_L, G, view)
    elif kind == "static_solid":
        wrapper(RADIUS, DENSITY, BULK, SHEAR, DEGREE_L, G, view)
    elif kind == "liquid":
        wrapper(FREQUENCY, RADIUS, DENSITY, BULK, DEGREE_L, G, view)
    elif kind == "solid_incomp":
        # An incompressible layer has no bulk modulus to take.
        wrapper(FREQUENCY, RADIUS, DENSITY, SHEAR, DEGREE_L, G, view)
    elif kind == "liquid_incomp":
        wrapper(FREQUENCY, RADIUS, DENSITY, DEGREE_L, G, view)
    else:
        raise AssertionError(kind)


@pytest.mark.parametrize("wrapper, kind, layer_type, is_static, is_incompressible, use_kamata, num_sols",
                         CASES, ids=[case[0].__name__ for case in CASES])
def test_wrapper_matches_the_driver(wrapper, kind, layer_type, is_static, is_incompressible,
                                    use_kamata, num_sols):
    """Calling a wrapper directly gives what the driver gives for the combination that selects it."""
    direct = _empty(layer_type, is_static, is_incompressible)
    _call(wrapper, kind, direct)

    through_driver = _empty(layer_type, is_static, is_incompressible)
    find_starting_conditions(
        layer_type, is_static, is_incompressible, use_kamata,
        FREQUENCY, RADIUS, DENSITY, BULK, SHEAR, DEGREE_L, G, through_driver)

    assert np.all(np.isfinite(direct.view(np.float64))), direct
    assert not np.all(direct == 0.0), "the wrapper wrote nothing"
    np.testing.assert_allclose(direct, through_driver, rtol=1.0e-12, atol=0.0)


def test_saito_matches_the_driver():
    """The static incompressible liquid wrapper takes only radius and degree."""
    direct = _empty(1, True, True)
    saito_liquid_static_incompressible(RADIUS, DEGREE_L, direct)

    through_driver = _empty(1, True, True)
    find_starting_conditions(
        1, True, True, False, FREQUENCY, RADIUS, DENSITY, BULK, SHEAR, DEGREE_L, G, through_driver)

    assert np.all(np.isfinite(direct.view(np.float64)))
    np.testing.assert_allclose(direct, through_driver, rtol=1.0e-12, atol=0.0)


# =====================================================================================================================
# The scalar helpers the Takeuchi conditions are built from
# =====================================================================================================================
@pytest.mark.parametrize("degree_l", (2, 3, 4))
def test_z_calc_is_finite_and_small_argument_limit(degree_l):
    """z -> x^2 / (2l + 3) as x^2 -> 0 (KMN15 Eq. B14)."""
    tiny = 1.0e-12 + 0.0j
    assert z_calc(tiny, degree_l) == pytest.approx(tiny / (2 * degree_l + 3), rel=1.0e-6)
    value = z_calc(0.5 + 0.25j, degree_l)
    assert np.isfinite(value.real) and np.isfinite(value.imag)


@pytest.mark.parametrize("degree_l", (2, 3))
def test_takeuchi_phi_psi_returns_three_finite_values(degree_l):
    """phi, phi_{l+1}, psi (TS72 Eq. 103); all three are finite and phi -> 1 at small argument."""
    result = takeuchi_phi_psi(1.0e-14 + 0.0j, degree_l)
    assert len(result) == 3
    assert all(np.isfinite(complex(value)) for value in result)
    assert complex(result[0]) == pytest.approx(1.0 + 0.0j, abs=1.0e-9)


# =====================================================================================================================
# Surface and interface helpers
# =====================================================================================================================
def test_apply_surface_bc_writes_a_finite_constant_vector():
    """A three-solution solid surface solve returns finite constants for a tidal boundary condition."""
    constants = np.zeros(3, dtype=np.complex128)
    bc = np.asarray([0.0, 0.0, (2.0 * DEGREE_L + 1.0) / RADIUS], dtype=np.float64)
    uppermost = np.zeros((3, 6), dtype=np.complex128, order="C")
    for index in range(3):
        uppermost[index, index] = 1.0 + 0.0j
        uppermost[index, 5] = 1.0e-6 * (index + 1)

    apply_surface_bc(constants, bc, uppermost, 1.62, G, 0, 0, False, False)

    assert np.all(np.isfinite(constants.view(np.float64))), constants


def test_top_to_bottom_interface_bc_runs_for_a_solid_under_solid():
    """The reversed interface helper fills the lower layer's constants from the layer above."""
    constants = np.zeros(3, dtype=np.complex128)
    above = np.asarray([1.0 + 0.0j, 0.5 + 0.0j, 0.25 + 0.0j], dtype=np.complex128)
    uppermost = np.zeros((3, 6), dtype=np.complex128, order="C")
    for index in range(3):
        uppermost[index, index] = 1.0 + 0.0j

    top_to_bottom_interface_bc(
        constants, above, uppermost, 1.62, 1.60, 3500.0, 3300.0, 0, 0, False, False, False, False)

    assert np.all(np.isfinite(constants.view(np.float64))), constants


# =====================================================================================================================
# The shared radius partition
# =====================================================================================================================
def test_partition_radius_by_layer_splits_at_the_repeated_interface():
    """Each interface radius appears twice; the first copy belongs to the layer below."""
    radius = np.asarray([0.0, 500.0, 1000.0, 1000.0, 1500.0, 2000.0], dtype=np.float64)
    uppers = np.asarray([1000.0, 2000.0], dtype=np.float64)

    first_index, num_slices = partition_radius_by_layer(radius, uppers)

    assert list(first_index) == [0, 3]
    assert list(num_slices) == [3, 3]
    assert sum(num_slices) == radius.size
