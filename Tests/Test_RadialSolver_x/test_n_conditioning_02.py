"""Rank, starting-condition, and starting-radius checks of the shooting method's surface solve.

- ``surface_solve_rcond`` is the reciprocal condition number of the surface boundary-condition system after
  equilibration. A healthy solve reports about 1e-3 to 1e-1; a system singular to working precision (below
  ``[numerical] minimum_surface_rcond``) fails with error code -13 instead of returning arbitrary constants, which
  ``surface_solve_amplification`` alone cannot flag (it reads about 1 for a singular system).
- The Kamata et al. (2015) starting conditions for a dynamic incompressible solid use a first solution that stays
  independent at long forcing periods, so a homogeneous body matches the closed-form incompressible Love number.
- A manual starting radius outside (0, max_start_radius_fraction * R] is refused on the world path as well as the
  standalone one.

The conditioning warning is logged through the C++ (spdlog) logger, which these tests route to a temporary file.
"""
import math

import numpy as np
import pytest

import TidalPy
from TidalPy.constants import G, update_constants_x
from TidalPy.exceptions import SolutionFailedError
from TidalPy.RadialSolver_x import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.RadialSolver_x.rs_solution import check_surface_solve_conditioning
from TidalPy.RadialSolver_x.starting.common import z_calc
from TidalPy.RadialSolver_x.starting.kamata import kamata_solid_dynamic_incompressible
from TidalPy.rheology_x import Elastic, Maxwell
from TidalPy.structures_x import build_world
from TidalPy.Utilities_x.logging_x.logger import flush_logger, init_logger
from TidalPy.initialize import build_logging_x_config

WARNING_TEXT = "poorly conditioned"

# Homogeneous Moon-like body for the closed-form comparison.
MOON_RADIUS = 1737.4e3    # [m]
MOON_DENSITY = 3344.0     # [kg m-3]
MOON_SHEAR = 6.0e10       # [Pa]
SECONDS_PER_DAY = 86400.0

IO_FREQUENCY = 4.11e-5    # [rad s-1]


# =====================================================================================================================
# Fixtures and helpers
# =====================================================================================================================
@pytest.fixture
def spdlog_text(tmp_path):
    """Route the C++ logger to a temporary file for the test and hand back a reader for its text."""
    log_path = tmp_path / "tidalpy_x.log"
    init_logger({"console_level": "off", "file_level": "warning", "log_to_file": True,
                 "log_file_path": str(log_path)})

    def read():
        flush_logger()
        return log_path.read_text(encoding="utf-8") if log_path.exists() else ""

    yield read
    init_logger(build_logging_x_config())


@pytest.fixture
def numerical_setter():
    """Set any ``[numerical]`` key for one test and restore every change afterwards."""
    originals = {}

    def set_value(key, value):
        originals.setdefault(key, TidalPy.config_x["numerical"][key])
        TidalPy.config_x["numerical"][key] = value
        update_constants_x()

    yield set_value
    for key, value in originals.items():
        TidalPy.config_x["numerical"][key] = value
    update_constants_x()


@pytest.fixture(scope="module")
def io():
    world = build_world("io")
    world.solve_eos()
    return world


def _static_one_layer_inputs():
    """A homogeneous static compressible solid (the same body the degree-1 validation test uses)."""
    return build_rs_input_homogeneous_layers(
        1.0e6, IO_FREQUENCY, (3000.0,), (1.0e11,), (5.0e10,), (1.0e30,), (1.0e20,),
        ("solid",), (True,), (False,), Maxwell(), Elastic(),
        radius_fraction_tuple=(1.0,), slice_per_layer=50)


def _moon_inputs(period_days):
    """A homogeneous elastic, dynamic, incompressible Moon forced at the given period."""
    frequency = 2.0 * np.pi / (period_days * SECONDS_PER_DAY)
    return build_rs_input_homogeneous_layers(
        MOON_RADIUS, frequency, (MOON_DENSITY,), (1.0e20,), (MOON_SHEAR,), (1.0e30,), (1.0e30,),
        ("solid",), (False,), (True,), Elastic(), Elastic(),
        radius_fraction_tuple=(1.0,), slice_per_layer=20)


def _closed_form_incompressible_k(degree_l):
    """Static Love number k_l of a homogeneous incompressible elastic sphere; k2 = (3/2) / (1 + 19 mu / (2 rho g R))."""
    gravity = (4.0 / 3.0) * np.pi * G * MOON_DENSITY * MOON_RADIUS
    rigidity = (2 * degree_l**2 + 4 * degree_l + 3) * MOON_SHEAR / (degree_l * MOON_DENSITY * gravity * MOON_RADIUS)
    return (3.0 / (2.0 * (degree_l - 1))) / (1.0 + rigidity)


# =====================================================================================================================
# The rank measure and its threshold
# =====================================================================================================================
def test_minimum_surface_rcond_is_wired_through():
    from TidalPy import constants
    assert TidalPy.config_x["numerical"]["minimum_surface_rcond"] == 1.0e-14
    assert constants.minimum_surface_rcond == TidalPy.config_x["numerical"]["minimum_surface_rcond"]


@pytest.mark.parametrize("use_kamata", (False, True))
@pytest.mark.parametrize("degree_l", (2, 3))
def test_regular_solve_reports_a_healthy_rcond(use_kamata, degree_l):
    """A regular solve is well conditioned: finite rcond, far above the threshold, at most 1."""
    solution = radial_solver(*_static_one_layer_inputs(), degree_l=degree_l, use_kamata=use_kamata)
    assert solution.success, solution.message
    rcond = solution.surface_solve_rcond
    assert np.isfinite(rcond)
    assert TidalPy.config_x["numerical"]["minimum_surface_rcond"] < 1.0e-3 < rcond <= 1.0
    diagnostics = solution.print_diagnostics(print_diagnostics=False, log_diagnostics=False)
    assert "rcond" in diagnostics


def test_singular_surface_system_fails():
    """A degree-1 load on a static body: a rigid translation meets every surface condition, so the constants are
    undetermined. The amplification reads about 1 here; the rank measure is what catches it."""
    solution = radial_solver(*_static_one_layer_inputs(), degree_l=1, solve_for=("loading",))
    assert not solution.success
    assert solution.error_code == -13
    assert "singular" in solution.message
    assert solution.surface_solve_rcond < TidalPy.config_x["numerical"]["minimum_surface_rcond"]
    assert np.isnan(solution.k)

    with pytest.raises(SolutionFailedError, match="singular"):
        radial_solver(*_static_one_layer_inputs(), degree_l=1, solve_for=("loading",), raise_on_fail=True)


def test_the_threshold_is_read_from_the_config(numerical_setter):
    """Raising minimum_surface_rcond above a healthy solve's value fails it, so C++ reads the configured value."""
    numerical_setter("minimum_surface_rcond", 0.5)
    solution = radial_solver(*_static_one_layer_inputs(), degree_l=2)
    assert not solution.success
    assert solution.error_code == -13
    assert "5.000e-01" in solution.message


def test_rcond_warning_and_backward_compatible_call(spdlog_text):
    """The warning fires when rcond is below the integration rtol; the two-argument form still works."""
    assert not check_surface_solve_conditioning(1.0, 1.0e-6)
    assert not check_surface_solve_conditioning(1.0, 1.0e-6, 0.1)
    assert WARNING_TEXT not in spdlog_text()
    assert check_surface_solve_conditioning(1.0, 1.0e-6, 1.0e-8)
    text = spdlog_text()
    assert WARNING_TEXT in text
    assert "reciprocal condition number" in text


# =====================================================================================================================
# Kamata starting conditions for a dynamic incompressible solid
# =====================================================================================================================
@pytest.mark.parametrize("degree_l", (2, 3))
def test_first_kamata_solution_is_the_published_difference(degree_l):
    """Solution 1 equals (KMN15 sn1 - sn2) gamma / w^2; the published sn1 (Eqs. B17-B28) is rebuilt here."""
    frequency, radius, density, shear = 3.0e-3, 5.0e5, 3000.0, 4.0e10 + 1.0e8j
    gamma = 4.0 * np.pi * G * density / 3.0
    starting = np.zeros((3, 6), dtype=np.complex128, order="C")
    kamata_solid_dynamic_incompressible(frequency, radius, density, shear, degree_l, G, starting)

    w2 = frequency**2
    llp1 = degree_l * (degree_l + 1.0)
    z_value = complex(z_calc(w2 * density / shear * radius**2 + 0.0j, degree_l))
    published_first = np.asarray((
        0.0,
        llp1 * (-density * gamma + 2.0 * shear * z_value / radius**2),
        z_value / radius,
        shear * (w2 * density / shear - 2.0 * z_value / radius**2),
        (degree_l + 1.0) * (degree_l * gamma - w2),
        (2.0 * degree_l + 1.0) * (degree_l + 1.0) * (degree_l * gamma - w2) / radius), dtype=np.complex128)
    expected_first = (published_first - starting[1]) * gamma / w2
    np.testing.assert_allclose(starting[0], expected_first, rtol=1.0e-9, atol=1.0e-12 * np.max(np.abs(expected_first)))


@pytest.mark.parametrize("frequency", (1.0e-14, 1.0e-9, 1.0e-5))
def test_kamata_solutions_stay_independent_at_low_frequency(frequency):
    """The three starting solutions stay finite and separated however low the frequency. After scaling each radial
    function and then each solution, their singular value ratio stays near 6e-4; the published pair gave 1e-8 at
    1e-5 rad/s and roundoff (1e-16) below that."""
    starting = np.zeros((3, 6), dtype=np.complex128, order="C")
    kamata_solid_dynamic_incompressible(frequency, 1.0e5, 3000.0, 4.0e10 + 0.0j, 2, G, starting)
    assert np.all(np.isfinite(starting.view(np.float64)))
    y_scale = np.max(np.abs(starting), axis=0)
    y_scale[y_scale == 0.0] = 1.0
    scaled = starting / y_scale
    scaled /= np.max(np.abs(scaled), axis=1)[:, None]
    singular_values = np.linalg.svd(scaled, compute_uv=False)
    assert singular_values[-1] / singular_values[0] > 1.0e-5


@pytest.mark.parametrize("degree_l", (2, 3))
@pytest.mark.parametrize("period_days", (16.0, 27.3, 80.0))
def test_kamata_low_frequency_matches_the_closed_form(period_days, degree_l):
    """At long periods a homogeneous incompressible Moon matches the static closed form to 1e-5 at the configured
    tolerance. The published basis missed by 3e-4 to 2% (k2) and up to 22% (k3) over these periods."""
    solution = radial_solver(*_moon_inputs(period_days), degree_l=degree_l, use_kamata=True)
    assert solution.success, solution.message
    expected = _closed_form_incompressible_k(degree_l)
    assert abs(solution.k - expected) / expected < 1.0e-5
    assert solution.surface_solve_amplification < 1.0e2
    assert solution.surface_solve_rcond > 1.0e-3


# =====================================================================================================================
# Starting radius
# =====================================================================================================================
@pytest.mark.parametrize("fraction", (1.0, 1.5))
def test_standalone_refuses_a_starting_radius_at_or_above_the_surface(fraction):
    with pytest.raises(ValueError, match="planet radius"):
        radial_solver(*_static_one_layer_inputs(), degree_l=2, starting_radius=fraction * 1.0e6)


@pytest.mark.parametrize("fraction", (0.95, 0.999, 1.0, 1.5))
def test_world_refuses_a_starting_radius_near_or_above_the_surface(io, fraction):
    """Above max_start_radius_fraction (0.9) of the radius, the world solve fails with a clear message instead of
    returning NaN (at or above the surface) or a wrong k2 marked as a success (just below it)."""
    result = io.solve_love_numbers(frequency=IO_FREQUENCY, starting_radius=fraction * io.radius)
    assert not result["success"]
    assert result["error_code"] == -5
    assert "starting radius" in result["message"]
    assert math.isnan(result["love_number_k"].real)


def test_world_accepts_a_starting_radius_inside_the_body(io):
    automatic = io.solve_love_numbers(frequency=IO_FREQUENCY)
    manual = io.solve_love_numbers(frequency=IO_FREQUENCY, starting_radius=0.5 * io.radius)
    assert automatic["success"] and manual["success"], manual["message"]
    assert np.isfinite(manual["love_number_k"].real)
