"""CyRK's explicit Tsit5, Vern7, and Vern8 integrators in the EOS and Love solves, by call, config, and world pin."""
import math

import numpy as np
import pytest

import TidalPy
from TidalPy.RadialSolver.helpers import homogeneous_love_numbers
from TidalPy.Structures import build_world

EXPLICIT_METHODS = ("Tsit5", "Vern7", "Vern8")
FREQUENCY_RAD_S = 1.0e-5
# Tight tolerances for the comparisons against DOP853, so each method's own error sits well below the check.
EOS_TOLERANCES = dict(rtol=1.0e-11, atol=1.0e-14)
LOVE_TOLERANCES = dict(rtol=1.0e-10, atol=1.0e-10)


def _relative(value, reference):
    return abs(value - reference) / abs(reference)


def _solved(world_name, integration_method, solve_temperature=None):
    world = build_world(world_name)
    result = world.solve_eos(
        integration_method=integration_method,
        solve_temperature=solve_temperature,
        **EOS_TOLERANCES,
    )
    assert result["success"], result["message"]
    return world


@pytest.mark.parametrize("integration_method", EXPLICIT_METHODS)
@pytest.mark.parametrize("world_name, solve_temperature", (
    ("io", None),
    ("earth_simple", None),
    ("earth_thermal", True),
))
def test_world_eos_solve_matches_dop853(world_name, solve_temperature, integration_method):
    """An isothermal and a thermal EOS solve reproduce the DOP853 mass, moment of inertia, and surface gravity."""
    reference = _solved(world_name, "DOP853", solve_temperature)
    world = _solved(world_name, integration_method, solve_temperature)
    assert _relative(world.planet_mass_eos, reference.planet_mass_eos) < 1.0e-8
    assert _relative(world.planet_moi_eos, reference.planet_moi_eos) < 1.0e-8
    assert _relative(world.calc_surface_gravity(), reference.calc_surface_gravity()) < 1.0e-8


@pytest.mark.parametrize("integration_method", EXPLICIT_METHODS)
@pytest.mark.parametrize("world_name", ("io", "earth_thermal"))
def test_world_love_solve_matches_dop853(world_name, integration_method):
    """The k2, h2, and l2 of each method agree with DOP853's at a tight tolerance on the same EOS solution."""
    world = _solved(world_name, "DOP853")
    reference = world.solve_love_numbers(frequency=FREQUENCY_RAD_S, integration_method="DOP853", **LOVE_TOLERANCES)
    assert reference["success"], reference["message"]
    reference_love = (world.love_number_k, world.love_number_h, world.love_number_l)

    result = world.solve_love_numbers(frequency=FREQUENCY_RAD_S, integration_method=integration_method,
                                      **LOVE_TOLERANCES)
    assert result["success"], result["message"]
    for love, love_reference in zip((world.love_number_k, world.love_number_h, world.love_number_l), reference_love):
        # Both parts carry about the same absolute error, so each is checked against the number's size. The
        # thermal Earth's differences reach about 1.5e-8 at these tolerances.
        assert abs(love.real - love_reference.real) < 1.0e-7 * abs(love_reference)
        assert abs(love.imag - love_reference.imag) < 1.0e-7 * abs(love_reference)


@pytest.mark.parametrize("integration_method", EXPLICIT_METHODS)
def test_standalone_radial_solver_matches_dop853(integration_method):
    """The standalone radial_solver accepts each method for both its radial and its EOS integrations."""
    def solve(method):
        solution = homogeneous_love_numbers(
            6.0e6,
            4000.0,
            5.0e10 + 1.0e9j,
            FREQUENCY_RAD_S,
            integration_method=method,
            eos_integration_method=method,
            integration_rtol=1.0e-10,
            integration_atol=1.0e-12,
        )
        assert solution.success, solution.message
        return complex(np.atleast_1d(solution.k)[0]), complex(np.atleast_1d(solution.h)[0])
    k_reference, h_reference = solve("DOP853")
    k_value, h_value = solve(integration_method.upper())
    assert _relative(k_value, k_reference) < 1.0e-8
    assert _relative(h_value, h_reference) < 1.0e-8


@pytest.mark.parametrize("integration_method", EXPLICIT_METHODS)
def test_configured_method_reaches_the_solves(integration_method, restore_config):
    """A method set in [radial_solver] and [eos_solver] is the one a call without an argument uses."""
    explicit_world = _solved("io", integration_method)
    explicit_world.solve_love_numbers(frequency=FREQUENCY_RAD_S, integration_method=integration_method)

    TidalPy.reinit(provided_config={
        "eos_solver": {"integration_method": integration_method.lower()},
        "radial_solver": {"integration_method": integration_method.upper()},
    })
    world = build_world("io")
    assert world.solve_eos(**EOS_TOLERANCES)["success"]
    assert world.solve_love_numbers(frequency=FREQUENCY_RAD_S)["success"]
    assert world.planet_mass_eos == explicit_world.planet_mass_eos
    assert world.love_number_k == explicit_world.love_number_k


@pytest.mark.parametrize("integration_method", EXPLICIT_METHODS)
def test_world_pin_round_trips(integration_method):
    """A world file may pin each method, in any case, and its config dict gives back the canonical name."""
    world = build_world("io", overrides={
        "eos_solver": {"integration_method": integration_method.upper()},
        "radial_solver": {"integration_method": integration_method.lower()},
    })
    pinned = world.get_solver_defaults()
    assert pinned["eos_solver"]["integration_method"] == integration_method
    assert pinned["radial_solver"]["integration_method"] == integration_method
    assert build_world(world.get_config_dict()).get_solver_defaults() == pinned
    assert world.solve_eos()["success"]
    assert world.solve_love_numbers(frequency=FREQUENCY_RAD_S)["success"]
    assert math.isfinite(world.love_number_k.real)


def test_unknown_method_raises():
    """A name close to a new method is still refused."""
    world = build_world("io")
    with pytest.raises(ValueError, match="Unsupported integration method"):
        world.solve_eos(integration_method="vern9")
