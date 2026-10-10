"""CyRK's explicit Tsit5, Vern7, and Vern8 integrators in the EOS and Love solves, by call, config, and world pin."""
import numpy as np
import pytest

import TidalPy
from TidalPy.RadialSolver.helpers import homogeneous_love_numbers
from TidalPy.Structures import build_world

EXPLICIT_METHODS = ("Tsit5", "Vern7", "Vern8")
FORCING_FREQUENCY = 1.0e-5  # [rad s-1]
# Tight tolerances for the comparisons against DOP853, so each method's own error sits well below the check.
EOS_TOLERANCES = dict(rtol=1.0e-11, atol=1.0e-14)
LOVE_TOLERANCES = dict(rtol=1.0e-10, atol=1.0e-10)
# Largest Love-number difference from DOP853, relative to the number's size, at LOVE_TOLERANCES. Measured worst cases
# are 6e-10 (Tsit5 and Vern7) and, for Vern8 on earth_thermal, 1.5e-8: Vern8's error estimate underestimates the
# error of a step across the kinks of a thermal profile.
LOVE_BOUND = 5.0e-9
LOVE_BOUND_VERN8_THERMAL = 1.0e-7


def relative_difference(value, reference):
    return abs(value - reference) / abs(reference)


def solved_world(world_name, integration_method, solve_temperature=None):
    world = build_world(world_name)
    result = world.solve_eos(
        integration_method=integration_method,
        solve_temperature=solve_temperature,
        **EOS_TOLERANCES,
    )
    assert result["success"], result["message"]
    return world, result


def love_numbers(world):
    return world.love_number_k, world.love_number_h, world.love_number_l


@pytest.mark.parametrize("integration_method", EXPLICIT_METHODS)
@pytest.mark.parametrize("world_name, solve_temperature", (
    ("io", None),
    ("earth_simple", None),
    ("earth_thermal", True),
))
def test_world_eos_solve_matches_dop853(world_name, solve_temperature, integration_method):
    """An isothermal and a thermal EOS solve reproduce the DOP853 mass, moment of inertia, and surface gravity."""
    reference, reference_result = solved_world(world_name, "DOP853", solve_temperature)
    world, result = solved_world(world_name, integration_method, solve_temperature)
    # Measured differences are at most 5e-12.
    assert relative_difference(world.planet_mass_eos, reference.planet_mass_eos) < 1.0e-10
    assert relative_difference(world.planet_moi_eos, reference.planet_moi_eos) < 1.0e-10
    assert relative_difference(result["surface_gravity"], reference_result["surface_gravity"]) < 1.0e-10
    # A different integration, not a fallback to DOP853.
    assert (world.planet_mass_eos, world.planet_moi_eos) != (reference.planet_mass_eos, reference.planet_moi_eos)


@pytest.mark.parametrize("integration_method", EXPLICIT_METHODS)
@pytest.mark.parametrize("world_name", ("io", "earth_thermal"))
def test_world_love_solve_matches_dop853(world_name, integration_method):
    """The k2, h2, and l2 of each method agree with DOP853's at a tight tolerance on the same EOS solution."""
    world, _ = solved_world(world_name, "DOP853")
    reference = world.solve_love_numbers(frequency=FORCING_FREQUENCY, integration_method="DOP853", **LOVE_TOLERANCES)
    assert reference["success"], reference["message"]
    reference_love = love_numbers(world)
    reference_steps = np.asarray(world.release_radial_solution().steps_taken).sum()

    result = world.solve_love_numbers(frequency=FORCING_FREQUENCY, integration_method=integration_method,
                                      **LOVE_TOLERANCES)
    assert result["success"], result["message"]
    bound = LOVE_BOUND_VERN8_THERMAL if (world_name, integration_method) == ("earth_thermal", "Vern8") else LOVE_BOUND
    for love, love_reference in zip(love_numbers(world), reference_love):
        # Both parts carry about the same absolute error, so each is checked against the number's size.
        assert abs(love.real - love_reference.real) < bound * abs(love_reference)
        assert abs(love.imag - love_reference.imag) < bound * abs(love_reference)
    assert np.asarray(world.release_radial_solution().steps_taken).sum() != reference_steps


@pytest.mark.parametrize("integration_method", EXPLICIT_METHODS)
def test_standalone_radial_solver_matches_dop853(integration_method):
    """The standalone radial_solver accepts each method for both its radial and its EOS integrations."""
    def solve(method):
        solution = homogeneous_love_numbers(
            6.0e6,
            4000.0,
            5.0e10 + 1.0e9j,
            FORCING_FREQUENCY,
            integration_method=method,
            eos_integration_method=method,
            integration_rtol=1.0e-10,
            integration_atol=1.0e-12,
        )
        assert solution.success, solution.message
        love = complex(np.atleast_1d(solution.k)[0]), complex(np.atleast_1d(solution.h)[0])
        return love, np.asarray(solution.steps_taken).sum()
    (k_reference, h_reference), reference_steps = solve("DOP853")
    (k_value, h_value), steps = solve(integration_method.upper())
    # Measured differences are at most 3.5e-11.
    assert relative_difference(k_value, k_reference) < 1.0e-9
    assert relative_difference(h_value, h_reference) < 1.0e-9
    # The step count shows the requested method ran (DOP853 83, Tsit5 209, Vern7 77, Vern8 44 per solution).
    assert steps != reference_steps


@pytest.mark.parametrize("integration_method", EXPLICIT_METHODS)
def test_configured_method_reaches_the_solves(integration_method, restore_config):
    """A method set in [radial_solver] and [eos_solver] is the one a call without an argument uses."""
    explicit_world, _ = solved_world("io", integration_method)
    explicit_world.solve_love_numbers(frequency=FORCING_FREQUENCY, integration_method=integration_method)
    dop853_world, _ = solved_world("io", "DOP853")
    dop853_world.solve_love_numbers(frequency=FORCING_FREQUENCY, integration_method="DOP853")

    TidalPy.reinit(provided_config={
        "eos_solver": {"integration_method": integration_method.lower()},
        "radial_solver": {"integration_method": integration_method.upper()},
    })
    world = build_world("io")
    assert world.solve_eos(**EOS_TOLERANCES)["success"]
    assert world.solve_love_numbers(frequency=FORCING_FREQUENCY)["success"]
    assert world.planet_mass_eos == explicit_world.planet_mass_eos
    assert world.love_number_k == explicit_world.love_number_k
    assert world.love_number_k != dop853_world.love_number_k


@pytest.mark.parametrize("integration_method", EXPLICIT_METHODS)
def test_world_pin_round_trips(integration_method, tmp_path):
    """A world file may pin each method, in any case: the config dict gives back the canonical name, the pin survives
    a binary round trip, and a pinned solve equals one given the method as an argument."""
    world = build_world("io", overrides={
        "eos_solver": {"integration_method": integration_method.upper()},
        "radial_solver": {"integration_method": integration_method.lower()},
    })
    pinned = world.get_solver_defaults()
    assert pinned["eos_solver"]["integration_method"] == integration_method
    assert pinned["radial_solver"]["integration_method"] == integration_method
    assert build_world(world.get_config_dict()).get_solver_defaults() == pinned

    path = str(tmp_path / "pinned.tpyb")
    world.save_binary(path)
    loaded = type(world)("placeholder", 1.0, 1.0)
    loaded.load_binary(path)
    assert loaded.get_solver_defaults() == pinned

    explicit = build_world("io")
    assert explicit.solve_eos(integration_method=integration_method)["success"]
    assert explicit.solve_love_numbers(frequency=FORCING_FREQUENCY, integration_method=integration_method)["success"]
    for pinned_world in (world, loaded):
        assert pinned_world.solve_eos()["success"]
        assert pinned_world.solve_love_numbers(frequency=FORCING_FREQUENCY)["success"]
        assert pinned_world.planet_mass_eos == explicit.planet_mass_eos
        assert love_numbers(pinned_world) == love_numbers(explicit)


def test_unknown_method_raises():
    """A name close to a new method is still refused."""
    world = build_world("io")
    with pytest.raises(ValueError, match="Unsupported integration method"):
        world.solve_eos(integration_method="vern9")
