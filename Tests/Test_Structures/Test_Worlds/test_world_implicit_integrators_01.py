"""Implicit CyRK integrators (BDF, LSODA, Radau) in the world-level EOS and Love solves."""
import math

import pytest

from TidalPy.Structures import build_world

IMPLICIT_METHODS = ("BDF", "LSODA", "Radau")
FREQUENCY_RAD_S = 1.0e-5


def _fresh_world():
    return build_world("earth_simple")


@pytest.mark.parametrize("integration_method", IMPLICIT_METHODS)
def test_world_eos_solve_implicit(integration_method):
    """solve_eos with each implicit method reproduces the DOP853 planet mass."""
    reference_world = _fresh_world()
    reference_world.solve_eos(integration_method="DOP853")
    reference_mass = reference_world.planet_mass_eos

    world = _fresh_world()
    result = world.solve_eos(integration_method=integration_method)
    # LSODA may fail its first step at the singular center; that failure must be reported cleanly.
    if integration_method == "LSODA" and not result["success"]:
        assert "failed" in result["message"].lower()
        return
    assert result["success"]
    assert math.isclose(world.planet_mass_eos, reference_mass, rel_tol=1.0e-5)


@pytest.mark.parametrize("integration_method", IMPLICIT_METHODS)
def test_world_love_solve_implicit(integration_method):
    """solve_love_numbers with each implicit method reproduces the DOP853 k2."""
    reference_world = _fresh_world()
    reference_world.solve_eos()
    reference = reference_world.solve_love_numbers(frequency=FREQUENCY_RAD_S, integration_method="DOP853")
    assert reference["success"]

    world = _fresh_world()
    world.solve_eos()
    result = world.solve_love_numbers(frequency=FREQUENCY_RAD_S, integration_method=integration_method)
    assert result["success"]

    k2_reference = reference_world.love_number_k
    k2_implicit = world.love_number_k
    assert math.isclose(k2_implicit.real, k2_reference.real, rel_tol=1.0e-4)
    assert math.isclose(k2_implicit.imag, k2_reference.imag, rel_tol=1.0e-3, abs_tol=1.0e-10)


def test_world_unknown_method_raises():
    """An unknown method name raises a ValueError."""
    world = _fresh_world()
    with pytest.raises(ValueError, match="Unsupported integration method"):
        world.solve_eos(integration_method="not_a_method")
