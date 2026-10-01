"""A world with pressure-dependent melting curves and an Anderson-Gruneisen expansivity: the thermal solve, the melt,
and the round trips through the world config and the binary file."""
import math
import os
import tempfile

import numpy as np
import pytest

from TidalPy.Structures import build_world
from TidalPy.Structures.worlds import TerrestrialWorld

# Monteux et al. (2016) peridotite melting curves, and mantle-silicate Anderson-Gruneisen constants.
_MONTEUX = {
    "solidus_k": 1661.2, "solidus_simon_a_pa": 1.336e9, "solidus_simon_c": 7.437,
    "solidus_transition_pressure_pa": 20.0e9, "solidus_high_k": 2081.8, "solidus_high_simon_a_pa": 1.0169e11,
    "solidus_high_simon_c": 1.226,
    "liquidus_k": 1982.1, "liquidus_simon_a_pa": 6.594e9, "liquidus_simon_c": 5.374,
    "liquidus_transition_pressure_pa": 20.0e9, "liquidus_high_k": 2006.8, "liquidus_high_simon_a_pa": 3.465e10,
    "liquidus_high_simon_c": 1.844,
}
_DELTA = 5.5
_KAPPA = 1.4
_SURFACE_TEMPERATURE = 300.0
_MANTLE_TEMPERATURE = 1600.0


def _earth(pressure_dependent):
    config = build_world("earth_simple").get_config_dict()
    mantle = config["layers"]["mantle"]
    mantle["temperature_k"] = _MANTLE_TEMPERATURE
    if pressure_dependent:
        mantle["material"]["anderson_gruneisen_parameter"] = _DELTA
        mantle["material"]["anderson_gruneisen_exponent"] = _KAPPA
        mantle["material"]["partial_melt"].update(_MONTEUX)
    return build_world(config)


def _solve(world):
    result = world.solve_eos(solve_temperature=True, surface_temperature=_SURFACE_TEMPERATURE)
    assert result["success"], result["message"]
    assert result["thermal_converged"]
    return result


def _mantle_melt(world, fraction):
    radius = world.mantle.radius_inner + fraction * (world.mantle.radius_outer - world.mantle.radius_inner)
    return world.get_state(radius)["melt_fraction"]


def test_pressure_dependent_curves_keep_the_lower_mantle_solid():
    """With constant melting curves the deep mantle of an Earth at a 1600 K upper mantle reads as molten; with the
    Monteux curves and a compressible expansivity it is solid throughout its interior. Melt, if any, is confined to
    the thermal boundary layer against the 4500 K outer core, as at Earth's core-mantle boundary."""
    constant = _earth(pressure_dependent=False)
    _solve(constant)
    assert _mantle_melt(constant, 0.05) > 0.5
    earth = _earth(pressure_dependent=True)
    result = _solve(earth)
    for name, molten_inner, molten_outer in earth.molten_regions:
        assert name == "mantle"
        assert molten_inner == pytest.approx(earth.mantle.radius_inner, rel=1e-9)
        assert molten_outer - molten_inner < result["layer_boundary_thickness"][2]
    for fraction in (0.05, 0.5, 0.95):
        assert _mantle_melt(earth, fraction) == 0.0
    # The adiabat warms the 1600 K upper mantle to well under the deep solidus.
    assert _MANTLE_TEMPERATURE < result["layer_base_temperature"][2] < 3000.0


def test_adiabat_base_follows_the_anderson_gruneisen_expansivity():
    """The reported base of the mantle's adiabat is T exp(int alpha(rho) g / c_p dr) over the solved structure."""
    earth = _earth(pressure_dependent=True)
    result = _solve(earth)
    material = earth.mantle.get_config_dict()["material"]
    alpha0 = material["thermal_expansion_1_k"]
    heat_capacity = material["heat_capacity_j_kgk"]
    reference_density = material["reference_density_kg_m3"]
    boundary = result["layer_boundary_thickness"][2]
    radii = np.linspace(earth.mantle.radius_inner + boundary, earth.mantle.radius_outer - boundary, 4001)
    density = earth.get_density(radii)
    alpha = alpha0 * np.exp((_DELTA / _KAPPA) * ((reference_density / density) ** _KAPPA - 1.0))
    exponent = np.trapezoid(alpha * earth.get_gravity(radii), radii) / heat_capacity
    assert result["layer_base_temperature"][2] == pytest.approx(_MANTLE_TEMPERATURE * math.exp(exponent), rel=1e-5)


def test_parameters_round_trip_through_the_world_config():
    """The world's config carries the curve and expansivity keys, and a world rebuilt from it solves identically."""
    earth = _earth(pressure_dependent=True)
    config = earth.get_config_dict()
    material = config["layers"]["mantle"]["material"]
    assert material["anderson_gruneisen_parameter"] == _DELTA
    assert material["anderson_gruneisen_exponent"] == _KAPPA
    for key, value in _MONTEUX.items():
        assert material["partial_melt"][key] == value, key
    rebuilt = build_world(config)
    assert rebuilt.get_config_dict() == config
    assert _solve(rebuilt)["layer_base_temperature"] == _solve(earth)["layer_base_temperature"]


def test_parameters_round_trip_through_the_binary_file():
    """A world saved to binary and loaded back keeps the curves and the expansivity law."""
    earth = _earth(pressure_dependent=True)
    path = os.path.join(tempfile.mkdtemp(prefix="tidalpy_melt_"), "earth.tpyb")
    earth.save_binary(path)
    loaded = TerrestrialWorld("placeholder", 1.0, 1.0)
    loaded.load_binary(path)
    assert loaded.get_config_dict()["layers"] == earth.get_config_dict()["layers"]
    assert _solve(loaded)["layer_base_temperature"] == _solve(earth)["layer_base_temperature"]
