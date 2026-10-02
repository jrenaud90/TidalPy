"""A world with pressure-dependent melting curves and an expansivity that falls with compression: the thermal solve,
the melt, and the round trips through the world config and the binary file."""
import math
import os
import tempfile

import numpy as np
import pytest

from TidalPy.Material.laws import make_eos
from TidalPy.Structures import build_world
from TidalPy.Structures.worlds import TerrestrialWorld
from numpy_compat import trapezoid

# Monteux et al. (2016) peridotite melting curves.
_MONTEUX_SOLIDUS = {
    "model": "simon_glatzel_2", "temperature_k": 1661.2, "simon_a_pa": 1.336e9, "simon_c": 7.437,
    "transition_pressure_pa": 20.0e9, "high_temperature_k": 2081.8, "high_simon_a_pa": 1.0169e11,
    "high_simon_c": 1.226}
_MONTEUX_LIQUIDUS = {
    "model": "simon_glatzel_2", "temperature_k": 1982.1, "simon_a_pa": 6.594e9, "simon_c": 5.374,
    "transition_pressure_pa": 20.0e9, "high_temperature_k": 2006.8, "high_simon_a_pa": 3.465e10,
    "high_simon_c": 1.844}
_SURFACE_TEMPERATURE = 300.0
_MANTLE_TEMPERATURE = 1600.0
# Silicate and iron thermal constants: conductivity [W m-1 K-1], heat capacity [J kg-1 K-1], expansivity [1/K].
_ROCK_THERMAL = {"thermal_conductivity_w_mk": 3.75, "heat_capacity_j_kgk": 1200.0}
_ROCK_EXPANSION = 5.2e-5
_IRON_THERMAL = {"thermal_conductivity_w_mk": 7.95, "heat_capacity_j_kgk": 840.0}
_IRON_EXPANSION = 1.2e-5


def _thermal_earth_config(pressure_dependent):
    """Bundled earth_simple, ready for a thermal solve: iron cores, and a convecting mantle with silicate thermal
    constants that melts into a Murnaghan melt of 0.2 Pa s through Henning weakening. Its melting curves are a constant
    1600 K solidus and 2000 K liquidus, or with pressure_dependent the Monteux curves, which follow the pressure."""
    config = build_world("earth_simple").get_config_dict()
    for name in ("inner_core", "outer_core"):
        layer = config["layers"][name]
        for phase in layer["material"].values():
            if isinstance(phase, dict):
                phase.update(_IRON_THERMAL)
                phase["eos"]["thermal_expansion_1_k"] = _IRON_EXPANSION
        layer["cooling"] = {"model": "off"}
        layer["radiogenics"] = {"model": "off"}
    mantle = config["layers"]["mantle"]
    material = mantle["material"]
    material["solid"].update(_ROCK_THERMAL)
    material["solid"]["eos"]["thermal_expansion_1_k"] = _ROCK_EXPANSION
    material["liquid"] = {
        **_ROCK_THERMAL,
        "eos": {"model": "murnaghan", "reference_density_kg_m3": 2750.0, "reference_bulk_modulus_pa": 2.0e10,
                "bulk_modulus_derivative": 5.0, "thermal_expansion_1_k": _ROCK_EXPANSION},
        "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 0.2}}
    material["melting"] = {"solidus": {"model": "constant", "temperature_k": 1600.0},
                           "liquidus": {"model": "constant", "temperature_k": 2000.0},
                           "weakening": {"model": "henning"}}
    mantle["use_melting"] = True
    mantle["cooling"] = {"model": "convection", "convection_alpha": 1.0, "convection_beta": 1.0 / 3.0,
                         "critical_rayleigh": 1100.0}
    mantle["radiogenics"] = {"model": "isotope", "isotopes": "modern_day_chondritic"}
    if pressure_dependent:
        material["melting"]["solidus"] = dict(_MONTEUX_SOLIDUS)
        material["melting"]["liquidus"] = dict(_MONTEUX_LIQUIDUS)
        mantle["use_pressure_melting"] = True
    return config


def _earth(pressure_dependent):
    config = _thermal_earth_config(pressure_dependent)
    config["layers"]["mantle"]["temperature_k"] = _MANTLE_TEMPERATURE
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
    Monteux curves it is solid throughout its interior, since its expansivity falls with compression. Melt, if any, is confined to
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


def test_adiabat_base_follows_the_compressed_expansivity():
    """The reported base of the mantle's adiabat is T exp(int alpha g / c_p dr) over the solved structure, with the
    expansivity alpha0 K0 / K_T that the law's thermal pressure gives its density."""
    earth = _earth(pressure_dependent=True)
    result = _solve(earth)
    solid = earth.mantle.get_config_dict()["material"]["solid"]
    alpha0 = solid["eos"]["thermal_expansion_1_k"]
    heat_capacity = solid["heat_capacity_j_kgk"]
    law = make_eos(solid["eos"]["model"], solid["eos"])
    boundary = result["layer_boundary_thickness"][2]
    radii = np.linspace(earth.mantle.radius_inner + boundary, earth.mantle.radius_outer - boundary, 4001)
    temperature = earth.get_temperature(radii)
    bulk_modulus = law.calc_eos(
        earth.get_pressure(radii), temperature, thermal=earth.mantle.use_thermal_expansion)["bulk_modulus"]
    alpha = alpha0 * solid["eos"]["reference_bulk_modulus_pa"] / bulk_modulus
    exponent = trapezoid(alpha * earth.get_gravity(radii), radii) / heat_capacity
    assert result["layer_base_temperature"][2] == pytest.approx(_MANTLE_TEMPERATURE * math.exp(exponent), rel=1e-5)


def test_parameters_round_trip_through_the_world_config():
    """The world's config carries the curve and expansivity keys, and a world rebuilt from it solves identically."""
    earth = _earth(pressure_dependent=True)
    config = earth.get_config_dict()
    material = config["layers"]["mantle"]["material"]
    assert config["layers"]["mantle"]["use_pressure_melting"] is True
    for curve, expected in (("solidus", _MONTEUX_SOLIDUS), ("liquidus", _MONTEUX_LIQUIDUS)):
        for key, value in expected.items():
            assert material["melting"][curve][key] == value, (curve, key)
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
