"""End-to-end check of the rheology path: configuration -> material state -> complex moduli -> Love numbers.

A world built from a configuration carries a viscosity law, a rheology, and a temperature. What the radial solver
sees at each radius is the rheology applied to the viscosity the law gives at that radius's pressure and the layer's
temperature. Each link is checked here against the standalone models, so a break anywhere between the builder and
the radial solver (a model that is attached but not evaluated, a temperature or a pressure that does not reach it)
fails a test, which the per-module unit tests cannot see.
"""
import copy
import math

import numpy as np
import pytest

from TidalPy.rheology_x.rheology import make_rheology
from TidalPy.structures_x import build_world, build_world_from_dict
from TidalPy.viscosity_x import make_viscosity


_FREQUENCY = 2.0 * math.pi / (2.0 * 86400.0)
_TEMPERATURE = 1500.0

# An activation volume large enough that the pressure across this small mantle visibly changes the viscosity.
_VISCOSITY = {
    "model": "reference",
    "reference_viscosity_pas": 1.0e19,
    "reference_temperature_k": 1600.0,
    "molar_activation_energy_j_mol": 3.0e5,
    "molar_activation_volume_m3_mol": 1.0e-5,
}
_RHEOLOGY = {"model": "andrade", "alpha": 0.3, "zeta": 1.0}

_CONFIG = {
    "name": "rheology_path", "type": "terrestrial", "radius_m": 3.0e6, "mass_kg": 6.0e23,
    "tides": {"global_tidal_model": "rheology", "max_degree_l": 2, "eccentricity_trunc_lvl": 2},
    "layers": {
        "core": {
            "class": "physics", "type": "none", "radius_fraction": 0.45, "is_tidal": False,
            "material": {
                "model": "constant", "reference_density_kg_m3": 8000.0,
                "shear_modulus_static_pa": 1.0e11, "bulk_modulus_static_pa": 2.0e11},
            "shear_rheology": {"model": "elastic"}, "bulk_rheology": {"model": "elastic"},
        },
        "mantle": {
            "class": "solidliquid", "type": "none", "radius_fraction": 1.0, "temperature_k": _TEMPERATURE,
            "material": {
                "model": "constant", "reference_density_kg_m3": 3400.0,
                "shear_modulus_static_pa": 6.0e10, "bulk_modulus_static_pa": 1.5e11,
                "shear_viscosity": dict(_VISCOSITY),
                "bulk_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e24}},
            "shear_rheology": dict(_RHEOLOGY), "bulk_rheology": {"model": "elastic"},
        },
    },
}


def _solved_world(temperature=_TEMPERATURE):
    config = copy.deepcopy(_CONFIG)
    config["layers"]["mantle"]["temperature_k"] = temperature
    world = build_world(config)
    assert world.solve_eos()["success"]
    return world


def _standalone(model_config, factory):
    params = {key: value for key, value in model_config.items() if key != "model"}
    return factory(model_config["model"], params)


def _mantle_radii(world):
    mantle = world.mantle
    return np.linspace(mantle.radius_inner, mantle.radius_outer, 7)[1:-1]


def test_viscosity_law_sees_the_pressure_and_the_temperature():
    world = _solved_world()
    law = _standalone(_VISCOSITY, make_viscosity)
    radii = _mantle_radii(world)
    pressures = np.asarray(world.get_pressure(radii))
    expected = np.array([law.calc_viscosity(_TEMPERATURE, pressure) for pressure in pressures])
    np.testing.assert_allclose(np.asarray(world.get_shear_viscosity(radii)), expected, rtol=1.0e-10)
    # The law is doing real work here: the deep mantle is stiffer than the shallow mantle.
    assert expected[0] / expected[-1] > 1.5


def test_complex_modulus_is_the_rheology_of_the_material_state():
    world = _solved_world()
    rheology = _standalone(_RHEOLOGY, make_rheology)
    for radius in _mantle_radii(world):
        expected = rheology.calc_complex_modulus(
            float(world.get_shear_modulus(radius)), float(world.get_shear_viscosity(radius)), _FREQUENCY)
        found = world.calc_complex_shear_modulus(radius, _FREQUENCY)
        assert found == pytest.approx(expected, rel=1.0e-12)
        assert found.imag > 0.0


def test_love_number_responds_to_the_viscosity_law():
    """A colder mantle is stiffer in the law, and the radial solver has to see it."""
    hot = _solved_world(1700.0)
    cold = _solved_world(1300.0)
    k2 = {}
    for name, world in (("hot", hot), ("cold", cold)):
        result = world.solve_love_numbers(frequency=_FREQUENCY)
        assert result["success"], result["message"]
        k2[name] = result["love_number_k"]
    assert -k2["hot"].imag > 0.0 and -k2["cold"].imag > 0.0
    assert -k2["hot"].imag > 5.0 * -k2["cold"].imag
    assert k2["hot"].real > k2["cold"].real


def test_rebuilt_world_takes_the_same_path(tmp_path):
    """The file a world writes rebuilds a world whose radial solver sees the same material."""
    world = _solved_world()
    k2 = world.solve_love_numbers(frequency=_FREQUENCY)["love_number_k"]

    file_path = str(tmp_path / "rheology_path.toml")
    world.save_to_toml(file_path)
    for rebuilt in (build_world(file_path), build_world_from_dict(world.get_config_dict())):
        assert rebuilt.solve_eos()["success"]
        radius = float(_mantle_radii(world)[2])
        assert rebuilt.get_shear_viscosity(radius) == pytest.approx(world.get_shear_viscosity(radius), rel=1.0e-12)
        assert rebuilt.solve_love_numbers(frequency=_FREQUENCY)["love_number_k"] == pytest.approx(k2, rel=1.0e-9)


# =====================================================================================================================
# A liquid layer, against a limit that can fail
# =====================================================================================================================
def _core_world(core_is_liquid, core_shear=1.0e7, core_is_static=True):
    config = copy.deepcopy(_CONFIG)
    core = config["layers"]["core"]
    core["material"]["shear_modulus_static_pa"] = 0.0 if core_is_liquid else core_shear
    if core_is_liquid:
        core["is_solid"] = False
        core["is_static"] = core_is_static
    world = build_world(config)
    assert world.solve_eos()["success"]
    return world


# A dynamic liquid layer is only well conditioned at short forcing periods (at two days this world's surface solve
# amplifies errors by 1e8 and k2 is wrong in its second digit), so it is held to the limit where the method works.
@pytest.mark.parametrize("core_is_static, period_days", [(True, 2.0), (True, 10.0), (False, 0.5)])
def test_liquid_core_is_the_limit_of_a_vanishing_core_rigidity(core_is_static, period_days):
    """A solid core whose rigidity is small beside rho g R responds like a liquid one.

    The two solves share nothing past the interior model: the solid core is integrated with six radial functions
    and the liquid one with two (static) or four (dynamic) and their own interface conditions, so agreement is a
    check of the liquid-layer machinery that a wrong interface or starting condition breaks. The solid core's
    rigidity is 1e7 Pa against rho g R of order 1e10 Pa, so the two differ at the 1e-3 level at most; a rigid
    core (1e11 Pa) differs by tens of percent, which is what the second assertion holds the test to.
    """
    frequency = 2.0 * math.pi / (period_days * 86400.0)
    liquid_world = _core_world(True, core_is_static=core_is_static)
    liquid = liquid_world.solve_love_numbers(frequency=frequency)
    assert liquid_world.love_surface_amplification < 1.0e3, "the comparison needs a well conditioned liquid solve"
    soft = _core_world(False, core_shear=1.0e7).solve_love_numbers(frequency=frequency)
    rigid = _core_world(False, core_shear=1.0e11).solve_love_numbers(frequency=frequency)
    for result in (liquid, soft, rigid):
        assert result["success"], result["message"]

    k2_liquid, k2_soft, k2_rigid = (result["love_number_k"].real for result in (liquid, soft, rigid))
    assert k2_liquid == pytest.approx(k2_soft, rel=2.0e-3)
    assert abs(k2_rigid - k2_liquid) > 0.05 * k2_liquid
    h2_liquid, h2_soft = liquid["love_number_h"].real, soft["love_number_h"].real
    assert h2_liquid == pytest.approx(h2_soft, rel=2.0e-3)
