"""The world tide interface: optional orbital state, the calc_tides result, the near-synchronous warning, the tide
configuration round trip and its context manager, and tide model access."""
import math

import numpy as np
import pytest

from TidalPy.Structures import build_world
from TidalPy.Structures.system import System
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Tides.classes.tide import FixedQTide, make_tide

# Io-Jupiter-like orbital state.
HOST_MASS = 1.898e27
SEMI_MAJOR_AXIS = 4.2e8
ECCENTRICITY = 0.0041
ORBITAL_FREQUENCY = 2.05e-5
ORBIT = {
    "orbital_frequency": ORBITAL_FREQUENCY,
    "spin_frequency": ORBITAL_FREQUENCY,
    "eccentricity": ECCENTRICITY,
    "obliquity": 0.0,
    "semi_major_axis": SEMI_MAJOR_AXIS,
    "host_mass": HOST_MASS,
}
LAYER_MATERIAL = {"solid": {
    "eos": {"model": "constant", "reference_density_kg_m3": 3500.0, "bulk_modulus_pa": 2.0e11},
    "shear_modulus": {"model": "constant", "shear_modulus_pa": 6.0e10}}}


def cpl_world(name="moon"):
    """A two-layer terrestrial world with a constant-phase-lag tide model."""
    return build_world({
        "name": name,
        "type": "terrestrial",
        "radius_m": 1.6e6,
        "mass_kg": 8.9e22,
        "tides": {
            "global_tidal_model": "cpl",
            "max_degree_l": 2,
            "eccentricity_trunc_lvl": 2,
            "obliquity_trunc_lvl": "off",
            "fixed_k": [0.3],
            "fixed_q": [50.0],
        },
        "layers": {
            "core": {"radius_fraction": 0.5, "material": LAYER_MATERIAL},
            "mantle": {"radius_fraction": 1.0, "material": LAYER_MATERIAL},
        },
    })


def system_with(world, spin_frequency=None):
    """``world`` orbiting a Jupiter-mass host in a System, synchronous unless a spin is given."""
    system = System("pair")
    system.add_world(BaseWorld("host", 7.0e7, HOST_MASS))
    system.add_world(world, tidal_host="host", semi_major_axis=SEMI_MAJOR_AXIS, eccentricity=ECCENTRICITY)
    world.set_spin_frequency(system.calc_orbital_frequency(world) if spin_frequency is None else spin_frequency)
    return system


# =====================================================================================================================
# Optional orbital state (W1) and the calc_tides result (W2)
# =====================================================================================================================
def test_calc_tides_returns_the_named_results():
    world = cpl_world()
    result = world.calc_tides(**ORBIT)
    assert result["tidal_heating"] == world.get_tidal_heating()
    assert {key: result[key] for key in ("dU_dM", "dU_dw", "dU_dO")} == world.get_tidal_potential_derivatives()
    assert result["num_tidal_modes"] == world.get_num_tidal_modes()
    assert result["layer_tidal_heating"] == {"core": world.get_layer_tidal_heating(0),
                                             "mantle": world.get_layer_tidal_heating(1)}
    assert sum(result["layer_tidal_heating"].values()) == pytest.approx(result["tidal_heating"], rel=1.0e-12)


def test_potential_derivatives_are_named():
    world = cpl_world()
    assert all(math.isnan(value) for value in world.get_tidal_potential_derivatives().values())
    world.calc_tides(**ORBIT)
    derivatives = world.get_tidal_potential_derivatives()
    assert list(derivatives) == ["dU_dM", "dU_dw", "dU_dO"]
    assert derivatives["dU_dM"] != 0.0


def test_missing_state_outside_a_system_names_the_arguments():
    world = cpl_world()
    with pytest.raises(ValueError, match="calc_tides needs semi_major_axis, host_mass.*System"):
        world.calc_tides(ORBITAL_FREQUENCY, ORBITAL_FREQUENCY, ECCENTRICITY, 0.0)
    with pytest.raises(ValueError, match="needs orbital_frequency, spin_frequency"):
        world.calc_tides()


def test_missing_state_comes_from_the_system():
    world = cpl_world()
    # The system supplies the state only while it exists.
    system = system_with(world)
    state = world.get_tide_state()
    from_system = world.calc_tides()
    assert from_system == world.calc_tides(**state)
    # One value replaced, the rest from the system.
    changed = world.calc_tides(eccentricity=2.0 * ECCENTRICITY)
    assert changed["tidal_heating"] == pytest.approx(4.0 * from_system["tidal_heating"], rel=1.0e-9)


def test_3d_methods_take_the_state_from_the_system(io):
    world = io.copy()
    world.solve_eos()
    system = system_with(world)
    state = world.get_tide_state()
    radii = np.array([0.5, 0.8, 0.95]) * world.radius
    colatitudes = np.full(3, 0.7)
    np.testing.assert_array_equal(
        world.get_3d_tidal_heating_array(radii=radii, colatitudes=colatitudes),
        world.get_3d_tidal_heating_array(**state, radii=radii, colatitudes=colatitudes))
    assert world.get_3d_tidal_heating(radius=radii[1], colatitude=0.7) == \
        world.get_3d_tidal_heating(**state, radius=radii[1], colatitude=0.7)
    with pytest.raises(TypeError, match="needs radius"):
        world.get_3d_tidal_heating(colatitude=0.7)


# =====================================================================================================================
# Near-synchronous warning (W3)
# =====================================================================================================================
def test_near_synchronous_spin_warns_once(spdlog_text):
    world = cpl_world("near_sync")
    world.calc_tides(**{**ORBIT, "spin_frequency": ORBITAL_FREQUENCY * (1.0 + 2.0e-4)})
    world.calc_tides(**{**ORBIT, "spin_frequency": ORBITAL_FREQUENCY * (1.0 + 3.0e-4)})
    text = spdlog_text()
    assert text.count("world 'near_sync' spins at") == 1
    assert "set_synchronous_rotation" in text


@pytest.mark.parametrize("spin_ratio", [1.0, 1.01, 0.5])
def test_synchronous_or_clearly_asynchronous_spin_does_not_warn(spdlog_text, spin_ratio):
    world = cpl_world("quiet")
    world.calc_tides(**{**ORBIT, "spin_frequency": ORBITAL_FREQUENCY * spin_ratio})
    assert "spins at" not in spdlog_text()


def test_system_evolution_warns_too(spdlog_text):
    world = cpl_world("evolving")
    system = system_with(world)
    world.set_spin_frequency(world.spin_frequency * (1.0 + 5.0e-4))
    system.calc_world_evolution(world)
    assert "world 'evolving' spins at" in spdlog_text()


# =====================================================================================================================
# Tide configuration round trip and context manager (W7)
# =====================================================================================================================
def test_set_tide_config_takes_what_get_tide_config_returns():
    world = cpl_world()
    world.set_tide_config(love_fixed_dt=600.0, love_fixed_q=80.0, eccentricity_truncation="exact")
    saved = world.get_tide_config()
    assert saved["eccentricity_trunc_lvl"] == "exact"
    world.set_tide_config(eccentricity_truncation=4, love_fixed_dt=10.0)
    world.set_tide_config(**saved)
    assert world.get_tide_config() == saved
    world.set_tide_config(eccentricity_trunc_lvl=6, obliquity_trunc_lvl=2, love_fixed_dt_s=30.0)
    config = world.get_tide_config()
    assert (config["eccentricity_trunc_lvl"], config["obliquity_trunc_lvl"], config["love_fixed_dt_s"]) == (6, 2, 30.0)


@pytest.mark.parametrize("both", [
    {"eccentricity_truncation": 4, "eccentricity_trunc_lvl": 4},
    {"obliquity_truncation": 2, "obliquity_trunc_lvl": 2},
    {"love_fixed_dt": 1.0, "love_fixed_dt_s": 1.0},
])
def test_both_spellings_of_one_setting_raise(both):
    with pytest.raises(ValueError, match="two names for one setting"):
        cpl_world().set_tide_config(**both)


def test_temporary_tide_config_restores_the_previous_one():
    world = cpl_world()
    before = world.get_tide_config()
    heating_before = world.calc_tides(**ORBIT)["tidal_heating"]
    with world.temporary_tide_config(eccentricity_truncation=10, love_fixed_q=20.0) as inside:
        assert inside is world
        assert world.get_tide_config()["eccentricity_trunc_lvl"] == 10
        heating_at_10 = world.calc_tides(**ORBIT)["tidal_heating"]
    assert world.get_tide_config() == before
    assert world.calc_tides(**ORBIT)["tidal_heating"] == heating_before
    level_10_world = cpl_world()
    level_10_world.set_tide_config(eccentricity_truncation=10)
    assert heating_at_10 == level_10_world.calc_tides(**ORBIT)["tidal_heating"]


def test_temporary_tide_config_restores_after_an_exception():
    world = cpl_world()
    before = world.get_tide_config()
    with pytest.raises(RuntimeError, match="inside"):
        with world.temporary_tide_config(max_degree_l=3):
            raise RuntimeError("inside")
    assert world.get_tide_config() == before


# =====================================================================================================================
# Tide model access (W8)
# =====================================================================================================================
def test_tide_model_property():
    world = cpl_world()
    model = world.tide_model
    assert isinstance(model, FixedQTide)
    assert model.get_config_dict() == make_tide("cpl", fixed_k=[0.3], fixed_q=[50.0]).get_config_dict()
    world.tide_model = None
    assert world.tide_model is None
    assert not world.tide_model_set


def test_set_tide_model_takes_a_name_or_a_table():
    world = cpl_world()
    world.set_tide_model("fixed_dt")
    assert world.tide_model.model_name == "fixed_dt"
    world.tide_model = {"model": "ctl", "fixed_dt_s": [600.0], "fixed_k": [0.25]}
    assert world.tide_model.get_parameter("fixed_dt_s") == [600.0]
    # A world file's [tides] table, its settings included, round-trips through set_tide_model.
    tides = cpl_world().get_config_dict()["tides"]
    tides["eccentricity_trunc_lvl"] = 4
    world.set_tide_model(tides)
    assert world.tide_model.model_name == "fixed_q"
    assert world.get_tide_config()["eccentricity_trunc_lvl"] == 4


def test_set_tide_model_refuses_bad_tables():
    world = cpl_world()
    with pytest.raises(ValueError, match="names its model"):
        world.set_tide_model({"fixed_q": [10.0]})
    with pytest.raises(ValueError, match="no key 'fixed_qq'.*fixed_q"):
        world.set_tide_model({"model": "cpl", "fixed_qq": [10.0]})
    with pytest.raises(TypeError, match="tide model"):
        world.set_tide_model(3.0)
    assert world.tide_model.model_name == "fixed_q"


def test_make_tide_takes_keyword_parameters():
    by_config = make_tide("fixed_q", {"fixed_k": [0.2], "fixed_q": [40.0]})
    by_keyword = make_tide("cpl", fixed_k=[0.2], fixed_q=[40.0])
    merged = make_tide("cpl", {"fixed_k": [0.2], "fixed_q": [10.0]}, fixed_q=[40.0])
    assert by_config.get_config_dict() == by_keyword.get_config_dict() == merged.get_config_dict()
