"""A world's binary record carries every setting that changes a result: the tide model and its ``[tides]``
configuration on every world type, the spin model's moment-of-inertia factor and the pinned ``[eos_solver]`` and
``[radial_solver]`` settings on a layered world, and the luminosity model on a star. Solved state is not saved: a
loaded world re-solves its EOS and tides.
"""
import math

import pytest

from TidalPy.dynamics_x import Spin
from TidalPy.stellar_x.luminosity import make_luminosity
from TidalPy.structures_x import build_world
from TidalPy.structures_x.worlds.gasgiant import GasGiantWorld
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.worlds.stellar import StarWorld
from TidalPy.Tides_x.classes.tide import make_tide

# Io about Jupiter.
_IO_ORBIT = dict(orbital_frequency=4.1106e-5, spin_frequency=4.1106e-5, eccentricity=0.0041, obliquity=0.0,
                 semi_major_axis=4.217e8, host_mass=1.898e27)


def _round_trip(world, target, tmp_path):
    path = str(tmp_path / "world.tpyb")
    world.save_binary(path)
    target.load_binary(path)
    return target


def _tide_model_config(world):
    return world.get_config_dict()["tides"]


# =====================================================================================================================
# Settings
# =====================================================================================================================
def test_layered_world_keeps_its_tide_spin_and_solver_settings(tmp_path):
    world = build_world("io")
    world.set_tide_model(make_tide("fixed_q", {"fixed_k": [0.3, 0.1], "fixed_q": [50.0, 80.0]}))
    world.set_tide_config(
        min_degree_l=2,
        max_degree_l=3,
        eccentricity_truncation=5,
        obliquity_truncation=2,
        layer_tidal_heating=False)
    world.set_spin_model(Spin(moment_of_inertia_factor=0.33))
    world.set_solver_defaults(
        eos_solver={"rtol": 1.0e-9, "slices_per_layer": 70, "integration_method": "RK45"},
        radial_solver={"use_kamata": True, "max_num_steps": 12345})

    loaded = _round_trip(world, LayeredWorld("placeholder", 1.0, 1.0), tmp_path)

    assert loaded.tide_model_set
    assert loaded.get_tide_config() == world.get_tide_config()
    assert _tide_model_config(loaded) == _tide_model_config(world)
    assert loaded.get_solver_defaults() == world.get_solver_defaults()
    assert loaded.get_moment_of_inertia() == pytest.approx(0.33 * world.mass * world.radius**2, rel=1e-14)


def test_gas_giant_keeps_its_tide_model(tmp_path):
    world = build_world({
        "schema_version": "0.2.0", "name": "giant", "type": "gasgiant", "radius_m": 7.0e7, "mass_kg": 1.9e27,
        "layers": {"envelope": {"class": "gas", "type": "gas", "radius_fraction": 1.0}},
        "tides": {"global_tidal_model": "fixed_dt", "fixed_k": [0.4], "fixed_dt_s": [0.5],
                  "min_degree_l": 2, "max_degree_l": 2, "eccentricity_trunc_lvl": 10}})
    loaded = _round_trip(world, GasGiantWorld("placeholder", 1.0, 1.0), tmp_path)
    assert type(loaded) is GasGiantWorld
    assert _tide_model_config(loaded) == _tide_model_config(world)


def test_star_keeps_its_tide_and_luminosity_models(tmp_path):
    star = StarWorld("sun", 6.957e8, 1.988e30)
    star.set_tide_model(make_tide("fixed_q", {"fixed_k": [0.0289], "fixed_q": [1.0e6]}))
    star.set_luminosity_model(make_luminosity("power_law", {"power_law_coeff": 1.2, "power_law_exponent": 3.9}))
    loaded = _round_trip(star, StarWorld("placeholder", 1.0, 1.0), tmp_path)
    assert _tide_model_config(loaded) == _tide_model_config(star)
    assert loaded.luminosity_model_set
    assert loaded.calc_luminosity_from_mass() == pytest.approx(star.calc_luminosity_from_mass(), rel=1e-15)


def test_a_record_without_a_tide_model_clears_the_target_model(tmp_path):
    rigid = LayeredWorld("rigid", 1.0e6, 1.0e22)
    target = LayeredWorld("placeholder", 1.0, 1.0)
    target.set_tide_model(make_tide("fixed_q", {"fixed_k": [0.3], "fixed_q": [100.0]}))
    loaded = _round_trip(rigid, target, tmp_path)
    assert not loaded.tide_model_set


# =====================================================================================================================
# Results
# =====================================================================================================================
def test_loaded_world_reproduces_its_tidal_heating_without_reattaching(tmp_path):
    world = build_world("io")
    world.solve_eos()
    world.calc_tides(**_IO_ORBIT)
    reference = world.get_tidal_heating()
    assert math.isfinite(reference) and reference > 0.0

    loaded = _round_trip(world, LayeredWorld("placeholder", 1.0, 1.0), tmp_path)
    # Solved state is not saved: the loaded world has no tide result until it solves again.
    assert not loaded.tides_solved
    loaded.solve_eos()
    loaded.calc_tides(**_IO_ORBIT)
    assert loaded.get_tidal_heating() == pytest.approx(reference, rel=1e-12)
