"""A world with no layers saves and rebuilds when its tide model is analytic; one on the rheology model needs layers."""
import pytest

from TidalPy.Structures import build_world
from TidalPy.Structures.worlds.gasgiant import GasGiantWorld
from TidalPy.Tides.classes.tide import make_tide

JUPITER_RADIUS = 6.99e7
JUPITER_MASS = 1.898e27


def _layerless_gas_giant():
    world = GasGiantWorld("gg", JUPITER_RADIUS, JUPITER_MASS)
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.38], "fixed_q": [3.6e4]}))
    return world


def test_a_layerless_gas_giant_rebuilds_from_its_config():
    world = _layerless_gas_giant()
    config = world.get_config_dict()
    rebuilt = build_world(config)
    assert len(rebuilt) == 0
    assert rebuilt.radius == JUPITER_RADIUS
    assert rebuilt.get_config_dict()["tides"] == config["tides"]


def test_a_layerless_gas_giant_saves_and_reloads(tmp_path):
    path = tmp_path / "gg.toml"
    _layerless_gas_giant().save_to_toml(str(path))
    assert len(build_world(str(path))) == 0


def test_a_layerless_world_on_the_rheology_model_is_refused():
    config = {"schema_version": "0.2.0", "name": "rocky", "type": "terrestrial", "radius_m": 1.0e6, "mass_kg": 1.0e22}
    with pytest.raises(ValueError, match="'rheology' tide model requires at least one"):
        build_world(config)
