"""The passes of a thermal solve: their cap and tolerance come from [eos_solver] max_thermal_passes and thermal_tol,
a world may pin either, and a call may override both."""
from pathlib import Path

import pytest

import TidalPy
from TidalPy.Structures import build_world

_EARTH = str(Path(TidalPy.__file__).parent / "WorldPack" / "earth_simple.toml")
_SURFACE_TEMPERATURE = 300.0   # [K]


def _solve(world, **kwargs):
    result = world.solve_eos(solve_temperature=True, surface_temperature=_SURFACE_TEMPERATURE, **kwargs)
    assert result["success"], result["message"]
    return result


def test_the_defaults_come_from_the_configuration():
    eos_solver = TidalPy.config["eos_solver"]
    result = _solve(build_world(_EARTH))
    assert result["thermal_converged"]
    assert 0 < result["thermal_passes"] <= eos_solver["max_thermal_passes"]


def test_a_cap_the_passes_reach_leaves_the_solve_unconverged():
    result = _solve(build_world(_EARTH), max_thermal_passes=1)
    assert not result["thermal_converged"]
    assert result["thermal_passes"] == 1


def test_a_looser_tolerance_ends_the_passes_sooner():
    tight = _solve(build_world(_EARTH))
    loose = _solve(build_world(_EARTH), thermal_tol=1.0e-2)
    assert loose["thermal_converged"]
    assert loose["thermal_passes"] < tight["thermal_passes"]


def test_a_world_pins_its_pass_settings_and_keeps_them_through_a_binary_file(tmp_path):
    world = build_world(_EARTH)
    world.set_solver_defaults(eos_solver={"max_thermal_passes": 1, "thermal_tol": 1.0e-3})
    assert world.get_solver_defaults()["eos_solver"] == {"max_thermal_passes": 1, "thermal_tol": 1.0e-3}
    assert _solve(world)["thermal_passes"] == 1
    # A call's own argument wins over the pinned value.
    assert _solve(world, max_thermal_passes=20)["thermal_passes"] > 1

    path = str(tmp_path / "earth.tpyb")
    world.save_binary(path)
    loaded = build_world(_EARTH)
    loaded.load_binary(path)
    assert loaded.get_solver_defaults() == world.get_solver_defaults()


def test_a_world_file_pins_them():
    config = build_world(_EARTH).get_config_dict()
    config["eos_solver"] = {"max_thermal_passes": 3, "thermal_tol": 1.0e-4}
    world = build_world(config)
    assert world.get_solver_defaults()["eos_solver"] == {"max_thermal_passes": 3, "thermal_tol": 1.0e-4}
    assert world.get_config_dict()["eos_solver"] == {"max_thermal_passes": 3, "thermal_tol": 1.0e-4}


@pytest.mark.parametrize("table", [{"max_thermal_passes": 0}, {"thermal_tol": 0.0}])
def test_a_pinned_value_must_be_positive(table):
    with pytest.raises(ValueError, match="eos_solver"):
        build_world(_EARTH).set_solver_defaults(eos_solver=table)
