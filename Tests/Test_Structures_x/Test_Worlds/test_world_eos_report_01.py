"""The ``solve_eos`` result is copied out of the world under its call lock, in the same call as the solve.

Two threads solving one world each get a result dict that describes one completed solve: the profile arrays, the
per-layer lists, and the scalars all agree with one another, never a mix of two solves or a solution being freed.
"""
import math
import threading

import numpy as np

from TidalPy.constants import G
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS

_R          = 1.6e6               # [m]
_R_CORE     = 0.5 * _R            # [m]
_RHO_CORE   = 6000.0              # [kg m-3]
_RHO_MANTLE = 3300.0              # [kg m-3]
_SOLVES     = 25                  # per thread


def _two_layer_world():
    world = LayeredWorld("report", _R, (4.0 / 3.0) * math.pi * _R ** 3 * 4000.0)
    core = PhysicsLayer("core", 0, 0.0, _R_CORE, 0.0)
    core.set_eos(ConstantDensityEOS(reference_density=_RHO_CORE))
    mantle = PhysicsLayer("mantle", 1, _R_CORE, _R, 0.0)
    mantle.set_eos(ConstantDensityEOS(reference_density=_RHO_MANTLE))
    world.add_layer(core)
    world.add_layer(mantle)
    return world


def _expected_mass():
    return (4.0 / 3.0) * math.pi * (_RHO_CORE * _R_CORE ** 3 + _RHO_MANTLE * (_R ** 3 - _R_CORE ** 3))


def _check(result, slices):
    assert result["success"] is True, result["message"]
    n = len(result["radius"])
    assert n > 0
    for key in ("gravity", "pressure", "mass", "moi", "density", "temperature", "heat_flow"):
        assert len(result[key]) == n, key
    assert np.all(np.isfinite(result["density"]))
    assert math.isclose(result["planet_mass"], _expected_mass(), rel_tol=1.0e-8)
    assert math.isclose(result["mass"][-1], result["planet_mass"], rel_tol=1.0e-8)
    assert result["layer_radius_outer"] == [_R_CORE, _R]
    assert len(result["layer_temperature"]) == 2
    assert len(result["layer_temperature_rate"]) == 2
    assert n == 2 * slices


def test_concurrent_solves_each_return_one_consistent_result():
    world = _two_layer_world()
    errors = []

    def solve(slices):
        try:
            for _ in range(_SOLVES):
                _check(world.solve_eos(G_to_use=G, slices_per_layer=slices), slices)
        except BaseException as error:   # reported on the main thread
            errors.append(error)

    # The two threads ask for different grids, so a result that mixed the two solves would not line up.
    threads = [threading.Thread(target=solve, args=(slices,)) for slices in (40, 90)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert not errors, errors
    assert world.eos_solved


def test_the_result_can_be_rebuilt_from_the_world():
    world = _two_layer_world()
    result = world.solve_eos(G_to_use=G)
    rebuilt = world._build_eos_result()
    for key, value in result.items():
        if isinstance(value, np.ndarray):
            np.testing.assert_array_equal(rebuilt[key], value)
        else:
            assert rebuilt[key] == value, key


def test_an_unsolved_world_reports_empty_profiles():
    world = _two_layer_world()
    result = world._build_eos_result()
    assert result["success"] is False
    assert len(result["radius"]) == 0
    assert result["layer_temperature"] == []
