"""calc_tides spreads its Love-number solves over threads ([numerical] love_solve_threads) with identical results."""
import copy

import numpy as np
import pytest

import TidalPy
from TidalPy.constants import update_constants
from TidalPy.Structures import build_world

_IO_N = 4.11e-5
# Degrees 2 to 3, e^10, obliquity level 2, and a non-synchronous spin: 58 Love solves per calc_tides call.
_STATE = dict(orbital_frequency=_IO_N, spin_frequency=1.2 * _IO_N, eccentricity=0.05, obliquity=0.1,
              semi_major_axis=4.217e8, host_mass=1.898e27)
_TIDE_CONFIG = dict(min_degree_l=2, max_degree_l=3, eccentricity_truncation=10, obliquity_truncation=2)
_MODES = [(2, 0, 1, 1), (2, 2, 0, 0), (2, 1, 0, -1), (3, 1, 1, 0), (3, 3, 0, 2)]


@pytest.fixture
def restore_config():
    """Restore ``TidalPy.config`` and the C++ numerical settings after a test changes them."""
    original = copy.deepcopy(TidalPy.config)
    yield
    TidalPy.config = original
    update_constants()


def _set_love_threads(threads, min_parallel):
    TidalPy.config["numerical"]["love_solve_threads"] = threads
    TidalPy.config["numerical"]["love_solve_min_parallel"] = min_parallel
    update_constants()


def _tide_results(world):
    """Everything calc_tides reports: the heating, each layer's share, per-mode Love numbers, and the derivatives."""
    world.calc_tides(**_STATE)
    return (
        world.get_tidal_heating(),
        np.array([layer.get_tidal_heating() for layer in world]),
        np.array([world.get_tidal_love_k(*mode) for mode in _MODES]),
        np.array(world.get_tidal_potential_derivatives()))


@pytest.mark.parametrize("threads", [2, 16, 0])
@pytest.mark.parametrize("layer_tidal_heating", [True, False])
@pytest.mark.parametrize("love_method", ["radial_solver", "homogeneous"])
def test_threaded_love_solves_match_serial(restore_config, love_method, layer_tidal_heating, threads):
    """Any thread count gives bit-identical results to the calling thread alone."""
    io = build_world("io")
    io.solve_eos()
    io.set_tide_config(love_method=love_method, layer_tidal_heating=layer_tidal_heating, **_TIDE_CONFIG)
    _set_love_threads(1, 3)
    serial = _tide_results(io)
    _set_love_threads(threads, 1)
    threaded = _tide_results(io)
    assert threaded[0] == serial[0]
    for threaded_values, serial_values in zip(threaded[1:], serial[1:]):
        assert np.array_equal(threaded_values, serial_values, equal_nan=True)


def test_love_thread_settings_reach_the_cpp_config(restore_config):
    """The [numerical] values are what the C++ config singleton holds."""
    import TidalPy.constants as constants
    _set_love_threads(5, 12)
    assert constants.love_solve_threads == 5
    assert constants.love_solve_min_parallel == 12


def test_failed_love_solve_on_a_worker_thread_raises(restore_config):
    """A Love solve that fails on a worker thread raises on the calling thread and leaves the tides unsolved."""
    io = build_world("io")
    io.solve_eos()
    io.set_tide_config(**_TIDE_CONFIG)
    _set_love_threads(16, 1)
    # Too few integration steps to cross the first layer, so every solve fails.
    TidalPy.config["radial_solver"]["max_num_steps"] = 5
    update_constants()
    with pytest.raises(RuntimeError, match="Love-number solve failed during calc_tides"):
        io.calc_tides(**_STATE)
    assert io.tides_solved is False
