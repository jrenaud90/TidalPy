"""Threads sharing one world or system get the serial answers (per-world lock with the GIL released)."""
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from TidalPy.Structures import build_world, build_system

_IO_N = 4.11e-5
_STATE = dict(orbital_frequency=_IO_N, spin_frequency=_IO_N, eccentricity=0.0041, obliquity=0.0,
              semi_major_axis=4.217e8, host_mass=1.898e27)


def _heating(world):
    world.calc_tides(**_STATE)
    total = world.get_tidal_heating()
    grid = world.calc_3d_tides(**_STATE, latitude_summed=True, longitude_summed=True, radial_summed=True)["total"]
    return total, grid


def test_threads_sharing_one_world_get_the_serial_answers():
    """Concurrent calc_tides and calc_3d_tides on one world match the serial result."""
    io = build_world("io")
    io.solve_eos()
    serial = _heating(io)
    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(lambda _: _heating(io), range(16)))
    for total, grid in results:
        assert total == serial[0]
        assert grid == serial[1]


def test_system_evolution_releases_the_gil_and_agrees_with_serial():
    """Concurrent calc_system_evolution calls match the serial result."""
    system = build_system("sol_system")
    system.earth.solve_eos()       # its rheology tides need the solved interior
    serial = system.calc_system_evolution()
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(lambda _: system.calc_system_evolution(), range(4)))
    for result in results:
        for row, expected in zip(result, serial):
            assert row["evolved"] == expected["evolved"]
            np.testing.assert_equal(row["da_dt"], expected["da_dt"])
