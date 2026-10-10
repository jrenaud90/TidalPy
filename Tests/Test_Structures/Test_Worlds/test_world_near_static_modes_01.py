"""Near-static tidal modes in a viscously relaxed solid mantle (``[numerical] minimum_complex_rigidity``).

A Maxwell solid forced near zero frequency has a complex shear modulus near i omega eta, which in a hot
mantle is tens of pascals; the solid radial equations divide by it. A Love solve raises a solid zone's complex modulus
to ``minimum_complex_rigidity`` times rho g R through its real part, keeping its imaginary (dissipative) part. The hot
state below is earth_thermal 9 kyr into an evolution at 0.05 au from a Sun-like star (a hot upper mantle with 2.5
percent melt under a cold lid), near a spin within 4e-10 of 2:1. The packaged floor still binds there (it moves
-Im(k2) by about 3e-3 at 1e-11 rad/s), and raising only the elastic part keeps the dissipation.
"""
import numpy as np
import pytest

import TidalPy
from TidalPy.constants import au, update_constants, year
from TidalPy.Structures import System
from TidalPy.Structures.configs import build_world


STATE_TIME = 0.00899433058579124e6 * year                                    # [s]
STATE_SEMI_MAJOR_AXIS = 0.05 * au * 0.9999904975663465                       # [m]
STATE_ECCENTRICITY = 0.29997949212106984
STATE_TEMPERATURES = (5500.000085285945, 4499.989805463227, 1854.5679085692806)   # [K] inner core, outer core, mantle
SUN_SPIN = 2.0 * np.pi / (25.38 * 86400.0)                                   # [rad s-1]


@pytest.fixture(scope="module")
def hot_pair():
    """earth_thermal hot about the Sun, its EOS solved with its temperature profile; returns (system, world, n)."""
    world = build_world("earth_thermal", {"radial_solver": {"rtol": 1.0e-10, "atol": 1.0e-14}})
    for layer, temperature in zip(world, STATE_TEMPERATURES):
        layer.use_heating = True
        layer.temperature = temperature
    sun = build_world("sol")
    sun.set_spin_frequency(SUN_SPIN)
    system = System("hot_pair")
    system.add_world(
        sun,
        is_star=True)
    system.add_world(
        world,
        tidal_host=sun,
        semi_major_axis=STATE_SEMI_MAJOR_AXIS,
        eccentricity=STATE_ECCENTRICITY)
    system.set_tidal_host(sun, world)
    world.solve_eos(
        solve_temperature=True,
        surface_temperature=system.calc_equilibrium_temperature(world),
        time=STATE_TIME)
    return system, world, system.calc_orbital_frequency(world)


@pytest.fixture
def complex_rigidity_floor():
    """Sets [numerical] minimum_complex_rigidity for one test and restores it."""
    numerical = TidalPy.config["numerical"]
    default = numerical["minimum_complex_rigidity"]

    def set_floor(value):
        numerical["minimum_complex_rigidity"] = value
        update_constants()

    yield set_floor
    set_floor(default)


@pytest.mark.parametrize("offset", [3.0e-10, 3.5e-10, 4.0e-10])
def test_near_static_mode_in_a_relaxed_mantle_solves(hot_pair, offset):
    system, world, orbital_frequency = hot_pair
    world.set_spin_frequency((2.0 + offset) * orbital_frequency)
    pair = system.calc_pair_evolution(world)
    for key in ("da_dt", "de_dt", "dn_dt"):
        assert np.isfinite(pair[key])
    world_part = pair["worlds"][world.name]
    assert np.isfinite(world_part["dspin_dt"])
    assert world_part["tidal_heating"] > 0.0


def test_floor_leaves_an_ordinary_solve_unchanged(complex_rigidity_floor):
    """At the M2 period no zone of earth_thermal comes near the floor, so the Love numbers are the floorless ones."""
    world = build_world("earth_thermal")
    world.solve_eos(solve_temperature=True, surface_temperature=288.0)
    frequency = 2.0 * np.pi / (12.4206 * 3600.0)
    love = []
    for floor in (0.0, TidalPy.config["numerical"]["minimum_complex_rigidity"]):
        complex_rigidity_floor(floor)
        world.solve_love_numbers(frequency=frequency)
        love.append(world.get_love_number_k())
    assert love[1] == love[0]


def test_floor_keeps_the_dissipation_of_a_relaxed_mantle(hot_pair, complex_rigidity_floor):
    """Toward a static frequency the hot mantle's upper part relaxes (|mu| near omega eta, below the floor): the floor
    raises only its elastic part, so k2 matches the floorless solve and -Im(k2) still falls toward zero with omega.

    The floorless solve is the reference only where it is well posed: at 1e-13 rad/s its Kamata start fails and its
    Takeuchi start gives Re(k2) differing by 1e-5 between platforms, while the floored solve agreed to 2e-7."""
    _, world, _ = hot_pair
    frequencies = (1.0e-11, 1.0e-12, 1.0e-13)   # [rad s-1]
    love = {}
    for floor in (0.0, TidalPy.config["numerical"]["minimum_complex_rigidity"]):
        complex_rigidity_floor(floor)
        love[floor] = []
        for frequency in frequencies:
            world.solve_love_numbers(frequency=frequency)
            love[floor].append(world.get_love_number_k())
    floorless, floored = np.array(love[0.0]), np.array(list(love.values())[-1])
    np.testing.assert_allclose(floored.imag[:2], floorless.imag[:2], rtol=1.0e-2)
    np.testing.assert_allclose(floored.real[:2], floorless.real[:2], rtol=1.0e-5)
    dissipation = -floored.imag
    assert dissipation[0] > dissipation[1] > dissipation[2] > 0.0
    assert dissipation[2] < 0.2 * dissipation[0]
