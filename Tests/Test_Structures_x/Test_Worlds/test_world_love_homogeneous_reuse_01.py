"""The homogeneous Love methods inside calc_tides reuse work within one call and never across calls.

A homogeneous-sphere Love solve averages the complex shear modulus over the tidal layers. Only the rheology in
that average depends on the forcing frequency: the static modulus and viscosity at each quadrature node do not,
and the average does not depend on the harmonic degree at all. So calc_tides reads the frequency-independent
node values once and forms the averaged modulus once per unique forcing frequency, instead of repeating the
whole integral for every unique (degree, frequency) pair. These tests pin the two properties that make that
safe: every mode still gets exactly the Love numbers a direct solve at its own frequency and degree gives, and
nothing is carried from one call to the next, so a changed interior is always seen.
"""
import importlib
import math

import pytest

from TidalPy.constants import G
from TidalPy.structures_x.configs import build_world
from TidalPy.Tides_x.classes import make_tide
from TidalPy.viscosity_x import make_viscosity

global_potential = importlib.import_module("TidalPy.Tides_x.potential.global").global_potential

_HOST_MASS = 1.898e27
_SEMI_MAJOR_AXIS = 4.217e8
_ECCENTRICITY = 0.0041
_ORBITAL_FREQUENCY = math.sqrt(6.674e-11 * _HOST_MASS / _SEMI_MAJOR_AXIS ** 3)
_ORBIT = (_ORBITAL_FREQUENCY, _ORBITAL_FREQUENCY, _ECCENTRICITY, 0.0, _SEMI_MAJOR_AXIS, _HOST_MASS)


def _io(love_method, eccentricity_truncation=10, max_degree_l=4):
    """The bundled Io: two tidal layers with different moduli and viscosities, so the volume average is real."""
    config = build_world("io").get_config_dict()
    config.setdefault("tides", {}).update({
        "love_method": love_method,
        "eccentricity_trunc_lvl": eccentricity_truncation,
        "max_degree_l": max_degree_l,
        "love_fixed_q": 50.0,
        "love_fixed_dt_s": 300.0,
    })
    world = build_world(config)
    world.solve_eos()
    world.set_tide_model(make_tide("rheology"))
    world.set_spin_frequency(_ORBITAL_FREQUENCY)
    return world


def _modes(world, eccentricity_truncation, max_degree_l):
    """Every active (l, m, p, q) mode of the orbit, paired with its forcing frequency magnitude."""
    mode_map = global_potential(
        world.radius,
        *_ORBIT,
        G,
        2,
        max_degree_l,
        eccentricity_truncation,
        0)[0]
    return [(key, abs(value[0])) for key, value in mode_map]


@pytest.mark.parametrize("love_method", ("homogeneous", "cpl", "ctl"))
def test_every_mode_gets_the_love_numbers_of_a_direct_solve(love_method):
    """Degrees 2 to 4 at e^10 on a two-tidal-layer interior: 123 modes sharing 9 forcing frequencies."""
    world = _io(love_method)
    world.calc_tides(*_ORBIT)
    modes = _modes(world, 10, 4)
    from_tides = {key: complex(world.get_tidal_love_k(*key)) for key, _ in modes}

    for key, frequency in modes:
        world.solve_love_numbers(frequency=frequency, degree_l=key[0], love_method=love_method)
        assert complex(world.love_number_k) == pytest.approx(from_tides[key], rel=1e-12), key
    assert len(modes) == world.get_num_tidal_modes()


def test_a_changed_interior_is_seen_by_the_next_call():
    """Nothing survives between calc_tides calls, so a stiffer tidal layer changes the answer at once."""
    stiffer = make_viscosity("constant", {"reference_viscosity_pas": 3.4239e14})

    world = _io("homogeneous", eccentricity_truncation=5, max_degree_l=2)
    world.calc_tides(*_ORBIT)
    before = world.get_tidal_heating()

    world.asthenosphere.set_shear_viscosity(stiffer)
    world.solve_eos()
    world.calc_tides(*_ORBIT)
    after = world.get_tidal_heating()

    reference = _io("homogeneous", eccentricity_truncation=5, max_degree_l=2)
    reference.asthenosphere.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 3.4239e14}))
    reference.solve_eos()
    reference.calc_tides(*_ORBIT)

    assert after != before
    assert after == pytest.approx(reference.get_tidal_heating(), rel=1e-12)


def test_calc_tides_is_repeatable_after_other_love_solves():
    """A direct Love solve in between must not leak into the next calc_tides."""
    world = _io("homogeneous", eccentricity_truncation=10, max_degree_l=3)
    world.calc_tides(*_ORBIT)
    first = (world.get_tidal_heating(), world.get_tidal_potential_derivatives())

    world.solve_love_numbers(frequency=7.0 * _ORBITAL_FREQUENCY, degree_l=3, love_method="homogeneous")
    world.calc_tides(*_ORBIT)
    assert (world.get_tidal_heating(), world.get_tidal_potential_derivatives()) == first
