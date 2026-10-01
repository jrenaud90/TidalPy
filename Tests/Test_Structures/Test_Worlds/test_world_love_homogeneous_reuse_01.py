"""Homogeneous Love methods in calc_tides reuse work within one call but never across calls."""
import importlib
import math

import pytest

from TidalPy.constants import G
from TidalPy.Material import Phase
from TidalPy.Structures.configs import build_world
from TidalPy.Tides.classes import make_tide
from TidalPy.Viscosity import make_viscosity

global_potential = importlib.import_module("TidalPy.Tides.potential.global").global_potential

_HOST_MASS = 1.898e27
_SEMI_MAJOR_AXIS = 4.217e8
_ECCENTRICITY = 0.0041
_ORBITAL_FREQUENCY = math.sqrt(6.674e-11 * _HOST_MASS / _SEMI_MAJOR_AXIS ** 3)
_ORBIT = (_ORBITAL_FREQUENCY, _ORBITAL_FREQUENCY, _ECCENTRICITY, 0.0, _SEMI_MAJOR_AXIS, _HOST_MASS)


def _io(love_method, eccentricity_truncation=10, max_degree_l=4):
    """Bundled Io, whose two tidal layers differ so the volume-averaged modulus is nontrivial."""
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


def _set_shear_viscosity(layer, viscosity):
    """Give a layer's solid phase a new shear viscosity law, keeping the rest of its material."""
    solid = layer.material.solid
    layer.material = layer.material.replace(solid=Phase(
        eos=solid.eos,
        shear_modulus=solid.shear_modulus,
        shear_viscosity=viscosity,
        bulk_viscosity=solid.bulk_viscosity,
        shear_rheology=solid.shear_rheology,
        bulk_rheology=solid.bulk_rheology,
        **solid.parameters))


def _modes(world, eccentricity_truncation, max_degree_l):
    """Every active (l, m, p, q) mode paired with its forcing frequency magnitude."""
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
    """Each calc_tides mode's k matches a direct solve at its own frequency and degree."""
    world = _io(love_method)
    world.calc_tides(*_ORBIT)
    modes = _modes(world, 10, 4)
    from_tides = {key: complex(world.get_tidal_love_k(*key)) for key, _ in modes}

    for key, frequency in modes:
        world.solve_love_numbers(frequency=frequency, degree_l=key[0], love_method=love_method)
        assert complex(world.love_number_k) == pytest.approx(from_tides[key], rel=1e-12), key
    assert len(modes) == world.get_num_tidal_modes()


def test_a_changed_interior_is_seen_by_the_next_call():
    """A changed tidal-layer viscosity is seen by the next calc_tides call."""
    stiffer = make_viscosity("constant", {"reference_viscosity_pas": 3.4239e14})

    world = _io("homogeneous", eccentricity_truncation=10, max_degree_l=2)
    world.calc_tides(*_ORBIT)
    before = world.get_tidal_heating()

    _set_shear_viscosity(world.asthenosphere, stiffer)
    world.solve_eos()
    world.calc_tides(*_ORBIT)
    after = world.get_tidal_heating()

    reference = _io("homogeneous", eccentricity_truncation=10, max_degree_l=2)
    _set_shear_viscosity(reference.asthenosphere, make_viscosity("constant", {"reference_viscosity_pas": 3.4239e14}))
    reference.solve_eos()
    reference.calc_tides(*_ORBIT)

    assert after != before
    assert after == pytest.approx(reference.get_tidal_heating(), rel=1e-12)


def test_calc_tides_is_repeatable_after_other_love_solves():
    """A direct Love solve between calc_tides calls does not change the next result."""
    world = _io("homogeneous", eccentricity_truncation=10, max_degree_l=3)
    world.calc_tides(*_ORBIT)
    first = (world.get_tidal_heating(), world.get_tidal_potential_derivatives())

    world.solve_love_numbers(frequency=7.0 * _ORBITAL_FREQUENCY, degree_l=3, love_method="homogeneous")
    world.calc_tides(*_ORBIT)
    assert (world.get_tidal_heating(), world.get_tidal_potential_derivatives()) == first
