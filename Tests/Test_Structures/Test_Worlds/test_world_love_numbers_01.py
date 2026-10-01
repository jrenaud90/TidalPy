"""BaseWorld.solve_love_numbers end to end (EOS then radial solver) on one- and two-layer Maxwell worlds."""

import cmath
import math

import pytest

from TidalPy.constants import G
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.layers import Layer
from TidalPy.Material import Material, Phase
from TidalPy.Rheology.rheology import Maxwell


_PLANET_RADIUS = 6.0e6     # [m]
_DENSITY       = 4000.0    # [kg/m^3]
_STATIC_SHEAR  = 6.0e10    # [Pa]
_STATIC_BULK   = 1.3e11    # [Pa]
_SHEAR_VISC    = 1.0e21    # [Pa s]; high viscosity keeps the solve near-elastic
_FREQ          = 1.0e-5    # [rad/s]


def _material(density, shear, bulk, shear_viscosity, bulk_viscosity=None):
    """Constant-density solid with constant viscosities; no bulk viscosity law when ``bulk_viscosity`` is None."""
    return Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": density, "bulk_modulus_pa": bulk},
        shear_modulus={"model": "constant", "shear_modulus_pa": shear},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": shear_viscosity},
        bulk_viscosity=None if bulk_viscosity is None else {
            "model": "constant", "reference_viscosity_pas": bulk_viscosity}))


def _solid_world():
    """Single solid uniform Maxwell sphere."""
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = BaseWorld("solid_planet", _PLANET_RADIUS, mass)
    world.add_layer(Layer(
        "mantle",
        0,
        0.0,
        _PLANET_RADIUS,
        mass,
        _material(_DENSITY, _STATIC_SHEAR, _STATIC_BULK, _SHEAR_VISC, 1.0e30),
        shear_rheology=Maxwell(),
        bulk_rheology=Maxwell()))
    return world


def _two_layer_solid_world():
    """Stiff iron core and rocky mantle, both solid."""
    r_core = 3.0e6    # [m]
    rho_c  = 8000.0   # [kg/m^3]
    rho_m  = 3300.0   # [kg/m^3]
    mu_c   = 1.5e11   # [Pa]
    mu_m   = _STATIC_SHEAR
    K_c    = 3.0e11   # [Pa]
    K_m    = _STATIC_BULK
    mass = (4.0 / 3.0) * math.pi * (
        rho_c * r_core ** 3 + rho_m * (_PLANET_RADIUS ** 3 - r_core ** 3)
    )
    world = BaseWorld("two_layer", _PLANET_RADIUS, mass)

    core = Layer(
        "core",
        0,
        0.0,
        r_core,
        0.0,
        _material(rho_c, mu_c, K_c, 1.0e21),
        shear_rheology=Maxwell(),
    )

    mantle = Layer(
        "mantle",
        1,
        r_core,
        _PLANET_RADIUS,
        0.0,
        _material(rho_m, mu_m, K_m, _SHEAR_VISC),
        shear_rheology=Maxwell(),
    )

    world.add_layer(core)
    world.add_layer(mantle)
    return world


def _solved_single_layer():
    world = _solid_world()
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    world.solve_love_numbers(frequency=_FREQ, verbose=False)
    return world


def _solved_two_layer():
    world = _two_layer_solid_world()
    world.solve_eos(G_to_use=G, verbose=False)
    world.solve_love_numbers(frequency=_FREQ, verbose=False)
    return world


def test_love_unsolved_before_any_solve():
    """A new world reports love_solved False."""
    world = _solid_world()
    assert world.love_solved is False


def test_solve_love_requires_eos_first():
    """A Love solve before the EOS solve raises."""
    world = _solid_world()
    with pytest.raises((ValueError, RuntimeError)):
        world.solve_love_numbers(frequency=_FREQ, verbose=False)


@pytest.mark.parametrize(
    "make_solved_world, attribute, is_in_range",
    [
        pytest.param(_solved_single_layer, "love_number_k", lambda value: 0.0 < value < 1.5, id="single_k2"),
        pytest.param(_solved_single_layer, "love_number_h", lambda value: 0.0 < value < 2.5, id="single_h2"),
        pytest.param(_solved_single_layer, "love_number_l", lambda value: value >= 0.0, id="single_l2"),
        # Any solid body sits below the fluid limit k2 = 1.5.
        pytest.param(_solved_two_layer, "love_number_k", lambda value: 0.0 < value < 1.5, id="two_layer_k2"),
    ],
)
def test_love_number_is_reasonable(make_solved_world, attribute, is_in_range):
    """The solve succeeds and the real part of each Love number lies in its physical range."""
    world = make_solved_world()
    assert world.love_solved is True
    love_number = getattr(world, attribute)
    assert not cmath.isnan(love_number)
    assert is_in_range(love_number.real)


def test_love_elastic_limit_small_imag():
    """In the elastic limit |Im(k2)| is much smaller than |Re(k2)|."""
    world = _solved_single_layer()
    k2 = world.love_number_k
    # tau = eta / mu = 1.67e10 s, so omega * tau ~ 1.67e5.
    assert abs(k2.imag) < 0.1 * abs(k2.real)
