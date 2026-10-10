"""Propagation-matrix Love solve on BaseWorld: analytic homogeneous k2 and clean rejection of unsupported worlds."""

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
_SHEAR_VISC    = 1.0e23    # [Pa s]; very high so the layer is in the elastic limit
_FREQ          = 1.0e-5    # [rad/s]


def _material(density, bulk_viscosity=None):
    """Constant-density elastic-limit solid; no bulk viscosity law when ``bulk_viscosity`` is None."""
    return Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": density, "bulk_modulus_pa": _STATIC_BULK},
        shear_modulus={"model": "constant", "shear_modulus_pa": _STATIC_SHEAR},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": _SHEAR_VISC},
        bulk_viscosity=None if bulk_viscosity is None else {
            "model": "constant", "reference_viscosity_pas": bulk_viscosity}))


def _maxwell_layer(
        name,
        layer_index,
        radius_bounds,
        density,
        mass=0.0,
):
    radius_inner, radius_outer = radius_bounds
    return Layer(
        name,
        layer_index,
        radius_inner,
        radius_outer,
        mass,
        _material(density),
        shear_rheology=Maxwell(),
    )


def _incompressible_solid_world():
    """Single solid, static, incompressible uniform sphere."""
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = BaseWorld("incompressible_planet", _PLANET_RADIUS, mass)
    layer = Layer(
        "mantle",
        0,
        0.0,
        _PLANET_RADIUS,
        mass,
        _material(_DENSITY, bulk_viscosity=1.0e30),
        shear_rheology=Maxwell(),
        bulk_rheology=Maxwell(),
    )
    layer.is_incompressible = True
    assert layer.is_incompressible is True
    world.add_layer(layer)
    return world


def _two_layer_incompressible_world():
    r_core = 3.0e6
    mass = (4.0 / 3.0) * math.pi * (
        8000.0 * r_core ** 3 + 3300.0 * (_PLANET_RADIUS ** 3 - r_core ** 3))
    world = BaseWorld("two_layer", _PLANET_RADIUS, mass)
    for name, layer_index, radius_bounds, density in (("core", 0, (0.0, r_core), 8000.0),
                                                      ("mantle", 1, (r_core, _PLANET_RADIUS), 3300.0)):
        layer = _maxwell_layer(name, layer_index, radius_bounds, density)
        layer.is_incompressible = True
        world.add_layer(layer)
    return world


def _compressible_world():
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = BaseWorld("compressible_planet", _PLANET_RADIUS, mass)
    layer = _maxwell_layer(
        "mantle",
        0,
        (0.0, _PLANET_RADIUS),
        _DENSITY,
        mass=mass,
    )
    assert layer.is_incompressible is False
    world.add_layer(layer)
    return world


def _analytic_k2():
    """Homogeneous incompressible sphere: k2 = 1.5 / (1 + 19 mu / (2 rho g R))."""
    g_surface = (4.0 / 3.0) * math.pi * G * _DENSITY * _PLANET_RADIUS
    mu_tilde  = 19.0 * _STATIC_SHEAR / (2.0 * _DENSITY * g_surface * _PLANET_RADIUS)
    return 1.5 / (1.0 + mu_tilde)


def test_prop_matrix_k2_matches_analytic():
    """The propagation matrix solves and reproduces the analytic k2 with a negligible imaginary part."""
    world = _incompressible_solid_world()
    world.solve_eos(G_to_use=G, verbose=False)
    world.solve_love_numbers(frequency=_FREQ, love_method='propagation_matrix', verbose=False)
    assert world.love_solved is True
    k2 = world.love_number_k
    analytic = _analytic_k2()
    assert not cmath.isnan(k2)
    assert k2.real == pytest.approx(analytic, rel=0.02)
    assert abs(k2.imag) < 0.05 * abs(k2.real)


def test_prop_matrix_close_to_shooting_for_same_world():
    """The propagation matrix and shooting methods agree on the same sphere."""
    world = _incompressible_solid_world()
    world.solve_eos(G_to_use=G, verbose=False)
    world.solve_love_numbers(frequency=_FREQ, love_method='propagation_matrix', verbose=False)
    assert world.love_solved
    k2_matrix = world.love_number_k
    # Shooting cannot start a static incompressible solid, so it solves the dynamic form (equal at this frequency).
    world.mantle.is_static = False
    world.solve_love_numbers(frequency=_FREQ, love_method='radial_solver', starting_method="kamata", verbose=False)
    assert world.love_solved
    k2_shoot = world.love_number_k
    assert k2_matrix.real == pytest.approx(k2_shoot.real, rel=1.0e-3)


@pytest.mark.parametrize(
    "make_world",
    [_two_layer_incompressible_world, _compressible_world],
    ids=["two_layers", "compressible_layer"],
)
def test_prop_matrix_rejects_unsupported_world(make_world):
    """An unsupported world fails the propagation-matrix solve cleanly."""
    world = make_world()
    world.solve_eos(G_to_use=G, verbose=False)
    world.solve_love_numbers(frequency=_FREQ, love_method='propagation_matrix', verbose=False)
    assert world.love_success is False
    assert world.love_error_code != 0
