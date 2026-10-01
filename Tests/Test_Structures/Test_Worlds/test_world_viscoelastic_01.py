"""solve_eos fills each layer's pre- and post-melt moduli and viscosities, read through the world's getters and
radius-resolved complex moduli.
"""

import math

import pytest

from TidalPy.constants import G
from TidalPy.Material import Material, Phase
from TidalPy.Rheology.rheology import Maxwell
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld


_PLANET_RADIUS  = 6.0e6     # [m]
_DENSITY        = 4000.0    # [kg/m^3]
_STATIC_SHEAR   = 6.0e10    # [Pa]
_STATIC_BULK    = 1.3e11    # [Pa]
_SHEAR_VISC     = 1.0e21    # [Pa s]
_BULK_VISC      = 1.0e30    # [Pa s]
_MASS           = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
_MID_RADIUS     = 0.5 * _PLANET_RADIUS
_SOLIDUS        = 1600.0    # [K]
_LIQUIDUS       = 2000.0    # [K]


def _whole_planet_layer(material, **layer_kwargs):
    return Layer(
        "mantle",
        0,
        0.0,
        _PLANET_RADIUS,
        _MASS,
        material,
        **layer_kwargs,
    )


def _material(shear_viscosity=_SHEAR_VISC, bulk_viscosity=_BULK_VISC, with_melt=False):
    """A constant-density solid with static moduli and constant viscosities (None: no law), and an optional melt
    phase that melts between the solidus and liquidus with Henning weakening."""
    solid = Phase(
        eos={"model": "constant", "reference_density_kg_m3": _DENSITY, "bulk_modulus_pa": _STATIC_BULK},
        shear_modulus={"model": "constant", "shear_modulus_pa": _STATIC_SHEAR},
        shear_viscosity=None if shear_viscosity is None else {
            "model": "constant", "reference_viscosity_pas": shear_viscosity},
        bulk_viscosity=None if bulk_viscosity is None else {
            "model": "constant", "reference_viscosity_pas": bulk_viscosity})
    if not with_melt:
        return Material(solid=solid)
    melt = Phase(
        eos={"model": "constant", "reference_density_kg_m3": _DENSITY, "bulk_modulus_pa": _STATIC_BULK},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0})
    return Material(
        solid=solid,
        liquid=melt,
        solidus={"model": "constant", "temperature_k": _SOLIDUS},
        liquidus={"model": "constant", "temperature_k": _LIQUIDUS},
        weakening="henning")


def _uniform_physics_world(with_viscosity=True, with_melt=False, with_rheology=False):
    world = BaseWorld("rocky", _PLANET_RADIUS, _MASS)
    if with_viscosity:
        material = _material(with_melt=with_melt)
    else:
        material = _material(shear_viscosity=None, bulk_viscosity=None, with_melt=with_melt)
    rheology = {"shear_rheology": Maxwell(), "bulk_rheology": Maxwell()} if with_rheology else {}
    world.add_layer(_whole_planet_layer(material, use_melting=with_melt, **rheology))
    return world


def test_nan_before_solve():
    """Moduli and viscosities are NaN before solve_eos."""
    world = _uniform_physics_world()
    assert math.isnan(world.get_shear_modulus(_MID_RADIUS))
    assert math.isnan(world.get_shear_viscosity(_MID_RADIUS))


def test_static_moduli_recovered_no_melt():
    """Without a melt model the solved moduli are the static values and nothing melts."""
    world = _uniform_physics_world(with_viscosity=True, with_melt=False)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    assert world.get_shear_modulus(_MID_RADIUS) == pytest.approx(_STATIC_SHEAR)
    assert world.get_bulk_modulus(_MID_RADIUS) == pytest.approx(_STATIC_BULK)
    assert world.get_melt_fraction(_MID_RADIUS) == 0.0


def test_viscosity_from_models():
    """Viscosities come from the attached viscosity models."""
    world = _uniform_physics_world(with_viscosity=True)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    assert world.get_shear_viscosity(_MID_RADIUS) == pytest.approx(_SHEAR_VISC)
    assert world.get_bulk_viscosity(_MID_RADIUS) == pytest.approx(_BULK_VISC)


def test_viscosity_nan_when_model_unset():
    """Without a viscosity model the viscosity is NaN while the moduli still come from the statics."""
    world = _uniform_physics_world(with_viscosity=False)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    assert math.isnan(world.get_shear_viscosity(_MID_RADIUS))
    assert world.get_shear_modulus(_MID_RADIUS) == pytest.approx(_STATIC_SHEAR)


def test_melt_weakens_shear_at_high_temperature():
    """Between solidus and liquidus the melt fraction is set and Henning weakens the shear modulus."""
    world = _uniform_physics_world(with_viscosity=True, with_melt=True)
    # 1800 K is halfway between the 1600 K solidus and the 2000 K liquidus.
    world.solve_eos(G_to_use=G, temperature=1800.0, verbose=False)
    assert world.get_shear_modulus(_MID_RADIUS) < _STATIC_SHEAR
    assert world.get_melt_fraction(_MID_RADIUS) == pytest.approx(0.5)


def test_no_melt_below_solidus():
    """Below the solidus the shear modulus is unweakened."""
    world = _uniform_physics_world(with_viscosity=True, with_melt=True)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    assert world.get_shear_modulus(_MID_RADIUS) == pytest.approx(_STATIC_SHEAR)


def test_complex_without_rheology_is_real_static():
    """Without a rheology the complex shear modulus is the real static value."""
    world = _uniform_physics_world(with_viscosity=True, with_rheology=False)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    value = world.calc_complex_shear_modulus(_MID_RADIUS, 1.0e-5)
    assert value.real == pytest.approx(_STATIC_SHEAR)
    assert value.imag == pytest.approx(0.0)


def test_complex_shear_matches_maxwell():
    """With a Maxwell rheology the complex shear modulus is Maxwell's, with a nonzero imaginary part."""
    world = _uniform_physics_world(with_viscosity=True, with_rheology=True)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    frequency = 1.0e-5
    expected = Maxwell().calc_complex_modulus(_STATIC_SHEAR, _SHEAR_VISC, frequency)
    value = world.calc_complex_shear_modulus(_MID_RADIUS, frequency)
    assert value.real == pytest.approx(expected.real, rel=1e-9)
    assert value.imag == pytest.approx(expected.imag, rel=1e-9)
    assert abs(value.imag) > 0.0


def test_complex_modulus_without_a_rheology_is_static():
    """A layer with no rheology reports its solved static shear modulus as a real complex modulus."""
    world = BaseWorld("geom", _PLANET_RADIUS, _MASS)
    world.add_layer(_whole_planet_layer(Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": _DENSITY}))))
    world.solve_eos(G_to_use=G, verbose=False)
    value = world.calc_complex_shear_modulus(_MID_RADIUS, 1.0e-5)
    assert value.real == world.get_shear_modulus(_MID_RADIUS)
    assert value.imag == 0.0


def test_layer_getters_match_world():
    """A layer with only a shear viscosity model reports its static shear modulus and viscosity."""
    world = BaseWorld("rocky", _PLANET_RADIUS, _MASS)
    world.add_layer(_whole_planet_layer(_material(bulk_viscosity=None)))
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    assert world.get_shear_modulus(_MID_RADIUS) == pytest.approx(_STATIC_SHEAR)
    assert world.get_shear_viscosity(_MID_RADIUS) == pytest.approx(_SHEAR_VISC)
