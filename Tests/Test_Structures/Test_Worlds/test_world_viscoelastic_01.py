"""solve_eos fills each layer's pre- and post-melt moduli and viscosities, read through the world's getters and
radius-resolved complex moduli.
"""

import math

import pytest

from TidalPy.constants import G
from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.PartialMelt import make_partial_melt
from TidalPy.Rheology.rheology import Maxwell
from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Structures.layers.physics import PhysicsLayer
from TidalPy.Structures.worlds.layered import LayeredWorld
from TidalPy.Viscosity import make_viscosity


_PLANET_RADIUS  = 6.0e6     # [m]
_DENSITY        = 4000.0    # [kg/m^3]
_STATIC_SHEAR   = 6.0e10    # [Pa]
_STATIC_BULK    = 1.3e11    # [Pa]
_SHEAR_VISC     = 1.0e21    # [Pa s]
_BULK_VISC      = 1.0e30    # [Pa s]
_MASS           = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
_MID_RADIUS     = 0.5 * _PLANET_RADIUS


def _whole_planet_layer(layer_class=PhysicsLayer):
    return layer_class(
        "mantle",
        0,
        0.0,
        _PLANET_RADIUS,
        _MASS,
    )


def _uniform_physics_world(with_viscosity=True, with_melt=False, with_rheology=False):
    world = LayeredWorld("rocky", _PLANET_RADIUS, _MASS)
    layer = _whole_planet_layer()
    layer.set_eos(ConstantDensityEOS(
        reference_density=_DENSITY, shear_modulus_static=_STATIC_SHEAR, bulk_modulus_static=_STATIC_BULK))
    if with_viscosity:
        layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _SHEAR_VISC}))
        layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _BULK_VISC}))
    if with_melt:
        layer.set_partial_melt(make_partial_melt("henning", {"solidus_k": 1600.0, "liquidus_k": 2000.0}))
    if with_rheology:
        layer.set_shear_rheology(Maxwell())
        layer.set_bulk_rheology(Maxwell())
    world.add_layer(layer)
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
    # 1800 K is halfway between the 1600 K solidus and 2000 K liquidus.
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


def test_complex_nan_on_geometry_layer():
    """A geometry-only BaseLayer world gives NaN complex moduli."""
    world = LayeredWorld("geom", _PLANET_RADIUS, _MASS)
    layer = _whole_planet_layer(BaseLayer)
    layer.set_eos(ConstantDensityEOS(reference_density=_DENSITY))
    world.add_layer(layer)
    world.solve_eos(G_to_use=G, verbose=False)
    value = world.calc_complex_shear_modulus(_MID_RADIUS, 1.0e-5)
    assert math.isnan(value.real)


def test_layer_getters_match_world():
    """A layer with only a shear viscosity model reports its static shear modulus and viscosity."""
    world = LayeredWorld("rocky", _PLANET_RADIUS, _MASS)
    layer = _whole_planet_layer()
    layer.set_eos(ConstantDensityEOS(
        reference_density=_DENSITY, shear_modulus_static=_STATIC_SHEAR, bulk_modulus_static=_STATIC_BULK))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _SHEAR_VISC}))
    world.add_layer(layer)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    assert world.get_shear_modulus(_MID_RADIUS) == pytest.approx(_STATIC_SHEAR)
    assert world.get_shear_viscosity(_MID_RADIUS) == pytest.approx(_SHEAR_VISC)
