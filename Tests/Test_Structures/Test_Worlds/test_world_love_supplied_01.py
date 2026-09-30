"""LayeredWorld.solve_love_numbers_supplied reproduces the rheology-driven solve when fed the rheology moduli."""

import cmath
import math

import numpy as np
import pytest

from TidalPy.constants import G


_PLANET_RADIUS = 6.0e6
_DENSITY       = 4000.0
_STATIC_SHEAR  = 6.0e10
_STATIC_BULK   = 1.3e11
_SHEAR_VISC    = 1.0e21
_FREQ          = 1.0e-5


def _maxwell_world():
    from TidalPy.Structures.worlds.layered import LayeredWorld
    from TidalPy.Structures.layers.physics import PhysicsLayer
    from TidalPy.Material.eos.material_eos import ConstantDensityEOS
    from TidalPy.Viscosity import make_viscosity
    from TidalPy.Rheology.rheology import Maxwell

    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("supplied_planet", _PLANET_RADIUS, mass)
    layer = PhysicsLayer(
        "mantle",
        0,
        0.0,
        _PLANET_RADIUS,
        mass,
    )
    layer.set_eos(ConstantDensityEOS(
        reference_density=_DENSITY, shear_modulus_static=_STATIC_SHEAR, bulk_modulus_static=_STATIC_BULK))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _SHEAR_VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e30}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Maxwell())
    world.add_layer(layer)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    return world


def _solve_with_rheology_moduli(world, **solver_kwargs):
    """Sample the rheology moduli on the EOS grid and solve through the supplied path."""
    eos = world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    radius = np.ascontiguousarray(eos["radius"], dtype=np.float64)
    shear = np.ascontiguousarray(world.calc_complex_shear_modulus(radius, _FREQ), dtype=np.complex128)
    bulk = np.ascontiguousarray(world.calc_complex_bulk_modulus(radius, _FREQ), dtype=np.complex128)
    return world.solve_love_numbers_supplied(shear, bulk, radius, frequency=_FREQ, **solver_kwargs)


def test_supplied_matches_rheology():
    """The supplied path fed the rheology moduli matches the rheology-driven k, h, l."""
    world = _maxwell_world()
    world.solve_love_numbers(frequency=_FREQ, verbose=False)
    k_ref = world.love_number_k
    h_ref = world.love_number_h
    l_ref = world.love_number_l

    res = _solve_with_rheology_moduli(world)
    assert res["success"] is True
    assert cmath.isclose(res["love_number_k"], k_ref, rel_tol=1e-6, abs_tol=1e-9)
    assert cmath.isclose(res["love_number_h"], h_ref, rel_tol=1e-6, abs_tol=1e-9)
    assert cmath.isclose(res["love_number_l"], l_ref, rel_tol=1e-6, abs_tol=1e-9)


def test_supplied_takes_the_pinned_radial_solver_settings():
    """The supplied path starts from the world's pinned [radial_solver] keys, as solve_love_numbers does."""
    loose = {"rtol": 1.0e-3, "atol": 1.0e-5}
    pinned = _maxwell_world()
    pinned.set_solver_defaults(radial_solver=loose)
    k_pinned = _solve_with_rheology_moduli(pinned)["love_number_k"]
    k_explicit = _solve_with_rheology_moduli(_maxwell_world(), **loose)["love_number_k"]
    k_default = _solve_with_rheology_moduli(_maxwell_world())["love_number_k"]
    assert k_pinned == k_explicit
    assert k_pinned != k_default


def test_supplied_k2_reasonable():
    """A supplied solve without a prior rheology solve gives 0 < Re(k2) < 1.5."""
    world = _maxwell_world()
    res = _solve_with_rheology_moduli(world)
    k2 = res["love_number_k"]
    assert not cmath.isnan(k2)
    assert 0.0 < k2.real < 1.5


def test_supplied_length_mismatch_raises():
    """Supplied arrays of different lengths raise a ValueError."""
    world = _maxwell_world()
    radius = np.linspace(0.0, _PLANET_RADIUS, 50, dtype=np.float64)
    shear = np.ones(50, dtype=np.complex128) * _STATIC_SHEAR
    bulk = np.ones(49, dtype=np.complex128) * _STATIC_BULK
    with pytest.raises(ValueError):
        world.solve_love_numbers_supplied(shear, bulk, radius, frequency=_FREQ)
