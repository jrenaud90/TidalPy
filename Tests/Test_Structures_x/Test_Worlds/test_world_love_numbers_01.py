"""
Tests for LayeredWorld.solve_love_numbers: end-to-end pipeline from EOS solve
through the radial (Love-number) solver.

A PhysicsLayer with Maxwell rheology and shear viscosity is used so that the
complex Love numbers are non-trivial.  Two world configurations are tested:

  * single solid-layer uniform sphere
  * two solid-layer planet (stiff iron core + softer rocky mantle)

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""

import cmath
import math

import pytest

from TidalPy.constants import G


_PLANET_RADIUS = 6.0e6     # [m]
_DENSITY       = 4000.0    # [kg/m^3]
_STATIC_SHEAR  = 6.0e10    # [Pa]
_STATIC_BULK   = 1.3e11    # [Pa]
_SHEAR_VISC    = 1.0e21    # [Pa s] — high viscosity => near-elastic regime
_FREQ          = 1.0e-5    # [rad/s]


def _imports():
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.rheology_x.rheology import Maxwell
    return LayeredWorld, PhysicsLayer, ConstantDensityEOS, make_viscosity, Maxwell


def _solid_world(with_rheology: bool = True):
    """Single solid uniform sphere fully equipped for love-number solving."""
    LayeredWorld, PhysicsLayer, ConstantDensityEOS, make_viscosity, Maxwell = _imports()
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("solid_planet", _PLANET_RADIUS, mass)
    layer = PhysicsLayer("mantle", 0, 0.0, _PLANET_RADIUS, mass,
                         shear_modulus_static_pa=_STATIC_SHEAR,
                         bulk_modulus_static_pa=_STATIC_BULK)
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _SHEAR_VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e30}))
    if with_rheology:
        layer.set_shear_rheology(Maxwell())
        layer.set_bulk_rheology(Maxwell())
    world.add_layer(layer)
    return world


def _two_layer_solid_world():
    """Stiff iron core + rocky mantle; both layers solid."""
    LayeredWorld, PhysicsLayer, ConstantDensityEOS, make_viscosity, Maxwell = _imports()
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
    world = LayeredWorld("two_layer", _PLANET_RADIUS, mass)

    core = PhysicsLayer("core", 0, 0.0, r_core, 0.0,
                        shear_modulus_static_pa=mu_c,
                        bulk_modulus_static_pa=K_c)
    core.set_eos(ConstantDensityEOS(reference_density_kg_m3=rho_c))
    core.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e21}))
    core.set_shear_rheology(Maxwell())

    mantle = PhysicsLayer("mantle", 1, r_core, _PLANET_RADIUS, 0.0,
                          shear_modulus_static_pa=mu_m,
                          bulk_modulus_static_pa=K_m)
    mantle.set_eos(ConstantDensityEOS(reference_density_kg_m3=rho_m))
    mantle.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _SHEAR_VISC}))
    mantle.set_shear_rheology(Maxwell())

    world.add_layer(core)
    world.add_layer(mantle)
    return world


# =====================================================================================================================
# Pre-condition checks
# =====================================================================================================================
def test_love_unsolved_before_any_solve():
    world = _solid_world()
    assert world.love_solved is False


def test_solve_love_requires_eos_first():
    world = _solid_world()
    with pytest.raises((ValueError, RuntimeError)):
        world.solve_love_numbers(frequency_rad_s=_FREQ, verbose=False)


# =====================================================================================================================
# Single solid layer — full pipeline
# =====================================================================================================================
def test_single_solid_layer_love_succeeds():
    world = _solid_world()
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    world.solve_love_numbers(frequency_rad_s=_FREQ, verbose=False)
    assert world.love_solved is True


def test_love_k2_is_reasonable():
    """k2 real part for a uniform solid sphere must be in (0, 1.5)."""
    world = _solid_world()
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    world.solve_love_numbers(frequency_rad_s=_FREQ, verbose=False)
    k2 = world.love_number_k
    assert not cmath.isnan(k2)
    assert 0.0 < k2.real < 1.5


def test_love_h2_is_reasonable():
    """h2 real part for a uniform solid sphere must be in (0, 2.5)."""
    world = _solid_world()
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    world.solve_love_numbers(frequency_rad_s=_FREQ, verbose=False)
    h2 = world.love_number_h
    assert not cmath.isnan(h2)
    assert 0.0 < h2.real < 2.5


def test_love_l2_is_nonnegative():
    """l2 (Shida number) real part for a solid sphere must be non-negative."""
    world = _solid_world()
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    world.solve_love_numbers(frequency_rad_s=_FREQ, verbose=False)
    l2 = world.love_number_l
    assert not cmath.isnan(l2)
    assert l2.real >= 0.0


def test_love_elastic_limit_small_imag():
    """At high viscosity (elastic limit) |Im(k2)| must be much smaller than |Re(k2)|."""
    world = _solid_world(with_rheology=True)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    world.solve_love_numbers(frequency_rad_s=_FREQ, verbose=False)
    k2 = world.love_number_k
    # eta = 1e21, mu = 6e10 => tau = 1.67e10 s, omega*tau ~ 1.67e5 >> 1 => elastic
    assert abs(k2.imag) < 0.1 * abs(k2.real)


# =====================================================================================================================
# Two-layer solid planet (stiff iron core + rocky mantle)
# =====================================================================================================================
def test_two_layer_solid_love_succeeds():
    world = _two_layer_solid_world()
    world.solve_eos(G_to_use=G, verbose=False)
    world.solve_love_numbers(frequency_rad_s=_FREQ, verbose=False)
    assert world.love_solved is True


def test_two_layer_love_k2_positive_real():
    """k2 real part must be positive for a two-layer solid planet."""
    world = _two_layer_solid_world()
    world.solve_eos(G_to_use=G, verbose=False)
    world.solve_love_numbers(frequency_rad_s=_FREQ, verbose=False)
    k2 = world.love_number_k
    assert not cmath.isnan(k2)
    assert k2.real > 0.0


def test_two_layer_love_k2_less_than_fluid_limit():
    """k2 for any solid body must be below the fluid limit 1.5."""
    world = _two_layer_solid_world()
    world.solve_eos(G_to_use=G, verbose=False)
    world.solve_love_numbers(frequency_rad_s=_FREQ, verbose=False)
    k2 = world.love_number_k
    assert k2.real < 1.5
