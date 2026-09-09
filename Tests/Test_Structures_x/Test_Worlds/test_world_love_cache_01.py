"""
Tests for the cached, reusable radial (Love-number) solver on LayeredWorld.

solve_love_numbers separates the frequency-independent setup (built once and
cached) from the frequency-dependent work (recomputed per call).  These tests
verify that caching does not change results: repeated and interleaved calls must
reproduce single-call answers exactly, and a re-solve of the EOS must transparently
rebuild the cache.

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
_SHEAR_VISC    = 1.0e21    # [Pa s]


def _solid_world():
    """Single solid uniform sphere fully equipped for love-number solving."""
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.rheology_x.rheology import Maxwell

    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("solid_planet", _PLANET_RADIUS, mass)
    layer = PhysicsLayer("mantle", 0, 0.0, _PLANET_RADIUS, mass,
                         shear_modulus_static_pa=_STATIC_SHEAR,
                         bulk_modulus_static_pa=_STATIC_BULK)
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _SHEAR_VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e30}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Maxwell())
    world.add_layer(layer)
    return world


def _solved_world():
    world = _solid_world()
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    return world


def _klh(world):
    return (world.love_number_k, world.love_number_h, world.love_number_l)


def _assert_close(a, b, rel=1.0e-12):
    assert cmath.isclose(a, b, rel_tol=rel, abs_tol=1.0e-300), f"{a} != {b}"


# =====================================================================================================================
# Repeatability: the cache must not change results
# =====================================================================================================================
def test_repeated_same_frequency_is_identical():
    """Two back-to-back solves at the same frequency must agree bit-for-bit."""
    world = _solved_world()
    world.solve_love_numbers(frequency_rad_s=1.0e-5, verbose=False)
    first = _klh(world)
    world.solve_love_numbers(frequency_rad_s=1.0e-5, verbose=False)
    second = _klh(world)
    for a, b in zip(first, second):
        _assert_close(a, b)


def test_interleaved_frequency_sweep_returns_same_value():
    """Solving A, then B, then A again must reproduce the first A result.

    This is the key cache-corruption guard: the reused storage and re-applied
    non-dim arrays must leave no state behind that perturbs a repeated solve.
    """
    world = _solved_world()
    world.solve_love_numbers(frequency_rad_s=1.0e-5, verbose=False)
    a_first = _klh(world)
    world.solve_love_numbers(frequency_rad_s=3.0e-6, verbose=False)
    _ = _klh(world)
    world.solve_love_numbers(frequency_rad_s=1.0e-5, verbose=False)
    a_again = _klh(world)
    for a, b in zip(a_first, a_again):
        _assert_close(a, b)


def test_cached_matches_fresh_world():
    """A swept-then-returned solve must match a brand-new world solved once."""
    swept = _solved_world()
    swept.solve_love_numbers(frequency_rad_s=2.0e-6, verbose=False)
    swept.solve_love_numbers(frequency_rad_s=1.0e-5, verbose=False)
    swept_result = _klh(swept)

    fresh = _solved_world()
    fresh.solve_love_numbers(frequency_rad_s=1.0e-5, verbose=False)
    fresh_result = _klh(fresh)

    for a, b in zip(swept_result, fresh_result):
        _assert_close(a, b)


# =====================================================================================================================
# Frequency dependence still works through the cached path
# =====================================================================================================================
def test_frequency_changes_dissipation():
    """Im(k2) must depend on frequency (the per-call rheology refill is live).

    With Maxwell tau = eta/mu ~ 1.7e10 s, dissipation peaks near omega ~ 1/tau.
    A near-resonant frequency must give a larger |Im(k2)| than a near-elastic one.
    """
    world = _solved_world()
    world.solve_love_numbers(frequency_rad_s=1.0e-5, verbose=False)   # omega*tau >> 1, near-elastic
    im_elastic = abs(world.love_number_k.imag)
    world.solve_love_numbers(frequency_rad_s=6.0e-11, verbose=False)  # omega*tau ~ 1, dissipative
    im_resonant = abs(world.love_number_k.imag)
    assert im_resonant > im_elastic


# =====================================================================================================================
# Re-solving the EOS rebuilds the cache transparently
# =====================================================================================================================
def test_eos_resolve_invalidates_cache():
    """After re-solving the EOS, a love solve must still succeed and be valid."""
    world = _solved_world()
    world.solve_love_numbers(frequency_rad_s=1.0e-5, verbose=False)
    k_before = world.love_number_k

    # Re-solve the EOS (same inputs) — must invalidate and rebuild the cache.
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    world.solve_love_numbers(frequency_rad_s=1.0e-5, verbose=False)
    assert world.love_solved is True
    k_after = world.love_number_k
    # Identical EOS inputs => identical Love numbers after the rebuild.
    _assert_close(k_before, k_after, rel=1.0e-9)
