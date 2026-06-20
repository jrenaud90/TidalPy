"""
Tests for the propagation-matrix radial-solver method on LayeredWorld.

``solve_love_numbers(use_prop_matrix=True)`` selects the propagation-matrix
technique instead of the shooting method. It is only valid for a single solid,
static, incompressible layer, for which the degree-2 tidal Love number of a
homogeneous sphere has the closed form::

    k2 = 1.5 / (1 + 19 mu / (2 rho g R))

with g = (4/3) pi G rho R the surface gravity. These tests check that the
world-integrated propagation-matrix path reproduces that analytic value and that
incompatible worlds are rejected gracefully.

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
_SHEAR_VISC    = 1.0e23    # [Pa s] — very high => elastic limit (mu ~ mu_static)
_FREQ          = 1.0e-5    # [rad/s]


def _imports():
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.rheology_x.rheology import Maxwell
    return LayeredWorld, PhysicsLayer, ConstantDensityEOS, make_viscosity, Maxwell


def _incompressible_solid_world():
    """Single solid, static, incompressible uniform sphere."""
    LayeredWorld, PhysicsLayer, ConstantDensityEOS, make_viscosity, Maxwell = _imports()
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("incompressible_planet", _PLANET_RADIUS, mass)
    layer = PhysicsLayer("mantle", 0, 0.0, _PLANET_RADIUS, mass,
                         shear_modulus_static_pa=_STATIC_SHEAR,
                         bulk_modulus_static_pa=_STATIC_BULK)
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _SHEAR_VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e30}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Maxwell())
    # Propagation matrix requires solid (default) + static (default) + incompressible.
    layer.is_incompressible = True
    assert layer.is_incompressible is True
    world.add_layer(layer)
    return world


def _analytic_k2():
    g_surface = (4.0 / 3.0) * math.pi * G * _DENSITY * _PLANET_RADIUS
    mu_tilde  = 19.0 * _STATIC_SHEAR / (2.0 * _DENSITY * g_surface * _PLANET_RADIUS)
    return 1.5 / (1.0 + mu_tilde)


# =====================================================================================================================
# Propagation matrix reproduces the analytic homogeneous-sphere Love number
# =====================================================================================================================
def test_prop_matrix_succeeds():
    world = _incompressible_solid_world()
    world.solve_eos(G_to_use=G, verbose=False)
    world.solve_love_numbers(frequency_rad_s=_FREQ, use_prop_matrix=True, verbose=False)
    assert world.love_solved is True


def test_prop_matrix_k2_matches_analytic():
    world = _incompressible_solid_world()
    world.solve_eos(G_to_use=G, verbose=False)
    world.solve_love_numbers(frequency_rad_s=_FREQ, use_prop_matrix=True, verbose=False)
    k2 = world.love_number_k
    analytic = _analytic_k2()
    assert not cmath.isnan(k2)
    assert k2.real == pytest.approx(analytic, rel=0.02)
    # Elastic limit => negligible imaginary part.
    assert abs(k2.imag) < 0.05 * abs(k2.real)


def test_prop_matrix_close_to_shooting_for_same_world():
    """For a near-incompressible sphere, the two methods should roughly agree."""
    world = _incompressible_solid_world()
    world.solve_eos(G_to_use=G, verbose=False)
    world.solve_love_numbers(frequency_rad_s=_FREQ, use_prop_matrix=True, verbose=False)
    k2_matrix = world.love_number_k
    world.solve_love_numbers(frequency_rad_s=_FREQ, use_prop_matrix=False, verbose=False)
    k2_shoot = world.love_number_k
    assert k2_matrix.real == pytest.approx(k2_shoot.real, rel=0.10)


# =====================================================================================================================
# Incompatible worlds are rejected gracefully (no crash)
# =====================================================================================================================
def test_prop_matrix_rejects_two_layers():
    LayeredWorld, PhysicsLayer, ConstantDensityEOS, make_viscosity, Maxwell = _imports()
    r_core = 3.0e6
    mass = (4.0 / 3.0) * math.pi * (
        8000.0 * r_core ** 3 + 3300.0 * (_PLANET_RADIUS ** 3 - r_core ** 3))
    world = LayeredWorld("two_layer", _PLANET_RADIUS, mass)
    for name, idx, ri, ro, rho in (("core", 0, 0.0, r_core, 8000.0),
                                   ("mantle", 1, r_core, _PLANET_RADIUS, 3300.0)):
        layer = PhysicsLayer(name, idx, ri, ro, 0.0,
                             shear_modulus_static_pa=_STATIC_SHEAR,
                             bulk_modulus_static_pa=_STATIC_BULK)
        layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=rho))
        layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _SHEAR_VISC}))
        layer.set_shear_rheology(Maxwell())
        layer.is_incompressible = True
        world.add_layer(layer)
    world.solve_eos(G_to_use=G, verbose=False)
    world.solve_love_numbers(frequency_rad_s=_FREQ, use_prop_matrix=True, verbose=False)
    # Must fail cleanly, not crash.
    assert world.love_success is False
    assert world.love_error_code != 0


def test_prop_matrix_rejects_compressible_layer():
    LayeredWorld, PhysicsLayer, ConstantDensityEOS, make_viscosity, Maxwell = _imports()
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("compressible_planet", _PLANET_RADIUS, mass)
    layer = PhysicsLayer("mantle", 0, 0.0, _PLANET_RADIUS, mass,
                         shear_modulus_static_pa=_STATIC_SHEAR,
                         bulk_modulus_static_pa=_STATIC_BULK)
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _SHEAR_VISC}))
    layer.set_shear_rheology(Maxwell())
    # Leave is_incompressible at its default (False) => prop matrix must reject.
    assert layer.is_incompressible is False
    world.add_layer(layer)
    world.solve_eos(G_to_use=G, verbose=False)
    world.solve_love_numbers(frequency_rad_s=_FREQ, use_prop_matrix=True, verbose=False)
    assert world.love_success is False
    assert world.love_error_code != 0
