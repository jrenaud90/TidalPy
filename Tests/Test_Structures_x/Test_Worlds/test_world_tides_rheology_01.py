"""Tests for the rheology tide model coupled to the radial solver in LayeredWorld.calc_tides.

The rheology global-tidal-model path differs from the analytic models (cpl/ctl/ctl_q): instead
of a fixed -Im[k_l], it runs the world radial solver at each unique tidal frequency to obtain the
complex Love numbers, then collapses the global modes with the per-mode k. These tests build a
single-layer Maxwell world (the same setup the radial-solver tests use, but with a viscosity in
the dissipative regime) and check the full calc_tides path:

  * the rheology model requires the EOS to be solved first,
  * calc_tides produces positive tidal heating and active modes,
  * the per-mode Love number is dissipative (Im[k] < 0) and sub-fluid (0 < Re[k] < 1.5),
  * per-layer heating is the world heating scaled by the layer's tidal_scale and is reported on
    the layer object itself,
  * analytic models leave the per-mode Love numbers unset (NaN).

Requires the Cython extensions to be compiled first (``uv pip install -v <repo_root>``).
"""
import cmath
import math

import pytest

from TidalPy.constants import G


# Single uniform Maxwell sphere in a clearly dissipative regime (omega*tau ~ a few).
_PLANET_RADIUS = 6.0e6     # [m]
_DENSITY       = 4000.0    # [kg/m^3]
_STATIC_SHEAR  = 6.0e10    # [Pa]
_STATIC_BULK   = 1.3e11    # [Pa]
_SHEAR_VISC    = 1.0e16    # [Pa s] -> tau = eta/mu ~ 1.7e5 s, omega*tau ~ a few

# Non-synchronous orbital/spin state so the semidiurnal mode (l=2, m=2, p=0, q=0) is active.
_HOST_MASS = 1.9e27        # [kg]
_SMA       = 4.0e8         # [m]
_N         = 2.0e-5        # mean motion [rad s-1]
_SPIN      = 0.5 * _N      # slow rotator -> |mode_2200| = |2n - 2*spin| = n != 0
_ECC       = 0.01
_TIDAL_SCALE = 0.8


def _imports():
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.rheology_x.rheology import Maxwell
    from TidalPy.Tides_x.classes.tide import make_tide
    return LayeredWorld, PhysicsLayer, ConstantDensityEOS, make_viscosity, Maxwell, make_tide


def _rheology_world(tidal_scale: float = _TIDAL_SCALE):
    """Single solid Maxwell sphere with a rheology tide model attached."""
    (LayeredWorld, PhysicsLayer, ConstantDensityEOS,
     make_viscosity, Maxwell, make_tide) = _imports()
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("rheo_planet", _PLANET_RADIUS, mass)
    layer = PhysicsLayer("mantle", 0, 0.0, _PLANET_RADIUS, mass,
                         shear_modulus_static_pa=_STATIC_SHEAR,
                         bulk_modulus_static_pa=_STATIC_BULK,
                         tidal_scale=tidal_scale)
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _SHEAR_VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e30}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Maxwell())
    world.add_layer(layer)
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=2,
                          eccentricity_truncation=2, obliquity_truncation=0)
    return world


def _solve(world):
    world.calc_tides(orbital_frequency=_N, spin_frequency=_SPIN, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=_SMA, host_mass=_HOST_MASS)


# =====================================================================================================================
# Pre-conditions
# =====================================================================================================================
def test_rheology_calc_tides_requires_eos_solved():
    """The rheology tide model needs the radial solver, so the EOS must be solved first."""
    world = _rheology_world()
    assert world.tide_model_set
    with pytest.raises(RuntimeError):
        _solve(world)


# =====================================================================================================================
# Full rheology pipeline
# =====================================================================================================================
def test_rheology_calc_tides_positive_heating():
    world = _rheology_world()
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    _solve(world)
    assert world.tides_solved
    assert world.get_num_tidal_modes() > 0
    assert world.get_tidal_heating() > 0.0


def test_rheology_per_mode_love_is_dissipative():
    """The semidiurnal (2,2,0,0) Love number must dissipate (Im[k] < 0) and stay sub-fluid."""
    world = _rheology_world()
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    _solve(world)
    k = world.get_tidal_love_k(2, 2, 0, 0)
    assert not cmath.isnan(k)
    assert 0.0 < k.real < 1.5
    assert k.imag < 0.0


def test_rheology_potential_derivatives_present():
    world = _rheology_world()
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    _solve(world)
    dUdM, dUdw, dUdO = world.get_tidal_potential_derivatives()
    assert abs(dUdM) > 0.0


# =====================================================================================================================
# Layer-side heating
# =====================================================================================================================
def test_rheology_layer_heating_scaled_by_tidal_scale():
    """Per-layer heating is the world heating scaled by the layer's tidal_scale.

    calc_tides also writes this value onto the C++ layer object (layer.get_tidal_heating);
    from Python it is read back through the world-side accessor.
    """
    world = _rheology_world(tidal_scale=0.8)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    _solve(world)
    total = world.get_tidal_heating()
    assert math.isclose(world.get_layer_tidal_heating(0), total * 0.8, rel_tol=1.0e-9)


def test_layer_heating_nan_before_solve():
    world = _rheology_world()
    assert math.isnan(world.get_layer_tidal_heating(0))


# =====================================================================================================================
# Analytic models carry no per-mode Love numbers
# =====================================================================================================================
def test_analytic_model_love_k_is_nan():
    """A cpl world's get_tidal_love_k is NaN: analytic models have no displacement Love numbers."""
    (LayeredWorld, PhysicsLayer, ConstantDensityEOS,
     make_viscosity, Maxwell, make_tide) = _imports()
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("cpl_planet", _PLANET_RADIUS, mass)
    layer = PhysicsLayer("mantle", 0, 0.0, _PLANET_RADIUS, mass,
                         shear_modulus_static_pa=_STATIC_SHEAR,
                         bulk_modulus_static_pa=_STATIC_BULK)
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    world.add_layer(layer)
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]}))
    world.set_tide_config(min_degree_l=2, max_degree_l=2,
                          eccentricity_truncation=2, obliquity_truncation=0)
    _solve(world)
    assert world.get_tidal_heating() > 0.0
    assert cmath.isnan(world.get_tidal_love_k(2, 2, 0, 0))
