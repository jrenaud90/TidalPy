"""Rheology tide model in LayeredWorld.calc_tides: per-frequency radial Love solves on a dissipative Maxwell sphere."""
import cmath
import math

import pytest

from TidalPy.constants import G
from TidalPy.Structures.worlds.layered import LayeredWorld
from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.Viscosity import make_viscosity
from TidalPy.Rheology.rheology import Maxwell
from TidalPy.Tides.classes.tide import make_tide


_PLANET_RADIUS = 6.0e6     # [m]
_DENSITY       = 4000.0    # [kg/m^3]
_STATIC_SHEAR  = 6.0e10    # [Pa]
_STATIC_BULK   = 1.3e11    # [Pa]
_SHEAR_VISC    = 1.0e16    # [Pa s]; tau = eta / mu ~ 1.7e5 s puts omega * tau at a few (dissipative)

_HOST_MASS = 1.9e27        # [kg]
_SMA       = 4.0e8         # [m]
_N         = 2.0e-5        # mean motion [rad s-1]
_SPIN      = 0.5 * _N      # non-synchronous so the (2, 2, 0, 0) mode is active
_ECC       = 0.01
_TIDAL_SCALE = 0.8


def _uniform_layer(tidal_scale=None):
    """A whole-planet layer and its mass."""
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    return BaseLayer(
        "mantle",
        0,
        0.0,
        _PLANET_RADIUS,
        mass,
        tidal_scale=tidal_scale,
    ), mass


def _rheology_world():
    """Single solid Maxwell sphere with a rheology tide model."""
    layer, mass = _uniform_layer(tidal_scale=_TIDAL_SCALE)
    world = LayeredWorld("rheo_planet", _PLANET_RADIUS, mass)
    layer.set_eos(ConstantDensityEOS(
        reference_density=_DENSITY, shear_modulus_static=_STATIC_SHEAR, bulk_modulus_static=_STATIC_BULK))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _SHEAR_VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e30}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Maxwell())
    world.add_layer(layer)
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
    return world


def _solve(world):
    world.calc_tides(
        orbital_frequency=_N,
        spin_frequency=_SPIN,
        eccentricity=_ECC,
        obliquity=0.0,
        semi_major_axis=_SMA,
        host_mass=_HOST_MASS,
    )


def _solved_rheology_world():
    world = _rheology_world()
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    _solve(world)
    return world


def test_rheology_calc_tides_requires_eos_solved():
    """The rheology tide model raises before the EOS is solved."""
    world = _rheology_world()
    assert world.tide_model_set
    with pytest.raises(RuntimeError):
        _solve(world)


def test_rheology_calc_tides_positive_heating():
    """calc_tides gives positive heating and active modes."""
    world = _solved_rheology_world()
    assert world.tides_solved
    assert world.get_num_tidal_modes() > 0
    assert world.get_tidal_heating() > 0.0


def test_rheology_per_mode_love_is_dissipative():
    """The (2, 2, 0, 0) Love number is dissipative and sub-fluid."""
    world = _solved_rheology_world()
    k = world.get_tidal_love_k(2, 2, 0, 0)
    assert not cmath.isnan(k)
    assert 0.0 < k.real < 1.5
    assert k.imag < 0.0


def test_rheology_potential_derivatives_present():
    """calc_tides fills a nonzero dU/dM."""
    world = _solved_rheology_world()
    dUdM, _, _ = world.get_tidal_potential_derivatives()
    assert abs(dUdM) > 0.0


def test_rheology_layer_heating_comes_from_the_radial_solution():
    """With the radial solver a lone layer takes all the heating whatever its tidal_scale."""
    world = _solved_rheology_world()
    total = world.get_tidal_heating()
    assert math.isclose(world.get_layer_tidal_heating(0), total, rel_tol=1.0e-12)


def test_layer_heating_nan_before_solve():
    """Layer heating is NaN before calc_tides."""
    world = _rheology_world()
    assert math.isnan(world.get_layer_tidal_heating(0))


def test_analytic_model_love_k_is_nan():
    """An analytic (cpl) model leaves the per-mode Love numbers NaN."""
    layer, mass = _uniform_layer()
    world = LayeredWorld("cpl_planet", _PLANET_RADIUS, mass)
    layer.set_eos(ConstantDensityEOS(
        reference_density=_DENSITY, shear_modulus_static=_STATIC_SHEAR, bulk_modulus_static=_STATIC_BULK))
    world.add_layer(layer)
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]}))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
    _solve(world)
    assert world.get_tidal_heating() > 0.0
    assert cmath.isnan(world.get_tidal_love_k(2, 2, 0, 0))
