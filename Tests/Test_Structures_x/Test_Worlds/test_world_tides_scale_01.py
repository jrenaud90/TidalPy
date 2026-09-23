"""Tidal scales and the tidal heating of each layer.

A layer's ``tidal_scale`` is its share of the planet in the quasi-homogeneous Love methods (``homogeneous``, ``cpl``,
``ctl``) and of an analytic tide model's heating; unset, it is the layer's volume over the planet's, and a layer that
is not tidal has none. How ``calc_tides`` resolves each layer's heating depends on where the Love numbers come from:

* an analytic tide model: the total times the layer's tidal scale;
* a quasi-homogeneous Love method: the heating of the layer's own scaled Love numbers, and the layers sum to the
  total, which is the sum of each layer as a homogeneous planet scaled by its tidal scale;
* the radial solver: the volume integral of the radial solution over the layer, scaled to the total (NaN when the
  tides config's ``layer_tidal_heating`` is off).
"""
import math

import pytest

from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
from TidalPy.rheology_x import Maxwell
from TidalPy.structures_x import build_world
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.Tides_x.classes.tide import make_tide
from TidalPy.Tides_x.love import calc_homogeneous_love_numbers
from TidalPy.viscosity_x import make_viscosity


_R         = 1.6e6
_R_CORE    = 0.5 * _R   # core volume fraction (0.5)^3 = 0.125
_HOST_MASS = 1.898e27
_SMA       = 4.2e8
_N         = 2.05e-5
_ECC       = 0.0041
_DENSITY   = 4000.0


def _layer(name, index, radius_inner, radius_outer, shear, viscosity, **kwargs):
    layer = PhysicsLayer(name, index, radius_inner, radius_outer, 0.0, **kwargs)
    layer.set_eos(ConstantDensityEOS(
        reference_density=_DENSITY, shear_modulus_static=shear, bulk_modulus_static=2.0e11))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": viscosity}))
    layer.set_shear_rheology(Maxwell())
    return layer


def _two_layer(core_scale=None, mantle_scale=None, core_tidal=True):
    mass = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY
    world = LayeredWorld("scaled", _R, mass)
    world.add_layer(_layer("core", 0, 0.0, _R_CORE, 8.0e10, 1.0e22, tidal_scale=core_scale, is_tidal=core_tidal))
    world.add_layer(_layer("mantle", 1, _R_CORE, _R, 5.0e10, 1.0e14, tidal_scale=mantle_scale))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
    return world


def _solve(world):
    world.calc_tides(orbital_frequency=_N, spin_frequency=_N, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=_SMA, host_mass=_HOST_MASS)


# =====================================================================================================================
# The tidal scale of a layer
# =====================================================================================================================
def test_an_unset_scale_is_the_volume_fraction():
    world = _two_layer()
    assert world.core.tidal_scale is None
    assert world.get_layer_tidal_scale(0) == pytest.approx(0.125, rel=1e-12)
    assert world.get_layer_tidal_scale(1) == pytest.approx(0.875, rel=1e-12)


def test_a_set_scale_wins_and_a_non_tidal_layer_has_none():
    world = _two_layer(mantle_scale=0.6, core_tidal=False)
    assert world.get_layer_tidal_scale(1) == 0.6
    assert world.get_layer_tidal_scale(0) == 0.0
    world.mantle.tidal_scale = None
    assert world.get_layer_tidal_scale(1) == pytest.approx(0.875, rel=1e-12)


def test_a_set_scale_is_written_and_an_unset_one_is_not():
    world = _two_layer(mantle_scale=0.6)
    layers = world.get_config_dict()["layers"]
    assert layers["mantle"]["tidal_scale"] == 0.6
    assert "tidal_scale" not in layers["core"]


# =====================================================================================================================
# An analytic tide model
# =====================================================================================================================
def test_an_analytic_model_shares_the_total_by_tidal_scale():
    world = _two_layer()
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]}))
    _solve(world)
    total = world.get_tidal_heating()
    assert world.get_layer_tidal_heating(0) == pytest.approx(0.125 * total, rel=1e-12)
    assert world.get_layer_tidal_heating(1) == pytest.approx(0.875 * total, rel=1e-12)


# =====================================================================================================================
# The quasi-homogeneous Love methods
# =====================================================================================================================
def _homogeneous_world(**kwargs):
    world = _two_layer(**kwargs)
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(love_method="homogeneous")
    world.solve_eos()
    return world


def test_the_world_love_number_is_the_scaled_sum_of_the_layers():
    world = _homogeneous_world()
    k2 = world.solve_love_numbers(frequency=_N, degree_l=2, love_method="homogeneous")["love_number_k"]
    parts = world.love_layer_parts
    assert [part["layer"] for part in parts] == ["core", "mantle"]
    gravity = world.surface_gravity_eos
    density = world.planet_mass_eos / ((4.0 / 3.0) * math.pi * _R ** 3)
    expected = 0.0
    for part in parts:
        mu = Maxwell().calc_complex_modulus(
            8.0e10 if part["layer"] == "core" else 5.0e10, 1.0e22 if part["layer"] == "core" else 1.0e14, _N)
        k_layer = calc_homogeneous_love_numbers(mu, density, gravity, _R, 2).k
        assert part["love_number_k"] == pytest.approx(k_layer, rel=1e-10)
        expected += part["tidal_scale"] * k_layer
    assert k2 == pytest.approx(expected, rel=1e-12)


def test_a_one_layer_planet_is_the_homogeneous_sphere():
    mass = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY
    world = LayeredWorld("single", _R, mass)
    world.add_layer(_layer("mantle", 0, 0.0, _R, 5.0e10, 1.0e14))
    world.solve_eos()
    k2 = world.solve_love_numbers(frequency=_N, degree_l=2, love_method="homogeneous")["love_number_k"]
    mu = Maxwell().calc_complex_modulus(5.0e10, 1.0e14, _N)
    expected = calc_homogeneous_love_numbers(mu, _DENSITY, world.surface_gravity_eos, _R, 2).k
    assert k2 == pytest.approx(expected, rel=1e-10)


def test_quasi_homogeneous_layer_heating_sums_to_the_total():
    """Each layer takes the heating of its own scaled Love numbers; a small, weak layer adds its share only."""
    world = _homogeneous_world()
    _solve(world)
    total = world.get_tidal_heating()
    core, mantle = world.get_layer_tidal_heating(0), world.get_layer_tidal_heating(1)
    assert core + mantle == pytest.approx(total, rel=1e-12)
    assert 0.0 < core < mantle
    assert world.core.get_tidal_heating() == core


# =====================================================================================================================
# The radial solver
# =====================================================================================================================
@pytest.fixture(scope="module")
def io_tides():
    io = build_world("io")
    io.solve_eos()
    n = 4.11e-5
    state = dict(orbital_frequency=n, spin_frequency=n, eccentricity=0.0041, obliquity=0.0,
                 semi_major_axis=4.217e8, host_mass=1.898e27)
    return io, state


def test_radial_layer_heating_is_the_volume_integral(io_tides):
    """The layers sum to the total and split as the 3D integral does."""
    io, state = io_tides
    io.calc_tides(**state)
    total = io.get_tidal_heating()
    heating = [io.get_layer_tidal_heating(i) for i in range(io.num_layers)]
    assert sum(heating) == pytest.approx(total, rel=1e-12)
    integral = io.calc_3d_tides(**state, latitude_summed=True, longitude_summed=True, radial_summed=True)["per_layer"]
    for layer_heating, layer_integral in zip(heating, integral):
        assert layer_heating / total == pytest.approx(layer_integral / sum(integral), rel=1e-9)
    # Most of Io's heat is in its asthenosphere.
    assert heating[2] / total > 0.9


def test_radial_layer_heating_can_be_switched_off(io_tides):
    io, state = io_tides
    io.set_tide_config(layer_tidal_heating=False)
    try:
        io.calc_tides(**state)
        assert math.isfinite(io.get_tidal_heating())
        assert all(math.isnan(io.get_layer_tidal_heating(i)) for i in range(io.num_layers))
    finally:
        io.set_tide_config(layer_tidal_heating=True)


def test_calc_tides_leaves_the_world_love_solve_alone(io_tides):
    """The tide paths solve into workspaces of their own; the world's last solve_love_numbers result stays."""
    io, state = io_tides
    k2 = io.solve_love_numbers(frequency=1.0e-6, degree_l=2)["love_number_k"]
    io.set_tide_config(max_degree_l=3)
    try:
        io.calc_tides(**state)
        io.calc_3d_tides(**state, latitude_summed=True, longitude_summed=True, radial_summed=True)
    finally:
        io.set_tide_config(max_degree_l=2)
    assert io.love_number_k == k2
