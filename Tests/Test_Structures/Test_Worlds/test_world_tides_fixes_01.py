"""Layered-world tides: Love solves past 16-bit frequency indices, orbit validation on the 3D paths, NaN where a
radius has no radial solution, and a failed calc_tides leaving no stale results."""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.Rheology.rheology import Elastic, Maxwell
from TidalPy.Structures import build_world
from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Tides.potential import global_potential
from TidalPy.Utilities.logging.logger import flush_logger, init_logger
from TidalPy.Viscosity import make_viscosity
from TidalPy.initialize import build_logging_config


_R = 1.8216e6
_MASS = 8.9319e22
_HOST = 1.898e27
_SMA = 4.217e8
_N = math.sqrt(G * (_HOST + _MASS) / _SMA ** 3)
_DENSITY = _MASS / ((4.0 / 3.0) * math.pi * _R ** 3)
_STATE = dict(orbital_frequency=_N, spin_frequency=_N, eccentricity=0.05, obliquity=0.0, semi_major_axis=_SMA,
              host_mass=_HOST)
# The automatic starting radius at l = 2 is R sqrt(1e-5), about 0.003 R: the first radius has no radial solution.
_RADII = np.array([1.0e-4, 0.5, 0.9]) * _R
_SUMMED = dict(latitude_summed=True, longitude_summed=True, radial_summed=True, num_threads=1)
_RANGE_WARNING = "can underestimate the tides by 10% or more"
_MISSING_NODES_WARNING = "have no radial solution"


@pytest.fixture
def spdlog_text(tmp_path):
    """Route the C++ logger to a temporary file for the test and hand back a reader for its text."""
    log_path = tmp_path / "tidalpy.log"
    init_logger({"console_level": "off", "file_level": "warning", "log_to_file": True,
                 "log_file_path": str(log_path)})

    def read():
        flush_logger()
        return log_path.read_text(encoding="utf-8") if log_path.exists() else ""

    yield read
    init_logger(build_logging_config())


def _world(eccentricity_truncation=2):
    """One compressible, static Maxwell layer solved with the shooting method."""
    world = BaseWorld("homogeneous", _R, _MASS)
    layer = BaseLayer("mantle", 0, 0.0, _R, _MASS)
    layer.set_eos(ConstantDensityEOS(
        reference_density=_DENSITY, shear_modulus_static=6.0e10, bulk_modulus_static=2.0e11))
    layer.is_static = True
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e15}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e30}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world.add_layer(layer)
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(
        min_degree_l=2, max_degree_l=2, eccentricity_truncation=eccentricity_truncation, obliquity_truncation=0)
    world.solve_eos(G_to_use=G)
    return world


# ======================================================================================================================
# Love solves past 16-bit frequency indices
# ======================================================================================================================
def test_love_numbers_past_65536_unique_frequencies():
    """Every mode gets the Love numbers of its own frequency, however many unique frequencies there are.

    With 'exact' eccentricity at e = 0.97, degrees 2 to 10 hold about 98,600 unique frequencies. A 16-bit solve key
    once folded the indices past 65,535 onto earlier ones, so those modes took another frequency's Love numbers.
    """
    orbital_frequency = 4.1e-5
    spin_frequency = 1.37 * orbital_frequency
    eccentricity = 0.97
    world = build_world("io")
    world.solve_eos()
    world.set_tide_config(
        min_degree_l=2,
        max_degree_l=10,
        eccentricity_truncation="exact",
        obliquity_truncation=0,
        love_method="homogeneous",
        layer_tidal_heating=False)
    world.calc_tides(orbital_frequency, spin_frequency, eccentricity, 0.0, _SMA, _HOST)

    tolerance = world.get_tide_config()["eccentricity_exact_tolerance"]
    potential = global_potential(
        world.radius, orbital_frequency, spin_frequency, eccentricity, 0.0, _SMA, _HOST, G, 2, 10, "exact", 0,
        tolerance)
    frequency_index = potential[1]
    assert len(potential[2]) > 65536
    # The modes whose frequency index is past the 16-bit range, spread over the degrees that reach it.
    late_modes = [key for key in potential[3] if frequency_index[key] >= 65536]
    assert len(late_modes) > 100
    rng = np.random.default_rng(3)
    checked = [late_modes[i] for i in rng.choice(len(late_modes), 12, replace=False)]
    love_k = [world.get_tidal_love_k(*mode) for mode in checked]
    for (degree_l, order_m, p, q), k_mode in zip(checked, love_k):
        frequency = abs((degree_l - 2 * p + q) * orbital_frequency - order_m * spin_frequency)
        k_reference = world.solve_love_numbers(
            frequency=frequency, degree_l=degree_l, love_method="homogeneous")["love_number_k"]
        k_reference = complex(np.ravel(k_reference)[0])
        assert abs(k_mode - k_reference) <= 1e-9 * abs(k_reference), (degree_l, order_m, p, q)


# ======================================================================================================================
# Orbit validation on the 3D paths
# ======================================================================================================================
_CALLS_3D = {
    "get_3d_tidal_heating": lambda world, state: world.get_3d_tidal_heating(
        **state, radius=0.5 * _R, colatitude=1.0),
    "get_3d_tidal_heating_array": lambda world, state: world.get_3d_tidal_heating_array(
        **state, radii=np.array([0.5 * _R]), colatitudes=np.array([1.0]), num_threads=1),
    "calc_3d_tides": lambda world, state: world.calc_3d_tides(**state, **_SUMMED),
    "calc_3d_stress_strain": lambda world, state: world.calc_3d_stress_strain(
        **state, radii=[0.5 * _R], colatitudes=[1.0], longitudes=[0.3], times=[0.0], num_threads=1),
    "calc_3d_displacements": lambda world, state: world.calc_3d_displacements(
        **state, radii=[0.5 * _R], colatitudes=[1.0], longitudes=[0.3], times=[0.0], num_threads=1),
}


@pytest.mark.parametrize("method", list(_CALLS_3D))
@pytest.mark.parametrize("key, value, message", [
    ("eccentricity", 1.5, "eccentricity"),
    ("eccentricity", -0.1, "eccentricity"),
    ("semi_major_axis", -_SMA, "semi-major axis"),
])
def test_3d_paths_reject_what_calc_tides_rejects(method, key, value, message):
    world = _world()
    state = dict(_STATE, **{key: value})
    with pytest.raises(ValueError, match=message):
        world.calc_tides(**state)
    with pytest.raises(ValueError, match=message):
        _CALLS_3D[method](world, state)


@pytest.mark.parametrize("method", list(_CALLS_3D))
def test_3d_paths_reject_a_non_finite_orbital_frequency(method):
    world = _world()
    state = dict(_STATE, orbital_frequency=math.nan)
    with pytest.raises(ValueError, match="orbital frequency"):
        world.calc_tides(**state)
    with pytest.raises(ValueError, match="orbital frequency"):
        _CALLS_3D[method](world, state)


@pytest.mark.parametrize("method", list(_CALLS_3D))
def test_3d_paths_warn_past_the_eccentricity_truncation_range(spdlog_text, method):
    """Level 6 holds to 10% up to e of about 0.3 (test_world_eccentricity_range_01.py)."""
    world = _world(eccentricity_truncation=6)
    _CALLS_3D[method](world, dict(_STATE, eccentricity=0.2))
    assert _RANGE_WARNING not in spdlog_text()
    _CALLS_3D[method](world, dict(_STATE, eccentricity=0.35))
    _CALLS_3D[method](world, dict(_STATE, eccentricity=0.35))
    assert spdlog_text().count(_RANGE_WARNING) == 1


# ======================================================================================================================
# NaN where a radius has no radial solution
# ======================================================================================================================
@pytest.mark.parametrize("arguments", [
    dict(latitude_summed=True, longitude_summed=True),
    dict(latitude_summed=True, longitude_summed=True, latitude_analytic=False),
    dict(colatitudes=np.array([0.4, 1.2]), longitude_summed=True),
    dict(longitudes=np.array([0.0, 1.0]), latitude_summed=True),
    dict(latitude_summed=True, longitude_summed=True, orbit_averaged=False, times=np.array([0.0, 1.0e4])),
    dict(colatitudes=np.array([1.2]), longitudes=np.array([0.3])),
], ids=["analytic_profile", "quadrature_profile", "longitude_mean_map", "latitude_summed_map",
        "instantaneous_profile", "raw_grid"])
def test_radius_without_solution_is_nan(arguments):
    """A radius below the solver start is NaN on every output with a radius axis, not 0 (a liquid's value)."""
    heating = _world().calc_3d_tides(**_STATE, radii=_RADII, num_threads=1, **arguments)["heating"]
    assert np.all(np.isnan(heating[0]))
    assert np.all(np.isfinite(heating[1:])) and np.all(heating[1:] != 0.0)


def test_volume_integral_warns_when_nodes_are_left_out(spdlog_text):
    """Nodes below the starting radius are left out of the volume integral with a warning when they matter."""
    world = _world()
    world.calc_tides(**_STATE)
    heating_1d = world.get_tidal_heating()
    total = world.calc_3d_tides(**_STATE, **_SUMMED)["total"]
    assert math.isclose(total, heating_1d, rel_tol=5e-8)
    assert _MISSING_NODES_WARNING not in spdlog_text()

    # A starting radius near 0.32 R leaves about 3% of the volume without a radial solution.
    world.set_solver_defaults(radial_solver={"start_radius_tolerance": 0.1})
    total = world.calc_3d_tides(**_STATE, **_SUMMED)["total"]
    assert np.isfinite(total) and total < heating_1d
    text = spdlog_text()
    assert _MISSING_NODES_WARNING in text and "homogeneous" in text


# ======================================================================================================================
# A failed calc_tides leaves no stale results
# ======================================================================================================================
@pytest.mark.parametrize("bad_state", [dict(eccentricity=1.5), dict(semi_major_axis=-_SMA)])
def test_failed_calc_tides_clears_the_previous_results(bad_state):
    world = _world()
    world.calc_tides(**_STATE)
    assert np.isfinite(world.get_tidal_heating())
    assert all(np.isfinite(layer.get_tidal_heating()) for layer in world)
    assert np.isfinite(world.get_layer_tidal_heating(0))

    with pytest.raises(ValueError):
        world.calc_tides(**dict(_STATE, **bad_state))
    assert math.isnan(world.get_tidal_heating())
    assert all(math.isnan(layer.get_tidal_heating()) for layer in world)
    assert math.isnan(world.get_layer_tidal_heating(0))

    # A good call after the failure solves again.
    world.calc_tides(**_STATE)
    assert np.isfinite(world.get_tidal_heating())
    assert all(np.isfinite(layer.get_tidal_heating()) for layer in world)
