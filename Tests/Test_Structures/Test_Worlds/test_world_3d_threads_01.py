"""3D grid methods reproduce the one-thread result exactly (NaN cells included) for any ``num_threads``."""

import functools
import math

import numpy as np
import pytest

import TidalPy
from TidalPy.constants import G, update_constants
from TidalPy.Material import Material
from TidalPy.Rheology.rheology import Maxwell
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Tides.classes.tide import make_tide

from shared_materials import constant_solid


_RADIUS          = 1.8216e6
_MASS            = 8.9319e22
_HOST_MASS       = 1.898e27
_SEMI_MAJOR_AXIS = 4.217e8
_THREAD_COUNTS   = [0, 3, 64]    # 0: the automatic count

# The center has no depth-resolved solution, so cells at radius 0 are NaN (or skipped when an axis is summed).
_RADII       = np.array([0.0, 0.2, 0.45, 0.7, 0.95]) * _RADIUS
_COLATITUDES = np.linspace(0.05, np.pi - 0.05, 7)
_LONGITUDES  = np.linspace(0.0, 2.0 * np.pi, 9, endpoint=False)
_TIMES       = np.linspace(0.0, 1.5e5, 3)

_COLLAPSE_CASES = {
    "secular_grid": dict(radii=_RADII, colatitudes=_COLATITUDES, longitudes=_LONGITUDES),
    "instantaneous_grid": dict(radii=_RADII, colatitudes=_COLATITUDES, longitudes=_LONGITUDES, times=_TIMES,
                               orbit_averaged=False),
    "secular_radial_profile_numeric": dict(radii=_RADII, latitude_summed=True, longitude_summed=True,
                                           latitude_analytic=False),
    # The analytic colatitude collapse, which spreads its radii rather than colatitude rows over the threads.
    "secular_radial_profile_analytic": dict(radii=_RADII, latitude_summed=True, longitude_summed=True),
    "secular_total_analytic": dict(latitude_summed=True, longitude_summed=True, radial_summed=True),
    "secular_colatitude_profile": dict(colatitudes=_COLATITUDES, longitude_summed=True, radial_summed=True),
    "secular_total_numeric": dict(latitude_summed=True, longitude_summed=True, radial_summed=True,
                                  latitude_analytic=False),
    "secular_total_band": dict(latitude_summed=True, longitude_summed=True, radial_summed=True,
                               colatitude_min=0.3, colatitude_max=1.2),
    "instantaneous_total": dict(times=_TIMES, orbit_averaged=False, latitude_summed=True, longitude_summed=True,
                                radial_summed=True, latitude_nodes=6, longitude_nodes=8),
    "instantaneous_longitude_map": dict(radii=_RADII, longitudes=_LONGITUDES, times=_TIMES, orbit_averaged=False,
                                        latitude_summed=True, latitude_nodes=6),
}


def _material(density, shear, viscosity):
    return constant_solid(
        density, bulk_modulus=2.0e11, shear_modulus=shear, shear_viscosity=viscosity, bulk_viscosity=1.0e30)


@functools.lru_cache(maxsize=1)
def _world():
    """Io-sized two-layer Maxwell world at degrees 2 to 3, so many waves share frequencies and (l, m) pairs."""
    density = _MASS / ((4.0 / 3.0) * math.pi * _RADIUS ** 3)
    world = BaseWorld("threads", _RADIUS, _MASS)
    layers = (("core", 0.0, 0.45 * _RADIUS, 1.0e11, 1.0e22), ("mantle", 0.45 * _RADIUS, _RADIUS, 6.0e10, 1.0e15))
    for index, (name, r_inner, r_outer, shear, viscosity) in enumerate(layers):
        mass = (4.0 / 3.0) * math.pi * (r_outer ** 3 - r_inner ** 3) * density
        world.add_layer(Layer(name, index, r_inner, r_outer, mass, _material(density, shear, viscosity),
                              shear_rheology=Maxwell(), bulk_rheology=Maxwell()))
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=3, eccentricity_truncation=10, obliquity_truncation=2)
    world.solve_eos(G_to_use=G)
    return world


def _orbit():
    mean_motion = math.sqrt(G * (_HOST_MASS + _MASS) / _SEMI_MAJOR_AXIS ** 3)
    return (mean_motion, 1.2 * mean_motion, 0.1, 0.2, _SEMI_MAJOR_AXIS, _HOST_MASS)


@functools.lru_cache(maxsize=None)
def _one_thread_collapse(case):
    return _world().calc_3d_tides(*_orbit(), num_threads=1, **_COLLAPSE_CASES[case])


def _assert_same_outputs(threaded, serial):
    assert threaded.keys() == serial.keys()
    for key in serial:
        if key == "components":
            assert threaded[key] == serial[key]
            continue
        np.testing.assert_array_equal(np.asarray(threaded[key]), np.asarray(serial[key]), err_msg=key)


@pytest.mark.parametrize("num_threads", _THREAD_COUNTS)
@pytest.mark.parametrize("case", list(_COLLAPSE_CASES))
def test_calc_3d_tides_is_identical_for_any_thread_count(case, num_threads):
    serial = _one_thread_collapse(case)
    threaded = _world().calc_3d_tides(*_orbit(), num_threads=num_threads, **_COLLAPSE_CASES[case])
    _assert_same_outputs(threaded, serial)


def test_secular_grid_keeps_nan_at_the_center():
    """The compared grid really holds NaN cells: NaN at radius 0 and finite elsewhere."""
    heating = _one_thread_collapse("secular_grid")["heating"]
    assert np.isnan(heating[0]).all()
    assert np.isfinite(heating[1:]).all()


@pytest.mark.parametrize("num_threads", _THREAD_COUNTS)
@pytest.mark.parametrize("method", ["calc_3d_stress_strain", "calc_3d_displacements"])
def test_time_grid_is_identical_for_any_thread_count(method, num_threads):
    arguments = dict(radii=_RADII, colatitudes=_COLATITUDES, longitudes=_LONGITUDES, times=_TIMES)
    serial = getattr(_world(), method)(*_orbit(), num_threads=1, **arguments)
    threaded = getattr(_world(), method)(*_orbit(), num_threads=num_threads, **arguments)
    _assert_same_outputs(threaded, serial)


@pytest.mark.parametrize("num_threads", _THREAD_COUNTS)
def test_batch_heating_is_identical_for_any_thread_count(num_threads):
    """Paired points with repeated colatitudes, a NaN colatitude, and a radius without a solution."""
    radius_grid, colatitude_grid = np.meshgrid(_RADII, _COLATITUDES, indexing="ij")
    radii = np.append(radius_grid.ravel(), 0.5 * _RADIUS)
    colatitudes = np.append(colatitude_grid.ravel(), np.nan)
    serial = _world().get_3d_tidal_heating_array(*_orbit(), radii, colatitudes, num_threads=1)
    threaded = _world().get_3d_tidal_heating_array(*_orbit(), radii, colatitudes, num_threads=num_threads)
    np.testing.assert_array_equal(threaded, serial)
    assert np.isnan(serial[-1])


_METHOD_CALLS = {
    "calc_3d_tides": lambda world, count: world.calc_3d_tides(
        *_orbit(), radii=_RADII, colatitudes=_COLATITUDES, longitudes=_LONGITUDES, num_threads=count),
    "calc_3d_stress_strain": lambda world, count: world.calc_3d_stress_strain(
        *_orbit(), radii=_RADII, colatitudes=_COLATITUDES, longitudes=_LONGITUDES, times=_TIMES, num_threads=count),
    "calc_3d_displacements": lambda world, count: world.calc_3d_displacements(
        *_orbit(), radii=_RADII, colatitudes=_COLATITUDES, longitudes=_LONGITUDES, times=_TIMES, num_threads=count),
    "get_3d_tidal_heating_array": lambda world, count: world.get_3d_tidal_heating_array(
        *_orbit(), _RADII, np.full(_RADII.shape, 1.0), num_threads=count),
}


@pytest.mark.parametrize("num_threads", [-1, -2])
@pytest.mark.parametrize("method", list(_METHOD_CALLS))
def test_negative_thread_count_raises(method, num_threads):
    with pytest.raises(ValueError, match="num_threads"):
        _METHOD_CALLS[method](_world(), num_threads)


@pytest.mark.parametrize("method", list(_METHOD_CALLS))
def test_failed_radial_solve_on_a_worker_thread_raises(restore_config, method):
    """A radial solve that fails on a worker thread raises on the calling thread."""
    # Too few integration steps to cross the first layer, so every radial solve fails.
    TidalPy.config["radial_solver"]["max_num_steps"] = 5
    update_constants()
    with pytest.raises(RuntimeError, match="radial solve failed"):
        _METHOD_CALLS[method](_world(), 16)
