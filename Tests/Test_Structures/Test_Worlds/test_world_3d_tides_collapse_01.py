"""The 3D tidal heating grid and its collapsed, secular, and instantaneous forms (``BaseWorld.calc_3d_tides``)."""
import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.Material import Material
from TidalPy.Rheology.rheology import Elastic, Maxwell
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Utilities.conversions import orbital_motion2semi_a

from shared_materials import constant_solid
from numpy_compat import trapezoid


_R = 1.0e6
_DENSITY = 5000.0
_SHEAR = 5.0e10
_BULK = 1.0e11
_VISC = 1.0e19
_N = 2.0 * np.pi / 86400.0
_ECC = 0.05
_HOST = mass_trap1
_MASS = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY
_SOFT_SHEAR = 1.0e9
_SOFT_VISC = 1.0e13
_SUMMED = dict(latitude_summed=True, longitude_summed=True, radial_summed=True)


def _material(shear=_SHEAR, shear_visc=_VISC):
    return constant_solid(
        _DENSITY, bulk_modulus=_BULK, shear_modulus=shear, shear_viscosity=shear_visc, bulk_viscosity=_VISC)


def _build_world(two_layer=False, soft_shell=False):
    world = BaseWorld("w", _R, _MASS)

    def _mk(name, idx, r_in, r_out, shear=_SHEAR, shear_visc=_VISC):
        mass = (4.0 / 3.0) * math.pi * (r_out ** 3 - r_in ** 3) * _DENSITY
        return Layer(name, idx, r_in, r_out, mass, _material(shear, shear_visc), is_static=False,
                     shear_rheology=Maxwell(), bulk_rheology=Elastic())

    if soft_shell:
        # A thin dissipating shell whose modulus jumps 50x at its base.
        world.add_layer(_mk("core", 0, 0.0, 0.5 * _R))
        world.add_layer(_mk("mantle", 1, 0.5 * _R, 0.95 * _R))
        world.add_layer(_mk("shell", 2, 0.95 * _R, _R, shear=_SOFT_SHEAR, shear_visc=_SOFT_VISC))
    elif two_layer:
        world.add_layer(_mk("core", 0, 0.0, 0.5 * _R))
        world.add_layer(_mk("mantle", 1, 0.5 * _R, _R))
    else:
        world.add_layer(_mk("mantle", 0, 0.0, _R))

    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=2,
                          eccentricity_truncation=6, obliquity_truncation=0)
    world.solve_eos(G_to_use=G)
    return world


def _args(spin, sma):
    return (_N, spin, _ECC, 0.0, sma, _HOST)


def _default_grid():
    return dict(radii=np.linspace(0.05e6, 0.99e6, 6),
                colatitudes=np.linspace(0.2, np.pi - 0.2, 5),
                longitudes=np.linspace(0.0, 2.0 * np.pi, 4, endpoint=False))


def _heating_1d(world, spin, sma):
    world.calc_tides(orbital_frequency=_N, spin_frequency=spin, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=sma, host_mass=_HOST)
    return world.get_tidal_heating()


# =====================================================================================================================
# Default full 3D grid
# =====================================================================================================================
def test_default_grid_shape_and_longitude_independence():
    """The default secular grid is (nr, ncolat, nlon) and constant in longitude when no waves share a frequency."""
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    world = _build_world()
    grid = _default_grid()
    # At spin 1.37 n no two active waves share a frequency.
    res = world.calc_3d_tides(*_args(1.37 * _N, sma), **grid)
    heat = res['heating']
    assert heat.shape == (6, 5, 4)
    for k in range(1, 4):
        np.testing.assert_allclose(heat[:, :, k], heat[:, :, 0], rtol=1e-12, equal_nan=True)


def test_grid_longitude_mean_matches_scalar():
    """The scalar get_3d_tidal_heating is the longitude mean of a longitude-dependent secular grid."""
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    world = _build_world()
    grid = _default_grid()
    grid['longitudes'] = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    # At spin 1.5 n several degree-2 waves share a frequency, so the density varies with longitude.
    res = world.calc_3d_tides(*_args(1.5 * _N, sma), **grid)
    mean = res['heating'].mean(axis=-1)
    assert np.ptp(res['heating'][3, 2, :]) > 0.0
    for i, r in enumerate(grid['radii']):
        for j, c in enumerate(grid['colatitudes']):
            scal = world.get_3d_tidal_heating(*_args(1.5 * _N, sma), r, c)
            assert math.isclose(mean[i, j], scal, rel_tol=1e-10)


# =====================================================================================================================
# Collapsed totals + profiles
# =====================================================================================================================
@pytest.mark.parametrize("world_kwargs, extra_kwargs, num_layers, rtol", [
    pytest.param({}, {}, 1, 1e-2, id="homogeneous"),
    # Radial nodes stay inside each layer, so a few per layer already resolve the 50x stiffness jump.
    pytest.param(dict(soft_shell=True), dict(radial_slices=4), 3, 1e-4, id="soft_shell_4_slices"),
    pytest.param(dict(soft_shell=True), dict(radial_slices=16), 3, 1e-4, id="soft_shell_16_slices"),
])
def test_total_matches_1d(world_kwargs, extra_kwargs, num_layers, rtol):
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.37 * _N
    world = _build_world(**world_kwargs)
    h_1d = _heating_1d(world, spin, sma)

    res = world.calc_3d_tides(*_args(spin, sma), **extra_kwargs, **_SUMMED)
    assert res['per_layer'].shape == (num_layers,)
    assert math.isclose(res['total'], h_1d, rel_tol=rtol), \
        f"collapsed total {res['total']:.6e} != 1D {h_1d:.6e} (ratio {res['total'] / h_1d:.6f})"


def test_per_layer_sums_to_total():
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    world = _build_world(two_layer=True)
    res = world.calc_3d_tides(*_args(1.5 * _N, sma), **_SUMMED)
    assert res['per_layer'].shape == (2,)
    assert np.all(res['per_layer'] > 0.0)
    assert math.isclose(res['per_layer'].sum(), res['total'], rel_tol=1e-10)


@pytest.mark.parametrize("axis_name, axis, summed", [
    pytest.param("radii", np.linspace(1.0e4, _R, 400), dict(latitude_summed=True, longitude_summed=True),
                 id="radial"),
    pytest.param("colatitudes", np.linspace(1e-4, np.pi - 1e-4, 400), dict(longitude_summed=True, radial_summed=True),
                 id="colatitude"),
])
def test_profile_integrates_to_total(axis_name, axis, summed):
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.37 * _N
    world = _build_world()
    total = world.calc_3d_tides(*_args(spin, sma), **_SUMMED)['total']

    prof = world.calc_3d_tides(*_args(spin, sma), **{axis_name: axis}, **summed)
    assert prof['heating'].shape == axis.shape
    integ = trapezoid(prof['heating'], axis)
    assert math.isclose(integ, total, rel_tol=2e-2), f"integrated profile {integ:.4e} != total {total:.4e}"


# =====================================================================================================================
# Instantaneous (orbit_averaged=False)
# =====================================================================================================================
def test_instantaneous_time_average_matches_secular():
    """The instantaneous power averaged over the common mode period is the secular density at that point."""
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.5 * _N
    world = _build_world()
    r, colat, lon = 0.6 * _R, 1.1, 0.7

    h_bar = world.calc_3d_tides(*_args(spin, sma), radii=np.array([r]), colatitudes=np.array([colat]),
                                longitudes=np.array([lon]))['heating'][0, 0, 0]

    # Two orbits are an exact common period for spin 1.5 n, so the uniform trapezoid is exact.
    period = 2.0 * (2.0 * np.pi / _N)
    times = np.linspace(0.0, period, 4001)
    res = world.calc_3d_tides(*_args(spin, sma),
                              radii=np.array([r]), colatitudes=np.array([colat]),
                              longitudes=np.array([lon]), times=times, orbit_averaged=False)
    p = res['heating'][0, 0, 0, :]
    avg = trapezoid(p, times) / period
    assert math.isclose(avg, h_bar, rel_tol=1e-8), f"time-avg {avg:.4e} != secular {h_bar:.4e}"


def test_instantaneous_varies_with_longitude_and_time():
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    world = _build_world()
    times = np.linspace(0.0, 5.0e4, 7)
    res = world.calc_3d_tides(*_args(1.5 * _N, sma),
                              radii=np.array([0.6e6]), colatitudes=np.array([1.1]),
                              longitudes=np.array([0.0, 1.0, 2.0]), times=times, orbit_averaged=False)
    heat = res['heating']
    assert heat.shape == (1, 1, 3, 7)
    assert 'times' in res
    assert not np.allclose(heat[0, 0, 0, :], heat[0, 0, 1, :])
    assert np.ptp(heat[0, 0, 0, :]) > 0.0


def test_instantaneous_total_time_average_matches_1d():
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.5 * _N
    world = _build_world()
    h_1d = _heating_1d(world, spin, sma)

    period = 2.0 * (2.0 * np.pi / _N)   # Exact common period for spin 1.5 n.
    times = np.linspace(0.0, period, 801)
    res = world.calc_3d_tides(*_args(spin, sma), times=times, orbit_averaged=False, **_SUMMED)
    total_t = res['total']
    assert total_t.shape == times.shape
    avg = trapezoid(total_t, times) / period
    assert math.isclose(avg, h_1d, rel_tol=3e-2), f"time-avg total {avg:.4e} != 1D {h_1d:.4e}"


# =====================================================================================================================
# Guards
# =====================================================================================================================
@pytest.mark.parametrize("extra_kwargs", [
    pytest.param({}, id="longitudes_not_summed_but_missing"),
    pytest.param(dict(orbit_averaged=False, longitudes=np.array([0.0])), id="instantaneous_without_times"),
])
def test_missing_axis_raises(extra_kwargs):
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    world = _build_world()
    with pytest.raises(ValueError):
        world.calc_3d_tides(*_args(1.5 * _N, sma), radii=np.array([0.5e6]), colatitudes=np.array([1.0]),
                            **extra_kwargs)
