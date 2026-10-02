"""Coherent summation of same-frequency tidal waves in the secular 3D heating, at synchronous rotation."""
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
_SMA = orbital_motion2semi_a(_N, _HOST, _MASS)
_SYNC = (_N, _N, _ECC, 0.0, _SMA, _HOST)   # Synchronous rotation: every active mode sits at a multiple of n.


def _material():
    return constant_solid(
        _DENSITY, bulk_modulus=_BULK, shear_modulus=_SHEAR, shear_viscosity=_VISC, bulk_viscosity=_VISC)


def _build_world(max_degree_l=2):
    world = BaseWorld("w", _R, _MASS)
    world.add_layer(Layer("mantle", 0, 0.0, _R, _MASS, _material(), is_static=False, shear_rheology=Maxwell(),
                          bulk_rheology=Elastic()))
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=max_degree_l,
                          eccentricity_truncation=6, obliquity_truncation=0)
    world.solve_eos(G_to_use=G)
    return world


def _heating_1d(world):
    world.calc_tides(orbital_frequency=_N, spin_frequency=_N, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=_SMA, host_mass=_HOST)
    return world.get_tidal_heating()


@pytest.mark.parametrize("max_degree_l", [2, 4])
def test_synchronous_total_matches_1d(max_degree_l):
    """At synchronous rotation the collapsed secular total equals the 1D heating, zonal pairs included."""
    world = _build_world(max_degree_l=max_degree_l)
    h_1d = _heating_1d(world)
    total = world.calc_3d_tides(*_SYNC, latitude_summed=True, longitude_summed=True, radial_summed=True,
                                radial_slices=64)['total']
    assert math.isclose(total, h_1d, rel_tol=1.0e-3), \
        f"synchronous 3D total {total:.6e} != 1D {h_1d:.6e} (ratio {total / h_1d:.5f})"
    # Summing the m = 0 pair incoherently would lose 4.5/84 of a degree-2 total; stay well inside that gap.
    assert abs(total / h_1d - 1.0) < 0.2 * (4.5 / 84.0)


def test_analytic_matches_quadrature_with_shared_frequencies():
    """With degrees 2 to 4 sharing frequencies, analytic and quadrature collapses agree for total and profile."""
    world = _build_world(max_degree_l=4)
    kw = dict(latitude_summed=True, longitude_summed=True, radial_summed=True)
    analytic = world.calc_3d_tides(*_SYNC, latitude_analytic=True, **kw)['total']
    numeric = world.calc_3d_tides(*_SYNC, latitude_analytic=False, latitude_nodes=64, **kw)['total']
    assert math.isclose(analytic, numeric, rel_tol=1.0e-10), f"analytic {analytic:.6e} != numeric {numeric:.6e}"

    radii = np.linspace(1.0e3, _R, 120)
    kw = dict(radii=radii, latitude_summed=True, longitude_summed=True)
    analytic = world.calc_3d_tides(*_SYNC, latitude_analytic=True, **kw)['heating']
    numeric = world.calc_3d_tides(*_SYNC, latitude_analytic=False, latitude_nodes=64, **kw)['heating']
    np.testing.assert_allclose(analytic, numeric, rtol=1.0e-9)


@pytest.mark.parametrize("max_degree_l", [2, 4])
def test_secular_grid_is_time_average_of_instantaneous(max_degree_l):
    """The secular grid equals the one-period time average of the instantaneous power, longitude included."""
    world = _build_world(max_degree_l=max_degree_l)
    # Secular and instantaneous differ by terms past e^N; level 20 puts those far below the tolerance.
    world.set_tide_config(eccentricity_truncation=20)
    r, colat = 0.9 * _R, 1.0
    lons = np.array([0.0, 0.5, 1.0, 0.5 * np.pi, 2.5])
    secular = world.calc_3d_tides(*_SYNC, radii=np.array([r]), colatitudes=np.array([colat]),
                                  longitudes=lons)['heating'][0, 0, :]
    # One orbit is the exact common period of every synchronous mode, so the uniform trapezoid is exact.
    period = 2.0 * np.pi / _N
    times = np.linspace(0.0, period, 2001)
    inst = world.calc_3d_tides(*_SYNC, radii=np.array([r]), colatitudes=np.array([colat]),
                               longitudes=lons, times=times, orbit_averaged=False)['heating'][0, 0, :, :]
    averaged = trapezoid(inst, times, axis=-1) / period
    np.testing.assert_allclose(averaged, secular, rtol=1.0e-8)


def test_synchronous_heating_varies_with_longitude():
    """The synchronous secular heating varies with longitude, symmetric about the sub-host meridian."""
    world = _build_world()
    lons = np.linspace(0.0, 2.0 * np.pi, 9)
    secular = world.calc_3d_tides(*_SYNC, radii=np.array([0.9 * _R]), colatitudes=np.array([1.0]),
                                  longitudes=lons)['heating'][0, 0, :]
    assert secular.max() > 1.3 * secular.min()
    assert math.isclose(secular[0], secular[-1], rel_tol=1.0e-10)
    np.testing.assert_allclose(secular[1:], secular[-2::-1], rtol=1.0e-10)


def test_scalar_and_batch_are_the_longitude_mean():
    """The scalar and batch paths return the longitude mean of the secular grid."""
    world = _build_world()
    radii = np.array([0.5 * _R, 0.9 * _R])
    colats = np.array([0.7, 1.9])
    lons = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
    grid = world.calc_3d_tides(*_SYNC, radii=radii, colatitudes=colats, longitudes=lons)['heating']
    mean = grid.mean(axis=-1)
    for i, r in enumerate(radii):
        for j, c in enumerate(colats):
            scalar = world.get_3d_tidal_heating(*_SYNC, r, c)
            assert math.isclose(scalar, mean[i, j], rel_tol=1.0e-10)
    batch = world.get_3d_tidal_heating_array(*_SYNC, np.repeat(radii, 2), np.tile(colats, 2))
    np.testing.assert_allclose(batch, mean.ravel(), rtol=1.0e-12)
