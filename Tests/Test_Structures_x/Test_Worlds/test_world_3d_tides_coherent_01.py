"""Coherent summation of same-frequency tidal waves in the secular 3D heating.

The secular (cycle-averaged) heating at one frequency is ``(|omega|/2) Im(sigma_c : conj(eps_c))`` of the TOTAL
complex amplitude at that frequency, so waves that share a frequency must be summed before the bilinear form.
The ``m = 0`` modes always come in such pairs (``(l, 0, p, q)`` at ``+omega`` and ``(l, 0, l-p, -q)`` at ``-omega``
are the same real sinusoid), and for a synchronously rotating body every active mode sits at a multiple of the
mean motion. Summing the modes' powers separately undercounts the zonal terms by 2x, which for a homogeneous
degree-2 body at zero obliquity is 4.5/84 = 5.36% of the total at synchronous rotation.

These tests pin the coherent behaviour: the collapsed 3D total equals the 1D heating at synchronous rotation,
the analytic and Gauss-Legendre colatitude collapses agree when waves of different degree share a frequency,
the secular grid is the exact time average of the instantaneous power point by point (including its longitude
dependence), and its longitude mean is what the scalar and batch paths return.
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.Utilities_x.conversions import orbital_motion2semi_a


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
_SYNC = (_N, _N, _ECC, 0.0, _SMA, _HOST)   # synchronous rotation: every active mode at a multiple of n


def _build_world(max_degree_l=2):
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.rheology_x.rheology import Maxwell, Elastic
    from TidalPy.Tides_x.classes.tide import make_tide

    world = LayeredWorld("w", _R, _MASS)
    layer = PhysicsLayer("mantle", 0, 0.0, _R, _MASS,
                         shear_modulus_static=_SHEAR, bulk_modulus_static=_BULK)
    layer.is_static = False
    layer.set_eos(ConstantDensityEOS(reference_density=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world.add_layer(layer)
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=max_degree_l,
                          eccentricity_truncation=3, obliquity_truncation=0)
    world.solve_eos(G_to_use=G)
    return world


def _heating_1d(world):
    world.calc_tides(orbital_frequency=_N, spin_frequency=_N, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=_SMA, host_mass=_HOST)
    return world.get_tidal_heating()


@pytest.mark.parametrize("max_degree_l", [2, 4])
def test_synchronous_total_matches_1d(max_degree_l):
    """At synchronous rotation the collapsed secular total equals the 1D heating (radial trapezoid accuracy).

    Degree 4 adds waves that share frequencies with the degree-2 ones (cross-degree coherent terms, which
    integrate to zero over the volume).
    """
    world = _build_world(max_degree_l=max_degree_l)
    h_1d = _heating_1d(world)
    total = world.calc_3d_tides(*_SYNC, latitude_summed=True, longitude_summed=True, radial_summed=True,
                                radial_slices=64)['total']
    assert math.isclose(total, h_1d, rel_tol=1.0e-3), \
        f"synchronous 3D total {total:.6e} != 1D {h_1d:.6e} (ratio {total / h_1d:.5f})"


def test_analytic_matches_quadrature_with_shared_frequencies():
    """The analytic collapse (cross-degree Gram integrals) equals the Gauss-Legendre collapse at synchronous
    rotation with degrees 2 to 4 active, for the total and the radial profile."""
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


def test_secular_grid_is_time_average_of_instantaneous():
    """The secular grid equals the time average of the instantaneous power point by point, longitude included.

    At synchronous rotation one orbital period is the exact common period of every active mode, and the
    instantaneous power is a trigonometric polynomial in time, so a uniform trapezoid over the period is exact.
    """
    world = _build_world()
    r, colat = 0.9 * _R, 1.0
    lons = np.array([0.0, 0.5, 1.0, 0.5 * np.pi, 2.5])
    secular = world.calc_3d_tides(*_SYNC, radii=np.array([r]), colatitudes=np.array([colat]),
                                  longitudes=lons)['heating'][0, 0, :]
    period = 2.0 * np.pi / _N
    times = np.linspace(0.0, period, 2001)
    inst = world.calc_3d_tides(*_SYNC, radii=np.array([r]), colatitudes=np.array([colat]),
                               longitudes=lons, times=times, orbit_averaged=False)['heating'][0, 0, :, :]
    averaged = np.trapezoid(inst, times, axis=-1) / period
    np.testing.assert_allclose(averaged, secular, rtol=1.0e-8)


def test_synchronous_heating_varies_with_longitude():
    """Waves at one frequency with different longitude structure make the secular heating longitude-dependent:
    at synchronous rotation the pattern is symmetric about the sub-host meridian and far from flat."""
    world = _build_world()
    lons = np.linspace(0.0, 2.0 * np.pi, 9)
    secular = world.calc_3d_tides(*_SYNC, radii=np.array([0.9 * _R]), colatitudes=np.array([1.0]),
                                  longitudes=lons)['heating'][0, 0, :]
    assert secular.max() > 1.3 * secular.min()
    # periodic, and symmetric under phi -> 2 pi - phi (index k <-> 8 - k)
    assert math.isclose(secular[0], secular[-1], rel_tol=1.0e-10)
    np.testing.assert_allclose(secular[1:], secular[-2::-1], rtol=1.0e-10)


def test_scalar_and_batch_are_the_longitude_mean():
    """The scalar and batch paths return the longitude mean of the secular grid (a uniform longitude grid
    integrates its trigonometric polynomial exactly)."""
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


def test_zonal_pairs_are_one_wave():
    """The m = 0 modes at +q n and -q n are one real sinusoid: the 1D formula counts the pair through its
    (2 - delta_m0) weighting, and the coherent 3D sum reproduces it. Dropping the coherence loses 2x on the
    zonal terms, i.e. 4.5/84 of a degree-2, zero-obliquity total at leading order in e; this pins the 3D/1D
    ratio well inside that gap."""
    world = _build_world()
    h_1d = _heating_1d(world)
    total = world.calc_3d_tides(*_SYNC, latitude_summed=True, longitude_summed=True, radial_summed=True,
                                radial_slices=64)['total']
    assert abs(total / h_1d - 1.0) < 0.2 * (4.5 / 84.0)
