"""3D tidal heating grid + collapse flavors (``LayeredWorld.calc_3d_tides``).

Covers the full (radius, colatitude, longitude[, time]) grid, the summed/averaged reductions, the
secular vs instantaneous quantities, and the physical consistency checks: the fully-collapsed total
equals the 1D global heating, per-layer totals sum to it, profiles integrate to it, the density grid
matches the scalar path, the secular grid is longitude-independent, and the instantaneous power
orbit-averages back to the secular density.
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.utilities.conversions import orbital_motion2semi_a


_R = 1.0e6
_DENSITY = 5000.0
_SHEAR = 5.0e10
_BULK = 1.0e11
_VISC = 1.0e19
_N = 2.0 * np.pi / 86400.0
_ECC = 0.05
_HOST = mass_trap1
_MASS = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY


def _build_world(two_layer=False):
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.rheology_x.rheology import Maxwell, Elastic
    from TidalPy.Tides_x.classes.tide import make_tide

    world = LayeredWorld("w", _R, _MASS)

    def _mk(name, idx, r_in, r_out):
        mass = (4.0 / 3.0) * math.pi * (r_out ** 3 - r_in ** 3) * _DENSITY
        layer = PhysicsLayer(name, idx, r_in, r_out, mass,
                             shear_modulus_static_pa=_SHEAR, bulk_modulus_static_pa=_BULK)
        layer.is_static = False
        layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
        layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
        layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
        layer.set_shear_rheology(Maxwell())
        layer.set_bulk_rheology(Elastic())
        return layer

    if two_layer:
        world.add_layer(_mk("core", 0, 0.0, 0.5 * _R))
        world.add_layer(_mk("mantle", 1, 0.5 * _R, _R))
    else:
        world.add_layer(_mk("mantle", 0, 0.0, _R))

    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=2,
                          eccentricity_truncation=3, obliquity_truncation=0)
    world.solve_eos(G_to_use=G)
    return world


def _args(spin, sma):
    return (_N, spin, _ECC, 0.0, sma, _HOST)


def _default_grid():
    return dict(radii=np.linspace(0.05e6, 0.99e6, 6),
                colatitudes=np.linspace(0.2, np.pi - 0.2, 5),
                longitudes=np.linspace(0.0, 2.0 * np.pi, 4, endpoint=False))


# =====================================================================================================================
# Default full 3D grid
# =====================================================================================================================
def test_default_grid_shape_and_longitude_independence():
    """Default output is the secular density grid (nr, ncolat, nlon); constant along longitude."""
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    world = _build_world()
    grid = _default_grid()
    res = world.calc_3d_tides(*_args(1.37 * _N, sma), **grid)
    heat = res['heating']
    assert heat.shape == (6, 5, 4)
    # Secular density is longitude-independent: every longitude slice is identical.
    for k in range(1, 4):
        np.testing.assert_allclose(heat[:, :, k], heat[:, :, 0], rtol=1e-12, equal_nan=True)


def test_grid_matches_scalar():
    """The uncollapsed secular grid equals the scalar get_3d_tidal_heating point-for-point."""
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    world = _build_world()
    grid = _default_grid()
    res = world.calc_3d_tides(*_args(1.5 * _N, sma), **grid)
    for i, r in enumerate(grid['radii']):
        for j, c in enumerate(grid['colatitudes']):
            scal = world.get_3d_tidal_heating(*_args(1.5 * _N, sma), r, c)
            assert math.isclose(res['heating'][i, j, 0], scal, rel_tol=1e-12)


# =====================================================================================================================
# Collapsed totals + profiles
# =====================================================================================================================
def test_total_matches_1d():
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.37 * _N
    world = _build_world()
    world.calc_tides(orbital_frequency=_N, spin_frequency=spin, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=sma, host_mass=_HOST)
    h_1d = world.get_tidal_heating()

    res = world.calc_3d_tides(*_args(spin, sma),
                              latitude_summed=True, longitude_summed=True, radial_summed=True)
    assert math.isclose(res['total'], h_1d, rel_tol=1e-2), \
        f"collapsed total {res['total']:.4e} != 1D {h_1d:.4e}"


def test_per_layer_sums_to_total():
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    world = _build_world(two_layer=True)
    res = world.calc_3d_tides(*_args(1.5 * _N, sma),
                              latitude_summed=True, longitude_summed=True, radial_summed=True)
    assert res['per_layer'].shape == (2,)
    assert np.all(res['per_layer'] > 0.0)
    assert math.isclose(res['per_layer'].sum(), res['total'], rel_tol=1e-10)


def test_radial_profile_integrates_to_total():
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.37 * _N
    world = _build_world()
    total = world.calc_3d_tides(*_args(spin, sma),
                                latitude_summed=True, longitude_summed=True, radial_summed=True)['total']

    radii = np.linspace(1.0e3, _R, 400)
    prof = world.calc_3d_tides(*_args(spin, sma), radii=radii,
                               latitude_summed=True, longitude_summed=True)
    assert prof['heating'].shape == radii.shape
    integ = np.trapezoid(prof['heating'], radii)
    assert math.isclose(integ, total, rel_tol=2e-2), f"int dP/dr {integ:.4e} != total {total:.4e}"


def test_colatitude_profile_integrates_to_total():
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.37 * _N
    world = _build_world()
    total = world.calc_3d_tides(*_args(spin, sma),
                                latitude_summed=True, longitude_summed=True, radial_summed=True)['total']

    colat = np.linspace(1e-4, np.pi - 1e-4, 400)
    prof = world.calc_3d_tides(*_args(spin, sma), colatitudes=colat,
                               longitude_summed=True, radial_summed=True)
    assert prof['heating'].shape == colat.shape
    integ = np.trapezoid(prof['heating'], colat)
    assert math.isclose(integ, total, rel_tol=2e-2), f"int dP/dtheta {integ:.4e} != total {total:.4e}"


# =====================================================================================================================
# Instantaneous (orbit_averaged=False): sigma:eps_dot, orbit-averages to the secular density
# =====================================================================================================================
def test_instantaneous_time_average_matches_secular():
    """Averaging the instantaneous power over the mode common period recovers the secular density.

    For spin = 1.5 n the mode frequencies are all multiples of 0.5 n, so two orbital periods is an exact
    common period; averaging sigma:eps_dot over it must return the secular density h_bar.
    """
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.5 * _N
    world = _build_world()
    r, colat, lon = 0.6 * _R, 1.1, 0.7

    h_bar = world.get_3d_tidal_heating(*_args(spin, sma), r, colat)

    period = 2.0 * (2.0 * np.pi / _N)   # two orbital periods = exact common period for spin = 1.5 n
    times = np.linspace(0.0, period, 4001)
    res = world.calc_3d_tides(*_args(spin, sma),
                              radii=np.array([r]), colatitudes=np.array([colat]),
                              longitudes=np.array([lon]), times=times, orbit_averaged=False)
    p = res['heating'][0, 0, 0, :]
    avg = np.trapezoid(p, times) / period
    assert math.isclose(avg, h_bar, rel_tol=1e-2), f"time-avg {avg:.4e} != secular {h_bar:.4e}"


def test_instantaneous_varies_with_longitude_and_time():
    """Unlike the secular density, the instantaneous power depends on longitude and time."""
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    world = _build_world()
    times = np.linspace(0.0, 5.0e4, 7)
    res = world.calc_3d_tides(*_args(1.5 * _N, sma),
                              radii=np.array([0.6e6]), colatitudes=np.array([1.1]),
                              longitudes=np.array([0.0, 1.0, 2.0]), times=times, orbit_averaged=False)
    heat = res['heating']
    assert heat.shape == (1, 1, 3, 7)
    assert 'times' in res
    # Varies across longitude and across time.
    assert not np.allclose(heat[0, 0, 0, :], heat[0, 0, 1, :])
    assert np.ptp(heat[0, 0, 0, :]) > 0.0


def test_instantaneous_total_time_average_matches_1d():
    """The instantaneous whole-planet power, time-averaged over a period, equals the 1D heating."""
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.5 * _N
    world = _build_world()
    world.calc_tides(orbital_frequency=_N, spin_frequency=spin, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=sma, host_mass=_HOST)
    h_1d = world.get_tidal_heating()

    period = 2.0 * (2.0 * np.pi / _N)   # two orbital periods = exact common period for spin = 1.5 n
    times = np.linspace(0.0, period, 801)
    res = world.calc_3d_tides(*_args(spin, sma), times=times, orbit_averaged=False,
                              latitude_summed=True, longitude_summed=True, radial_summed=True)
    total_t = res['total']
    assert total_t.shape == times.shape
    avg = np.trapezoid(total_t, times) / period
    assert math.isclose(avg, h_1d, rel_tol=3e-2), f"time-avg total {avg:.4e} != 1D {h_1d:.4e}"


# =====================================================================================================================
# Guards
# =====================================================================================================================
def test_missing_axis_raises():
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    world = _build_world()
    with pytest.raises(ValueError):  # longitude not summed but no longitudes
        world.calc_3d_tides(*_args(1.5 * _N, sma),
                            radii=np.array([0.5e6]), colatitudes=np.array([1.0]))
    with pytest.raises(ValueError):  # instantaneous but no times
        world.calc_3d_tides(*_args(1.5 * _N, sma), orbit_averaged=False,
                            radii=np.array([0.5e6]), colatitudes=np.array([1.0]),
                            longitudes=np.array([0.0]))
