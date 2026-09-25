"""The 3D heating path agrees with itself and with the 1D path in the corner cases.

* The analytic colatitude integral is of the longitude mean, so a call that keeps its longitudes must use the
  quadrature and keep the longitude dependence.
* The 3D path keeps the same modes as the 1D path, including a slow mode near a Maxwell peak.
* A liquid point carries no shear dissipation: its heating is 0, not NaN.
* The two poles are treated alike.
"""
import math

import numpy as np
import pytest

from TidalPy.structures_x import build_world

_IO_N = 4.11e-5
_IO_STATE = dict(orbital_frequency=_IO_N, spin_frequency=_IO_N, eccentricity=0.0041, obliquity=0.0,
                 semi_major_axis=4.217e8, host_mass=1.898e27)


@pytest.fixture(scope="module")
def io():
    world = build_world("io")
    world.solve_eos()
    return world


def test_kept_longitudes_keep_their_dependence(io):
    """Latitude summed, longitudes kept: the analytic default and the quadrature agree, and vary with longitude."""
    radii = np.array([1.75e6])
    longitudes = np.linspace(0.0, math.pi, 7)
    analytic = io.calc_3d_tides(**_IO_STATE, radii=radii, longitudes=longitudes, latitude_summed=True)["heating"]
    quadrature = io.calc_3d_tides(**_IO_STATE, radii=radii, longitudes=longitudes, latitude_summed=True,
                                  latitude_analytic=False)["heating"]
    np.testing.assert_allclose(analytic, quadrature, rtol=1e-12)
    assert np.ptp(analytic) > 0.1 * np.mean(analytic)


def test_a_slow_mode_near_synchronous_is_kept_by_both_paths(io):
    """At a spin just off synchronous the slow semi-diurnal mode is in both the 1D and the 3D totals."""
    state = dict(_IO_STATE, spin_frequency=_IO_N + 3.0e-11)
    io.calc_tides(**state)
    total_1d = io.get_tidal_heating()
    total_3d = io.calc_3d_tides(**state, latitude_summed=True, longitude_summed=True, radial_summed=True)["total"]
    assert total_3d == pytest.approx(total_1d, rel=1e-3)


def test_liquid_points_have_no_heating():
    """A point in a liquid layer takes no shear dissipation: zero heating, while its stress and strain are NaN."""
    mercury = build_world("mercury")
    mercury.solve_eos()
    n = 2.0 * math.pi / (87.9691 * 86400.0)
    state = dict(orbital_frequency=n, spin_frequency=1.5 * n, eccentricity=0.2056, obliquity=0.0,
                 semi_major_axis=5.7909e10, host_mass=1.989e30)
    core = mercury.outer_core
    radius = 0.5 * (core.radius_inner + core.radius_outer)
    assert not core.is_solid
    assert mercury.get_3d_tidal_heating(radius=radius, colatitude=1.0, **state) == 0.0
    mantle_radius = 0.5 * (mercury.mantle.radius_inner + mercury.mantle.radius_outer)
    assert mercury.get_3d_tidal_heating(radius=mantle_radius, colatitude=1.0, **state) > 0.0


def test_both_poles_are_treated_alike(io):
    heating = io.get_3d_tidal_heating_array(
        **_IO_STATE, radii=np.array([1.75e6, 1.75e6, 1.75e6]), colatitudes=np.array([0.0, math.pi, 1.0]))
    assert math.isnan(heating[0]) == math.isnan(heating[1])
    assert math.isfinite(heating[2])
