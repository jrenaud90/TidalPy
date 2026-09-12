"""Tests for the instantaneous 3D tidal displacement grid (LayeredWorld.calc_3d_displacements).

A uniform, effectively incompressible sphere on a circular, non-synchronous, zero-obliquity orbit has a single
active tidal mode, so the displacement field can be checked against the mode's potential amplitude and the
radial functions (y1 = h/g and y3 = l/g at the surface).
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
from TidalPy.rheology_x import Elastic, Maxwell
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.Tides_x.classes import make_tide
from TidalPy.Tides_x.potential import tidal_potential_3d_modes
from TidalPy.viscosity_x import make_viscosity

RADIUS = 1.8e6
DENSITY = 3500.0
MASS = (4.0 / 3.0) * math.pi * RADIUS**3 * DENSITY
HOST_MASS = 1.9e27
SEMI_MAJOR_AXIS = 4.2e8
ORBITAL_FREQ = 4.1e-5
SPIN_FREQ = 3.0e-5          # non-synchronous: one active mode at omega = 2 (n - Omega)
PERIOD = 2.0 * math.pi / (2.0 * abs(ORBITAL_FREQ - SPIN_FREQ))
ORBIT = dict(orbital_frequency=ORBITAL_FREQ, spin_frequency=SPIN_FREQ, eccentricity=0.0, obliquity=0.0,
             semi_major_axis=SEMI_MAJOR_AXIS, host_mass=HOST_MASS)


@pytest.fixture(scope="module")
def world():
    layer = PhysicsLayer("mantle", 0, 0.0, RADIUS, MASS, shear_modulus_static=6.0e10,
                         bulk_modulus_static=1.0e15)
    layer.set_eos(ConstantDensityEOS(reference_density=DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e19}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e30}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world = LayeredWorld("io_like", RADIUS, MASS)
    world.add_layer(layer)
    world.solve_eos()
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
    return world


def _single_mode(colatitude, longitude):
    """The one active (l=2) mode's signed frequency and complex potential row at a point."""
    degrees, freqs, pots = tidal_potential_3d_modes(
        RADIUS, ORBITAL_FREQ, SPIN_FREQ, 0.0, 0.0, SEMI_MAJOR_AXIS, HOST_MASS, G, colatitude, longitude,
        max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
    active = np.abs(freqs) > 1e-12
    assert np.count_nonzero(active) == 1
    index = int(np.flatnonzero(active)[0])
    assert degrees[index] == 2
    return float(freqs[index]), pots[index]


def _expected(y1, y3, colatitude, longitude, times):
    omega, pot = _single_mode(colatitude, longitude)
    phase = np.exp(1j * omega * np.asarray(times))
    return (np.real(y1 * pot[0] * phase),
            np.real(y3 * pot[1] * phase),
            np.real(y3 * pot[2] / math.sin(colatitude) * phase))


def test_grid_shapes_and_axes(world):
    out = world.calc_3d_displacements(**ORBIT, radii=[0.5 * RADIUS, RADIUS], colatitudes=[0.7, 1.2, 2.0],
                                      longitudes=np.linspace(0.0, 2.0 * math.pi, 5), times=[0.0, 100.0])
    for key in ("radial", "polar", "azimuthal"):
        assert out[key].shape == (2, 3, 5, 2) and out[key].dtype == np.float64
        assert np.all(np.isfinite(out[key]))
    assert out["radii"].shape == (2,) and out["times"].shape == (2,)
    scalar = world.calc_3d_displacements(**ORBIT, radii=RADIUS, colatitudes=1.2, longitudes=0.3, times=0.0)
    assert scalar["radial"].shape == (1, 1, 1, 1)


def test_interior_point_matches_radial_functions(world):
    radius, colatitude, longitude = 0.8 * RADIUS, 1.1, 0.4
    times = np.array([0.0, 0.13 * PERIOD, 0.61 * PERIOD])
    omega, _ = _single_mode(colatitude, longitude)
    world.solve_love_numbers(frequency=abs(omega), degree_l=2)
    y1 = world.get_love_radial_y(radius, 0, 0)
    y3 = world.get_love_radial_y(radius, 0, 2)
    out = world.calc_3d_displacements(**ORBIT, radii=radius, colatitudes=colatitude, longitudes=longitude,
                                      times=times)
    expected = _expected(y1, y3, colatitude, longitude, times)
    for key, reference in zip(("radial", "polar", "azimuthal"), expected):
        np.testing.assert_allclose(out[key][0, 0, 0, :], reference, rtol=1e-8, atol=1e-12)


def test_surface_displacement_follows_love_numbers(world):
    colatitude, longitude = 1.3, 2.2
    omega, _ = _single_mode(colatitude, longitude)
    result = world.solve_love_numbers(frequency=abs(omega), degree_l=2)
    gravity = float(world.get_gravity(RADIUS))
    y1 = result["love_number_h"] / gravity
    y3 = result["love_number_l"] / gravity
    times = np.linspace(0.0, PERIOD, 9)
    out = world.calc_3d_displacements(**ORBIT, radii=RADIUS, colatitudes=colatitude, longitudes=longitude,
                                      times=times)
    expected = _expected(y1, y3, colatitude, longitude, times)
    np.testing.assert_allclose(out["radial"][0, 0, 0, :], expected[0], rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(out["polar"][0, 0, 0, :], expected[1], rtol=1e-6, atol=1e-9)
    # Metre-scale tides on an Io-like body, averaging to zero over a cycle.
    assert 0.1 < np.max(np.abs(out["radial"])) < 100.0
    assert abs(np.mean(out["radial"][0, 0, 0, :-1])) < 1e-6 * np.max(np.abs(out["radial"]))


def test_center_is_nan_and_preconditions(world):
    out = world.calc_3d_displacements(**ORBIT, radii=[0.0, 0.5 * RADIUS], colatitudes=1.0, longitudes=0.0,
                                      times=0.0)
    assert np.all(np.isnan(out["radial"][0])) and np.all(np.isfinite(out["radial"][1]))
    with pytest.raises(ValueError, match="at least one"):
        world.calc_3d_displacements(**ORBIT, radii=[], colatitudes=1.0, longitudes=0.0, times=0.0)
    world.set_tide_config(max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0, love_method="homogeneous")
    try:
        with pytest.raises(RuntimeError, match="homogeneous"):
            world.calc_3d_displacements(**ORBIT, radii=RADIUS, colatitudes=1.0, longitudes=0.0, times=0.0)
    finally:
        world.set_tide_config(max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
    bare = LayeredWorld("bare", RADIUS, MASS)
    with pytest.raises(RuntimeError, match="no tide model"):
        bare.calc_3d_displacements(**ORBIT, radii=RADIUS, colatitudes=1.0, longitudes=0.0, times=0.0)
