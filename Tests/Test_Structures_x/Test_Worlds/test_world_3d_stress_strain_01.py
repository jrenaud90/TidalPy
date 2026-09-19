"""Instantaneous tidal stress and strain grids on a layered world (``LayeredWorld.calc_3d_stress_strain``).

The grids are checked two independent ways. Every value equals the point-wise kernel helpers of
``Tides_x.multilayer.stress_strain`` assembled mode by mode, and the strain equals the symmetric gradient of the
world's own displacement grid by central differences, which checks dy1/dr and the corrected angular strain forms
against the displacement field. Also covered: the output layout, skipping one tensor, NaN at a radius without a
depth-resolved solution, and input validation.
"""
import math

import numpy as np
import pytest

import TidalPy.constants as tidalpy_constants
from TidalPy.constants import G, mass_trap1
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
from TidalPy.rheology_x.rheology import Elastic, Maxwell
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.structures_x.worlds.layered import STRESS_STRAIN_COMPONENTS, LayeredWorld
from TidalPy.Tides_x.classes.tide import make_tide
from TidalPy.Tides_x.multilayer.stress_strain import strain_stress_heating_point
from TidalPy.Tides_x.potential import tidal_potential_3d_modes
from TidalPy.Utilities_x.conversions import orbital_motion2semi_a
from TidalPy.viscosity_x import make_viscosity


_R = 1.0e6
_DENSITY = 5000.0
_SHEAR = 5.0e10
_BULK = 1.0e11
_VISC = 1.0e19
_N = 2.0 * np.pi / 86400.0
_ECC = 0.05
_ECC_TRUNCATION = 3
_HOST = mass_trap1
_MASS = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY
_SMA = orbital_motion2semi_a(_N, _HOST, _MASS)

# (spin / mean motion, maximum degree, obliquity [rad], obliquity truncation): synchronous with frequency-sharing
# waves, non-synchronous, two degrees, and an obliquity case with m = 1 waves.
_CASES = [(1.0, 2, 0.0, 0), (1.5, 2, 0.0, 0), (1.5, 3, 0.0, 0), (1.2, 2, 0.2, 2)]


def _build_world(max_degree_l=2, obliquity_truncation=0, tide_model="rheology", solve_eos=True):
    world = LayeredWorld("w", _R, _MASS)
    layer = PhysicsLayer(
        "mantle",
        0,
        0.0,
        _R,
        _MASS,
        shear_modulus_static=_SHEAR,
        bulk_modulus_static=_BULK)
    layer.is_static = False
    layer.set_eos(ConstantDensityEOS(reference_density=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _VISC}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world.add_layer(layer)
    world.set_tide_model(make_tide(tide_model))
    world.set_tide_config(
        min_degree_l=2,
        max_degree_l=max_degree_l,
        eccentricity_truncation=_ECC_TRUNCATION,
        obliquity_truncation=obliquity_truncation)
    if solve_eos:
        world.solve_eos(G_to_use=G)
    return world


def _state(spin_ratio, obliquity):
    return (_N, spin_ratio * _N, _ECC, obliquity, _SMA, _HOST)


def _helper_tensors(world, spin_ratio, obliquity, max_degree_l, obliquity_truncation, point, times):
    """Stress and strain, shaped (ntime, 6), at one (radius, colatitude, longitude) point from the point helpers."""
    radius, colatitude, longitude = point
    degrees, frequencies, rows = tidal_potential_3d_modes(
        _R,
        _N,
        spin_ratio * _N,
        _ECC,
        obliquity,
        _SMA,
        _HOST,
        G,
        colatitude,
        longitude,
        min_degree_l=2,
        max_degree_l=max_degree_l,
        eccentricity_truncation=_ECC_TRUNCATION,
        obliquity_truncation=obliquity_truncation)
    stress = np.zeros((times.size, 6))
    strain = np.zeros((times.size, 6))
    for index in range(degrees.size):
        frequency = float(frequencies[index])
        magnitude = abs(frequency)
        if magnitude <= tidalpy_constants.min_spin_orbit_diff:
            continue
        world.solve_love_numbers(frequency=magnitude, degree_l=int(degrees[index]))
        y = np.array([world.get_love_radial_y(radius, 0, y_index) for y_index in range(6)], dtype=np.complex128)
        strain_amplitude, stress_amplitude, _ = strain_stress_heating_point(
            y,
            complex(world.calc_complex_shear_modulus(radius, magnitude)),
            complex(world.calc_complex_bulk_modulus(radius, magnitude)),
            radius,
            float(degrees[index]),
            True,
            False,
            rows[index] if frequency > 0.0 else np.conj(rows[index]),
            colatitude)
        phasor = np.exp(1j * magnitude * times)[:, np.newaxis]
        stress += np.real(stress_amplitude[np.newaxis, :] * phasor)
        strain += np.real(strain_amplitude[np.newaxis, :] * phasor)
    return stress, strain


# =====================================================================================================================
# Layout and validation
# =====================================================================================================================
def test_output_layout_and_center():
    world = _build_world()
    radii = np.array([0.0, 0.5 * _R, _R])
    colatitudes = np.array([0.4, 1.2])
    longitudes = np.array([0.0, 1.0, 2.0])
    times = np.array([0.0, 1.0e4])
    result = world.calc_3d_stress_strain(
        *_state(1.0, 0.0),
        radii=radii,
        colatitudes=colatitudes,
        longitudes=longitudes,
        times=times)
    assert result["components"] == ("rr", "theta_theta", "phi_phi", "r_theta", "r_phi", "theta_phi")
    assert result["components"] == STRESS_STRAIN_COMPONENTS
    for axis_name, axis in (("radii", radii), ("colatitudes", colatitudes), ("longitudes", longitudes),
                            ("times", times)):
        np.testing.assert_array_equal(result[axis_name], axis)
    for name in ("stress", "strain"):
        assert result[name].shape == (3, 2, 3, 2, 6)
        assert result[name].dtype == np.float64
        assert np.all(np.isnan(result[name][0]))       # The center has no depth-resolved solution
        assert np.all(np.isfinite(result[name][1:]))
        assert np.any(result[name][1:] != 0.0)


def test_skipping_a_tensor_leaves_the_other_unchanged():
    world = _build_world()
    grid = dict(radii=[0.0, 0.6 * _R], colatitudes=[0.5, 1.5], longitudes=[0.3, 2.0], times=[0.0, 3.0e4])
    both = world.calc_3d_stress_strain(*_state(1.5, 0.0), **grid)
    stress_only = world.calc_3d_stress_strain(*_state(1.5, 0.0), return_strain=False, **grid)
    strain_only = world.calc_3d_stress_strain(*_state(1.5, 0.0), return_stress=False, **grid)
    assert "strain" not in stress_only and "stress" not in strain_only
    np.testing.assert_array_equal(stress_only["stress"], both["stress"])
    np.testing.assert_array_equal(strain_only["strain"], both["strain"])


def test_invalid_requests_raise():
    grid = dict(radii=[0.5 * _R], colatitudes=[1.0], longitudes=[0.0], times=[0.0])
    world = _build_world()
    with pytest.raises(ValueError):
        world.calc_3d_stress_strain(*_state(1.0, 0.0), return_stress=False, return_strain=False, **grid)
    with pytest.raises(ValueError):
        world.calc_3d_stress_strain(*_state(1.0, 0.0), **dict(grid, times=[]))
    with pytest.raises(RuntimeError):
        _build_world(solve_eos=False).calc_3d_stress_strain(*_state(1.0, 0.0), **grid)
    with pytest.raises(RuntimeError):
        _build_world(tide_model="cpl").calc_3d_stress_strain(*_state(1.0, 0.0), **grid)


# =====================================================================================================================
# Physics
# =====================================================================================================================
@pytest.mark.parametrize("spin_ratio, max_degree_l, obliquity, obliquity_truncation", _CASES)
def test_grid_matches_point_helpers(spin_ratio, max_degree_l, obliquity, obliquity_truncation):
    """Every grid value equals the point-wise kernel helpers assembled mode by mode at that point and time."""
    world = _build_world(max_degree_l, obliquity_truncation)
    point = (0.8 * _R, 1.1, 0.7)
    times = np.array([0.0, 1.3e4, 5.1e4])
    grid = world.calc_3d_stress_strain(
        *_state(spin_ratio, obliquity),
        radii=[point[0]],
        colatitudes=[point[1]],
        longitudes=[point[2]],
        times=times)
    stress, strain = _helper_tensors(world, spin_ratio, obliquity, max_degree_l, obliquity_truncation, point, times)
    for name, expected in (("stress", stress), ("strain", strain)):
        np.testing.assert_allclose(
            grid[name][0, 0, 0],
            expected,
            rtol=1.0e-9,
            atol=1.0e-12 * np.max(np.abs(expected)),
            err_msg=name)


@pytest.mark.parametrize("spin_ratio, obliquity, obliquity_truncation", [(1.0, 0.0, 0), (1.2, 0.2, 2)])
def test_strain_is_symmetric_gradient_of_displacement(spin_ratio, obliquity, obliquity_truncation):
    """The strain grid equals the symmetric gradient of the displacement grid in spherical coordinates, taken by
    central differences in radius, colatitude, and longitude."""
    world = _build_world(2, obliquity_truncation)
    state = _state(spin_ratio, obliquity)
    radius, colatitude, longitude, time = 0.7 * _R, 1.0, 0.6, 2.0e4
    radius_step, angle_step = 1.0e-3 * _R, 1.0e-4
    offsets = np.array([-1.0, 0.0, 1.0])
    displacement = world.calc_3d_displacements(
        *state,
        radii=radius + radius_step * offsets,
        colatitudes=colatitude + angle_step * offsets,
        longitudes=longitude + angle_step * offsets,
        times=[time])
    u_r, u_theta, u_phi = (displacement[name][..., 0] for name in ("radial", "polar", "azimuthal"))

    def d_radius(u):
        return (u[2, 1, 1] - u[0, 1, 1]) / (2.0 * radius_step)

    def d_colatitude(u):
        return (u[1, 2, 1] - u[1, 0, 1]) / (2.0 * angle_step)

    def d_longitude(u):
        return (u[1, 1, 2] - u[1, 1, 0]) / (2.0 * angle_step)

    sin_theta = np.sin(colatitude)
    cot_theta = np.cos(colatitude) / sin_theta
    ur, ut, up = u_r[1, 1, 1], u_theta[1, 1, 1], u_phi[1, 1, 1]
    expected = np.array([
        d_radius(u_r),
        (d_colatitude(u_theta) + ur) / radius,
        (d_longitude(u_phi) / sin_theta + ur + ut * cot_theta) / radius,
        0.5 * (d_colatitude(u_r) / radius + d_radius(u_theta) - ut / radius),
        0.5 * (d_longitude(u_r) / (radius * sin_theta) + d_radius(u_phi) - up / radius),
        0.5 * (d_longitude(u_theta) / (radius * sin_theta) + d_colatitude(u_phi) / radius - up * cot_theta / radius),
    ])

    strain = world.calc_3d_stress_strain(
        *state,
        radii=[radius],
        colatitudes=[colatitude],
        longitudes=[longitude],
        times=[time],
        return_stress=False)["strain"][0, 0, 0, 0]
    np.testing.assert_allclose(strain, expected, rtol=1.0e-4, atol=1.0e-6 * np.max(np.abs(expected)))
