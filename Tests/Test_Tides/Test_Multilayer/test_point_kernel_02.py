"""Point-wise 3D kernels with complex potential rows, and rebuilding the world's 3D heating and displacements."""
import math

import numpy as np
import pytest

import TidalPy.constants as tidalpy_constants
from TidalPy.constants import G, mass_trap1
from TidalPy.Material import Material
from TidalPy.Rheology.rheology import Elastic, Maxwell
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Tides.multilayer.stress_strain import displacement_point, strain_stress_heating_point, volumetric_heating
from TidalPy.Tides.potential import tidal_potential_3d_modes
from TidalPy.Utilities.conversions import orbital_motion2semi_a

from shared_materials import constant_solid


_R = 1.0e6
_DENSITY = 5000.0
_SHEAR = 5.0e10
_BULK = 1.0e11
_VISC = 1.0e19
_N = 2.0 * np.pi / 86400.0
_ECC = 0.05
# The helpers use the unsquared potential while the world cuts products at e^N; at level 20 the gap is negligible.
_ECC_TRUNCATION = 20
_HOST = mass_trap1
_MASS = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY
_SMA = orbital_motion2semi_a(_N, _HOST, _MASS)

# Every potential entry has a distinct real and imaginary part.
_Y = np.array([1.5 - 0.2j, 3.0 + 0.4j, 0.7 + 0.1j, 2.0 - 0.3j, 0.5, 0.1], dtype=np.complex128)
_ROW = np.array([10.0 + 2.0j, -3.0 + 1.0j, 4.0j, 1.0 - 2.0j, -2.0 + 0.0j, 0.5 + 0.5j], dtype=np.complex128)
_COLATITUDE = 1.1

# (spin / mean motion, maximum degree): synchronous with frequency-sharing waves, non-synchronous, and two degrees.
_WORLD_CASES = [(1.0, 2), (1.5, 2), (1.5, 3)]


def _kernel(row):
    return strain_stress_heating_point(
        _Y,
        complex(_SHEAR, 1.0e8),
        complex(_BULK, 0.0),
        0.8 * _R,
        2.0,
        _N,
        True,
        False,
        row,
        _COLATITUDE)


def _material():
    return constant_solid(
        _DENSITY, bulk_modulus=_BULK, shear_modulus=_SHEAR, shear_viscosity=_VISC, bulk_viscosity=_VISC)


def _build_world(max_degree_l):
    world = BaseWorld("w", _R, _MASS)
    layer = Layer(
        "mantle",
        0,
        0.0,
        _R,
        _MASS,
        _material(),
        is_static=False,
        shear_rheology=Maxwell(),
        bulk_rheology=Elastic())
    world.add_layer(layer)
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(
        min_degree_l=2,
        max_degree_l=max_degree_l,
        eccentricity_truncation=_ECC_TRUNCATION,
        obliquity_truncation=0)
    world.solve_eos(G_to_use=G)
    return world


def _mode_amplitudes(
        world,
        spin,
        radius,
        colatitude,
        longitude,
        max_degree_l):
    """Each active mode's (|frequency|, strain, stress, displacement) from the helpers, evaluated at +|frequency|."""
    degrees, frequencies, rows = tidal_potential_3d_modes(
        _R,
        _N,
        spin,
        _ECC,
        0.0,
        _SMA,
        _HOST,
        G,
        colatitude,
        longitude,
        min_degree_l=2,
        max_degree_l=max_degree_l,
        eccentricity_truncation=_ECC_TRUNCATION,
        obliquity_truncation=0)
    amplitudes = []
    for index in range(degrees.size):
        degree = int(degrees[index])
        frequency = float(frequencies[index])
        magnitude = abs(frequency)
        if magnitude <= tidalpy_constants.min_frequency:
            # The world drops static modes: they neither dissipate nor move.
            continue
        world.solve_love_numbers(frequency=magnitude, degree_l=degree)
        y = np.array([world.get_love_radial_y(radius, 0, y_index) for y_index in range(6)], dtype=np.complex128)
        # A negative-frequency mode is evaluated at +|frequency| with its conjugate row.
        row = rows[index] if frequency > 0.0 else np.conj(rows[index])
        strain, stress, _ = strain_stress_heating_point(
            y,
            complex(world.calc_complex_shear_modulus(radius, magnitude)),
            complex(world.calc_complex_bulk_modulus(radius, magnitude)),
            radius,
            float(degree),
            magnitude,
            True,
            False,
            row,
            colatitude)
        amplitudes.append((magnitude, strain, stress, displacement_point(y, row, colatitude)))
    return amplitudes


def test_strain_stress_keep_imaginary_potential():
    """The kernel is linear in a complex potential row: R(a + i b) = R(a) + i R(b) and R(i a) = i R(a)."""
    strain_real, stress_real, heating_real = _kernel(_ROW.real.copy())
    strain_imag, stress_imag, _ = _kernel(_ROW.imag.copy())
    strain_full, stress_full, _ = _kernel(_ROW)
    np.testing.assert_allclose(strain_full, strain_real + 1j * strain_imag, rtol=1.0e-13)
    # The isotropic stress is a difference of much larger terms, so compare on the largest component's scale.
    stress_scale = 1.0e-13 * float(np.max(np.abs(stress_full)))
    np.testing.assert_allclose(stress_full, stress_real + 1j * stress_imag, rtol=1.0e-13, atol=stress_scale)

    strain_rotated, stress_rotated, heating_rotated = _kernel(1j * _ROW.real)
    np.testing.assert_allclose(strain_rotated, 1j * strain_real, rtol=1.0e-13)
    np.testing.assert_allclose(stress_rotated, 1j * stress_real, rtol=1.0e-13, atol=stress_scale)
    assert math.isclose(heating_rotated, heating_real, rel_tol=1.0e-13)


def test_real_row_forms_agree():
    """A tuple of floats and a complex array with zero imaginary parts are the same potential row."""
    from_tuple = _kernel(tuple(_ROW.real))
    from_complex = _kernel(_ROW.real.astype(np.complex128))
    for tuple_value, complex_value in zip(from_tuple, from_complex):
        assert np.array_equal(tuple_value, complex_value)


def test_displacement_keeps_imaginary_potential():
    """The displacement closed form holds for a complex row whose dU/dphi is purely imaginary."""
    displacement = displacement_point(_Y, _ROW, _COLATITUDE)
    expected = np.array([_Y[0] * _ROW[0], _Y[2] * _ROW[1], _Y[2] * _ROW[2] / np.sin(_COLATITUDE)])
    np.testing.assert_allclose(displacement, expected, rtol=1.0e-14)
    assert displacement[2] != 0.0


@pytest.mark.parametrize("bad_row", [np.ones(5), np.ones(7), np.ones((1, 6)), ()])
def test_potential_row_must_hold_six_values(bad_row):
    """The kernels reject a potential row without exactly six values."""
    with pytest.raises(ValueError):
        _kernel(bad_row)
    with pytest.raises(ValueError):
        displacement_point(_Y, bad_row, _COLATITUDE)


@pytest.mark.parametrize("num_stress, num_strain", [(5, 6), (6, 7)])
def test_volumetric_heating_requires_six_components(num_stress, num_strain):
    with pytest.raises(ValueError):
        volumetric_heating(np.ones(num_stress, dtype=np.complex128), np.ones(num_strain, dtype=np.complex128), _N)


@pytest.mark.parametrize("spin_ratio, max_degree_l", _WORLD_CASES)
def test_helpers_reproduce_world_secular_heating(spin_ratio, max_degree_l):
    """Summing amplitudes per frequency and forming their heating reproduces the world's pointwise secular heating."""
    world = _build_world(max_degree_l)
    spin = spin_ratio * _N
    radius, colatitude, longitude = 0.8 * _R, _COLATITUDE, 0.7

    # Modes sharing a frequency interfere, so their amplitudes are summed before forming the heating.
    frequency_groups = []   # [|frequency|, summed strain, summed stress]
    for magnitude, strain, stress, _ in _mode_amplitudes(world, spin, radius, colatitude, longitude, max_degree_l):
        for group in frequency_groups:
            if math.isclose(group[0], magnitude, rel_tol=1.0e-9):
                group[1] = group[1] + strain
                group[2] = group[2] + stress
                break
        else:
            frequency_groups.append([magnitude, strain.copy(), stress.copy()])
    from_helpers = sum(
        volumetric_heating(stress, strain, magnitude) for magnitude, strain, stress in frequency_groups)

    expected = world.calc_3d_tides(
        _N,
        spin,
        _ECC,
        0.0,
        _SMA,
        _HOST,
        radii=np.array([radius]),
        colatitudes=np.array([colatitude]),
        longitudes=np.array([longitude]))["heating"][0, 0, 0]
    assert expected > 0.0
    assert math.isclose(from_helpers, expected, rel_tol=1.0e-9), f"helpers {from_helpers:.10e} != world {expected:.10e}"


@pytest.mark.parametrize("spin_ratio, max_degree_l", _WORLD_CASES)
def test_helpers_reproduce_world_displacements(spin_ratio, max_degree_l):
    """Summing Re[u e^{i |frequency| t}] over the modes reproduces the world's instantaneous displacements."""
    world = _build_world(max_degree_l)
    spin = spin_ratio * _N
    radius, colatitude, longitude = 0.8 * _R, _COLATITUDE, 0.7
    times = np.array([0.0, 1.3e4, 5.1e4])

    from_helpers = np.zeros((3, times.size))
    for magnitude, _, _, displacement in _mode_amplitudes(world, spin, radius, colatitude, longitude, max_degree_l):
        from_helpers += np.real(displacement[:, np.newaxis] * np.exp(1j * magnitude * times)[np.newaxis, :])

    expected = world.calc_3d_displacements(
        _N,
        spin,
        _ECC,
        0.0,
        _SMA,
        _HOST,
        radii=np.array([radius]),
        colatitudes=np.array([colatitude]),
        longitudes=np.array([longitude]),
        times=times)
    for component_index, component in enumerate(("radial", "polar", "azimuthal")):
        reference = expected[component][0, 0, 0, :]
        np.testing.assert_allclose(
            from_helpers[component_index],
            reference,
            rtol=1.0e-9,
            atol=1.0e-12 * np.max(np.abs(reference)))
