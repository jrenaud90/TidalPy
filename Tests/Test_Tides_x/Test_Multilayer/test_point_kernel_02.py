"""Point-wise 3D kernels with complex potential rows (``TidalPy.Tides_x.multilayer.stress_strain``).

``tidal_potential_3d_modes`` returns complex phasor rows, so the point helpers must keep their imaginary parts. These
tests pin that, check the input validation, and rebuild the world's pointwise secular heating density and its
instantaneous displacements from the helpers alone, following the assembly rules in the module docstring: evaluate at
``|frequency|``, conjugate the rows of negative-frequency modes, and sum the amplitudes that share a frequency before
forming the heating.
"""
import math

import numpy as np
import pytest

import TidalPy.constants as tidalpy_constants
from TidalPy.constants import G, mass_trap1
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
from TidalPy.rheology_x.rheology import Elastic, Maxwell
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.Tides_x.classes.tide import make_tide
from TidalPy.Tides_x.multilayer.stress_strain import displacement_point, strain_stress_heating_point, volumetric_heating
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

# Generic inputs for the kernel-only checks: every potential entry has a distinct real and imaginary part.
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
        True,
        False,
        row,
        _COLATITUDE)


def _build_world(max_degree_l):
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
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(
        min_degree_l=2,
        max_degree_l=max_degree_l,
        eccentricity_truncation=_ECC_TRUNCATION,
        obliquity_truncation=0)
    world.solve_eos(G_to_use=G)
    return world


def _mode_amplitudes(world, spin, radius, colatitude, longitude, max_degree_l):
    """Each active mode's (|frequency|, strain, stress, displacement) amplitudes at +|frequency| from the helpers."""
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
        if magnitude <= tidalpy_constants.min_spin_orbit_diff:
            # The world drops static modes: they neither dissipate nor move.
            continue
        world.solve_love_numbers(frequency=magnitude, degree_l=degree)
        y = np.array([world.get_love_radial_y(radius, 0, y_index) for y_index in range(6)], dtype=np.complex128)
        row = rows[index] if frequency > 0.0 else np.conj(rows[index])
        strain, stress, _ = strain_stress_heating_point(
            y,
            complex(world.calc_complex_shear_modulus(radius, magnitude)),
            complex(world.calc_complex_bulk_modulus(radius, magnitude)),
            radius,
            float(degree),
            True,
            False,
            row,
            colatitude)
        amplitudes.append((magnitude, strain, stress, displacement_point(y, row, colatitude)))
    return amplitudes


# =====================================================================================================================
# Complex potential rows
# =====================================================================================================================
def test_strain_stress_keep_imaginary_potential():
    """The kernel is linear in the potential row, so the full complex row gives the response to its real part plus i
    times the response to its imaginary part, and i times a real row gives i times the response."""
    strain_real, stress_real, heating_real = _kernel(_ROW.real.copy())
    strain_imag, stress_imag, _ = _kernel(_ROW.imag.copy())
    strain_full, stress_full, _ = _kernel(_ROW)
    np.testing.assert_allclose(strain_full, strain_real + 1j * strain_imag, rtol=1.0e-13)
    np.testing.assert_allclose(stress_full, stress_real + 1j * stress_imag, rtol=1.0e-13)

    strain_rotated, stress_rotated, heating_rotated = _kernel(1j * _ROW.real)
    np.testing.assert_allclose(strain_rotated, 1j * strain_real, rtol=1.0e-13)
    np.testing.assert_allclose(stress_rotated, 1j * stress_real, rtol=1.0e-13)
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
    with pytest.raises(ValueError):
        _kernel(bad_row)
    with pytest.raises(ValueError):
        displacement_point(_Y, bad_row, _COLATITUDE)


def test_volumetric_heating_requires_six_components():
    with pytest.raises(ValueError):
        volumetric_heating(np.ones(5, dtype=np.complex128), np.ones(6, dtype=np.complex128))
    with pytest.raises(ValueError):
        volumetric_heating(np.ones(6, dtype=np.complex128), np.ones(7, dtype=np.complex128))


# =====================================================================================================================
# Assembling modes reproduces the world
# =====================================================================================================================
@pytest.mark.parametrize("spin_ratio, max_degree_l", _WORLD_CASES)
def test_helpers_reproduce_world_secular_heating(spin_ratio, max_degree_l):
    """Summing each frequency's amplitudes and applying |frequency| / 2 times the heating form reproduces the world's
    pointwise secular density, including the longitude dependence of the synchronous case."""
    world = _build_world(max_degree_l)
    spin = spin_ratio * _N
    radius, colatitude, longitude = 0.8 * _R, _COLATITUDE, 0.7

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
        0.5 * magnitude * volumetric_heating(stress, strain) for magnitude, strain, stress in frequency_groups)

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
