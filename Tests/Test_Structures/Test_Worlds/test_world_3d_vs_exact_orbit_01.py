"""The 1D and 3D degree-2 heating against a direct calculation from the host's exact Kepler position.

The direct calculation forms the degree-2 potential coefficients C_m(t) = (4 pi / 5) G M R^2 / r^3 conj(Y_2m) of the
host in the body frame, Fourier-transforms them over one period of the body-frame motion, and sums
(5 R / (4 pi G)) |chi| (-Im k2(|chi|)) |C_m(chi)|^2 over every signed frequency chi. It keeps every same-frequency cross
term exactly, at a chosen argument of periapsis. The 3D total is the heating with the periapsis at the ascending node;
the 1D heating is its average over the argument of periapsis.
"""
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


_R = 1.0e6
_DENSITY = 5000.0
_MASS = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY
_N = 2.0 * np.pi / 86400.0
_HOST = mass_trap1
_NUM_SAMPLES = 2048
_NUM_PERIAPSES = 8


def _material():
    return constant_solid(
        _DENSITY, bulk_modulus=1.0e11, shear_modulus=5.0e10, shear_viscosity=1.0e19, bulk_viscosity=1.0e19)


def _build_world():
    world = BaseWorld("w", _R, _MASS)
    world.add_layer(Layer("mantle", 0, 0.0, _R, _MASS, _material(), is_static=False, shear_rheology=Maxwell(),
                          bulk_rheology=Elastic()))
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=2,
                          eccentricity_truncation=50, obliquity_truncation="gen")
    world.solve_eos(G_to_use=G)
    return world


def _degree_2_harmonic_magnitudes(colatitude, longitude):
    """Y_2m for m = 0, 1, 2; |Y_2,-m| = |Y_2m|, and only magnitudes of the Fourier amplitudes enter the heating."""
    cos_theta = np.cos(colatitude)
    sin_theta = np.sin(colatitude)
    return (
        math.sqrt(5.0 / (16.0 * math.pi)) * (3.0 * cos_theta ** 2 - 1.0) + 0.0j,
        math.sqrt(15.0 / (8.0 * math.pi)) * sin_theta * cos_theta * np.exp(1.0j * longitude),
        math.sqrt(15.0 / (32.0 * math.pi)) * sin_theta ** 2 * np.exp(2.0j * longitude))


def _direct_heating(world, spin_factor, eccentricity, obliquity, semi_major_axis, periapsis, love_cache):
    """Degree-2 heating [W] from the exact orbit at one argument of periapsis [rad]."""
    # The body-frame motion repeats after two orbits for the half-integer spin factors used here.
    period = 2.0 * (2.0 * math.pi / _N)
    times = np.arange(_NUM_SAMPLES) * period / _NUM_SAMPLES
    mean_anomaly = _N * times
    eccentric_anomaly = mean_anomaly + eccentricity * np.sin(mean_anomaly)
    for _ in range(50):
        eccentric_anomaly -= (eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly) - mean_anomaly) / \
                             (1.0 - eccentricity * np.cos(eccentric_anomaly))
    true_anomaly = 2.0 * np.arctan2(math.sqrt(1.0 + eccentricity) * np.sin(eccentric_anomaly / 2.0),
                                    math.sqrt(1.0 - eccentricity) * np.cos(eccentric_anomaly / 2.0))
    distance = semi_major_axis * (1.0 - eccentricity * np.cos(eccentric_anomaly))

    # Host direction in the body frame: the orbit is inclined by the obliquity with its node on the x axis at t = 0.
    argument_of_latitude = periapsis + true_anomaly
    colatitude = np.arccos(np.clip(np.sin(argument_of_latitude) * math.sin(obliquity), -1.0, 1.0))
    longitude = np.arctan2(np.sin(argument_of_latitude) * math.cos(obliquity), np.cos(argument_of_latitude)) - \
        spin_factor * _N * times

    frequencies = 2.0 * math.pi * np.fft.fftfreq(_NUM_SAMPLES, d=period / _NUM_SAMPLES)
    radial_factor = (4.0 * math.pi / 5.0) * G * _HOST * _R ** 2 / distance ** 3
    total = 0.0
    for order_m, harmonic in enumerate(_degree_2_harmonic_magnitudes(colatitude, longitude)):
        power = np.abs(np.fft.fft(radial_factor * np.conj(harmonic)) / _NUM_SAMPLES) ** 2
        significant = power > 1.0e-20 * power.max()
        order_total = 0.0
        for frequency, amplitude_squared in zip(frequencies[significant], power[significant]):
            magnitude = abs(frequency)
            if magnitude < 1.0e-12 * _N:
                continue
            key = round(magnitude / _N, 8)
            if key not in love_cache:
                world.solve_love_numbers(frequency=magnitude, degree_l=2)
                love_cache[key] = -world.love_number_k.imag
            order_total += magnitude * love_cache[key] * amplitude_squared
        # m and -m carry the same power.
        total += order_total if order_m == 0 else 2.0 * order_total
    return (5.0 * _R / (4.0 * math.pi * G)) * total


@pytest.mark.parametrize("spin_factor, eccentricity, obliquity", [
    (1.0, 0.2, 0.3),
    (1.5, 0.2, 0.3),
    (2.5, 0.05, 0.1),
    (1.5, 0.2, 0.0),
])
def test_1d_is_periapsis_average_and_3d_is_fixed_periapsis(spin_factor, eccentricity, obliquity):
    """The 1D heating is the exact-orbit heating averaged over the periapsis; the 3D total is its value at zero."""
    world = _build_world()
    semi_major_axis = orbital_motion2semi_a(_N, _HOST, _MASS)
    state = (_N, spin_factor * _N, eccentricity, obliquity, semi_major_axis, _HOST)
    world.calc_tides(*state)
    heating_1d = world.get_tidal_heating()
    heating_3d = world.calc_3d_tides(
        *state,
        latitude_summed=True,
        longitude_summed=True,
        radial_summed=True)['total']

    love_cache = {}
    direct_by_periapsis = [
        _direct_heating(world, spin_factor, eccentricity, obliquity, semi_major_axis, periapsis, love_cache)
        for periapsis in np.arange(_NUM_PERIAPSES) * math.pi / _NUM_PERIAPSES]

    # The cross terms vary with cos(2 omega) and cos(4 omega), so eight periapses over pi average them out exactly.
    assert math.isclose(heating_1d, np.mean(direct_by_periapsis), rel_tol=1.0e-10), \
        (heating_1d, np.mean(direct_by_periapsis))
    # Tolerance is the radial quadrature of the 3D total (1e-9 here with the default 16 nodes).
    assert math.isclose(heating_3d, direct_by_periapsis[0], rel_tol=1.0e-7), (heating_3d, direct_by_periapsis[0])
