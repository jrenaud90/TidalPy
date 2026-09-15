"""Tests for the classic homogeneous-sphere Love-number helpers (TidalPy.tides.love1d).

`effective_rigidity_general` once multiplied mu / (rho g R) by 2 l^2 + 4 l + 3 / l instead of (2 l^2 + 4 l + 3) / l.
That gave 17.5 instead of 19/2 at degree 2, was wrong at every other degree, and made the classic quick tidal
dissipation return about 55 percent of the correct heating. These tests pin the formula at every degree, its agreement
with the degree-2 function and with the new backend, and the end-to-end heating of the quick path against the
closed-form leading-order result.
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.rheology_x import Maxwell
from TidalPy.tides.love1d import effective_rigidity, effective_rigidity_general
from TidalPy.Tides_x.love import calc_effective_rigidity, calc_homogeneous_love_numbers
from TidalPy.toolbox import quick_tidal_dissipation

SHEAR = 5.0e10
GRAVITY = 9.81
RADIUS = 6.371e6
DENSITY = 3500.0


def _prefactor(degree_l):
    return (2.0 * degree_l ** 2 + 4.0 * degree_l + 3.0) / degree_l


def test_general_matches_the_degree_two_function():
    """The two classic functions must agree where they overlap; the misplaced parenthesis broke exactly this."""
    general = effective_rigidity_general(SHEAR, GRAVITY, RADIUS, DENSITY, 2)
    assert general == pytest.approx(effective_rigidity(SHEAR, GRAVITY, RADIUS, DENSITY), rel=1e-15)
    assert general == pytest.approx(9.5 * SHEAR / (GRAVITY * RADIUS * DENSITY), rel=1e-15)


@pytest.mark.parametrize("degree_l", (2, 3, 4, 7, 10))
def test_general_prefactor(degree_l):
    expected = _prefactor(degree_l) * SHEAR / (GRAVITY * RADIUS * DENSITY)
    assert effective_rigidity_general(SHEAR, GRAVITY, RADIUS, DENSITY, degree_l) == pytest.approx(expected, rel=1e-15)


def test_general_vectorizes_over_shear_modulus():
    shear = np.linspace(1.0e9, 1.0e11, 5)
    result = effective_rigidity_general(shear, GRAVITY, RADIUS, DENSITY, 3)
    np.testing.assert_allclose(result, _prefactor(3) * shear / (GRAVITY * RADIUS * DENSITY), rtol=1e-15)


@pytest.mark.parametrize("degree_l", tuple(range(2, 11)))
def test_general_matches_the_new_backend(degree_l):
    """The new backend's calc_effective_rigidity takes its arguments as (mu, rho, g, R, l)."""
    classic = effective_rigidity_general(SHEAR, GRAVITY, RADIUS, DENSITY, degree_l)
    new = calc_effective_rigidity(SHEAR, DENSITY, GRAVITY, RADIUS, degree_l)
    assert classic == pytest.approx(new, rel=1e-14)


def test_quick_tidal_dissipation_matches_the_closed_form():
    """End to end: a homogeneous Maxwell body at degree 2, leading order in eccentricity.

    The closed form is (21/2) (-Im k2) G M_host^2 R^5 e^2 n / a^6, with k2 from the homogeneous-sphere formula and
    the semi-major axis from Kepler's law with both masses, which is how the quick path derives it. Before the fix
    the quick path returned 0.553 of this.
    """
    host_mass, radius, mass, gravity = 1.898e27, 1.82149e6, 8.9298e22, 1.796
    density = mass / ((4.0 / 3.0) * math.pi * radius ** 3)
    eccentricity, viscosity, shear = 0.0041, 1.0e19, 6.0e10
    orbital_frequency = math.sqrt(G * host_mass / 4.217e8 ** 3)
    semi_major_axis = (G * (host_mass + mass) / orbital_frequency ** 2) ** (1.0 / 3.0)

    result = quick_tidal_dissipation(
        host_mass,
        radius,
        mass,
        gravity,
        density,
        0.37685,
        viscosity=viscosity,
        shear_modulus=shear,
        rheology="Maxwell",
        eccentricity=eccentricity,
        obliquity=0.0,
        orbital_frequency=orbital_frequency,
        spin_frequency=orbital_frequency,
        max_tidal_order_l=2,
        eccentricity_truncation_lvl=2,
        use_obliquity=False)

    complex_shear = Maxwell().calc_complex_modulus(shear, viscosity, orbital_frequency)
    k2 = complex(calc_homogeneous_love_numbers(complex_shear, density, gravity, radius, 2).k)
    closed_form = ((21.0 / 2.0) * (-k2.imag) * G * host_mass ** 2 * radius ** 5 * eccentricity ** 2
                   * orbital_frequency / semi_major_axis ** 6)
    assert float(np.asarray(result["tidal_heating"])) == pytest.approx(closed_form, rel=1e-5)
