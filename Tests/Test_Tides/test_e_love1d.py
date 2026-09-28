"""Homogeneous-sphere Love helpers against TidalPy 0.7's frozen ``tides.love1d`` results.

Pins the fix to 0.7's misplaced parenthesis, which used 2 l^2 + 4 l + 3 / l as the effective-rigidity prefactor.
"""
import math
from pathlib import Path

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Rheology import Maxwell
from TidalPy.Tides.love import calc_effective_rigidity, calc_homogeneous_love_numbers

SHEAR = 5.0e10
GRAVITY = 9.81
RADIUS = 6.371e6
DENSITY = 3500.0

# Classic backend results (corrected formula), generated from the 0.8.0 development snapshot 8b8e0b12.
FROZEN_PATH = Path(__file__).parent / "frozen" / "test_e_love1d.npz"
with np.load(FROZEN_PATH, allow_pickle=False) as frozen_file:
    CLASSIC_REFERENCE = {key: frozen_file[key] for key in frozen_file.files}


def _prefactor(degree_l):
    return (2.0 * degree_l ** 2 + 4.0 * degree_l + 3.0) / degree_l


def test_general_matches_the_degree_two_function():
    """Degree 2 matches the classic degree-2 function and 19/2 mu / (rho g R)."""
    new = calc_effective_rigidity(SHEAR, DENSITY, GRAVITY, RADIUS, 2)
    assert new == pytest.approx(float(CLASSIC_REFERENCE["effective_rigidity_degree_2"]), rel=1e-15)
    assert new == pytest.approx(9.5 * SHEAR / (GRAVITY * RADIUS * DENSITY), rel=1e-15)


@pytest.mark.parametrize("degree_l", tuple(range(2, 11)))
def test_general_prefactor(degree_l):
    """calc_effective_rigidity(mu, rho, g, R, l) matches the formula and the classic general function."""
    expected = _prefactor(degree_l) * SHEAR / (GRAVITY * RADIUS * DENSITY)
    classic = float(CLASSIC_REFERENCE[f"effective_rigidity_general__degree_l_{degree_l}"])
    new = calc_effective_rigidity(SHEAR, DENSITY, GRAVITY, RADIUS, degree_l)
    assert new == pytest.approx(expected, rel=1e-15)
    assert new == pytest.approx(classic, rel=1e-15)


def test_general_vectorizes_over_shear_modulus():
    """Per-element evaluation matches the classic function's array of shear moduli."""
    shear = np.linspace(1.0e9, 1.0e11, 5)
    classic = CLASSIC_REFERENCE["effective_rigidity_general_vector__degree_l_3"]
    new = np.asarray([calc_effective_rigidity(shear_modulus, DENSITY, GRAVITY, RADIUS, 3) for shear_modulus in shear])
    np.testing.assert_allclose(new, classic, rtol=1e-15)
    np.testing.assert_allclose(new, _prefactor(3) * shear / (GRAVITY * RADIUS * DENSITY), rtol=1e-15)


def test_quick_tidal_dissipation_matches_the_closed_form():
    """Classic quick heating of a homogeneous Maxwell body matches (21/2) (-Im k2) G M_host^2 R^5 e^2 n / a^6."""
    host_mass, radius, mass, gravity = 1.898e27, 1.82149e6, 8.9298e22, 1.796
    density = mass / ((4.0 / 3.0) * math.pi * radius ** 3)
    eccentricity, viscosity, shear = 0.0041, 1.0e19, 6.0e10
    orbital_frequency = math.sqrt(G * host_mass / 4.217e8 ** 3)
    # The quick path derives the semi-major axis from Kepler's law with both masses.
    semi_major_axis = (G * (host_mass + mass) / orbital_frequency ** 2) ** (1.0 / 3.0)

    classic_heating = float(CLASSIC_REFERENCE["quick_tidal_dissipation_heating"])

    complex_shear = Maxwell().calc_complex_modulus(shear, viscosity, orbital_frequency)
    k2 = complex(calc_homogeneous_love_numbers(complex_shear, density, gravity, radius, 2).k)
    closed_form = ((21.0 / 2.0) * (-k2.imag) * G * host_mass ** 2 * radius ** 5 * eccentricity ** 2
                   * orbital_frequency / semi_major_axis ** 6)
    assert classic_heating == pytest.approx(closed_form, rel=1e-5)
