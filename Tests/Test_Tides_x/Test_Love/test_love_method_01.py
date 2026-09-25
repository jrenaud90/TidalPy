"""Tests for the Love-number method names and the homogeneous-sphere Love-number functions
(TidalPy.Tides_x.love: love_method_name, calc_effective_rigidity, calc_homogeneous_love_numbers,
apply_fixed_q, apply_fixed_dt).
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.RadialSolver_x import homogeneous_love_numbers
from TidalPy.Tides_x.love import (
    LoveNumbers,
    apply_fixed_dt,
    apply_fixed_q,
    calc_effective_rigidity,
    calc_homogeneous_love_numbers,
    love_method_name,
)

RADIUS = 6.0e6
DENSITY = 4000.0
GRAVITY = (4.0 / 3.0) * math.pi * G * DENSITY * RADIUS
SHEAR = 6.0e10


@pytest.mark.parametrize("alias, canonical", (
    ("radial_solver", "radial_solver"), ("shooting", "radial_solver"), ("RS", "radial_solver"),
    ("propagation_matrix", "propagation_matrix"), ("prop_matrix", "propagation_matrix"), ("pm", "propagation_matrix"),
    ("prop", "propagation_matrix"),
    ("homogeneous", "homogeneous"), ("Homogen", "homogeneous"),
    ("cpl", "cpl"), ("CTL", "ctl"),
    ("laterally_inhomogeneous", "laterally_inhomogeneous"), ("3d", "laterally_inhomogeneous"),
    ("lat_inhom", "laterally_inhomogeneous"),
))
def test_love_method_name_aliases(alias, canonical):
    assert love_method_name(alias) == canonical


def test_love_method_name_unknown():
    with pytest.raises(ValueError, match="unknown Love-number method"):
        love_method_name("shoot")


def test_effective_rigidity_degree_two_prefactor():
    """At degree 2 the effective rigidity is (19/2) mu / (rho g R)."""
    expected = 9.5 * SHEAR / (DENSITY * GRAVITY * RADIUS)
    assert calc_effective_rigidity(SHEAR, DENSITY, GRAVITY, RADIUS) == pytest.approx(expected, rel=1e-14)
    assert calc_effective_rigidity(SHEAR, DENSITY, GRAVITY, RADIUS, 2) == pytest.approx(expected, rel=1e-14)


@pytest.mark.parametrize("degree_l", (2, 3, 4, 7))
def test_effective_rigidity_general_prefactor(degree_l):
    prefactor = (2.0 * degree_l**2 + 4.0 * degree_l + 3.0) / degree_l
    expected = prefactor * SHEAR / (DENSITY * GRAVITY * RADIUS)
    assert calc_effective_rigidity(SHEAR, DENSITY, GRAVITY, RADIUS, degree_l) == pytest.approx(expected, rel=1e-14)


def test_effective_rigidity_complex_input():
    mu = SHEAR * (1.0 + 0.1j)
    out = calc_effective_rigidity(mu, DENSITY, GRAVITY, RADIUS)
    assert isinstance(out, complex)
    assert out == pytest.approx(9.5 * mu / (DENSITY * GRAVITY * RADIUS), rel=1e-14)


@pytest.mark.parametrize("bad", (
    dict(degree_l=1), dict(density=0.0), dict(gravity=-1.0), dict(radius=float("nan")),
))
def test_effective_rigidity_invalid_inputs(bad):
    kwargs = dict(shear_modulus=SHEAR, density=DENSITY, gravity=GRAVITY, radius=RADIUS, degree_l=2)
    kwargs.update(bad)
    with pytest.raises(ValueError):
        calc_effective_rigidity(**kwargs)


def test_homogeneous_love_fluid_limit():
    """mu -> 0 gives the fluid Love numbers k2 = 3/2, h2 = 5/2, l2 = 3/4."""
    love = calc_homogeneous_love_numbers(0.0, DENSITY, GRAVITY, RADIUS)
    assert isinstance(love, LoveNumbers)
    assert love.k == pytest.approx(1.5)
    assert love.h == pytest.approx(2.5)
    assert love.l == pytest.approx(0.75)


def test_homogeneous_love_rigid_limit():
    love = calc_homogeneous_love_numbers(1.0e30, DENSITY, GRAVITY, RADIUS)
    assert abs(love.k) < 1e-15 and abs(love.h) < 1e-15 and abs(love.l) < 1e-15


@pytest.mark.parametrize("degree_l", (2, 3, 5))
def test_homogeneous_love_closed_form(degree_l):
    mu_eff = calc_effective_rigidity(SHEAR, DENSITY, GRAVITY, RADIUS, degree_l)
    love = calc_homogeneous_love_numbers(SHEAR, DENSITY, GRAVITY, RADIUS, degree_l)
    response = 1.0 / (1.0 + mu_eff)
    assert love.k == pytest.approx(3.0 / (2.0 * (degree_l - 1.0)) * response, rel=1e-14)
    assert love.h == pytest.approx((2.0 * degree_l + 1.0) / (2.0 * (degree_l - 1.0)) * response, rel=1e-14)
    assert love.l == pytest.approx(3.0 / (2.0 * degree_l * (degree_l - 1.0)) * response, rel=1e-14)
    assert love.h / love.k == pytest.approx((2.0 * degree_l + 1.0) / 3.0, rel=1e-14)


def test_homogeneous_love_complex_modulus():
    mu = SHEAR * (1.0 - 0.05j)
    love = calc_homogeneous_love_numbers(mu, DENSITY, GRAVITY, RADIUS)
    assert love.k.imag > 0.0       # a lossy modulus (negative imaginary part) gives a lagging k
    assert love.k == pytest.approx(1.5 / (1.0 + calc_effective_rigidity(mu, DENSITY, GRAVITY, RADIUS)), rel=1e-14)


@pytest.mark.parametrize("degree_l", (2, 3))
def test_homogeneous_love_matches_radial_solver(degree_l):
    """The closed form reproduces the propagation-matrix solve of a static incompressible uniform sphere."""
    solution = homogeneous_love_numbers(
        RADIUS, DENSITY, SHEAR + 0.0j, 1.0e-5, num_slices=60, degree_l=degree_l,
        layer_is_static=True, layer_is_incompressible=True, love_method='propagation_matrix')
    assert solution.success
    love = calc_homogeneous_love_numbers(SHEAR, DENSITY, GRAVITY, RADIUS, degree_l)
    np.testing.assert_allclose(np.atleast_1d(solution.k)[0], love.k, rtol=1e-9)
    np.testing.assert_allclose(np.atleast_1d(solution.h)[0], love.h, rtol=1e-9)
    np.testing.assert_allclose(np.atleast_1d(solution.l)[0], love.l, rtol=1e-9)


def test_apply_fixed_q():
    static = calc_homogeneous_love_numbers(SHEAR, DENSITY, GRAVITY, RADIUS)
    lagged = apply_fixed_q(static, 50.0)
    assert lagged.k == pytest.approx(static.k * (1.0 - 1.0j / 50.0))
    assert lagged.h == pytest.approx(static.h * (1.0 - 1.0j / 50.0))
    assert lagged.l == pytest.approx(static.l * (1.0 - 1.0j / 50.0))
    assert -lagged.k.imag == pytest.approx(static.k.real / 50.0)
    for bad in (0.0, -10.0, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="fixed_q"):
            apply_fixed_q(static, bad)


def test_apply_fixed_dt():
    static = calc_homogeneous_love_numbers(SHEAR, DENSITY, GRAVITY, RADIUS)
    omega, dt = 2.0e-5, 600.0
    lagged = apply_fixed_dt(static, omega, dt)
    assert lagged.k == pytest.approx(static.k * (1.0 - 1.0j * omega * dt))
    assert apply_fixed_dt(static, -omega, dt).k == lagged.k          # the frequency magnitude is used
    assert apply_fixed_dt(static, omega, 0.0) == static
    with pytest.raises(ValueError, match="fixed_dt"):
        apply_fixed_dt(static, omega, -1.0)
