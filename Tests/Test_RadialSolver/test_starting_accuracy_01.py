"""Starting conditions against high-precision (mpmath) evaluations where double precision used to lose digits.

Takeuchi's phi, phi_{l+1}, and psi (TS72 Eq. 103) at every argument size, and the f(k2) term of the Kamata and
Takeuchi solid solutions when the forcing frequency squared is far above gamma = (4/3) pi G rho.
"""
import cmath
import math

import numpy as np
import pytest

from TidalPy.RadialSolver.starting.common import takeuchi_phi_psi, z_calc
from TidalPy.RadialSolver.starting.kamata import kamata_solid_dynamic_compressible
from TidalPy.RadialSolver.starting.takeuchi import takeuchi_solid_dynamic_compressible

mpmath = pytest.importorskip("mpmath")

G_TO_USE = 6.67430e-11
DENSITY = 7000.0
BULK = 100.0e9 + 0.0j
SHEAR = 50.0e9 + 5.0e7j


def _phi_psi_exact(z2, degree_l):
    """(phi_l, phi_{l+1}, psi_l) from their TS72 definitions at 40 digits."""
    with mpmath.workdps(40):
        x = mpmath.sqrt(mpmath.mpc(z2))

        def phi(n):
            return mpmath.fac2(2 * n + 1) * mpmath.sqrt(mpmath.pi / (2 * x)) * mpmath.besselj(n + 0.5, x) / x ** n

        return (complex(phi(degree_l)), complex(phi(degree_l + 1)),
                complex(2 * (2 * degree_l + 3) * (1 - phi(degree_l)) / x ** 2))


@pytest.mark.parametrize("degree_l", (2, 3, 5, 10))
@pytest.mark.parametrize("magnitude", (1e-4, 0.1, 0.4, 0.6, 1.0, 3.0, 10.0, 30.0, 100.0, 1000.0))
def test_takeuchi_phi_psi_matches_definition(degree_l, magnitude):
    """The series (small |z2|) and the definitions (large |z2|) both hold 1e-12 on each side of the switch."""
    # Real positive, complex, real negative (an imaginary argument), and complex in the lower half plane.
    for phase in (0.0, 0.5, math.pi, -2.5):
        z2 = magnitude * cmath.exp(1j * phase)
        got = takeuchi_phi_psi(z2, degree_l)
        ref = _phi_psi_exact(z2, degree_l)
        for name, g, r in zip(("phi", "phi_lp1", "psi"), got, ref):
            assert abs(g - r) <= 1e-12 * abs(r), f"{name} at z2={z2}: {g} vs {r}"


def _f_and_k2_pos_exact(frequency, degree_l):
    """KMN15 / TS72 f(k2+) = (beta2 k2+ - w^2) / gamma and k2+ at 60 digits, where the cancellation is harmless."""
    with mpmath.workdps(60):
        w2 = mpmath.mpf(frequency) ** 2
        rho = mpmath.mpf(DENSITY)
        mu = mpmath.mpc(SHEAR)
        lame = mpmath.mpc(BULK) - mpmath.mpf(2) / 3 * mu
        gamma = 4 * mpmath.pi * mpmath.mpf(G_TO_USE) * rho / 3
        alpha2 = (lame + 2 * mu) / rho
        beta2 = mu / rho
        quad_pos = w2 / beta2 + (w2 + 4 * gamma) / alpha2
        quad_neg = w2 / beta2 - (w2 + 4 * gamma) / alpha2
        k2_pos = (quad_pos + mpmath.sqrt(quad_neg ** 2 + 4 * degree_l * (degree_l + 1) * gamma ** 2 /
                                         (alpha2 * beta2))) / 2
        return complex((beta2 * k2_pos - w2) / gamma), complex(k2_pos)


# w^2 / gamma is 5e3 at 0.1 rad/s and 5e5 at 1 rad/s; the two radii put k2+ r^2 on each side of the phi/psi switch.
@pytest.mark.parametrize("degree_l", (2, 3))
@pytest.mark.parametrize("frequency", (0.1, 1.0))
@pytest.mark.parametrize("radius", (1.0e3, 1.0e5))
def test_kamata_solid_high_frequency(degree_l, frequency, radius):
    """Kamata solution 1's y1 = -f(k2+) z(k2+ r^2) / r keeps full precision when w^2 >> gamma."""
    f_pos, k2_pos = _f_and_k2_pos_exact(frequency, degree_l)
    new = np.full((3, 6), np.nan, dtype=np.complex128)
    kamata_solid_dynamic_compressible(frequency, radius, DENSITY, BULK, SHEAR, degree_l, G_TO_USE, new)
    expected = -f_pos * z_calc(k2_pos * radius ** 2, degree_l) / radius
    assert abs(new[0, 0] - expected) <= 1e-12 * abs(expected), (new[0, 0], expected)


@pytest.mark.parametrize("degree_l", (2, 3))
@pytest.mark.parametrize("frequency", (0.1, 1.0))
@pytest.mark.parametrize("radius", (1.0e3, 1.0e5))
def test_takeuchi_solid_high_frequency(degree_l, frequency, radius):
    """Takeuchi's k2+ solution (row 1) keeps full precision in y1 and y3 when w^2 >> gamma, at any k2+ r^2."""
    f_pos, k2_pos = _f_and_k2_pos_exact(frequency, degree_l)
    _, phi_lp1, psi = _phi_psi_exact(k2_pos * radius ** 2, degree_l)
    h_pos = f_pos - (degree_l + 1.0)
    prefactor = -radius ** (degree_l + 1) / (2.0 * degree_l + 3.0)
    expected_y1 = prefactor * (0.5 * degree_l * h_pos * psi + f_pos * phi_lp1)
    expected_y3 = prefactor * (0.5 * h_pos * psi - phi_lp1)
    new = np.full((3, 6), np.nan, dtype=np.complex128)
    takeuchi_solid_dynamic_compressible(frequency, radius, DENSITY, BULK, SHEAR, degree_l, G_TO_USE, new)
    assert abs(new[1, 0] - expected_y1) <= 1e-12 * abs(expected_y1), (new[1, 0], expected_y1)
    assert abs(new[1, 2] - expected_y3) <= 1e-12 * abs(expected_y3), (new[1, 2], expected_y3)
