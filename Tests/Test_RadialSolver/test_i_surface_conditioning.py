"""Surface boundary condition conditioning diagnostic and warning.

A deep manual starting radius at high harmonic degree makes the shooting method's surface collapse constants grow
enormous and cancel: roundoff and integration error are amplified into the surface solution and the Love numbers
derived from it. The radial solver records the amplification factor on the solution
(``surface_solve_amplification``) and logs a warning when the roundoff floor exceeds the requested tolerance (or the
amplification is severe).
"""
import logging

import numpy as np
import pytest

from TidalPy.rheology import Maxwell
from TidalPy.RadialSolver import radial_solver
from TidalPy.RadialSolver.solver import SEVERE_SURFACE_AMPLIFICATION, check_surface_solve_conditioning

# Homogeneous solid planet.
frequency = 2.0 * np.pi / (86400. * 1.0)
N = 10
radius_array = np.linspace(0.0, 6000.e3, N)
bulk_density = 3500.
density_array = bulk_density * np.ones_like(radius_array)
bulk_modulus_array = 1.0e11 * np.ones(N, dtype=np.complex128, order='C')
complex_shear_modulus_array = np.full(N, complex(Maxwell()(frequency, 5.0e10, 1.0e20)), dtype=np.complex128)
upper_radius_by_layer = np.asarray((radius_array[-1],))

WARNING_TEXT = "poorly conditioned"


def _run(starting_radius, warnings=True):
    """Degree-3 dynamic incompressible solve; a 0.1 m starting radius makes it pathological."""
    return radial_solver(
        radius_array, density_array, bulk_modulus_array, complex_shear_modulus_array,
        frequency, bulk_density,
        ('solid',), (False,), (True,), upper_radius_by_layer,
        degree_l=3, solve_for=('tidal',), use_kamata=True,
        integration_method='RK45', integration_rtol=1.0e-7, integration_atol=1.0e-10,
        scale_rtols_bylayer_type=False,
        max_num_steps=5_000_000, expected_size=250, max_step=0,
        verbose=False, nondimensionalize=False, starting_radius=starting_radius,
        raise_on_fail=True, warnings=warnings,
    )


def _warned(caplog):
    return any(WARNING_TEXT in record.getMessage() for record in caplog.records)


def test_pathological_solve_records_amplification_and_warns(caplog):
    """A 0.1 m starting radius at degree 3 must report severe amplification and log the warning."""
    with caplog.at_level(logging.WARNING, logger='TidalPy'):
        solution = _run(starting_radius=0.1)
    assert solution.success
    assert solution.surface_solve_amplification > SEVERE_SURFACE_AMPLIFICATION
    assert _warned(caplog)


def test_healthy_solve_is_silent(caplog):
    """The automatic starting radius must stay below the warning thresholds and log nothing."""
    with caplog.at_level(logging.WARNING, logger='TidalPy'):
        solution = _run(starting_radius=0.0)
    assert solution.success
    assert 0.0 < solution.surface_solve_amplification < SEVERE_SURFACE_AMPLIFICATION
    assert not _warned(caplog)


def test_diagnostic_is_skipped_without_warnings(caplog):
    """With warnings disabled the diagnostic is not computed (it stays 0) and nothing is logged."""
    with caplog.at_level(logging.WARNING, logger='TidalPy'):
        solution = _run(starting_radius=0.1, warnings=False)
    assert solution.success
    assert solution.surface_solve_amplification == 0.0
    assert not _warned(caplog)


@pytest.mark.parametrize(("amplification", "rtol", "expected"), (
    (1.0e3, 1.0e-7, False),
    (1.0e10, 1.0e-7, True),    # Past the severe threshold.
    (1.0e7, 1.0e-12, True),    # Roundoff floor (1e7 * eps ~ 2e-9) above the tolerance.
    (1.0e7, 1.0e-6, False),
))
def test_check_surface_solve_conditioning(amplification, rtol, expected):
    assert check_surface_solve_conditioning(amplification, rtol) is expected
