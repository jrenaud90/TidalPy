"""Surface boundary condition conditioning diagnostic and warning.

A deep manual starting radius at high harmonic degree makes the shooting method's surface collapse
constants grow enormous and cancel: roundoff and integration error are amplified into the surface
solution and the Love numbers derived from it. Both the original and `_x` radial solvers record the
amplification factor on the solution (``surface_solve_amplification``) and warn through the TidalPy
logger when the roundoff floor exceeds the requested tolerance (or the amplification is severe).
The layered world exposes the same diagnostic as ``love_surface_amplification``.
"""
import logging

import numpy as np
import pytest

from TidalPy.rheology.models import Maxwell
from TidalPy.RadialSolver import radial_solver as radial_solver_old
from TidalPy.RadialSolver_x.solver import radial_solver as radial_solver_new
from TidalPy.RadialSolver_x.rs_solution import SEVERE_SURFACE_AMPLIFICATION
from TidalPy.structures_x import build_world

# Homogeneous solid planet (matches the old-vs-new comparison test fixture).
frequency = 2.0 * np.pi / (86400. * 1.0)
N = 10
radius_array = np.linspace(0.0, 6000.e3, N)
bulk_density = 3500.
density_array = bulk_density * np.ones_like(radius_array)
bulk_modulus_array = 1.0e11 * np.ones(N, dtype=np.complex128, order='C')
viscosity_array = 1.0e20 * np.ones_like(radius_array)
shear_array = 5.0e10 * np.ones_like(radius_array)
_maxwell = Maxwell()
complex_shear_modulus_array = np.empty(N, dtype=np.complex128)
_maxwell.vectorize_modulus_viscosity(frequency, shear_array, viscosity_array, complex_shear_modulus_array)
upper_radius_by_layer = np.asarray((radius_array[-1],))

WARNING_TEXT = "poorly conditioned"


def _run(solver_func, starting_radius):
    """Degree-3 dynamic incompressible solve; a 0.1 m starting radius makes it pathological."""
    return solver_func(
        radius_array, density_array, bulk_modulus_array, complex_shear_modulus_array,
        frequency, bulk_density,
        ('solid',), (False,), (True,), upper_radius_by_layer,
        degree_l=3, solve_for=('tidal',), use_kamata=True,
        integration_method='RK45', integration_rtol=1.0e-7, integration_atol=1.0e-10,
        scale_rtols_bylayer_type=False,
        max_num_steps=5_000_000, expected_size=250, max_step=0,
        verbose=False, nondimensionalize=False, starting_radius=starting_radius,
        raise_on_fail=True, warnings=True,
    )


@pytest.mark.parametrize('solver_func', (radial_solver_old, radial_solver_new), ids=('old', 'new'))
def test_pathological_solve_records_amplification_and_warns(solver_func, caplog):
    """A 0.1 m starting radius at degree 3 must report severe amplification and log the warning."""
    with caplog.at_level(logging.WARNING, logger='TidalPy'):
        solution = _run(solver_func, starting_radius=0.1)
    assert solution.success
    assert solution.surface_solve_amplification > SEVERE_SURFACE_AMPLIFICATION
    assert any(WARNING_TEXT in record.message for record in caplog.records)


@pytest.mark.parametrize('solver_func', (radial_solver_old, radial_solver_new), ids=('old', 'new'))
def test_healthy_solve_is_silent(solver_func, caplog):
    """The automatic starting radius must stay below the warning thresholds and log nothing."""
    with caplog.at_level(logging.WARNING, logger='TidalPy'):
        solution = _run(solver_func, starting_radius=0.0)
    assert solution.success
    assert 0.0 < solution.surface_solve_amplification < SEVERE_SURFACE_AMPLIFICATION
    assert not any(WARNING_TEXT in record.message for record in caplog.records)


def test_world_level_amplification_property(caplog):
    """The layered world records the diagnostic after a solve and stays silent when healthy."""
    world = build_world("earth_simple")
    world.solve_eos()
    with caplog.at_level(logging.WARNING, logger='TidalPy'):
        result = world.solve_love_numbers(frequency_rad_s=1.0e-5)
    assert result['success']
    assert 0.0 < world.love_surface_amplification < SEVERE_SURFACE_AMPLIFICATION
    assert not any(WARNING_TEXT in record.message for record in caplog.records)
