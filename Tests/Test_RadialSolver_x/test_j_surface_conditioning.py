"""Surface boundary condition conditioning diagnostic and warning.

A deep manual starting radius at high harmonic degree makes the shooting method's surface collapse
constants grow enormous and cancel: roundoff and integration error are amplified into the surface
solution and the Love numbers derived from it. Both the original and `_x` radial solvers record the
amplification factor on the solution (``surface_solve_amplification``) and warn when the roundoff floor
exceeds the requested tolerance (or the amplification is severe). The check itself is shared `_x` code, so
both solvers warn through the C++ (spdlog) logger, which these tests route to a temporary file; the classic
Python logger is captured too in case a classic solver ever warns through it.
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
from TidalPy.Utilities_x.logging_x.logger import flush_logger, init_logger
from TidalPy.initialize import build_logging_x_config

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


@pytest.fixture
def spdlog_text(tmp_path):
    """Route the C++ logger to a temporary file for the test and hand back a reader for its text."""
    log_path = tmp_path / "tidalpy_x.log"
    init_logger({"console_level": "off", "file_level": "warning", "log_to_file": True,
                 "log_file_path": str(log_path)})

    def read():
        flush_logger()
        return log_path.read_text(encoding="utf-8") if log_path.exists() else ""

    yield read
    init_logger(build_logging_x_config())


def _warned(caplog, spdlog_text):
    """Whether the conditioning warning was emitted through either logger."""
    classic = any(WARNING_TEXT in record.message for record in caplog.records)
    return classic or (WARNING_TEXT in spdlog_text())


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
def test_pathological_solve_records_amplification_and_warns(solver_func, caplog, spdlog_text):
    """A 0.1 m starting radius at degree 3 must report severe amplification and log the warning."""
    with caplog.at_level(logging.WARNING, logger='TidalPy'):
        solution = _run(solver_func, starting_radius=0.1)
    assert solution.success
    assert solution.surface_solve_amplification > SEVERE_SURFACE_AMPLIFICATION
    assert _warned(caplog, spdlog_text)


@pytest.mark.parametrize('solver_func', (radial_solver_old, radial_solver_new), ids=('old', 'new'))
def test_healthy_solve_is_silent(solver_func, caplog, spdlog_text):
    """The automatic starting radius must stay below the warning thresholds and log nothing."""
    with caplog.at_level(logging.WARNING, logger='TidalPy'):
        solution = _run(solver_func, starting_radius=0.0)
    assert solution.success
    assert 0.0 < solution.surface_solve_amplification < SEVERE_SURFACE_AMPLIFICATION
    assert not _warned(caplog, spdlog_text)


def test_world_level_amplification_property(spdlog_text):
    """The layered world records the diagnostic after a solve and stays silent when healthy."""
    world = build_world("earth_simple")
    world.solve_eos()
    result = world.solve_love_numbers(frequency_rad_s=1.0e-5)
    assert result['success']
    assert 0.0 < world.love_surface_amplification < SEVERE_SURFACE_AMPLIFICATION
    assert WARNING_TEXT not in spdlog_text()
