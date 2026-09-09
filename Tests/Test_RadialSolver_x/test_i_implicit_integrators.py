"""Implicit (stiff) CyRK integrators through both radial solvers.

CyRK provides the implicit methods BDF, LSODA, and Radau alongside the explicit
Runge-Kutta family. These tests run a homogeneous Maxwell sphere through the
RadialSolver_x and original RadialSolver entry points with each implicit method and
require the tidal k2 to match the DOP853 reference, for both the radial shooting
integration and the equation-of-state solve.
"""
import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.rheology_x.rheology import Maxwell

from TidalPy.RadialSolver.solver import radial_solver as radial_solver_old
from TidalPy.RadialSolver_x.solver import radial_solver as radial_solver_new

IMPLICIT_METHODS = ("BDF", "LSODA", "Radau")

# Homogeneous sphere (Io-like scale) with a Maxwell viscoelastic shear response.
PLANET_RADIUS_M = 1.8e6
DENSITY_KG_M3 = 3500.0
SHEAR_MODULUS_PA = 6.0e10
VISCOSITY_PAS = 1.0e18
BULK_MODULUS_PA = 2.0e11
FREQUENCY_RAD_S = 4.1e-5
NUM_SLICES = 100

radius_array = np.linspace(0.0, PLANET_RADIUS_M, NUM_SLICES)
density_array = np.full(NUM_SLICES, DENSITY_KG_M3)
bulk_modulus_array = np.full(NUM_SLICES, BULK_MODULUS_PA, dtype=np.complex128)
complex_shear = Maxwell().calc_complex_modulus(SHEAR_MODULUS_PA, VISCOSITY_PAS, FREQUENCY_RAD_S)
complex_shear_array = np.full(NUM_SLICES, complex_shear, dtype=np.complex128)
upper_radius_array = np.asarray([PLANET_RADIUS_M])

COMMON_KWARGS = dict(
    degree_l=2,
    solve_for=("tidal",),
    use_kamata=True,
    integration_rtol=1.0e-7,
    integration_atol=1.0e-10,
    raise_on_fail=True,
)


def _solve_k2(solver_func, integration_method, eos_integration_method="DOP853"):
    solution = solver_func(
        np.copy(radius_array),
        np.copy(density_array),
        np.copy(bulk_modulus_array),
        np.copy(complex_shear_array),
        FREQUENCY_RAD_S,
        DENSITY_KG_M3,
        ("solid",),
        (False,),
        (False,),
        np.copy(upper_radius_array),
        integration_method=integration_method,
        eos_integration_method=eos_integration_method,
        **COMMON_KWARGS,
    )
    assert solution.success, solution.message
    return complex(np.atleast_1d(solution.k)[0])


@pytest.mark.parametrize("integration_method", IMPLICIT_METHODS)
def test_implicit_radial_integration_new_solver(integration_method):
    """RadialSolver_x accepts each implicit method and reproduces the DOP853 k2."""
    k2_reference = _solve_k2(radial_solver_new, "DOP853")
    k2_implicit = _solve_k2(radial_solver_new, integration_method)
    assert np.isclose(k2_implicit.real, k2_reference.real, rtol=1.0e-4)
    assert np.isclose(k2_implicit.imag, k2_reference.imag, rtol=1.0e-3, atol=1.0e-10)


@pytest.mark.parametrize("integration_method", IMPLICIT_METHODS)
def test_implicit_radial_integration_old_solver(integration_method):
    """The original RadialSolver accepts each implicit method and reproduces the DOP853 k2."""
    k2_reference = _solve_k2(radial_solver_old, "DOP853")
    k2_implicit = _solve_k2(radial_solver_old, integration_method)
    assert np.isclose(k2_implicit.real, k2_reference.real, rtol=1.0e-4)
    assert np.isclose(k2_implicit.imag, k2_reference.imag, rtol=1.0e-3, atol=1.0e-10)


@pytest.mark.parametrize("eos_method", IMPLICIT_METHODS)
@pytest.mark.parametrize("solver_func", [radial_solver_new, radial_solver_old],
                         ids=["radial_solver_x", "radial_solver"])
def test_implicit_eos_integration(solver_func, eos_method):
    """Both solvers accept each implicit method for the equation-of-state integration."""
    k2_reference = _solve_k2(solver_func, "DOP853", eos_integration_method="DOP853")
    k2_implicit = _solve_k2(solver_func, "DOP853", eos_integration_method=eos_method)
    assert np.isclose(k2_implicit.real, k2_reference.real, rtol=1.0e-4)
    assert np.isclose(k2_implicit.imag, k2_reference.imag, rtol=1.0e-3, atol=1.0e-10)


@pytest.mark.parametrize("solver_func", [radial_solver_new, radial_solver_old],
                         ids=["radial_solver_x", "radial_solver"])
def test_unknown_integration_method_raises(solver_func):
    """An unknown method name raises with a message listing the supported methods."""
    with pytest.raises(Exception, match="[Uu]nsupported"):
        _solve_k2(solver_func, "not_a_method")
