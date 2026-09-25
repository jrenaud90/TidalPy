"""Implicit (stiff) CyRK integrators through the radial solver.

CyRK provides the implicit methods BDF, LSODA, and Radau alongside the explicit Runge-Kutta family. These tests run a
homogeneous Maxwell sphere through `radial_solver` with each implicit method and require the tidal k2 to match the
DOP853 reference, for both the radial shooting integration and the equation-of-state solve.
"""
import numpy as np
import pytest

from TidalPy.rheology import Maxwell
from TidalPy.RadialSolver import radial_solver

IMPLICIT_METHODS = ("BDF", "LSODA", "Radau")

# Homogeneous sphere (Io-like scale) with a Maxwell viscoelastic shear response.
PLANET_RADIUS = 1.8e6
DENSITY = 3500.0
SHEAR_MODULUS = 6.0e10
VISCOSITY = 1.0e18
BULK_MODULUS = 2.0e11
FREQUENCY = 4.1e-5
NUM_SLICES = 100

radius_array = np.linspace(0.0, PLANET_RADIUS, NUM_SLICES)
density_array = np.full(NUM_SLICES, DENSITY)
bulk_modulus_array = np.full(NUM_SLICES, BULK_MODULUS, dtype=np.complex128)
complex_shear = complex(Maxwell()(FREQUENCY, SHEAR_MODULUS, VISCOSITY))
complex_shear_array = np.full(NUM_SLICES, complex_shear, dtype=np.complex128)
upper_radius_array = np.asarray([PLANET_RADIUS])

# The EOS tolerances are pinned loose: LSODA's startup at the singular center fails cleanly at tight
# tolerances (a documented limitation), and a homogeneous sphere's structure is exact at any tolerance.
COMMON_KWARGS = dict(
    degree_l=2,
    solve_for=("tidal",),
    use_kamata=True,
    integration_rtol=1.0e-7,
    integration_atol=1.0e-10,
    eos_rtol=1.0e-4,
    eos_atol=1.0e-6,
    raise_on_fail=True,
)


def _solve_k2(integration_method, eos_integration_method="DOP853"):
    solution = radial_solver(
        np.copy(radius_array),
        np.copy(density_array),
        np.copy(bulk_modulus_array),
        np.copy(complex_shear_array),
        FREQUENCY,
        DENSITY,
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
def test_implicit_radial_integration(integration_method):
    """Each implicit method reproduces the DOP853 k2."""
    k2_reference = _solve_k2("DOP853")
    k2_implicit = _solve_k2(integration_method)
    assert np.isclose(k2_implicit.real, k2_reference.real, rtol=1.0e-4)
    assert np.isclose(k2_implicit.imag, k2_reference.imag, rtol=1.0e-3, atol=1.0e-10)


@pytest.mark.parametrize("eos_method", IMPLICIT_METHODS)
def test_implicit_eos_integration(eos_method):
    """Each implicit method is accepted for the equation-of-state integration."""
    k2_reference = _solve_k2("DOP853", eos_integration_method="DOP853")
    k2_implicit = _solve_k2("DOP853", eos_integration_method=eos_method)
    assert np.isclose(k2_implicit.real, k2_reference.real, rtol=1.0e-4)
    assert np.isclose(k2_implicit.imag, k2_reference.imag, rtol=1.0e-3, atol=1.0e-10)


@pytest.mark.parametrize("method_name", ("rk45", "Dop853", "bdf", "RADAU"))
def test_method_names_are_case_insensitive(method_name):
    assert np.isfinite(_solve_k2(method_name).real)


@pytest.mark.parametrize("argument", ("integration_method", "eos_integration_method"))
def test_unknown_integration_method_raises(argument):
    """An unknown method name raises with a message listing the supported methods."""
    kwargs = {"integration_method": "DOP853", "eos_integration_method": "DOP853", argument: "not_a_method"}
    with pytest.raises(Exception, match="Supported: rk23, rk45, dop853, bdf, lsoda, radau"):
        _solve_k2(**kwargs)
