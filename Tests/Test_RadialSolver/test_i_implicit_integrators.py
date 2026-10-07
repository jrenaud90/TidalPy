"""Implicit CyRK integrators (BDF, LSODA, Radau) for the radial and EOS integrations of a homogeneous Maxwell sphere."""
from pathlib import Path

import numpy as np
import pytest

from TidalPy.Rheology.rheology import Maxwell
from TidalPy.RadialSolver.solver import radial_solver

IMPLICIT_METHODS = ("BDF", "LSODA", "Radau")

# Frozen classic k2 keyed 'radial__<method>' (implicit radial, DOP853 EOS) and 'eos__<method>' (the reverse).
FROZEN_PATH = Path(__file__).parent / "frozen" / "test_i_implicit_integrators.npz"
with np.load(FROZEN_PATH, allow_pickle=False) as frozen_file:
    CLASSIC_K2 = {key: complex(frozen_file[key]) for key in frozen_file.files}

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

# Loose EOS tolerances: LSODA's startup at the singular center fails at tight ones, and a homogeneous structure is
# exact at any tolerance.
COMMON_KWARGS = dict(
    degree_l=2,
    solve_for=("tidal",),
    starting_method="kamata",
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


@pytest.mark.parametrize("reference", ("DOP853", "classic"))
@pytest.mark.parametrize("stage", ("radial", "eos"))
@pytest.mark.parametrize("method", IMPLICIT_METHODS)
def test_implicit_method_reproduces_k2(method, stage, reference):
    """An implicit radial or EOS integration reproduces the all-DOP853 k2 and the classic k2 for that method."""
    if stage == "radial":
        k2_implicit = _solve_k2(method)
    else:
        k2_implicit = _solve_k2("DOP853", eos_integration_method=method)
    k2_reference = _solve_k2("DOP853") if reference == "DOP853" else CLASSIC_K2[f"{stage}__{method}"]
    assert np.isclose(k2_implicit.real, k2_reference.real, rtol=1.0e-4)
    assert np.isclose(k2_implicit.imag, k2_reference.imag, rtol=1.0e-3, atol=1.0e-10)


def test_unknown_integration_method_raises():
    """An unknown method name raises with a message naming it unsupported."""
    with pytest.raises(Exception, match="[Uu]nsupported"):
        _solve_k2("not_a_method")
