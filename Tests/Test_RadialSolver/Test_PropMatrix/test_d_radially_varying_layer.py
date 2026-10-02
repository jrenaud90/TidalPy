"""The propagation matrix on a radially varying layer converges to the shooting method as slices shrink."""
import numpy as np
import pytest

from TidalPy.RadialSolver.solver import radial_solver
from numpy_compat import trapezoid

PLANET_RADIUS = 6.0e6                      # [m]
FORCING_FREQUENCY = 2.0 * np.pi / 86400.0  # [rad s-1]


def _profiles(case, num_slices):
    radius = np.linspace(0.0, PLANET_RADIUS, num_slices)
    fraction = radius / PLANET_RADIUS
    if case == "shear":
        shear = 2.0e10 + 8.0e10 * fraction
        density = np.full(num_slices, 5000.0)
    else:
        shear = np.full(num_slices, 5.0e10)
        density = 7000.0 - 4000.0 * fraction
    bulk_density = 3.0 * trapezoid(density * radius ** 2, radius) / PLANET_RADIUS ** 3
    return radius, np.ascontiguousarray(density), np.ascontiguousarray(shear + 0j), bulk_density


def _k2(
        case,
        num_slices,
        bulk_modulus,
        is_incompressible,
        **kwargs):
    radius, density, shear, bulk_density = _profiles(case, num_slices)
    solution = radial_solver(
        radius,
        density,
        np.full(num_slices, bulk_modulus + 0j),
        shear,
        FORCING_FREQUENCY,
        bulk_density,
        ("solid",),
        (True,),
        (is_incompressible,),
        np.array([PLANET_RADIUS]),
        warnings=False,
        **kwargs)
    return solution.k


def _matrix_k2(case, num_slices):
    return _k2(
        case,
        num_slices,
        1.0e11,
        True,
        love_method="propagation_matrix")


def _shooting_k2(case, num_slices):
    # Nearly incompressible shooting is the closest counterpart of the incompressible matrix method.
    return _k2(
        case,
        num_slices,
        1.0e17,
        False,
        integration_rtol=1.0e-10,
        integration_atol=1.0e-14)


@pytest.mark.parametrize("case", ("shear", "density"))
def test_matrix_converges_to_shooting_as_slices_shrink(case):
    """The matrix k2 converges to shooting first order in slice width (a neighbor-slice propagator would not)."""
    reference = _shooting_k2(case, 800)
    errors = [abs(_matrix_k2(case, num_slices) - reference) / abs(reference) for num_slices in (50, 200, 800)]
    assert errors[0] > errors[1] > errors[2]
    assert errors[2] < 2.0e-3
    assert errors[0] / errors[2] > 8.0
