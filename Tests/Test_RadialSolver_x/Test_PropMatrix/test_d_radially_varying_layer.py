"""The propagation matrix on a layer whose properties vary with radius.

Each slice is a shell of its own material, so the product of shell propagators must converge to the shooting
method's answer as the slices shrink. A propagator built from neighbouring slices' matrices instead cancels
every interior factor and returns the uniform body of the surface material, whatever the slice count.
"""
import numpy as np
import pytest

from TidalPy.RadialSolver_x.solver import radial_solver

PLANET_RADIUS = 6.0e6                     # [m]
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
    bulk_density = 3.0 * np.trapezoid(density * radius ** 2, radius) / PLANET_RADIUS ** 3
    return radius, np.ascontiguousarray(density), np.ascontiguousarray(shear + 0j), bulk_density


def _matrix_k2(case, num_slices):
    radius, density, shear, bulk_density = _profiles(case, num_slices)
    solution = radial_solver(
        radius, density, np.full(num_slices, 1.0e11 + 0j), shear, FORCING_FREQUENCY, bulk_density,
        ("solid",), (True,), (True,), np.array([PLANET_RADIUS]),
        love_method="propagation_matrix",
        warnings=False)
    return solution.k


def _shooting_k2(case, num_slices):
    # Dynamic and nearly incompressible, the closest shooting counterpart of the static incompressible matrix.
    radius, density, shear, bulk_density = _profiles(case, num_slices)
    solution = radial_solver(
        radius, density, np.full(num_slices, 1.0e17 + 0j), shear, FORCING_FREQUENCY, bulk_density,
        ("solid",), (True,), (False,), np.array([PLANET_RADIUS]),
        integration_rtol=1.0e-10,
        integration_atol=1.0e-14,
        warnings=False)
    return solution.k


@pytest.mark.parametrize("case", ("shear", "density"))
def test_matrix_converges_to_shooting_as_slices_shrink(case):
    reference = _shooting_k2(case, 800)
    errors = [abs(_matrix_k2(case, num_slices) - reference) / abs(reference) for num_slices in (50, 200, 800)]
    # First order in the slice width: each fourfold refinement cuts the error about fourfold.
    assert errors[0] > errors[1] > errors[2]
    assert errors[2] < 2.0e-3
    assert errors[0] / errors[2] > 8.0
