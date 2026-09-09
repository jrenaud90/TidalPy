"""Dense, SI-aware EOS readout (``RadialSolverSolution.eos_call_si``) for the standalone radial solver.

These lock in the fix for the 3D-tides blocker: a post-solve dense EOS query at an arbitrary SI radius
must route through the solver's own dense interpolant and return correct SI values. The earlier code left
the cysolver's interpolation inputs as ``c_radial_solver`` stack locals (they dangled after return) and
re-dimensionalized the input arrays, so the dense extra-output re-invoke produced doubly-scaled garbage
(density off by the bulk-density factor, shear off by the pascal factor). The inputs are now persisted in
the solution storage in non-dim solve units, so ``call`` redimensionalizes them back to SI exactly.

Dense output layout (SI): [0] gravity [1] pressure [2] mass [3] moi [4] density [5,6] shear re/im
[7,8] bulk re/im.
"""
import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.rheology.models import Maxwell, Elastic
from TidalPy.RadialSolver_x.solver import radial_solver as rs_x


def _build_homogeneous(nondimensionalize=True):
    n_slices = 20
    planet_radius = 1.0e6
    frequency = 2.0 * np.pi / 86400.0
    radius_array = np.linspace(0.0, planet_radius, n_slices)
    density_array = 5000.0 * np.ones(n_slices)
    shear_array = 5.0e10 * np.ones(n_slices)
    bulk_array = 1.0e11 * np.ones(n_slices)
    viscosity_array = 1.0e19 * np.ones(n_slices)

    complex_shear = np.empty(n_slices, np.complex128)
    Maxwell().vectorize_modulus_viscosity(frequency, shear_array, viscosity_array, complex_shear)
    complex_bulk = np.empty(n_slices, np.complex128)
    Elastic().vectorize_modulus_viscosity(frequency, bulk_array, viscosity_array, complex_bulk)

    shell_volume = (4.0 / 3.0) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
    bulk_density = float(np.sum(shell_volume * density_array[1:]) / np.sum(shell_volume))

    solution = rs_x(
        radius_array.copy(), density_array.copy(), complex_bulk.copy(), complex_shear.copy(),
        frequency, bulk_density,
        ('solid',), (False,), (False,), np.asarray((planet_radius,)),
        degree_l=2, solve_for=('tidal',), nondimensionalize=nondimensionalize,
        integration_method='DOP853', integration_rtol=1e-8, integration_atol=1e-10,
        max_num_steps=5_000_000, raise_on_fail=True)
    assert solution.success
    return solution, complex_shear[0], complex_bulk[0], 5000.0, planet_radius


@pytest.mark.parametrize('nondimensionalize', (True, False))
def test_eos_call_si_matches_supplied_moduli_at_grid_radii(nondimensionalize):
    """At each grid radius the dense getter reproduces the gridded SI moduli (and the supplied constants)."""
    solution, fed_shear, fed_bulk, fed_density, _ = _build_homogeneous(nondimensionalize)
    radius_array = np.asarray(solution.radius_array)
    shear_grid = np.asarray(solution.shear_modulus_array)
    bulk_grid = np.asarray(solution.bulk_modulus_array)

    # Skip the exact center (r = 0): the structure ODE zeros its derivatives there.
    for index in range(1, radius_array.size):
        eos = solution.eos_call_si(float(radius_array[index]))
        dense_shear = complex(eos[5], eos[6])
        dense_bulk = complex(eos[7], eos[8])
        # Dense readout reproduces the gridded modulus arrays (which are the redimensionalized solve output).
        assert dense_shear == pytest.approx(shear_grid[index], rel=1e-6)
        assert dense_bulk == pytest.approx(bulk_grid[index], rel=1e-6)
        # The body is homogeneous, so both equal the supplied constants - NOT the bug's *bulk-density garbage.
        assert dense_shear == pytest.approx(fed_shear, rel=1e-6)
        assert dense_bulk == pytest.approx(fed_bulk, rel=1e-6)
        assert eos[4] == pytest.approx(fed_density, rel=1e-6)   # density, was ~rho^2 before the fix


def test_eos_call_si_structure_is_physical():
    """Dense gravity/density at an off-grid radius match the homogeneous-sphere analytic values (SI)."""
    solution, _, _, density, planet_radius = _build_homogeneous(nondimensionalize=True)
    radius = 0.5454 * planet_radius   # deliberately between grid slices
    eos = solution.eos_call_si(radius)
    # Homogeneous sphere: g(r) = (4/3) pi G rho r.
    analytic_gravity = (4.0 / 3.0) * np.pi * G * density * radius
    assert eos[0] == pytest.approx(analytic_gravity, rel=1e-4)
    assert eos[4] == pytest.approx(density, rel=1e-6)
    # Off-grid moduli stay at the homogeneous constants (finite, correctly scaled).
    assert np.isfinite(eos[5]) and np.isfinite(eos[7])
    assert complex(eos[5], eos[6]).real == pytest.approx(5.0e10, rel=1e-3)


def test_eos_call_si_two_layers_distinct_moduli():
    """In a two-layer body the dense getter returns each layer's own moduli (correct layer location)."""
    n_per = 14
    planet_radius = 1.0e6
    interface = 0.5e6
    frequency = 2.0 * np.pi / 86400.0
    lower = np.linspace(0.0, interface, n_per)
    upper = np.linspace(interface, planet_radius, n_per)
    radius_array = np.concatenate([lower, upper])           # duplicated interface radius (per-layer convention)
    density_array = np.where(radius_array <= interface, 6000.0, 4000.0)

    shear_si = np.where(radius_array <= interface, 8.0e10, 3.0e10)
    bulk_si = np.where(radius_array <= interface, 2.0e11, 1.0e11)
    viscosity_array = 1.0e19 * np.ones_like(radius_array)
    complex_shear = np.empty(radius_array.size, np.complex128)
    Maxwell().vectorize_modulus_viscosity(frequency, shear_si, viscosity_array, complex_shear)
    complex_bulk = np.empty(radius_array.size, np.complex128)
    Elastic().vectorize_modulus_viscosity(frequency, bulk_si, viscosity_array, complex_bulk)

    shell_volume = (4.0 / 3.0) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
    weights = np.clip(shell_volume, 0.0, None)
    bulk_density = float(np.sum(weights * density_array[1:]) / np.sum(weights))

    solution = rs_x(
        radius_array.copy(), density_array.copy(), complex_bulk.copy(), complex_shear.copy(),
        frequency, bulk_density,
        ('solid', 'solid'), (False, False), (False, False), np.asarray((interface, planet_radius)),
        degree_l=2, solve_for=('tidal',), nondimensionalize=True,
        integration_method='DOP853', integration_rtol=1e-8, integration_atol=1e-10,
        max_num_steps=5_000_000, raise_on_fail=True)
    assert solution.success

    deep = solution.eos_call_si(0.25e6)     # well inside the lower layer
    shallow = solution.eos_call_si(0.75e6)  # well inside the upper layer
    assert complex(deep[5], deep[6]).real == pytest.approx(8.0e10, rel=1e-3)
    assert complex(shallow[5], shallow[6]).real == pytest.approx(3.0e10, rel=1e-3)
    assert deep[4] == pytest.approx(6000.0, rel=1e-3)
    assert shallow[4] == pytest.approx(4000.0, rel=1e-3)
