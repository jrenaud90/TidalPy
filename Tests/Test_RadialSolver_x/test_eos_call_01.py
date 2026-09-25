"""Dense EOS readout (``RadialSolverSolution.eos_call``) for the standalone radial solver.

A post-solve dense EOS query at an arbitrary SI radius routes through the solver's own dense interpolant and
answers in SI with named fields: a float radius gives scalars, an array of radii gives arrays of the same shape.
The layout is frequency-independent, so its moduli are the unrelaxed ones; the ``complex_*`` fields come from
the same read ``get_complex_shear_modulus`` makes, which on this path reports the supplied complex arrays.
"""
import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.rheology_x import Maxwell, Elastic
from TidalPy.RadialSolver_x.rs_solution import EOS_CALL_FIELDS
from TidalPy.RadialSolver_x.solver import radial_solver as rs_x

COMPLEX_FIELDS = ("complex_shear_modulus", "complex_bulk_modulus")


def _build_homogeneous(nondimensionalize=True):
    n_slices = 20
    planet_radius = 1.0e6
    frequency = 2.0 * np.pi / 86400.0
    radius_array = np.linspace(0.0, planet_radius, n_slices)
    density_array = 5000.0 * np.ones(n_slices)
    shear_array = 5.0e10 * np.ones(n_slices)
    bulk_array = 1.0e11 * np.ones(n_slices)
    viscosity_array = 1.0e19 * np.ones(n_slices)

    complex_shear = Maxwell().calc_complex_modulus_vectorize_modulus(shear_array, viscosity_array, frequency)
    complex_bulk = Elastic().calc_complex_modulus_vectorize_modulus(bulk_array, viscosity_array, frequency)

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


def test_the_fields_are_the_layout_in_order():
    assert EOS_CALL_FIELDS == (
        "gravity", "pressure", "mass", "moi", "density", "shear_modulus", "bulk_modulus",
        "shear_viscosity", "bulk_viscosity", "temperature", "heat_flow", "melt_fraction")
    solution, _, _, _, planet_radius = _build_homogeneous()
    eos = solution.eos_call(0.5 * planet_radius)
    assert tuple(eos) == EOS_CALL_FIELDS + COMPLEX_FIELDS
    for name in EOS_CALL_FIELDS:
        assert type(eos[name]) is float, name
    for name in COMPLEX_FIELDS:
        assert type(eos[name]) is complex, name


@pytest.mark.parametrize('nondimensionalize', (True, False))
def test_eos_call_matches_the_radius_getters_and_the_supplied_moduli(nondimensionalize):
    """At each grid radius the dict reproduces the per-quantity getters and the supplied constants."""
    solution, fed_shear, fed_bulk, fed_density, _ = _build_homogeneous(nondimensionalize)
    radius_array = solution.sample_radii()

    # Skip the exact center (r = 0): the structure ODE zeros its derivatives there.
    for index in range(1, radius_array.size):
        radius = float(radius_array[index])
        eos = solution.eos_call(radius)
        assert eos["shear_modulus"] == pytest.approx(solution.get_shear_modulus(radius), rel=1e-6)
        assert eos["bulk_modulus"] == pytest.approx(solution.get_bulk_modulus(radius), rel=1e-6)
        assert eos["gravity"] == pytest.approx(solution.get_gravity(radius), rel=1e-12)
        assert eos["pressure"] == pytest.approx(solution.get_pressure(radius), rel=1e-12)
        # The body is homogeneous, so the moduli are the supplied constants.
        assert eos["complex_shear_modulus"] == pytest.approx(fed_shear, rel=1e-6)
        assert eos["complex_bulk_modulus"] == pytest.approx(fed_bulk, rel=1e-6)
        assert eos["complex_shear_modulus"] == solution.get_complex_shear_modulus(radius)
        assert eos["density"] == pytest.approx(fed_density, rel=1e-6)


def test_eos_call_structure_is_physical():
    """Dense gravity and density at an off-grid radius match the homogeneous-sphere analytic values."""
    solution, _, _, density, planet_radius = _build_homogeneous(nondimensionalize=True)
    radius = 0.5454 * planet_radius   # deliberately between grid slices
    eos = solution.eos_call(radius)
    # Homogeneous sphere: g(r) = (4/3) pi G rho r.
    analytic_gravity = (4.0 / 3.0) * np.pi * G * density * radius
    assert eos["gravity"] == pytest.approx(analytic_gravity, rel=1e-4)
    assert eos["density"] == pytest.approx(density, rel=1e-6)
    assert eos["mass"] == pytest.approx((4.0 / 3.0) * np.pi * density * radius**3, rel=1e-4)
    assert eos["shear_modulus"] == pytest.approx(5.0e10, rel=1e-3)
    assert np.isfinite(eos["bulk_modulus"])


def test_an_array_of_radii_gives_arrays_of_the_same_shape():
    solution, _, _, density, planet_radius = _build_homogeneous(nondimensionalize=True)
    radii = planet_radius * np.array([[0.1, 0.4], [0.7, 0.95]])
    profile = solution.eos_call(radii)
    for name in EOS_CALL_FIELDS + COMPLEX_FIELDS:
        assert profile[name].shape == radii.shape, name
    assert profile["density"].dtype == np.float64
    assert profile["complex_shear_modulus"].dtype == np.complex128
    # Element by element, the array answer is the scalar answer.
    for index in np.ndindex(radii.shape):
        single = solution.eos_call(float(radii[index]))
        for name in EOS_CALL_FIELDS + COMPLEX_FIELDS:
            assert np.array_equal(profile[name][index], single[name], equal_nan=True), name
    assert np.allclose(profile["density"], density)
    # A list works as an array does.
    from_list = solution.eos_call([0.3 * planet_radius, 0.6 * planet_radius])
    assert from_list["gravity"].shape == (2,)


def test_outside_the_body_is_nan():
    solution, _, _, _, planet_radius = _build_homogeneous(nondimensionalize=True)
    outside = solution.eos_call(1.5 * planet_radius)
    for name in EOS_CALL_FIELDS + COMPLEX_FIELDS:
        assert np.isnan(outside[name]), name
    mixed = solution.eos_call([0.5 * planet_radius, 2.0 * planet_radius])
    assert np.isfinite(mixed["density"][0]) and np.isnan(mixed["density"][1])
    assert np.isnan(mixed["complex_shear_modulus"][1])


def test_eos_call_is_independent_of_query_order():
    """The dense readout carries a running search seed between calls; a stale seed must not change it.

    The interpolated EOS pre-evaluation seeds each binary search with the slice the previous call used, so
    the value at a radius could in principle depend on which radii were queried before it. Ascending,
    descending, and shuffled sweeps over the same radii must agree exactly.
    """
    solution, _, _, _, planet_radius = _build_homogeneous(nondimensionalize=True)
    # Off-grid radii spanning the body, plus a couple of jumps to defeat a sequential seed.
    radii = [float(fraction * planet_radius)
             for fraction in (0.05, 0.17, 0.33, 0.49, 0.61, 0.78, 0.86, 0.97)]

    def row(radius):
        eos = solution.eos_call(radius)
        return np.array([eos[name] for name in EOS_CALL_FIELDS])

    ascending = [row(radius) for radius in radii]
    descending = [row(radius) for radius in reversed(radii)][::-1]
    shuffled_order = [3, 7, 0, 5, 1, 6, 2, 4]
    shuffled = [None] * len(radii)
    for index in shuffled_order:
        shuffled[index] = row(radii[index])

    # equal_nan: the viscosity entries are NaN for a supplied-moduli solve.
    for index, radius in enumerate(radii):
        assert np.array_equal(ascending[index], descending[index], equal_nan=True), radius
        assert np.array_equal(ascending[index], shuffled[index], equal_nan=True), radius


def test_eos_call_two_layers_distinct_moduli():
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
    complex_shear = Maxwell().calc_complex_modulus_vectorize_modulus(shear_si, viscosity_array, frequency)
    complex_bulk = Elastic().calc_complex_modulus_vectorize_modulus(bulk_si, viscosity_array, frequency)

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

    deep = solution.eos_call(0.25e6)     # well inside the lower layer
    shallow = solution.eos_call(0.75e6)  # well inside the upper layer
    assert deep["shear_modulus"] == pytest.approx(8.0e10, rel=1e-3)
    assert shallow["shear_modulus"] == pytest.approx(3.0e10, rel=1e-3)
    assert deep["density"] == pytest.approx(6000.0, rel=1e-3)
    assert shallow["density"] == pytest.approx(4000.0, rel=1e-3)
