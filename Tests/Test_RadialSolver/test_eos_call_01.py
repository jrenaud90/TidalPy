"""The dense EOS readout `RadialSolverSolution.eos_call` for the standalone radial solver."""
import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Rheology import Maxwell, Elastic
from TidalPy.RadialSolver.rs_solution import EOS_CALL_FIELDS
from TidalPy.RadialSolver.solver import radial_solver

COMPLEX_FIELDS = ("complex_shear_modulus", "complex_bulk_modulus")
ALL_FIELDS = EOS_CALL_FIELDS + COMPLEX_FIELDS
PLANET_RADIUS = 1.0e6
FREQUENCY = 2.0 * np.pi / 86400.0
DENSITY = 5000.0


def _solve(
        radius_array,
        density_array,
        shear_si,
        bulk_si,
        upper_radii,
        nondimensionalize=True):
    """Solve a solid body with Maxwell shear and elastic bulk moduli; returns the solution and the complex moduli."""
    viscosity_array = 1.0e19 * np.ones_like(radius_array)
    complex_shear = Maxwell().calc_complex_modulus_vectorize_modulus(shear_si, viscosity_array, FREQUENCY)
    complex_bulk = Elastic().calc_complex_modulus_vectorize_modulus(bulk_si, viscosity_array, FREQUENCY)

    shell_volume = (4.0 / 3.0) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
    bulk_density = float(np.sum(shell_volume * density_array[1:]) / np.sum(shell_volume))

    num_layers = len(upper_radii)
    solution = radial_solver(
        radius_array.copy(),
        density_array.copy(),
        complex_bulk.copy(),
        complex_shear.copy(),
        FREQUENCY,
        bulk_density,
        ('solid',) * num_layers,
        (False,) * num_layers,
        (False,) * num_layers,
        np.asarray(upper_radii),
        degree_l=2,
        solve_for=('tidal',),
        nondimensionalize=nondimensionalize,
        integration_method='DOP853',
        integration_rtol=1e-8,
        integration_atol=1e-10,
        max_num_steps=5_000_000,
        raise_on_fail=True)
    assert solution.success
    return solution, complex_shear, complex_bulk


def _build_homogeneous(nondimensionalize=True):
    """A homogeneous 20-slice solid; returns the solution and its supplied complex shear and bulk moduli."""
    num_slices = 20
    solution, complex_shear, complex_bulk = _solve(
        np.linspace(0.0, PLANET_RADIUS, num_slices),
        DENSITY * np.ones(num_slices),
        5.0e10 * np.ones(num_slices),
        1.0e11 * np.ones(num_slices),
        (PLANET_RADIUS,),
        nondimensionalize)
    return solution, complex_shear[0], complex_bulk[0]


def test_the_fields_are_the_layout_in_order():
    """The field names follow the EOS layout, then the complex moduli, with float and complex scalar values."""
    assert EOS_CALL_FIELDS == (
        "gravity", "pressure", "mass", "moi", "density", "shear_modulus", "bulk_modulus",
        "shear_viscosity", "bulk_viscosity", "temperature", "heat_flow", "melt_fraction")
    solution, _, _ = _build_homogeneous()
    eos = solution.eos_call(0.5 * PLANET_RADIUS)
    assert tuple(eos) == ALL_FIELDS
    for name in EOS_CALL_FIELDS:
        assert type(eos[name]) is float, name
    for name in COMPLEX_FIELDS:
        assert type(eos[name]) is complex, name


@pytest.mark.parametrize('nondimensionalize', (True, False))
def test_eos_call_matches_the_radius_getters_and_the_supplied_moduli(nondimensionalize):
    """At each grid radius the dict reproduces the per-quantity getters and the supplied constants."""
    solution, fed_shear, fed_bulk = _build_homogeneous(nondimensionalize)
    radius_array = solution.sample_radii()

    # Skip r = 0, where the structure ODE zeros its derivatives.
    for index in range(1, radius_array.size):
        radius = float(radius_array[index])
        eos = solution.eos_call(radius)
        assert eos["shear_modulus"] == pytest.approx(solution.get_shear_modulus(radius), rel=1e-6)
        assert eos["bulk_modulus"] == pytest.approx(solution.get_bulk_modulus(radius), rel=1e-6)
        assert eos["gravity"] == pytest.approx(solution.get_gravity(radius), rel=1e-12)
        assert eos["pressure"] == pytest.approx(solution.get_pressure(radius), rel=1e-12)
        assert eos["complex_shear_modulus"] == pytest.approx(fed_shear, rel=1e-6)
        assert eos["complex_bulk_modulus"] == pytest.approx(fed_bulk, rel=1e-6)
        assert eos["complex_shear_modulus"] == solution.get_complex_shear_modulus(radius)
        assert eos["density"] == pytest.approx(DENSITY, rel=1e-6)


def test_eos_call_structure_is_physical():
    """Gravity, density, and mass at an off-grid radius match the homogeneous-sphere values."""
    solution, _, _ = _build_homogeneous()
    radius = 0.5454 * PLANET_RADIUS
    eos = solution.eos_call(radius)
    assert eos["gravity"] == pytest.approx((4.0 / 3.0) * np.pi * G * DENSITY * radius, rel=1e-4)
    assert eos["density"] == pytest.approx(DENSITY, rel=1e-6)
    assert eos["mass"] == pytest.approx((4.0 / 3.0) * np.pi * DENSITY * radius**3, rel=1e-4)
    assert eos["shear_modulus"] == pytest.approx(5.0e10, rel=1e-3)
    assert np.isfinite(eos["bulk_modulus"])


def test_an_array_of_radii_gives_arrays_of_the_same_shape():
    """An array (or list) of radii returns same-shape arrays that equal the scalar answers element by element."""
    solution, _, _ = _build_homogeneous()
    radii = PLANET_RADIUS * np.array([[0.1, 0.4], [0.7, 0.95]])
    profile = solution.eos_call(radii)
    for name in ALL_FIELDS:
        assert profile[name].shape == radii.shape, name
    assert profile["density"].dtype == np.float64
    assert profile["complex_shear_modulus"].dtype == np.complex128
    for index in np.ndindex(radii.shape):
        single = solution.eos_call(float(radii[index]))
        for name in ALL_FIELDS:
            assert np.array_equal(profile[name][index], single[name], equal_nan=True), name
    assert np.allclose(profile["density"], DENSITY)
    from_list = solution.eos_call([0.3 * PLANET_RADIUS, 0.6 * PLANET_RADIUS])
    assert from_list["gravity"].shape == (2,)


def test_outside_the_body_is_nan():
    """Every field is NaN above the surface, for scalar and mixed array queries."""
    solution, _, _ = _build_homogeneous()
    outside = solution.eos_call(1.5 * PLANET_RADIUS)
    for name in ALL_FIELDS:
        assert np.isnan(outside[name]), name
    mixed = solution.eos_call([0.5 * PLANET_RADIUS, 2.0 * PLANET_RADIUS])
    assert np.isfinite(mixed["density"][0]) and np.isnan(mixed["density"][1])
    assert np.isnan(mixed["complex_shear_modulus"][1])


def test_eos_call_is_independent_of_query_order():
    """Ascending, descending, and shuffled sweeps agree exactly despite the readout's running search seed."""
    solution, _, _ = _build_homogeneous()
    # Off-grid radii with jumps between them to defeat a sequential seed.
    radii = [float(fraction * PLANET_RADIUS) for fraction in (0.05, 0.17, 0.33, 0.49, 0.61, 0.78, 0.86, 0.97)]

    def row(radius):
        eos = solution.eos_call(radius)
        return np.array([eos[name] for name in EOS_CALL_FIELDS])

    ascending = [row(radius) for radius in radii]
    descending = [row(radius) for radius in reversed(radii)][::-1]
    shuffled = [None] * len(radii)
    for index in [3, 7, 0, 5, 1, 6, 2, 4]:
        shuffled[index] = row(radii[index])

    # The viscosity entries are NaN for a supplied-moduli solve.
    for index, radius in enumerate(radii):
        assert np.array_equal(ascending[index], descending[index], equal_nan=True), radius
        assert np.array_equal(ascending[index], shuffled[index], equal_nan=True), radius


def test_eos_call_two_layers_distinct_moduli():
    """In a two-layer body each layer's interior reports its own moduli and density."""
    num_per_layer = 14
    interface = 0.5e6
    # The interface radius is duplicated, one copy per layer.
    radius_array = np.concatenate([np.linspace(0.0, interface, num_per_layer),
                                   np.linspace(interface, PLANET_RADIUS, num_per_layer)])
    is_lower = radius_array <= interface
    solution, _, _ = _solve(
        radius_array,
        np.where(is_lower, 6000.0, 4000.0),
        np.where(is_lower, 8.0e10, 3.0e10),
        np.where(is_lower, 2.0e11, 1.0e11),
        (interface, PLANET_RADIUS))

    deep = solution.eos_call(0.25e6)
    shallow = solution.eos_call(0.75e6)
    assert deep["shear_modulus"] == pytest.approx(8.0e10, rel=1e-3)
    assert shallow["shear_modulus"] == pytest.approx(3.0e10, rel=1e-3)
    assert deep["density"] == pytest.approx(6000.0, rel=1e-3)
    assert shallow["density"] == pytest.approx(4000.0, rel=1e-3)
