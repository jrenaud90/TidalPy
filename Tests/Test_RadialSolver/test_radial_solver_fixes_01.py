"""Degree-1 reference-frame singularity, Love-number Q and lag signs, static-ocean h, and non-positive bulk moduli."""
import cmath
import math
from pathlib import Path

import numpy as np
import pytest

from TidalPy.RadialSolver import radial_solver, build_rs_input_homogeneous_layers, build_rs_input_from_data
from TidalPy.RadialSolver.love import LoveNumbers
from TidalPy.Rheology import Elastic
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.layers import Layer
from TidalPy.Material import Material, Phase

_EARTH_RADIUS = 6.371e6               # [m]
_ONE_DAY      = 2.0 * np.pi / 86400.  # [rad s-1]
_GUO_DATA     = Path(__file__).resolve().parents[2] / 'Benchmarks' / 'RadialSolver' / 'Guo+2004.npy'


# ======================================================================================================================
# Degree-1 solves: a static body has a rigid-translation mode, a dynamic one does not
# ======================================================================================================================
def _core_mantle_inputs(core_type, is_static):
    """A 0.55 R core (11000 kg m-3) under a 4500 kg m-3 solid mantle (mu = 70 GPa)."""
    core_shear = 0.0 if core_type == 'liquid' else 1.2e11
    return build_rs_input_homogeneous_layers(
        _EARTH_RADIUS,
        _ONE_DAY,
        (11000., 4500.),
        (1.2e11, 2.0e11),
        (core_shear, 7.0e10),
        (1.0e30, 1.0e30),
        (1.0e3, 1.0e30),
        (core_type, 'solid'),
        is_static,
        (False, False),
        'elastic',
        'elastic',
        radius_fraction_tuple=(0.55, 1.0),
        slice_per_layer=40)


@pytest.mark.parametrize('rtol', (1.0e-6, 1.0e-8, 1.0e-10))
@pytest.mark.parametrize('core_type', ('liquid', 'solid'))
def test_degree_one_static_body_is_singular(core_type, rtol):
    """Every layer static: a degree-1 load fails as singular (-13) at every tolerance, liquid or solid core.

    Before the structural check the static liquid core reported success with k' = 0.098, 0.296, 0.996 at rtol 1e-6,
    1e-8, 1e-10, because integration error held its surface rcond at 1e-11 to 4.5e-14, above the 1e-14 threshold.
    """
    solution = radial_solver(
        *_core_mantle_inputs(core_type, (True, True)),
        degree_l=1,
        solve_for=('loading',),
        integration_rtol=rtol,
        integration_atol=rtol * 1.0e-4,
        warnings=False)
    assert not solution.success
    assert solution.error_code == -13
    assert 'singular' in solution.message
    assert 'reference frame' in solution.message
    assert np.isnan(solution.k)


def test_degree_two_static_body_still_solves():
    """The check is degree-1 only."""
    solution = radial_solver(
        *_core_mantle_inputs('liquid', (True, True)), degree_l=2, solve_for=('loading',), warnings=False)
    assert solution.success, solution.message


def _guo_earth_inputs(is_static):
    data = np.load(_GUO_DATA)
    radius = data[:, 0] * 1.0e3
    density = data[:, 1] * 1.0e3
    shear = (data[:, 2] * 1.0e3)**2 * density
    bulk = (data[:, 3] * 1.0e3)**2 * density - (4.0 / 3.0) * shear
    dummy_viscosity = np.full_like(radius, 1.0e30)
    return build_rs_input_from_data(
        _ONE_DAY,
        radius,
        density,
        bulk,
        shear,
        dummy_viscosity,
        dummy_viscosity,
        (1.2225e6, 3.4810e6, 6.3710e6),
        ('solid', 'liquid', 'solid'),
        is_static,
        (False, False, False),
        Elastic(),
        Elastic(),
        perform_checks=False,
        warnings=False)


@pytest.mark.skipif(not _GUO_DATA.exists(), reason='Guo+2004.npy benchmark data not found')
def test_degree_one_dynamic_earth_is_determined():
    """Dynamic solid layers remove the translation: the Guo et al. (2004) Earth solves and matches their h'.

    A frame change shifts h', l', and k' together, so h' - k' is frame independent; it equals h' in the frame of the
    Earth's center of mass (k' = 0), where Guo et al. report h' = -0.2856.
    """
    solution = radial_solver(
        *_guo_earth_inputs((False, True, False)),
        degree_l=1,
        solve_for=('loading',),
        integration_method='RK45',
        integration_rtol=1.0e-10,
        integration_atol=1.0e-14,
        warnings=False)
    assert solution.success, solution.message
    h_load = complex(solution.h).real
    k_load = complex(solution.k).real
    assert h_load == pytest.approx(-0.28732, abs=2.0e-5)
    assert k_load == pytest.approx(-0.00169, abs=2.0e-5)
    assert h_load - k_load == pytest.approx(-0.2856, abs=2.0e-4)

    # The same Earth with every layer static fails instead of returning a frame-dependent answer.
    static_solution = radial_solver(
        *_guo_earth_inputs((True, True, True)),
        degree_l=1,
        solve_for=('loading',),
        integration_method='RK45',
        warnings=False)
    assert static_solution.error_code == -13


def _propagation_matrix_body(num_slices=25, radius=6.0e6, density=5500.0):
    radius_array = np.linspace(0.0, radius, num_slices)
    return (
        radius_array,
        np.full(num_slices, density),
        np.full(num_slices, 1.0e14 + 0.0j),
        np.full(num_slices, 5.0e10 + 1.0e9j),
        _ONE_DAY,
        density,
        ('solid',),
        (True,),
        (True,),
        np.asarray((radius,)))


def test_degree_one_propagation_matrix_is_singular():
    """The propagation matrix only takes static layers, so degree 1 always fails as singular."""
    solution = radial_solver(
        *_propagation_matrix_body(), degree_l=1, solve_for=('loading',), love_method='propagation_matrix')
    assert not solution.success
    assert solution.error_code == -13
    assert 'singular' in solution.message


# ======================================================================================================================
# Quality factor and phase lag for Love numbers of either sign
# ======================================================================================================================
@pytest.mark.parametrize('value', (
    0.3 - 0.01j,              # tidal k with a lag
    -0.5029 + 0.00195j,       # loading k' with a lag (Maxwell body, 1-day period)
    -2.336 + 0.00695j,        # loading h' with a lag
    1.1 - 0.4j,
    -1.1 + 0.4j,
))
def test_dissipative_love_numbers_have_positive_q_and_lag(value):
    """A lagging response has Q > 0 and lag > 0 for either sign of its real part, and sin(lag) = 1/Q."""
    love = LoveNumbers(value, value, value)
    for component in ('k', 'h', 'l'):
        quality = getattr(love, f'Q_{component}')
        lag = getattr(love, f'lag_{component}')
        assert quality > 0.0
        assert lag > 0.0
        assert math.sin(lag) == pytest.approx(1.0 / quality, rel=1.0e-12)
        assert math.tan(lag) == pytest.approx(abs(value.imag) / abs(value.real), rel=1.0e-12)


def test_negative_love_number_quality_factor_value():
    """The audit's Maxwell loading k' = -0.5029 + 0.00195j: Q = |k'| / Im(k'), about 258, not -258."""
    value = -0.5029 + 0.00195j
    love = LoveNumbers(value, value, value)
    assert love.Q_k == pytest.approx(abs(value) / value.imag, rel=1.0e-14)
    assert love.lag_k == pytest.approx(math.atan2(value.imag, -value.real), rel=1.0e-14)


@pytest.mark.parametrize('value', (0.3 + 0.01j, -0.5 - 0.002j))
def test_leading_love_numbers_have_negative_q_and_lag(value):
    """A response that leads its forcing has Q < 0 and lag < 0, whatever the sign of its real part."""
    love = LoveNumbers(value, value, value)
    assert love.Q_k < 0.0
    assert love.lag_k < 0.0


@pytest.mark.parametrize('value', (0.3 + 0.0j, -0.5 + 0.0j))
def test_elastic_love_numbers_have_infinite_q(value):
    love = LoveNumbers(value, value, value)
    assert math.isinf(love.Q_k) and love.Q_k > 0.0
    assert love.lag_k == 0.0


def test_maxwell_tidal_and_loading_q_are_positive():
    """A Maxwell sphere: the tidal k and the negative loading k' and h' all report positive Q and lag."""
    inputs = build_rs_input_homogeneous_layers(
        _EARTH_RADIUS,
        _ONE_DAY,
        (5000.,),
        (1.0e11,),
        (5.0e10,),
        (1.0e30,),
        (1.0e17,),
        ('solid',),
        (True,),
        (False,),
        'maxwell',
        'elastic',
        radius_fraction_tuple=(1.0,),
        slice_per_layer=40)
    solution = radial_solver(*inputs, solve_for=('tidal', 'loading'), warnings=False)
    assert solution.success, solution.message
    assert solution.k[1].real < 0.0 and solution.h[1].real < 0.0
    assert np.all(solution.Q_k > 0.0) and np.all(solution.lag_k > 0.0)
    assert np.all(solution.Q_h > 0.0) and np.all(solution.lag_h > 0.0)
    # The tidal and loading k lag by about the same angle for this body.
    assert solution.lag_k[1] == pytest.approx(solution.lag_k[0], rel=1.0e-2)


# ======================================================================================================================
# Static liquid surface layer: h is defined at the free surface, l is not
# ======================================================================================================================
def _ocean_world_inputs(ocean_is_static):
    """A 2500 km solid body under a 5% (by radius) ocean."""
    return build_rs_input_homogeneous_layers(
        2.5e6,
        2.0 * np.pi / (3.5 * 86400.),
        (3500., 1000.),
        (1.0e11, 2.0e9),
        (5.0e10, 0.0),
        (1.0e30, 1.0e30),
        (1.0e20, 1.0e-3),
        ('solid', 'liquid'),
        (True, ocean_is_static),
        (False, False),
        'elastic',
        'elastic',
        radius_fraction_tuple=(0.95, 1.0),
        slice_per_layer=40)


def test_static_ocean_surface_h_is_finite():
    """h = y2/rho + y5 at a static liquid free surface: 1 + k for a tide, the isostatic depression for a load."""
    solution = radial_solver(
        *_ocean_world_inputs(True),
        solve_for=('tidal', 'loading'),
        integration_rtol=1.0e-10,
        integration_atol=1.0e-14,
        warnings=False)
    assert solution.success, solution.message
    degree_l = 2
    k_tidal, k_load = solution.k
    h_tidal, h_load = solution.h
    assert cmath.isfinite(h_tidal) and cmath.isfinite(h_load)
    assert h_tidal == pytest.approx(1.0 + k_tidal, rel=1.0e-12)
    surface_density = float(solution.get_density(solution.radius))
    loading_y2_over_rho = -(2.0 * degree_l + 1.0) * solution.density_bulk / (3.0 * surface_density)
    assert h_load == pytest.approx(loading_y2_over_rho + 1.0 + k_load, rel=1.0e-8)
    # l needs y3, which a static liquid does not define.
    assert np.all(np.isnan(solution.l))

    # The dense getter agrees at the surface: y1 g = h, y2 = the boundary condition (0 for a tide).
    y_surface = solution.get_radial_solution(solution.radius, 0)
    assert y_surface[0] * solution.surface_gravity == pytest.approx(h_tidal, rel=1.0e-10)
    assert abs(y_surface[1]) == 0.0
    # Inside the static ocean y1 stays undefined.
    y_inside = solution.get_radial_solution(0.975 * solution.radius, 0)
    assert np.isnan(y_inside[0])


def test_static_ocean_h_is_close_to_dynamic_ocean():
    """A dynamic ocean at a 3.5-day period gives nearly the same h (1.2756 against 1 + k = 1.2722)."""
    static = radial_solver(*_ocean_world_inputs(True), solve_for=('tidal',), warnings=False,
                           integration_rtol=1.0e-10, integration_atol=1.0e-14)
    dynamic = radial_solver(*_ocean_world_inputs(False), solve_for=('tidal',), warnings=False,
                            integration_rtol=1.0e-10, integration_atol=1.0e-14)
    assert static.success and dynamic.success
    assert complex(static.h) == pytest.approx(complex(dynamic.h), rel=5.0e-3)


# ======================================================================================================================
# Compressible layers need a positive bulk modulus
# ======================================================================================================================
def _uniform_world(bulk_modulus, is_incompressible=False):
    """A 1000 km uniform body; a bulk modulus of None gives the constant-density EOS a zero bulk modulus."""
    radius = 1.0e6
    density = 3000.0
    mass = (4.0 / 3.0) * math.pi * radius**3 * density
    world = BaseWorld("uniform", radius, mass)
    material = Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": density,
             "bulk_modulus_pa": 0.0 if bulk_modulus is None else bulk_modulus},
        shear_modulus={"model": "constant", "shear_modulus_pa": 5.0e10}))
    layer = Layer("mantle", 0, 0.0, radius, mass, material, is_incompressible=is_incompressible)
    world.add_layer(layer)
    world.solve_eos()
    return world


def test_world_compressible_layer_without_bulk_modulus_fails():
    """Before the check this solved with k2 = 0.008681 (Poisson ratio -1), against 0.007928 at K = 100 GPa."""
    result = _uniform_world(None).solve_love_numbers(frequency=1.0e-5)
    assert not result['success']
    assert result['error_code'] == -15
    assert 'bulk modulus' in result['message']
    assert 'incompressible' in result['message']


def test_world_compressible_layer_with_bulk_modulus_solves():
    result = _uniform_world(1.0e11).solve_love_numbers(frequency=1.0e-5)
    assert result['success'], result['message']
    assert result['love_number_k'].real == pytest.approx(0.0079276, rel=1.0e-4)


def test_world_incompressible_layer_needs_no_bulk_modulus():
    result = _uniform_world(None, is_incompressible=True).solve_love_numbers(
        frequency=1.0e-5, love_method='propagation_matrix')
    assert result['success'], result['message']
    assert result['love_number_k'].real == pytest.approx(0.0079039, rel=1.0e-4)


def test_standalone_compressible_solid_with_zero_bulk_fails():
    """A zero bulk-modulus array reaches the solve (a negative one is refused by the input checks)."""
    radius_array = np.linspace(0.0, 1.0e6, 30)
    solution = radial_solver(
        radius_array,
        np.full(30, 3000.0),
        np.full(30, 0.0 + 0.0j),
        np.full(30, 5.0e10 + 0.0j),
        1.0e-5,
        3000.0,
        ('solid',),
        (True,),
        (False,),
        np.asarray((1.0e6,)),
        warnings=False)
    assert not solution.success
    assert solution.error_code == -15


def test_static_liquid_layer_ignores_its_bulk_modulus():
    """A static liquid never reads the bulk modulus, so zero there is allowed."""
    inputs = build_rs_input_homogeneous_layers(
        2.5e6,
        2.0 * np.pi / (3.5 * 86400.),
        (3500., 1000.),
        (1.0e11, 0.0),
        (5.0e10, 0.0),
        (1.0e30, 1.0e30),
        (1.0e20, 1.0e-3),
        ('solid', 'liquid'),
        (True, True),
        (False, False),
        'elastic',
        'elastic',
        radius_fraction_tuple=(0.95, 1.0),
        slice_per_layer=40)
    solution = radial_solver(*inputs, warnings=False)
    assert solution.success, solution.message
