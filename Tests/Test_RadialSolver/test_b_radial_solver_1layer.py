"""The shooting method on a 1-layer planet for each layer flag, integrator, degree, and solve type."""
from functools import lru_cache

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.RadialSolver.solver import radial_solver
from TidalPy.Rheology import Maxwell

from starting_methods import STARTING_METHODS

frequency = 1.0 / (86400. * 0.2)
N = 100
radius_array = np.linspace(0.0, 6000.e3, N)
density_array = 3500. * np.ones_like(radius_array)
bulk_modulus_array = 1.0e11 * np.ones(radius_array.size, dtype=np.complex128, order='C')
viscosity_array = 1.0e20 * np.ones_like(radius_array)
shear_array = 5.0e10 * np.ones_like(radius_array)
complex_shear_modulus_array = Maxwell().calc_complex_modulus_vectorize_modulus(shear_array, viscosity_array, frequency)
planet_bulk_density = float(density_array[0])
upper_radius_by_layer = np.asarray((radius_array[-1],))
integration_settings = dict(
    integration_rtol=1.0e-7,
    integration_atol=1.0e-10,
    scale_rtols_bylayer_type=False,
    max_num_steps=5_000_000,
    expected_size=250,
    max_step=0,
    verbose=False,
    nondimensionalize=True,
    starting_radius=0.2 * radius_array[-1])

# The incompressible, elastic homogeneous sphere's k_l = (3 / (2 (l - 1))) / (1 + (2 l^2 + 4 l + 3) mu / (l rho g R)).
# This planet is compressible (K = 2 mu) and barely viscous (Maxwell time about 40 orbits), which moves k_l by well
# under this tolerance; a wrong boundary condition or start moves it by far more.
K_TOLERANCE = 0.25

# A liquid planet is checked against the homogeneous inviscid sphere (see _homogeneous_liquid_love). The regular starts
# match it to 2e-5 (RK23 at this rtol). Unity's unit vectors at 0.2 R carry singular content that decays only as
# (0.2)^(2l + 1) on the way out, which leaves 1.4e-3 at degree 2.
LIQUID_TOLERANCE = 1.0e-4
UNITY_LIQUID_TOLERANCE = 5.0e-3


def _homogeneous_k(degree_l):
    radius = float(radius_array[-1])
    surface_gravity = 4.0 / 3.0 * np.pi * 6.67430e-11 * planet_bulk_density * radius
    rigidity = (2 * degree_l**2 + 4 * degree_l + 3) * float(shear_array[0]) / (
        degree_l * planet_bulk_density * surface_gravity * radius)
    return (3.0 / (2.0 * (degree_l - 1))) / (1.0 + rigidity)


@lru_cache(maxsize=None)
def _dynamic_compressible_load_reference(degree_l):
    """(k, h, l) of a load on the dynamic compressible liquid planet, solved at rtol 1e-10."""
    solution = radial_solver(
        radius_array,
        density_array,
        bulk_modulus_array,
        complex_shear_modulus_array,
        frequency,
        planet_bulk_density,
        ('liquid',),
        (False,),
        (False,),
        upper_radius_by_layer,
        degree_l=degree_l,
        solve_for=('loading',),
        starting_method='takeuchi',
        integration_method='DOP853',
        integration_rtol=1.0e-10,
        integration_atol=1.0e-12,
        starting_radius=0.2 * radius_array[-1],
        raise_on_fail=True)
    return complex(solution.k), complex(solution.h), complex(solution.l)


def _homogeneous_liquid_love(degree_l, is_static, is_incompressible, solve_type):
    """(k, h, l) of a homogeneous inviscid liquid sphere of this planet's density and radius.

    Static, a tide gives k = 3 / (2 (l - 1)) and h = (2 l + 1) / (2 (l - 1)), and a load k' = -1 and h' = -(2 l + 1) / 3
    (the displaced liquid's mass cancels the load's). Dynamic, both scale by 1 / (1 - omega^2 / omega_l^2), with
    omega_l^2 = (2 l (l - 1) / (2 l + 1)) (4 pi G rho / 3) the sphere's fundamental mode; its flow is a potential flow,
    so the Shida number is h / l. A static liquid has no horizontal displacement, so l is NaN.

    Assumptions
    -----------
    - The static liquid equations do not read the bulk modulus, so a compressible static liquid gives the same numbers.
    - A tide leaves the Lagrangian pressure perturbation zero throughout a constant-density sphere, so the dynamic
      compressible liquid compresses nowhere and matches the incompressible forms.
    - A load's weight does compress it: k' (set by the conserved mass) keeps its form, but h' and l' do not, and are
      taken from a solve at rtol 1e-10.
    """
    surface_gravity_over_radius = 4.0 / 3.0 * np.pi * G * planet_bulk_density
    if is_static:
        resonance = 1.0
    else:
        mode_frequency2 = 2.0 * degree_l * (degree_l - 1) / (2.0 * degree_l + 1.0) * surface_gravity_over_radius
        resonance = 1.0 / (1.0 - frequency**2 / mode_frequency2)
    if solve_type == 'tidal':
        k_love = 3.0 / (2.0 * (degree_l - 1)) * resonance
        h_love = (2.0 * degree_l + 1.0) / (2.0 * (degree_l - 1)) * resonance
    else:
        k_love = -resonance
        h_love = -(2.0 * degree_l + 1.0) / 3.0 * resonance
    l_love = np.nan if is_static else h_love / degree_l
    if (solve_type == 'loading') and not (is_static or is_incompressible):
        _, h_love, l_love = _dynamic_compressible_load_reference(degree_l)
    return k_love, h_love, l_love


def _solve_1layer(
        layer_type,
        is_static,
        is_incompressible,
        method,
        degree_l,
        starting_method,
        solve_for,
        **kwargs):
    """Solve a 1-layer planet: every start succeeds, a solid with a tidal k_l near the homogeneous sphere's and a liquid
    with the homogeneous liquid sphere's tidal and load Love numbers."""
    out = radial_solver(
        radius_array,
        density_array,
        bulk_modulus_array,
        complex_shear_modulus_array,
        frequency,
        planet_bulk_density,
        (layer_type,),
        (is_static,),
        (is_incompressible,),
        upper_radius_by_layer,
        degree_l=degree_l,
        solve_for=solve_for,
        starting_method=starting_method,
        integration_method=method,
        raise_on_fail=True,
        **integration_settings,
        **kwargs)
    assert out.success
    assert type(out.message) is str
    assert type(out.result) is np.ndarray
    assert out.result.shape == (len(solve_for) * 6, N)
    for index, solve_type in enumerate(solve_for):
        k_value = np.atleast_1d(out.k)[index]
        if layer_type == 'solid':
            if solve_type == 'tidal':
                assert k_value.real == pytest.approx(_homogeneous_k(degree_l), rel=K_TOLERANCE)
            continue
        if solve_type == 'free':
            continue
        tolerance = UNITY_LIQUID_TOLERANCE if starting_method == 'unity' else LIQUID_TOLERANCE
        expected = _homogeneous_liquid_love(degree_l, is_static, is_incompressible, solve_type)
        found = (k_value, np.atleast_1d(out.h)[index], np.atleast_1d(out.l)[index])
        for name, value, target in zip('khl', found, expected):
            if np.isnan(target):
                assert np.isnan(value), (solve_type, name, value)
            else:
                assert value == pytest.approx(target, rel=tolerance), (solve_type, name, value, target)
    return out


@pytest.mark.parametrize('layer_type', ("solid", "liquid"))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('method', ("rk23", "rk45", "dop853"))
@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('starting_method', STARTING_METHODS)
@pytest.mark.parametrize('solve_for', (('free',), ('tidal',), ('loading',)))
def test_radial_solver_1layer(
        layer_type,
        is_static,
        is_incompressible,
        method,
        degree_l,
        starting_method,
        solve_for):
    """A single solve type succeeds with a (6, N) result; the log_info kwarg is accepted."""
    _solve_1layer(
        layer_type,
        is_static,
        is_incompressible,
        method,
        degree_l,
        starting_method,
        solve_for,
        log_info=True)


@pytest.mark.parametrize('layer_type', ("solid", "liquid"))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('method', ("rk23", "rk45", "dop853"))
@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('starting_method', STARTING_METHODS)
def test_radial_solver_1layer_solve_for_both(
        layer_type,
        is_static,
        is_incompressible,
        method,
        degree_l,
        starting_method):
    """Solving for tidal and loading together succeeds and every solution attribute has its expected type and shape."""
    solve_for = ('tidal', 'loading')
    out = _solve_1layer(
        layer_type,
        is_static,
        is_incompressible,
        method,
        degree_l,
        starting_method,
        solve_for)
    num_solve_for = len(solve_for)

    assert type(out.error_code) is int
    assert out.error_code == 0
    assert type(out.eos_error_code) is int
    assert out.eos_error_code == 0
    assert type(out.eos_message) is str
    assert out.eos_success == True
    assert type(out.eos_pressure_error) is float
    assert type(out.eos_iterations) is int
    assert out.eos_iterations >= 1
    assert type(out.eos_steps_taken) is np.ndarray
    mid_radius = 0.5 * out.radius
    for getter in (out.get_gravity, out.get_pressure, out.get_mass, out.get_moi, out.get_density,
                   out.get_shear_modulus, out.get_bulk_modulus):
        assert type(getter(mid_radius)) is float
        assert type(getter(out.sample_radii(8))) is np.ndarray
    assert out.get_density(mid_radius) > 0.0
    assert type(out.sample_radii()) is np.ndarray
    assert type(out.layer_upper_radius_array) is np.ndarray
    for name in ('radius', 'volume', 'mass', 'moi', 'moi_factor', 'moi_sphere_ratio', 'density_bulk',
                 'central_pressure', 'surface_pressure', 'surface_gravity'):
        assert type(getattr(out, name)) is float, name
    for name in ('num_ytypes', 'num_layers', 'degree_l'):
        assert type(getattr(out, name)) is int, name
    assert type(out['tidal']) is np.ndarray
    assert type(out['loading']) is np.ndarray
    assert type(out.love) is np.ndarray
    assert out.love.shape == (num_solve_for, 3)
    # Some builds expose Q / lag aliases, others only Q_k / lag_k.
    quality = out.Q if hasattr(out, 'Q') else out.Q_k
    lag = out.lag if hasattr(out, 'lag') else out.lag_k
    for array in (out.k, out.h, out.l, quality, lag):
        assert type(array) is np.ndarray
        assert array.shape == (num_solve_for,)
    assert type(out.steps_taken) is np.ndarray

    diagnostics = out.print_diagnostics(print_diagnostics=False, log_diagnostics=False)
    # Every dimensional value of the structure carries its unit.
    for label, unit in (("Pressure Error", " Pa"), ("Central Pressure", " Pa"), ("Mass", " kg"), ("MOI", " kg m2"),
                        ("Surface gravity", " m s-2")):
        line = next(line for line in diagnostics.splitlines() if line.strip().startswith(label + ":"))
        assert unit in line, line

    eos_result = out.eos_call(radius=1.5e6)
    assert type(eos_result) is dict
    assert np.isfinite(eos_result["density"])
