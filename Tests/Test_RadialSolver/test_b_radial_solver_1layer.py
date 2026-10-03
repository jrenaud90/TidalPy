"""The shooting method on a 1-layer planet for each layer flag, integrator, degree, and solve type."""
import numpy as np
import pytest

from TidalPy.RadialSolver.solver import radial_solver
from TidalPy.Rheology import Maxwell

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

# The incompressible, elastic homogeneous sphere's k_l = (3 / (2 (l - 1))) / (1 + (2 l^2 + 4 l + 3) mu / (l rho g R)).
# This planet is compressible (K = 2 mu) and barely viscous (Maxwell time about 40 orbits), which moves k_l by well
# under this tolerance; a wrong boundary condition or start moves it by far more.
K_TOLERANCE = 0.25


def _homogeneous_k(degree_l):
    radius = float(radius_array[-1])
    surface_gravity = 4.0 / 3.0 * np.pi * 6.67430e-11 * planet_bulk_density * radius
    rigidity = (2 * degree_l**2 + 4 * degree_l + 3) * float(shear_array[0]) / (
        degree_l * planet_bulk_density * surface_gravity * radius)
    return (3.0 / (2.0 * (degree_l - 1))) / (1.0 + rigidity)


def _solve_1layer(
        layer_type,
        is_static,
        is_incompressible,
        method,
        degree_l,
        use_kamata,
        solve_for,
        supported,
        **kwargs):
    """Solve a 1-layer planet: an unsupported start raises NotImplementedError, and a supported one succeeds with a
    k_l near the homogeneous sphere's. All-liquid planets are skipped."""
    # A 1-layer all-liquid planet is numerically unstable.
    if layer_type != 'solid':
        pytest.skip('Planets with 1-layer liquid are not currently very stable. Skipping tests.')

    def solve():
        return radial_solver(
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
            use_kamata=use_kamata,
            integration_method=method,
            integration_rtol=1.0e-7,
            integration_atol=1.0e-10,
            scale_rtols_bylayer_type=False,
            max_num_steps=5_000_000,
            expected_size=250,
            max_step=0,
            verbose=False,
            nondimensionalize=True,
            starting_radius=0.2 * radius_array[-1],
            raise_on_fail=True,
            **kwargs)

    if not supported(layer_type, is_static, is_incompressible, use_kamata):
        with pytest.raises(NotImplementedError):
            solve()
        return None
    out = solve()
    assert out.success
    assert type(out.message) is str
    assert type(out.result) is np.ndarray
    assert out.result.shape == (len(solve_for) * 6, N)
    if 'tidal' in solve_for:
        k_tidal = np.atleast_1d(out.k)[solve_for.index('tidal')]
        assert k_tidal.real == pytest.approx(_homogeneous_k(degree_l), rel=K_TOLERANCE)
    return out



@pytest.mark.parametrize('layer_type', ("solid", "liquid"))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('method', ("rk23", "rk45", "dop853"))
@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('use_kamata', (True, False))
@pytest.mark.parametrize('solve_for', (('free',), ('tidal',), ('loading',)))
def test_radial_solver_1layer(
        layer_type,
        is_static,
        is_incompressible,
        method,
        degree_l,
        use_kamata,
        solve_for,
        starting_conditions_supported):
    """A single solve type succeeds with a (6, N) result; the log_info kwarg is accepted."""
    _solve_1layer(
        layer_type,
        is_static,
        is_incompressible,
        method,
        degree_l,
        use_kamata,
        solve_for,
        starting_conditions_supported,
        log_info=True)


@pytest.mark.parametrize('layer_type', ("solid", "liquid"))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('method', ("rk23", "rk45", "dop853"))
@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('use_kamata', (True, False))
def test_radial_solver_1layer_solve_for_both(
        layer_type,
        is_static,
        is_incompressible,
        method,
        degree_l,
        use_kamata,
        starting_conditions_supported):
    """Solving for tidal and loading together succeeds and every solution attribute has its expected type and shape."""
    solve_for = ('tidal', 'loading')
    out = _solve_1layer(
        layer_type,
        is_static,
        is_incompressible,
        method,
        degree_l,
        use_kamata,
        solve_for,
        starting_conditions_supported)
    if out is None:
        return
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
