"""The shooting method on a solid-liquid-solid planet for each layer flag, integrator, degree, and starting radius."""
import numpy as np
import pytest

from TidalPy.exceptions import SolutionFailedError
from TidalPy.RadialSolver.solver import radial_solver
from TidalPy.Rheology import Maxwell

frequency = 1.0 / (86400. * 0.2)
N = 20
planet_r = 6000.0e3
icb_r = planet_r * (1. / 3.0)
cmb_r = planet_r * (2. / 3.0)
radius_array = np.concatenate((
    np.linspace(0.0, icb_r, N),
    np.linspace(icb_r, cmb_r, N),
    np.linspace(cmb_r, planet_r, N)
    ))
# Inner core, outer core, and mantle values, N slices each.
density_array = np.repeat((8500., 7000., 3500.), N)
bulk_modulus_array = 1.0e11 * np.ones(radius_array.size, dtype=np.complex128, order='C')
viscosity_array = np.repeat((1.0e26, 1.0e6, 1.0e20), N)
shear_array = np.repeat((1.0e11, 0.0, 5.0e10), N)
complex_shear_modulus_array = Maxwell().calc_complex_modulus_vectorize_modulus(shear_array, viscosity_array, frequency)

planet_bulk_density = np.average(density_array)
upper_radius_by_layer = np.asarray((icb_r, cmb_r, planet_r))
layer_types = ("solid", "liquid", "solid")
# Inner core, automatic, mid outer core, and mid mantle.
starting_radii = (0.2 * planet_r, 0.0, icb_r + (cmb_r - icb_r) / 2.0, cmb_r + (planet_r - cmb_r) / 2.0)


def _check_3layer(
        solid_is_static,
        liquid_is_static,
        solid_is_incompressible,
        liquid_is_incompressible,
        method,
        degree_l,
        use_kamata,
        solve_for,
        starting_radius):
    """Solve the 3-layer planet and check the output; skip unsupported inputs and the known unstable pairing."""
    # A compressible dynamic liquid over an incompressible solid can fail to integrate (step size underflow).
    known_unstable = (not liquid_is_incompressible) and solid_is_incompressible
    unstable_reason = 'Integration Failed. Compressible liquid with incompressible solid below is not very stable.'
    try:
        out = radial_solver(
            radius_array,
            density_array,
            bulk_modulus_array,
            complex_shear_modulus_array,
            frequency,
            planet_bulk_density,
            layer_types,
            (solid_is_static, liquid_is_static, solid_is_static),
            (solid_is_incompressible, liquid_is_incompressible, solid_is_incompressible),
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
            starting_radius=starting_radius,
            raise_on_fail=True,
            verbose=False,
            nondimensionalize=True)
    except NotImplementedError as error:
        pytest.skip(f'function does not currently support requested inputs. Skipping Test. Details: {error}')
    except SolutionFailedError as error:
        # Only the liquid layer's step-size collapse is expected; any other failure is a real one.
        if known_unstable and ('Integration problem at layer 2' in str(error)):
            pytest.skip(unstable_reason)
        raise

    assert out.success
    assert type(out.message) is str
    assert type(out.result) is np.ndarray
    assert out.result.shape == (len(solve_for) * 6, radius_array.size)


@pytest.mark.parametrize('solid_is_static', (True, False))
@pytest.mark.parametrize('liquid_is_static', (True, False))
@pytest.mark.parametrize('solid_is_incompressible', (True, False))
@pytest.mark.parametrize('liquid_is_incompressible', (True, False))
@pytest.mark.parametrize('method', ("rk23", "rk45", "dop853"))
@pytest.mark.parametrize('degree_l', (2, 3, 10))
@pytest.mark.parametrize('use_kamata', (True,))
@pytest.mark.parametrize('solve_for', (('free',), ('tidal',), ('loading',)))
@pytest.mark.parametrize('starting_radius', starting_radii)
def test_radial_solver_3layer(
        solid_is_static,
        liquid_is_static,
        solid_is_incompressible,
        liquid_is_incompressible,
        method,
        degree_l,
        use_kamata,
        solve_for,
        starting_radius):
    """A single solve type succeeds with a (6, N) result."""
    _check_3layer(
        solid_is_static,
        liquid_is_static,
        solid_is_incompressible,
        liquid_is_incompressible,
        method,
        degree_l,
        use_kamata,
        solve_for,
        starting_radius)


@pytest.mark.parametrize('solid_is_static', (True, False))
@pytest.mark.parametrize('liquid_is_static', (True, False))
@pytest.mark.parametrize('solid_is_incompressible', (True, False))
@pytest.mark.parametrize('liquid_is_incompressible', (True, False))
@pytest.mark.parametrize('method', ("rk23", "rk45", "dop853"))
@pytest.mark.parametrize('degree_l', (2,))
@pytest.mark.parametrize('use_kamata', (True,))
def test_radial_solver_3layer_solve_for_both(
        solid_is_static,
        liquid_is_static,
        solid_is_incompressible,
        liquid_is_incompressible,
        method,
        degree_l,
        use_kamata):
    """Solving for tidal and loading together succeeds with a (12, N) result."""
    _check_3layer(
        solid_is_static,
        liquid_is_static,
        solid_is_incompressible,
        liquid_is_incompressible,
        method,
        degree_l,
        use_kamata,
        ('tidal', 'loading'),
        0.2 * planet_r)
