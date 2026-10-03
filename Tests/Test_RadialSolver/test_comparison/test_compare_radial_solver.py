"""Full radial_solver Love numbers and radial solutions against the frozen classic solver's (1-layer solid).

The frozen `case_status` per case key is 0 (classic succeeded), 1 (NotImplementedError), or 2 (unsuccessful);
successful cases' `result` and `love` arrays are stacked in `compared_case_keys` order.
"""
from pathlib import Path

import numpy as np
import pytest

from TidalPy.Rheology import Maxwell
from TidalPy.RadialSolver.solver import radial_solver

CLASSIC_NOT_IMPLEMENTED = 1
CLASSIC_UNSUCCESSFUL = 2

FROZEN_PATH = Path(__file__).parent / "frozen" / "test_compare_radial_solver.npz"
with np.load(FROZEN_PATH, allow_pickle=False) as frozen_file:
    CLASSIC_STATUS = dict(zip(frozen_file["case_keys"].tolist(), frozen_file["case_status"].tolist()))
    CLASSIC_INDEX = {key: index for index, key in enumerate(frozen_file["compared_case_keys"].tolist())}
    CLASSIC_RESULT_NUM_ROWS = frozen_file["result_num_rows"]
    CLASSIC_RESULTS = frozen_file["results"]
    CLASSIC_LOVES = frozen_file["loves"]

frequency = 2.0 * np.pi / (86400. * 1.0)
N = 10
radius_array = np.linspace(0.0, 6000.e3, N)
bulk_density = 3500.
density_array = bulk_density * np.ones_like(radius_array)
bulk_modulus_array = 1.0e11 * np.ones(N, dtype=np.complex128, order='C')
viscosity_array = 1.0e20 * np.ones_like(radius_array)
shear_array = 5.0e10 * np.ones_like(radius_array)
complex_shear_modulus_array = Maxwell().calc_complex_modulus_vectorize_modulus(shear_array, viscosity_array, frequency)
upper_radius_by_layer = np.asarray((radius_array[-1],))

MANUAL_STARTING_RADIUS = 0.1 * radius_array[-1]


def _case_key(
        is_static,
        is_incompressible,
        degree_l,
        use_kamata,
        solve_for,
        starting_radius,
        integration_method,
        nondimensionalize):
    starting_label = "auto" if starting_radius == 0.0 else "manual"
    return (f"static_{is_static}__incompressible_{is_incompressible}__degree_l_{degree_l}__kamata_{use_kamata}"
            f"__solve_for_{'+'.join(solve_for)}__start_{starting_label}__method_{integration_method}"
            f"__nondimensionalize_{nondimensionalize}")


def _classic_solution(case_key):
    """The classic `result` and `love` arrays for a case it solved."""
    index = CLASSIC_INDEX[case_key]
    num_rows = CLASSIC_RESULT_NUM_ROWS[index]
    return CLASSIC_RESULTS[index, :num_rows, :], CLASSIC_LOVES[index, :num_rows // 6, :]


@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('use_kamata', (True, False))
@pytest.mark.parametrize('solve_for', (('tidal',), ('loading',), ('tidal', 'loading')))
@pytest.mark.parametrize('starting_radius', (0.0, MANUAL_STARTING_RADIUS))
@pytest.mark.parametrize('integration_method', ('RK45', 'DOP853'))
@pytest.mark.parametrize('nondimensionalize', (False, True))
def test_compare_radial_solver_1layer_solid(
        is_static,
        is_incompressible,
        degree_l,
        use_kamata,
        solve_for,
        starting_radius,
        integration_method,
        nondimensionalize):
    """Love numbers and radial solutions match the classic solver for 1-layer solid planets."""
    case_key = _case_key(
        is_static,
        is_incompressible,
        degree_l,
        use_kamata,
        solve_for,
        starting_radius,
        integration_method,
        nondimensionalize)
    classic_status = CLASSIC_STATUS[case_key]
    if classic_status == CLASSIC_NOT_IMPLEMENTED:
        pytest.skip('Not implemented in original RadialSolver.')
    if classic_status == CLASSIC_UNSUCCESSFUL:
        pytest.skip("Old solver was not successful, no need to compare.")

    old_result, old_love = _classic_solution(case_key)

    try:
        new_out = radial_solver(
            radius_array,
            density_array,
            bulk_modulus_array,
            complex_shear_modulus_array,
            frequency,
            bulk_density,
            ('solid',),
            (is_static,),
            (is_incompressible,),
            upper_radius_by_layer,
            degree_l=degree_l,
            solve_for=solve_for,
            use_kamata=use_kamata,
            integration_method=integration_method,
            integration_rtol=1.0e-7,
            integration_atol=1.0e-10,
            scale_rtols_bylayer_type=False,
            max_num_steps=5_000_000,
            expected_size=250,
            max_step=0,
            verbose=False,
            nondimensionalize=nondimensionalize,
            starting_radius=starting_radius,
            # The EOS defaults differ from the classic solver's; a deep starting radius amplifies any structure
            # difference, so both EOS solves were pinned tight.
            eos_rtol=1.0e-10,
            eos_atol=1.0e-14,
            raise_on_fail=True,
        )
    except NotImplementedError:
        pytest.skip('Not implemented in RadialSolver.')

    assert new_out.success, f"New solver failed: {new_out.message}"

    assert old_result.shape == new_out.result.shape
    # Different LU and ODE implementations agree to about 1e-4 in the interior. The classic Kamata start for a dynamic
    # incompressible solid has two nearly parallel solutions at low frequency and loses digits there.
    interior_rtol = 3.0e-4 if (use_kamata and (not is_static) and is_incompressible) else 1.0e-4
    np.testing.assert_allclose(
        new_out.result[:, :-1],
        old_result[:, :-1],
        rtol=interior_rtol,
        atol=1e-6,
        err_msg="Interior radial solution arrays differ.")
    # The surface rows pinned to homogeneous boundary conditions are cancellation residuals (constants near 1e11
    # cancel to 1e-4 for a deep start), pure roundoff that differs between LAPACK and Eigen, hence the loose atol.
    np.testing.assert_allclose(
        new_out.result[:, -1],
        old_result[:, -1],
        rtol=1e-2,
        atol=1e-2,
        err_msg="Surface radial solution values differ.")
    assert old_love.shape == new_out.love.shape
    np.testing.assert_allclose(new_out.love, old_love, rtol=1e-3, err_msg="Love numbers differ.")
