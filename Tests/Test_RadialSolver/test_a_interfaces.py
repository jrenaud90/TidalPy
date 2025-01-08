import pytest

import numpy as np

from TidalPy.RadialSolver.shooting import find_num_shooting_solutions
from TidalPy.RadialSolver.interfaces.interfaces import solve_upper_y_at_interface

tpy_0p4_results = {
    (True, True, True, True): np.asarray(
        [(0.1+0.1j), (0.2+0.2j), (0.3+0.3j), (0.4+0.4j), (0.5+0.5j), (0.6+0.6j),
         (-0.1-0.1j), (-0.2-0.2j), (-0.3-0.3j), (-0.4-0.4j), (-0.5-0.5j), (-0.6-0.6j),
         (0.16000000000000003+0.16000000000000003j), (0.34+0.34j), (0.54+0.54j), (0.76+0.76j), (1+1j), (12.06+12.06j)],
        dtype=np.complex128, order="C"),
    (True, True, True, False): np.asarray(
        [(0.1+0.1j), (0.2+0.2j), (0.3+0.3j), (0.4+0.4j), (0.5+0.5j), (0.6+0.6j),
         (-0.1-0.1j), (-0.2-0.2j), (-0.3-0.3j), (-0.4-0.4j), (-0.5-0.5j), (-0.6-0.6j),
         (0.16000000000000003+0.16000000000000003j), (0.34+0.34j), (0.54+0.54j), (0.76+0.76j), (1+1j), (12.06+12.06j)],
        dtype=np.complex128, order="C"),
    (True, True, False, True): np.asarray(
        [0j, 0j, np.nan, np.nan, np.nan, np.nan],
        dtype=np.complex128, order="C"),
    (True, True, False, False): np.asarray(
        [(0.01578947368421052+0.01578947368421052j), (0.02105263157894738+0.02105263157894738j), (-0.02631578947368418-0.02631578947368418j), (-5.747368421052632-5.747368421052632j), np.nan, np.nan,
         (-0.01578947368421052-0.01578947368421052j), (-0.02105263157894738-0.02105263157894738j), (0.02631578947368418+0.02631578947368418j), (5.747368421052632+5.747368421052632j), np.nan, np.nan],
        dtype=np.complex128, order="C"),
    (True, False, True, True): np.asarray(
        [(0.1+0.1j), (0.2+0.2j), (0.3+0.3j), (0.4+0.4j), (0.5+0.5j), (0.6+0.6j),
         (-0.1-0.1j), (-0.2-0.2j), (-0.3-0.3j), (-0.4-0.4j), (-0.5-0.5j), (-0.6-0.6j),
         (0.16000000000000003+0.16000000000000003j), (0.34+0.34j), (0.54+0.54j), (0.76+0.76j), (1+1j), (12.06+12.06j)],
        dtype=np.complex128, order="C"),
    (True, False, True, False): np.asarray(
        [(0.1+0.1j), (0.2+0.2j), (0.3+0.3j), (0.4+0.4j), (0.5+0.5j), (0.6+0.6j),
         (-0.1-0.1j), (-0.2-0.2j), (-0.3-0.3j), (-0.4-0.4j), (-0.5-0.5j), (-0.6-0.6j),
         (0.16000000000000003+0.16000000000000003j), (0.34+0.34j), (0.54+0.54j), (0.76+0.76j), (1+1j), (12.06+12.06j)],
        dtype=np.complex128, order="C"),
    (True, False, False, True): np.asarray(
        [0j, 0j, np.nan, np.nan, np.nan, np.nan],
        dtype=np.complex128, order="C"),
    (True, False, False, False): np.asarray(
        [(0.01578947368421052+0.01578947368421052j), (0.02105263157894738+0.02105263157894738j), (-0.02631578947368418-0.02631578947368418j), (-5.747368421052632-5.747368421052632j), np.nan, np.nan,
         (-0.01578947368421052-0.01578947368421052j), (-0.02105263157894738-0.02105263157894738j), (0.02631578947368418+0.02631578947368418j), (5.747368421052632+5.747368421052632j), np.nan, np.nan],
        dtype=np.complex128, order="C"),
    (False, True, True, True): np.asarray(
        [0j, (-3800-3800j), 0j, 0j, (0.5+0.5j), (0.600001180416904-9.599998819583096j),
         (1+0j), (20520+0j), 0j, 0j, 0j, (-6.374251281747724e-06+0j),
         0j, 0j, (1+0j), 0j, 0j, 0j],
        dtype=np.complex128, order="C"),
    (False, True, True, False): np.asarray(
        [0j, (-3800-3800j), 0j, 0j, (0.5+0.5j), (0.600001180416904-9.599998819583096j),
         (1+0j), (20520+0j), 0j, 0j, 0j, (-6.374251281747724e-06+0j),
         0j, 0j, (1+0j), 0j, 0j, 0j],
        dtype=np.complex128, order="C"),
    (False, True, False, True): np.asarray(
        [(0.5+0.5j), (0.6-9.6j), np.nan, np.nan, np.nan, np.nan],
        dtype=np.complex128, order="C"),
    (False, True, False, False): np.asarray(
        [0j, (-3800-3800j), (0.5+0.5j), (0.600001180416904-9.599998819583096j), np.nan, np.nan,
         (1+0j), (20520+0j), 0j, (-6.374251281747724e-06+0j), np.nan, np.nan],
        dtype=np.complex128, order="C"),
    (False, False, True, True): np.asarray(
        [(0.1+0.1j), (0.2+0.2j), 0j, 0j, (0.5+0.5j), (0.6+0.6j),
         (0.16000000000000003+0.16000000000000003j), (0.34+0.34j), 0j, 0j, (1+1j), (12.06+12.06j),
         0j, 0j, (1+0j), 0j, 0j, 0j],
        dtype=np.complex128, order="C"),
    (False, False, True, False): np.asarray(
        [(0.1+0.1j), (0.2+0.2j), 0j, 0j, (0.5+0.5j), (0.6+0.6j),
         (0.16000000000000003+0.16000000000000003j), (0.34+0.34j), 0j, 0j, (1+1j), (12.06+12.06j),
         0j, 0j, (1+0j), 0j, 0j, 0j],
        dtype=np.complex128, order="C"),
    (False, False, False, True): np.asarray(
        [(0.09505598613897154+0.09505598613897154j), (-4.283624807144646-4.283624807144646j), np.nan, np.nan, np.nan, np.nan],
        dtype=np.complex128, order="C"),
    (False, False, False, False): np.asarray(
        [(0.1+0.1j), (0.2+0.2j), (0.5+0.5j), (0.6+0.6j), np.nan, np.nan,
         (0.16000000000000003+0.16000000000000003j), (0.34+0.34j), (1+1j), (12.06+12.06j), np.nan, np.nan],
        dtype=np.complex128, order="C")
}

static_liquid_density=7600.
interface_gravity=2.7
G_to_use=6.67430e-11

y_lower_solid = np.asarray(
    ((0.1+0.1j, 0.2+0.2j, 0.3+0.3j, 0.4+0.4j, 0.5+0.5j, 0.6+0.6j),
    (-0.1-0.1j, -0.2-0.2j, -0.3-0.3j, -0.4-0.4j, -0.5-0.5j,-0.6-0.6j),
    (1.6*(0.1+0.1j), 1.7*(0.2+0.2j), 1.8*(0.3+0.3j), 1.9*(0.4+0.4j), 2.0*(0.5+0.5j), 20.1*(0.6+0.6j))), dtype=np.complex128
)
y_lower_liquid = np.asarray(
    ((0.1+0.1j, 0.2+0.2j, 0.5+0.5j, 0.6+0.6j, np.nan, np.nan),
    (1.6*(0.1+0.1j), 1.7*(0.2+0.2j), 2.0*(0.5+0.5j), 20.1*(0.6+0.6j), np.nan, np.nan)), dtype=np.complex128
)
y_lower_staticliq = np.asarray(
    ((0.5+0.5j, 0.6-9.6j, np.nan, np.nan, np.nan, np.nan),), dtype=np.complex128
)

MAX_NUM_Ys = 6

@pytest.mark.parametrize('lower_layer_type', (0, 1))
@pytest.mark.parametrize('lower_is_static', (True, False))
@pytest.mark.parametrize('upper_layer_type', (0, 1))
@pytest.mark.parametrize('upper_is_static', (True, False))
@pytest.mark.parametrize('lower_is_incompressible', (True, False))
@pytest.mark.parametrize('upper_is_incompressible', (True, False))
def test_interface_driver(lower_layer_type, lower_is_static, upper_layer_type, upper_is_static,
                          lower_is_incompressible, upper_is_incompressible):

    # TODO: Currently the RadialSolver.interfaces does not properly use incompressibility so there are no tests.
    if lower_is_incompressible or upper_is_incompressible:
        pytest.skip('Currently the RadialSolver.interfaces does not properly use incompressibility so there are no tests.')

    num_sols_lower = find_num_shooting_solutions(lower_layer_type, lower_is_static, lower_is_incompressible)
    num_ys_lower   = num_sols_lower * 2
    num_sols_upper = find_num_shooting_solutions(upper_layer_type, upper_is_static, upper_is_incompressible)
    num_ys_upper   = num_sols_upper * 2

    # Determine which lower layer y solution to use
    if (lower_layer_type == 0):
        lower_y = y_lower_solid
    elif lower_is_static:
        lower_y = y_lower_staticliq
    else:
        lower_y = y_lower_liquid
    
    # Create array full of nans for the "output" array. If things worked okay then they should no longer be nan.
    upper_y = np.nan * np.ones((3, 6), dtype=np.complex128, order='C')

    solve_upper_y_at_interface(
        lower_y, upper_y,
        lower_layer_type, lower_is_static, lower_is_incompressible,
        upper_layer_type, upper_is_static, upper_is_incompressible,
        interface_gravity, static_liquid_density, G_to_use
        )

    # Make sure the correct number of array elements were set.
    assert np.sum(~np.isnan(upper_y)) == (num_sols_upper * num_ys_upper)


@pytest.mark.parametrize('lower_layer_type', (0, 1))
@pytest.mark.parametrize('lower_is_static', (True, False))
@pytest.mark.parametrize('upper_layer_type', (0, 1))
@pytest.mark.parametrize('upper_is_static', (True, False))
@pytest.mark.parametrize('lower_is_incompressible', (True, False))
@pytest.mark.parametrize('upper_is_incompressible', (True, False))
def test_interface_accuracy(lower_layer_type, lower_is_static, upper_layer_type, upper_is_static,
                            lower_is_incompressible, upper_is_incompressible):
    """ Use results from TidalPy 0.4.0 to test the accuracy of the current interface solver. """

    # TODO: Currently the RadialSolver.interfaces does not properly use incompressibility so there are no tests.
    if lower_is_incompressible or upper_is_incompressible:
        pytest.skip('Currently the RadialSolver.interfaces does not properly use incompressibility so there are no tests.')

    lower_layer_is_solid = lower_layer_type == 0
    upper_layer_is_solid = upper_layer_type == 0
    if (lower_layer_is_solid, lower_is_static, upper_layer_is_solid, upper_is_static) not in tpy_0p4_results:
        pytest.skip(f'Combination {(lower_layer_type, lower_is_static, upper_layer_type, upper_is_static)} not found (or not implemented) in pre-calculated TidalPy v0.4 results.')

    else:
        num_sols_lower = find_num_shooting_solutions(lower_layer_type, lower_is_static, lower_is_incompressible)
        num_ys_lower   = num_sols_lower * 2
        num_sols_upper = find_num_shooting_solutions(upper_layer_type, upper_is_static, upper_is_incompressible)
        num_ys_upper   = num_sols_upper * 2

        # Determine which lower layer y solution to use
        if (lower_layer_type == 0):
            lower_y = y_lower_solid
        elif lower_is_static:
            lower_y = y_lower_staticliq
        else:
            lower_y = y_lower_liquid
        
        # Create array full of nans for the "output" array. If things worked okay then they should no longer be nan.
        upper_y = np.nan * np.ones((3, 6), dtype=np.complex128, order='C')
        
        solve_upper_y_at_interface(
            lower_y, upper_y,
            lower_layer_type, lower_is_static, lower_is_incompressible,
            upper_layer_type, upper_is_static, upper_is_incompressible,
            interface_gravity, static_liquid_density, G_to_use
            )

        # Compare results to TidalPy 0.4.0 pre-calculated interface solution.
        comparison_results = tpy_0p4_results[(lower_layer_is_solid, lower_is_static, upper_layer_is_solid, upper_is_static)]
        comparison_results_no_nan = comparison_results[~np.isnan(comparison_results)]
        upper_y_no_nan = upper_y[~np.isnan(upper_y)]
        np.allclose(upper_y_no_nan.flatten(), comparison_results_no_nan)
