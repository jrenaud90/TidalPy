"""Interface solver: filled-value counts and accuracy against TidalPy 0.4 results (compressible layers)."""
import numpy as np
import pytest

from TidalPy.RadialSolver.derivatives.odes import find_num_shooting_solutions
from TidalPy.RadialSolver.interfaces.interfaces import solve_upper_y_at_interface

# Keyed by (lower is solid, lower is static, upper is solid, upper is static).
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

STATIC_LIQUID_DENSITY = 7600.
INTERFACE_GRAVITY = 2.7
G_TO_USE = 6.67430e-11

Y_LOWER_SOLID = np.asarray(
    ((0.1+0.1j, 0.2+0.2j, 0.3+0.3j, 0.4+0.4j, 0.5+0.5j, 0.6+0.6j),
    (-0.1-0.1j, -0.2-0.2j, -0.3-0.3j, -0.4-0.4j, -0.5-0.5j, -0.6-0.6j),
    (1.6*(0.1+0.1j), 1.7*(0.2+0.2j), 1.8*(0.3+0.3j), 1.9*(0.4+0.4j), 2.0*(0.5+0.5j), 20.1*(0.6+0.6j))),
    dtype=np.complex128
)
Y_LOWER_LIQUID = np.asarray(
    ((0.1+0.1j, 0.2+0.2j, 0.5+0.5j, 0.6+0.6j, np.nan, np.nan),
    (1.6*(0.1+0.1j), 1.7*(0.2+0.2j), 2.0*(0.5+0.5j), 20.1*(0.6+0.6j), np.nan, np.nan)), dtype=np.complex128
)
Y_LOWER_STATIC_LIQUID = np.asarray(((0.5+0.5j, 0.6-9.6j, np.nan, np.nan, np.nan, np.nan),), dtype=np.complex128)


def _solve_interface(lower_layer_type, lower_is_static, upper_layer_type, upper_is_static):
    """Upper-layer y values (NaN where unused) for compressible layers."""
    if lower_layer_type == 0:
        lower_y = Y_LOWER_SOLID
    elif lower_is_static:
        lower_y = Y_LOWER_STATIC_LIQUID
    else:
        lower_y = Y_LOWER_LIQUID

    upper_y = np.nan * np.ones((3, 6), dtype=np.complex128, order='C')
    solve_upper_y_at_interface(
        lower_y,
        upper_y,
        lower_layer_type,
        lower_is_static,
        False,
        upper_layer_type,
        upper_is_static,
        False,
        INTERFACE_GRAVITY,
        STATIC_LIQUID_DENSITY,
        G_TO_USE
    )
    return upper_y


@pytest.mark.parametrize('lower_layer_type', (0, 1))
@pytest.mark.parametrize('lower_is_static', (True, False))
@pytest.mark.parametrize('upper_layer_type', (0, 1))
@pytest.mark.parametrize('upper_is_static', (True, False))
@pytest.mark.parametrize('lower_is_incompressible', (True, False))
@pytest.mark.parametrize('upper_is_incompressible', (True, False))
def test_interface_driver(
        lower_layer_type,
        lower_is_static,
        upper_layer_type,
        upper_is_static,
        lower_is_incompressible,
        upper_is_incompressible):
    """The interface solver fills exactly num_solutions * num_ys upper-layer values."""
    if lower_is_incompressible or upper_is_incompressible:
        pytest.skip('Incompressible interface tests not yet implemented.')

    num_sols_upper = find_num_shooting_solutions(upper_layer_type, upper_is_static, upper_is_incompressible)
    upper_y = _solve_interface(lower_layer_type, lower_is_static, upper_layer_type, upper_is_static)
    assert np.sum(~np.isnan(upper_y)) == (num_sols_upper * num_sols_upper * 2)


@pytest.mark.parametrize('lower_layer_type', (0, 1))
@pytest.mark.parametrize('lower_is_static', (True, False))
@pytest.mark.parametrize('upper_layer_type', (0, 1))
@pytest.mark.parametrize('upper_is_static', (True, False))
@pytest.mark.parametrize('lower_is_incompressible', (True, False))
@pytest.mark.parametrize('upper_is_incompressible', (True, False))
def test_interface_accuracy(
        lower_layer_type,
        lower_is_static,
        upper_layer_type,
        upper_is_static,
        lower_is_incompressible,
        upper_is_incompressible):
    """The finite upper-layer values match TidalPy 0.4's."""
    if lower_is_incompressible or upper_is_incompressible:
        pytest.skip('Incompressible interface tests not yet implemented.')

    key = (lower_layer_type == 0, lower_is_static, upper_layer_type == 0, upper_is_static)
    if key not in tpy_0p4_results:
        pytest.skip('Combination not found in pre-calculated TidalPy v0.4 results.')

    upper_y = _solve_interface(lower_layer_type, lower_is_static, upper_layer_type, upper_is_static)
    comparison_results = tpy_0p4_results[key]
    assert np.allclose(upper_y[~np.isnan(upper_y)].flatten(), comparison_results[~np.isnan(comparison_results)])
