"""The dense radial-solution interpolants against frozen TidalPy 0.7 grid solutions and Love numbers."""
import os

import numpy as np
import pytest

from TidalPy.RadialSolver.solver import radial_solver

_HERE = os.path.dirname(os.path.abspath(__file__))
# The dynamic-liquid case uses a short forcing period: at long periods that formulation is ill-conditioned and both
# solvers diverge.
_CONFIGS = ('1layer', '2solid', '3layer_dynliq')


def _load(name):
    path = os.path.join(_HERE, f'dense_benchmark_{name}.npz')
    if not os.path.exists(path):
        pytest.skip(f"Frozen benchmark {path} missing; run generate_dense_benchmark_data.py.")
    return np.load(path, allow_pickle=False)


def _run_new(data):
    return radial_solver(
        data['radius_array'].copy(),
        data['density_array'].copy(),
        data['bulk_modulus_array'].copy(),
        data['complex_shear_modulus_array'].copy(),
        float(data['frequency']),
        float(data['planet_bulk_density']),
        tuple(str(layer_type) for layer_type in data['layer_types']),
        tuple(bool(flag) for flag in data['is_static']),
        tuple(bool(flag) for flag in data['is_incompressible']),
        data['upper_radius_by_layer'].copy(),
        degree_l=int(data['degree_l']),
        solve_for=tuple(str(solve_type) for solve_type in data['solve_for']),
        nondimensionalize=True,
        integration_method='DOP853',
        integration_rtol=1.0e-8,
        integration_atol=1.0e-10,
        max_num_steps=5_000_000,
        expected_size=500,
        raise_on_fail=True,
    )


def _old_y_grid(data):
    """Original grid y-solution for ytype 0 as (num_slices, 6)."""
    return np.ascontiguousarray(data['old_result'][0:6, :].T)


@pytest.mark.parametrize('name', _CONFIGS)
def test_love_numbers_match_frozen(name):
    """Love numbers match the frozen original to 1e-3."""
    data = _load(name)
    new_out = _run_new(data)
    assert new_out.success, new_out.message
    np.testing.assert_allclose(
        new_out.love,
        data['old_love'],
        rtol=1.0e-3,
        atol=1.0e-9,
        err_msg=f"[{name}] dense-solver Love numbers differ from frozen original.")


@pytest.mark.parametrize('name', _CONFIGS)
def test_dense_matches_grid_at_surface_and_interfaces(name):
    """The dense solution on the original grid matches it at clean interior slices and at the surface."""
    data = _load(name)
    radius_array = data['radius_array']
    new_out = _run_new(data)
    assert new_out.success, new_out.message

    new_dense = new_out.get_radial_solution_array(radius_array, 0)
    old_y = _old_y_grid(data)
    assert new_dense.shape == old_y.shape

    # Both are NaN below the automatic starting radius.
    finite = np.isfinite(new_dense) & np.isfinite(old_y)
    # The dense getter resolves a duplicated interface radius to its lower layer, so interface slices are skipped.
    is_interface = np.zeros(radius_array.size, dtype=bool)
    for interface_radius in data['upper_radius_by_layer'][:-1]:
        is_interface |= np.isclose(radius_array, interface_radius, rtol=0.0, atol=1.0e-6)

    scale = np.nanmax(np.abs(old_y[finite])) if finite.any() else 1.0

    clean = finite & ~is_interface[:, None]
    np.testing.assert_allclose(
        new_dense[clean],
        old_y[clean],
        rtol=2.0e-3,
        atol=1.0e-3 * scale,
        err_msg=f"[{name}] dense vs original grid differ at clean slices.")

    surface_finite = np.isfinite(old_y[-1]) & np.isfinite(new_dense[-1])
    np.testing.assert_allclose(
        new_dense[-1][surface_finite],
        old_y[-1][surface_finite],
        rtol=2.0e-3,
        atol=1.0e-3 * scale,
        err_msg=f"[{name}] dense vs original grid differ at the surface.")


@pytest.mark.parametrize('name', _CONFIGS)
def test_dense_offgrid_is_finite_and_consistent(name):
    """Off-grid dense values are finite and within 5% (of the y1 scale) of a linear interpolation of the grid."""
    data = _load(name)
    radius_array = data['radius_array']
    new_out = _run_new(data)
    assert new_out.success, new_out.message

    old_y = _old_y_grid(data)
    finite_rows = np.all(np.isfinite(old_y), axis=1)
    radius_finite = radius_array[finite_rows]
    midpoints = np.asarray([0.5 * (lower + upper) for lower, upper in zip(radius_finite[:-1], radius_finite[1:])
                            if upper - lower > 1.0])
    assert midpoints.size > 0

    dense_mid = new_out.get_radial_solution_array(midpoints, 0)
    # y1, y2, y5, and y6 are defined in every layer type.
    assert np.all(np.isfinite(dense_mid[:, [0, 1, 4, 5]])), f"[{name}] off-grid dense has non-finite y."

    linear_y1 = np.interp(midpoints, radius_array[finite_rows], old_y[finite_rows, 0].real)
    scale = np.nanmax(np.abs(old_y[finite_rows, 0].real))
    np.testing.assert_allclose(
        dense_mid[:, 0].real,
        linear_y1,
        rtol=0.0,
        atol=0.05 * scale,
        err_msg=f"[{name}] off-grid dense y1 wildly inconsistent with grid.")
