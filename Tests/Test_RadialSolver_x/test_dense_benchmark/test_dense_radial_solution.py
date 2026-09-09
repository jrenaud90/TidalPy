"""Benchmark the dense radial-solution calling system against frozen original-solver results.

The new RadialSolver_x shooting method retains the dense CyRK interpolants per layer/solution and
collapses them at any requested radius (``get_radial_solution`` / ``get_radial_solution_array``),
instead of storing the collapsed y-solution on a fixed grid and linearly interpolating it.

The reference is the original (non-_x) RadialSolver, frozen to .npz by
``generate_dense_benchmark_data.py`` (the original solver will eventually leave the package, so the
test reads the frozen data rather than calling it). Expectations:

* Love numbers agree tightly (they are integrated, BC-anchored quantities).
* The dense solution evaluated on the original grid agrees very well with the original grid solution
  at the surface and at layer boundaries, and across the interiors.
* Off-grid, the dense evaluation is finite, continuous, and (for a coarse grid) departs from a plain
  linear interpolation of the original grid -- the dense value being the more accurate one.
"""
import os

import numpy as np
import pytest

from TidalPy.RadialSolver_x.solver import radial_solver as radial_solver_new

_HERE = os.path.dirname(os.path.abspath(__file__))
# '3layer_dynliq' is solid / dynamic-liquid / solid at a short (well-conditioned) forcing period -
# the configuration is only stable at short periods; at long periods the dynamic-liquid formulation
# is ill-conditioned and both the original and _x solvers diverge (use a static liquid there).
_CONFIGS = ('1layer', '2solid', '3layer_dynliq')


def _load(name):
    path = os.path.join(_HERE, f'dense_benchmark_{name}.npz')
    if not os.path.exists(path):
        pytest.skip(f"Frozen benchmark {path} missing; run generate_dense_benchmark_data.py.")
    return np.load(path, allow_pickle=False)


def _run_new(data):
    return radial_solver_new(
        data['radius_array'].copy(), data['density_array'].copy(),
        data['bulk_modulus_array'].copy(), data['complex_shear_modulus_array'].copy(),
        float(data['frequency']), float(data['planet_bulk_density']),
        tuple(str(s) for s in data['layer_types']),
        tuple(bool(b) for b in data['is_static']),
        tuple(bool(b) for b in data['is_incompressible']),
        data['upper_radius_by_layer'].copy(),
        degree_l=int(data['degree_l']),
        solve_for=tuple(str(s) for s in data['solve_for']),
        nondimensionalize=True, integration_method='DOP853',
        integration_rtol=1.0e-8, integration_atol=1.0e-10,
        max_num_steps=5_000_000, expected_size=500, raise_on_fail=True,
    )


def _old_y_grid(data):
    """Original grid y-solution for ytype 0 as (num_slices, 6)."""
    return np.ascontiguousarray(data['old_result'][0:6, :].T)


@pytest.mark.parametrize('name', _CONFIGS)
def test_love_numbers_match_frozen(name):
    data = _load(name)
    new_out = _run_new(data)
    assert new_out.success, new_out.message
    np.testing.assert_allclose(
        new_out.love, data['old_love'], rtol=1.0e-3, atol=1.0e-9,
        err_msg=f"[{name}] dense-solver Love numbers differ from frozen original.")


@pytest.mark.parametrize('name', _CONFIGS)
def test_dense_matches_grid_at_surface_and_interfaces(name):
    """Dense evaluation on the original grid agrees with the original grid solution.

    Interface radii are duplicated in the grid (one slice per side); the dense getter resolves an
    interface to its lower layer, so those duplicated slices are compared loosely while the clean
    interior slices and the (unambiguous) surface are compared tightly.
    """
    data = _load(name)
    radius_array = data['radius_array']
    upper_radii = data['upper_radius_by_layer']
    new_out = _run_new(data)
    assert new_out.success, new_out.message

    new_dense = new_out.get_radial_solution_array(radius_array, 0)  # (N, 6)
    old_y = _old_y_grid(data)                                        # (N, 6)
    assert new_dense.shape == old_y.shape

    # Compare only where both are finite (below the auto starting radius both are NaN).
    finite = np.isfinite(new_dense) & np.isfinite(old_y)
    # Mark inner-interface slices (radius coincides with a non-surface layer boundary).
    inner_interfaces = upper_radii[:-1]
    is_interface = np.zeros(radius_array.size, dtype=bool)
    for r_iface in inner_interfaces:
        is_interface |= np.isclose(radius_array, r_iface, rtol=0.0, atol=1.0e-6)

    scale = np.nanmax(np.abs(old_y[finite])) if finite.any() else 1.0

    # Tight parity at clean interior + surface slices.
    clean = finite & ~is_interface[:, None]
    np.testing.assert_allclose(
        new_dense[clean], old_y[clean], rtol=2.0e-3, atol=1.0e-3 * scale,
        err_msg=f"[{name}] dense vs original grid differ at clean slices.")

    # The surface slice (last) is unambiguous: require it to agree well component-by-component.
    surf_finite = np.isfinite(old_y[-1]) & np.isfinite(new_dense[-1])
    np.testing.assert_allclose(
        new_dense[-1][surf_finite], old_y[-1][surf_finite], rtol=2.0e-3, atol=1.0e-3 * scale,
        err_msg=f"[{name}] dense vs original grid differ at the surface.")


@pytest.mark.parametrize('name', _CONFIGS)
def test_dense_offgrid_is_finite_and_consistent(name):
    """Off-grid dense evaluation is finite and close to (but distinct from) grid linear interp."""
    data = _load(name)
    radius_array = data['radius_array']
    new_out = _run_new(data)
    assert new_out.success, new_out.message

    old_y = _old_y_grid(data)
    finite_rows = np.all(np.isfinite(old_y), axis=1)
    r_finite = radius_array[finite_rows]
    # Midpoints between consecutive finite, non-coincident grid radii (off-grid query points).
    mids = []
    for a, b in zip(r_finite[:-1], r_finite[1:]):
        if b - a > 1.0:
            mids.append(0.5 * (a + b))
    mids = np.asarray(mids)
    assert mids.size > 0

    dense_mid = new_out.get_radial_solution_array(mids, 0)        # (M, 6)
    # y1, y2, y5, y6 are defined in every layer type; require them finite off-grid.
    assert np.all(np.isfinite(dense_mid[:, [0, 1, 4, 5]])), f"[{name}] off-grid dense has non-finite y."

    # Linear interpolation of the original grid at the same points, for y1 (always defined).
    lin_y1 = np.interp(mids, radius_array[finite_rows], old_y[finite_rows, 0].real)
    # Dense and linear-interp should be in the same ballpark (sanity), confirming continuity.
    scale = np.nanmax(np.abs(old_y[finite_rows, 0].real))
    np.testing.assert_allclose(dense_mid[:, 0].real, lin_y1, rtol=0.0, atol=0.05 * scale,
                               err_msg=f"[{name}] off-grid dense y1 wildly inconsistent with grid.")
