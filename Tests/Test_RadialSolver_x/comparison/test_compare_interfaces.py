"""Compare RadialSolver_x interface solver against original RadialSolver."""
import pytest
import numpy as np

from TidalPy.RadialSolver.interfaces.interfaces import solve_upper_y_at_interface as solve_old
from TidalPy.RadialSolver_x.interfaces.interfaces import solve_upper_y_at_interface as solve_new

static_liquid_density = 7600.
interface_gravity = 2.7
G_to_use = 6.67430e-11

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


@pytest.mark.parametrize('lower_layer_type', (0, 1))
@pytest.mark.parametrize('lower_is_static', (True, False))
@pytest.mark.parametrize('upper_layer_type', (0, 1))
@pytest.mark.parametrize('upper_is_static', (True, False))
def test_compare_interfaces(lower_layer_type, lower_is_static, upper_layer_type, upper_is_static):
    """Interface conditions should match between old and new implementations."""

    if (lower_layer_type == 0):
        lower_y = y_lower_solid
    elif lower_is_static:
        lower_y = y_lower_staticliq
    else:
        lower_y = y_lower_liquid

    upper_y_old = np.nan * np.ones((3, 6), dtype=np.complex128, order='C')
    upper_y_new = np.nan * np.ones((3, 6), dtype=np.complex128, order='C')

    lower_is_incompressible = False
    upper_is_incompressible = False

    solve_old(
        lower_y, upper_y_old,
        lower_layer_type, lower_is_static, lower_is_incompressible,
        upper_layer_type, upper_is_static, upper_is_incompressible,
        interface_gravity, static_liquid_density, G_to_use
    )
    solve_new(
        lower_y, upper_y_new,
        lower_layer_type, lower_is_static, lower_is_incompressible,
        upper_layer_type, upper_is_static, upper_is_incompressible,
        interface_gravity, static_liquid_density, G_to_use
    )

    # Compare non-NaN values.
    old_no_nan = upper_y_old[~np.isnan(upper_y_old)]
    new_no_nan = upper_y_new[~np.isnan(upper_y_new)]
    assert len(old_no_nan) == len(new_no_nan)
    np.testing.assert_allclose(new_no_nan, old_no_nan, rtol=1e-12)
