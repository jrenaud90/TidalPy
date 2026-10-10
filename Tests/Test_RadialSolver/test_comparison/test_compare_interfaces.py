"""Interface solutions against the frozen classic solver's."""
from pathlib import Path

import numpy as np
import pytest

from TidalPy.RadialSolver.interfaces.interfaces import solve_upper_y_at_interface

FROZEN_PATH = Path(__file__).parent / "frozen" / "test_compare_interfaces.npz"
with np.load(FROZEN_PATH, allow_pickle=False) as frozen_file:
    CLASSIC_UPPER_Y = {key: frozen_file[key] for key in frozen_file.files}

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


@pytest.mark.parametrize('lower_layer_type', (0, 1))
@pytest.mark.parametrize('lower_is_static', (True, False))
@pytest.mark.parametrize('upper_layer_type', (0, 1))
@pytest.mark.parametrize('upper_is_static', (True, False))
def test_compare_interfaces(lower_layer_type, lower_is_static, upper_layer_type, upper_is_static):
    """The finite upper-layer values match the classic ones for every compressible layer pairing."""
    if lower_layer_type == 0:
        lower_y = Y_LOWER_SOLID
    elif lower_is_static:
        lower_y = Y_LOWER_STATIC_LIQUID
    else:
        lower_y = Y_LOWER_LIQUID

    upper_y_old = CLASSIC_UPPER_Y[
        f"lower_type_{lower_layer_type}__lower_static_{lower_is_static}"
        f"__upper_type_{upper_layer_type}__upper_static_{upper_is_static}"]
    upper_y_new = np.nan * np.ones((3, 6), dtype=np.complex128, order='C')

    solve_upper_y_at_interface(
        lower_y,
        upper_y_new,
        lower_layer_type,
        lower_is_static,
        upper_layer_type,
        upper_is_static,
        INTERFACE_GRAVITY,
        STATIC_LIQUID_DENSITY,
        G_TO_USE
    )

    old_no_nan = upper_y_old[~np.isnan(upper_y_old)]
    new_no_nan = upper_y_new[~np.isnan(upper_y_new)]
    assert len(old_no_nan) == len(new_no_nan)
    # Y_LOWER_SOLID's second solution is minus its first, so the solid to static-liquid solution is exactly zero.
    # Fused multiply-adds (clang's default on arm64) leave roundoff there on the scale of the lower-layer values.
    atol = 1e-14 * np.nanmax(np.abs(lower_y))
    np.testing.assert_allclose(new_no_nan, old_no_nan, rtol=1e-12, atol=atol)
