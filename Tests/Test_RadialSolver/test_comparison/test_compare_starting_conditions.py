"""Shooting-solution counts and starting conditions against the frozen classic solver's."""
from pathlib import Path

import numpy as np
import pytest

from TidalPy.Rheology import Maxwell
from TidalPy.RadialSolver.derivatives.odes import find_num_shooting_solutions
from TidalPy.RadialSolver.starting.driver import find_starting_conditions

FROZEN_PATH = Path(__file__).parent / "frozen" / "test_compare_starting_conditions.npz"
with np.load(FROZEN_PATH, allow_pickle=False) as frozen_file:
    CLASSIC_REFERENCE = {key: frozen_file[key] for key in frozen_file.files}

FREQUENCY = 0.1
RADIUS = 0.1
DENSITY = 7000.
BULK_MODULUS = 100.0e9
COMPLEX_SHEAR = Maxwell().calc_complex_modulus(50.0e9, 1.0e20, FREQUENCY)
G_TO_USE = 6.67430e-11


def _layer_key(layer_type, is_static, is_incompressible):
    return f"layer_type_{layer_type}__static_{is_static}__incompressible_{is_incompressible}"


@pytest.mark.parametrize('layer_type', (0, 1))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
def test_compare_num_shooting_solutions(layer_type, is_static, is_incompressible):
    """The number of shooting solutions matches the classic count."""
    old_num = int(CLASSIC_REFERENCE["num_solutions__" + _layer_key(layer_type, is_static, is_incompressible)])
    assert old_num == find_num_shooting_solutions(layer_type, is_static, is_incompressible)


@pytest.mark.parametrize('layer_type', (0, 1))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('use_kamata', (True, False))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_compare_starting_conditions(
        layer_type,
        is_static,
        is_incompressible,
        use_kamata,
        degree_l):
    """Starting conditions match the classic ones."""
    is_solid = layer_type == 0
    if ((not use_kamata) and is_incompressible and (is_solid or ((not is_solid) and (not is_static)))) \
            or (use_kamata and is_static and is_incompressible and is_solid):
        pytest.skip('Combination not yet implemented.')

    # Copy the frozen array: the Kamata branch below rewrites its first row.
    old_arr = np.array(CLASSIC_REFERENCE[
        f"starting__{_layer_key(layer_type, is_static, is_incompressible)}__kamata_{use_kamata}__degree_l_{degree_l}"])
    new_arr = np.nan * np.ones(old_arr.shape, dtype=np.complex128, order='C')

    find_starting_conditions(
        layer_type,
        is_static,
        is_incompressible,
        use_kamata,
        FREQUENCY,
        RADIUS,
        DENSITY,
        BULK_MODULUS,
        COMPLEX_SHEAR,
        degree_l,
        G_TO_USE,
        new_arr
    )

    if use_kamata and is_solid and (not is_static) and is_incompressible:
        # The new first solution is (s1 - s2) gamma / omega^2 of the classic basis, which stays independent as
        # omega -> 0; compare against that same combination of the classic solutions.
        gamma = 4.0 * np.pi * G_TO_USE * DENSITY / 3.0
        old_arr[0] = (old_arr[0] - old_arr[1]) * gamma / FREQUENCY ** 2

    np.testing.assert_allclose(new_arr, old_arr, rtol=1e-10)
