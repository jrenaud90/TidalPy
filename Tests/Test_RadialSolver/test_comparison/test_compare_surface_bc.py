"""Surface boundary conditions against the frozen classic solver's."""
from pathlib import Path

import numpy as np
import pytest

from TidalPy.RadialSolver.boundaries.surface_bc import get_surface_bc

FROZEN_PATH = Path(__file__).parent / "frozen" / "test_compare_surface_bc.npz"
with np.load(FROZEN_PATH, allow_pickle=False) as frozen_file:
    CLASSIC_SURFACE_BC = {key: frozen_file[key] for key in frozen_file.files}


@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('bc_models', ((0,), (1,), (2,), (0, 1), (1, 2), (0, 1, 2)))
def test_compare_surface_bc(degree_l, bc_models):
    """Surface boundary conditions match the classic ones to roundoff."""
    old_result = CLASSIC_SURFACE_BC[f"degree_l_{degree_l}__bc_models_" + "+".join(str(model) for model in bc_models)]
    new_result = get_surface_bc(np.asarray(bc_models, dtype=np.intc), 1000., 2000., degree_l)
    np.testing.assert_allclose(new_result, old_result, rtol=1e-15)
