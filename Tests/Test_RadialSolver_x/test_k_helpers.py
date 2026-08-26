"""Tests for the `_x` radial solver convenience helpers (TidalPy.RadialSolver_x.helpers)."""
import numpy as np

from TidalPy.rheology_x import Maxwell
from TidalPy.RadialSolver_x import radial_solver, homogeneous_love_numbers

_R = 1.8216e6
_RHO = 3529.0
_FREQ = 4.1e-5
_SLICES = 60
_MU = Maxwell().calc_complex_modulus(60.0e9, 1.0e15, _FREQ)


def test_homogeneous_helper_matches_direct_call():
    """The helper reproduces a direct radial_solver call on the same hand-built arrays."""
    solution = homogeneous_love_numbers(_R, _RHO, _MU, _FREQ, num_slices=_SLICES)
    assert solution.success

    radius_array = np.linspace(0.0, _R, _SLICES)
    direct = radial_solver(
        radius_array, np.full(_SLICES, _RHO), np.full(_SLICES, 200.0e9 + 0j),
        np.full(_SLICES, _MU, dtype=np.complex128),
        _FREQ, _RHO, ("solid",), (True,), (False,), np.asarray((_R,)),
        degree_l=2, solve_for=("tidal",))
    assert direct.success
    np.testing.assert_allclose(np.atleast_1d(solution.k), np.atleast_1d(direct.k), rtol=1e-14)
    np.testing.assert_allclose(np.atleast_1d(solution.h), np.atleast_1d(direct.h), rtol=1e-14)


def test_homogeneous_helper_kwarg_passthrough():
    """Extra keyword arguments reach radial_solver (here: the propagation matrix path)."""
    solution = homogeneous_love_numbers(
        _R, _RHO, 60.0e9 + 0.0j, _FREQ, num_slices=_SLICES,
        layer_is_incompressible=True, use_prop_matrix=True)
    assert solution.success
    k2 = complex(np.atleast_1d(solution.k)[0])
    # Static incompressible homogeneous sphere: 0 < k2 < 1.5 with no imaginary part (real modulus).
    assert 0.0 < k2.real < 1.5
    assert k2.imag == 0.0


def test_homogeneous_helper_degree_l():
    """Higher harmonic degrees produce smaller Love numbers, as expected."""
    k2 = complex(np.atleast_1d(homogeneous_love_numbers(_R, _RHO, _MU, _FREQ).k)[0])
    k3 = complex(np.atleast_1d(homogeneous_love_numbers(_R, _RHO, _MU, _FREQ, degree_l=3).k)[0])
    assert abs(k3) < abs(k2)
