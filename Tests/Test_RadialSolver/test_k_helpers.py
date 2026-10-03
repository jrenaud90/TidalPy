"""The homogeneous_love_numbers convenience helper."""
import numpy as np

from TidalPy.Rheology import Maxwell
from TidalPy.RadialSolver import radial_solver, homogeneous_love_numbers

RADIUS = 1.8216e6
DENSITY = 3529.0
FREQUENCY = 4.1e-5
NUM_SLICES = 60
COMPLEX_SHEAR = Maxwell().calc_complex_modulus(60.0e9, 1.0e15, FREQUENCY)


def _helper_k(**kwargs):
    """The helper's first k on the Maxwell body."""
    return complex(np.atleast_1d(homogeneous_love_numbers(RADIUS, DENSITY, COMPLEX_SHEAR, FREQUENCY, **kwargs).k)[0])


def test_homogeneous_helper_matches_direct_call():
    """The helper reproduces a direct radial_solver call on the same hand-built arrays."""
    solution = homogeneous_love_numbers(
        RADIUS,
        DENSITY,
        COMPLEX_SHEAR,
        FREQUENCY,
        num_slices=NUM_SLICES)
    assert solution.success

    direct = radial_solver(
        np.linspace(0.0, RADIUS, NUM_SLICES),
        np.full(NUM_SLICES, DENSITY),
        np.full(NUM_SLICES, 200.0e9 + 0j),
        np.full(NUM_SLICES, COMPLEX_SHEAR, dtype=np.complex128),
        FREQUENCY,
        DENSITY,
        ("solid",),
        (True,),
        (False,),
        np.asarray((RADIUS,)),
        degree_l=2,
        solve_for=("tidal",))
    assert direct.success
    np.testing.assert_allclose(np.atleast_1d(solution.k), np.atleast_1d(direct.k), rtol=1e-14)
    np.testing.assert_allclose(np.atleast_1d(solution.h), np.atleast_1d(direct.h), rtol=1e-14)


def test_homogeneous_helper_kwarg_passthrough():
    """Extra keyword arguments reach radial_solver (here: the propagation matrix path)."""
    solution = homogeneous_love_numbers(
        RADIUS,
        DENSITY,
        60.0e9 + 0.0j,
        FREQUENCY,
        num_slices=NUM_SLICES,
        layer_is_incompressible=True,
        love_method='propagation_matrix')
    assert solution.success
    k2 = complex(np.atleast_1d(solution.k)[0])
    # A real modulus gives a purely real k2.
    assert 0.0 < k2.real < 1.5
    assert k2.imag == 0.0


def test_homogeneous_helper_degree_l():
    """Degree 3 gives a smaller Love number than degree 2."""
    assert abs(_helper_k(degree_l=3)) < abs(_helper_k())
