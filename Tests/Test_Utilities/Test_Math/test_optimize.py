import pytest
import math

from TidalPy.utilities.math.optimize import test_bracket as test_bracket_func

def func(x, *args):
    return x**5 - x**3


@pytest.mark.parametrize('x0', (-1.1, -0.9, -0.1, 0.1, 0.9, 1.1, -10, 10)) # Test functions roots are -1, 0, 1. Start on either side of each.
@pytest.mark.parametrize('ratio', (0.5, 1.618))
@pytest.mark.parametrize('dx', (0.01, 0.1, 1.0))
@pytest.mark.parametrize('maxiter', (100, 1))
def test_bracket(x0, dx, ratio, maxiter):

    try:
        import burnman as bm
    except ImportError:
        pytest.skip("BurnMan not installed, can not test math.optimize.bracket.")

    xa, xb, fa, fb, iters_p1, iters_p2, error_code, message = \
        test_bracket_func(x0, dx, ratio=ratio, maxiter=maxiter)
    
    continue_test = True
    if ratio < 1.0:
        # Ratio must be greater than 1
        continue_test = False
        assert error_code == -1
    
    if continue_test:
        tpy_result = test_bracket_func(x0, dx, ratio=ratio, maxiter=maxiter)
        try:
            bm_result = bm.utils.math.bracket(func, x0, dx, ratio=ratio, maxiter=maxiter)
        except ValueError:
            # BurnMan could not find a zero. Let's make sure TidalPy didn't either.
            assert error_code != 0
        else:
            for tpy_subresult, bm_subresult in zip(tpy_result, bm_result):
                assert math.isclose(tpy_subresult, bm_subresult, rel_tol=1e-10, abs_tol=1e-10)
    