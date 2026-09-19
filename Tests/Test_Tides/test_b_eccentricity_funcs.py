from math import isclose

import pytest
import numpy as np

from TidalPy.Tides_x.eccentricity import eccentricity_func
from TidalPy.Utilities_x.lookups import IntMap3


@pytest.mark.parametrize('degree_l', (2, 3, 4, 5, 6, 7, 8, 9, 10))
@pytest.mark.parametrize('truncation', (1, 2, 3, 4, 5, 10, 0))
def test_eccentricity_funcs(degree_l, truncation):
    """Tests TidalPy's eccentricity functions for various degree_ls and various truncation levels."""

    eccentricity = 0.5

    if degree_l == 1:
        # Degree l is currently not supported, check that it raises an error.
        with pytest.raises(NotImplementedError):
            ecc_results_by_lpq, ecc_results_by_lp = eccentricity_func(eccentricity, degree_l, truncation)
    elif truncation == 0:
        # Truncation = 5 is not supported, check that it raises an error.
        with pytest.raises(NotImplementedError):
            ecc_results_by_lpq, ecc_results_by_lp = eccentricity_func(eccentricity, degree_l, truncation)
    else:
        ecc_results_by_lpq, ecc_results_by_lp = eccentricity_func(eccentricity, degree_l, truncation)

        # Check return types
        assert isinstance(ecc_results_by_lpq, IntMap3)
        assert isinstance(ecc_results_by_lp, dict)

        # No matter the assumptions, these should have some size to them.
        assert len(ecc_results_by_lpq) > 0
        for (l, p, q), ecc_result in ecc_results_by_lpq:
            assert isinstance(l, int)
            assert isinstance(p, int)
            assert isinstance(q, int)
            assert l == degree_l
            assert p <= l
            assert isinstance(ecc_result, float)
        
        for (l, m), ecc_results_by_p in ecc_results_by_lp.items():
            assert isinstance(l, int)
            assert isinstance(m, int)
            assert l == degree_l
            assert m <= l
            for (q,), ecc_result in ecc_results_by_p:
                assert isinstance(q, int)
                assert isinstance(ecc_result, float)
        
        # Spot checks
        if degree_l == 2:
            if truncation == 1:
                assert len(ecc_results_by_lpq) == 9
                assert isclose(ecc_results_by_lpq[(2, 0, 0)], 1.0)
                assert isclose(ecc_results_by_lpq[(2, 0, -1)], -0.5*eccentricity)
                assert isclose(ecc_results_by_lpq[(2, 1, 0)], (1.0 - eccentricity**2)**-1.5)
                assert len(ecc_results_by_lp) == 3
                assert len(ecc_results_by_lp[(2, 1)]) == 3
                assert isclose(ecc_results_by_lp[(2, 1)][(0,)], (1.0 - eccentricity**2)**-1.5)
                assert isclose(ecc_results_by_lp[(2, 1)][(1,)], 1.5 * eccentricity)
                assert len(ecc_results_by_lp[(2, 2)]) == 3
            elif truncation == 3:
                assert len(ecc_results_by_lpq) == 19
                assert isclose(ecc_results_by_lpq[(2, 0, 0)], 1.0 - 2.5*eccentricity**2)
                assert isclose(ecc_results_by_lpq[(2, 0, -3)], 0.020833333333333333*eccentricity**3)
                assert isclose(ecc_results_by_lpq[(2, 1, 1)], 1.6875*eccentricity**3 + 1.5*eccentricity)
                assert len(ecc_results_by_lp) == 3
                assert len(ecc_results_by_lp[(2, 1)]) == 7
                assert isclose(ecc_results_by_lp[(2, 1)][(-2,)], 2.25*eccentricity**2)
                assert isclose(ecc_results_by_lp[(2, 1)][(3,)], 3.3125*eccentricity**3)
                assert len(ecc_results_by_lp[(2, 2)]) == 6
            elif truncation == 4:
                # Kaula (1964); Murray and Dermott (1999): G_200 and G_202 through e^4.
                assert len(ecc_results_by_lpq) == 25
                assert isclose(ecc_results_by_lpq[(2, 0, 0)], 1.0 - 2.5*eccentricity**2 + 0.8125*eccentricity**4)
                assert isclose(ecc_results_by_lpq[(2, 0, 2)], 8.5*eccentricity**2 - (115.0 / 6.0)*eccentricity**4)
                assert len(ecc_results_by_lp[(2, 1)]) == 9
                assert len(ecc_results_by_lp[(2, 2)]) == 8


@pytest.mark.parametrize('degree_l', (2, 3, 4, 5, 6, 7, 8, 9, 10))
@pytest.mark.parametrize('truncation', (1, 2, 3, 4, 5))
def test_eccentricity_truncation_order(degree_l, truncation):
    """Truncation N keeps every term through e^N: halving e shrinks the gap to the e^20 table by at least 2^(N+1).

    Modes the truncation drops count with a value of zero. Levels 10 and above are not checked this way because their
    gaps at small e are below double precision.
    """

    # Small enough that the next-order term does not bend the ratio at high degree, where Hansen coefficients grow.
    large_eccentricity = 0.01
    small_eccentricity = 0.005
    gaps = dict()
    for eccentricity in (large_eccentricity, small_eccentricity):
        reference = dict(eccentricity_func(eccentricity, degree_l, 20)[0])
        truncated = dict(eccentricity_func(eccentricity, degree_l, truncation)[0])
        assert set(truncated) <= set(reference)
        gaps[eccentricity] = {key: abs(truncated.get(key, 0.0) - value) for key, value in reference.items()}

    num_checked = 0
    for key, gap_small in gaps[small_eccentricity].items():
        # Skip gaps at round-off, where the ratio carries no information.
        if gap_small < 1.0e-13:
            continue
        num_checked += 1
        assert gaps[large_eccentricity][key] / gap_small > 0.9 * 2**(truncation + 1), key
    assert num_checked > 0
