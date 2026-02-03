from math import isclose

import pytest
import numpy as np

from TidalPy.Tides_x.obliquity import obliquity_func
from TidalPy.utilities.lookups import IntMap3


@pytest.mark.parametrize('degree_l', (1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
@pytest.mark.parametrize('truncation', ('gen', 2, 4, 'off', 5))
def test_obliquity_funcs(degree_l, truncation):
    """Tests TidalPy's obliquity functions for various degree_ls and various truncation levels."""

    obliquity = 0.5

    if degree_l == 1:
        # Degree l is currently not supported, check that it raises an error.
        with pytest.raises(NotImplementedError):
            ob_results_by_lmp, ob_results_by_lm = obliquity_func(obliquity, degree_l, truncation)
    elif truncation == 5:
        # Truncation = 5 is not supported, check that it raises an error.
        with pytest.raises(NotImplementedError):
            ob_results_by_lmp, ob_results_by_lm = obliquity_func(obliquity, degree_l, truncation)
    else:
        ob_results_by_lmp, ob_results_by_lm = obliquity_func(obliquity, degree_l, truncation)

        # Check return types
        assert isinstance(ob_results_by_lmp, IntMap3)
        assert isinstance(ob_results_by_lm, dict)

        # No matter the assumptions, these should have some size to them.
        assert len(ob_results_by_lmp) > 0
        for (l, m, p), ob_result in ob_results_by_lmp:
            assert isinstance(l, int)
            assert isinstance(m, int)
            assert isinstance(p, int)
            assert l == degree_l
            assert m <= l
            assert p <= l
            assert isinstance(ob_result, float)
        
        for (l, m), ob_results_by_p in ob_results_by_lm.items():
            assert isinstance(l, int)
            assert isinstance(m, int)
            assert l == degree_l
            assert m <= l
            for (p,), ob_result in ob_results_by_p:
                assert isinstance(p, int)
                assert isinstance(ob_result, float)
                assert p <= degree_l
        
        # Spot checks
        if degree_l == 2:
            if truncation == 'off':
                assert len(ob_results_by_lmp) == 2
                assert isclose(ob_results_by_lmp[(2, 0, 1)], -0.5)
                assert isclose(ob_results_by_lmp[(2, 2, 0)], 3.0)
                assert len(ob_results_by_lm) == 2
                assert len(ob_results_by_lm[(2, 0)]) == 1
                assert isclose(ob_results_by_lm[(2, 0)][(1,)], -0.5)
                assert len(ob_results_by_lm[(2, 2)]) == 1
                assert isclose(ob_results_by_lm[(2, 2)][(0,)], 3.0)
            elif truncation == 2:
                assert len(ob_results_by_lmp) == 4
                assert isclose(ob_results_by_lmp[(2, 0, 1)], -0.5)
                assert isclose(ob_results_by_lmp[(2, 1, 0)], 1.5 * obliquity)
                assert len(ob_results_by_lm) == 3
                assert len(ob_results_by_lm[(2, 1)]) == 2
                assert isclose(ob_results_by_lm[(2, 1)][(0,)], 1.5 * obliquity)
                assert isclose(ob_results_by_lm[(2, 1)][(1,)], -1.5 * obliquity)
                assert len(ob_results_by_lm[(2, 2)]) == 1
            elif truncation == 4:
                assert len(ob_results_by_lmp) == 8
                assert isclose(ob_results_by_lmp[(2, 0, 1)], 0.75*obliquity**2 - 0.5)
                assert isclose(ob_results_by_lmp[(2, 1, 0)], -0.625*obliquity**3 + 1.5*obliquity)
                assert len(ob_results_by_lm) == 3
                assert len(ob_results_by_lm[(2, 1)]) == 3
                assert isclose(ob_results_by_lm[(2, 1)][(0,)], -0.625*obliquity**3 + 1.5*obliquity)
                assert isclose(ob_results_by_lm[(2, 1)][(1,)], obliquity**3 - 1.5*obliquity)
                assert len(ob_results_by_lm[(2, 2)]) == 2
            elif truncation == 'gen':
                assert len(ob_results_by_lmp) == 9
                assert isclose(ob_results_by_lmp[(2, 0, 1)], -np.sin(obliquity/2)**4 + np.sin(obliquity/2)**2 + 0.5*np.sin(obliquity)**2 - 0.5)
                assert isclose(ob_results_by_lmp[(2, 1, 0)], 3.0*np.sin(obliquity/2)*np.cos(obliquity/2)**3)
                assert len(ob_results_by_lm) == 3
                assert len(ob_results_by_lm[(2, 1)]) == 3
                assert len(ob_results_by_lm[(2, 2)]) == 3
                assert isclose(ob_results_by_lm[(2, 2)][(0,)], 3.0*np.cos(obliquity/2)**4)
                assert isclose(ob_results_by_lm[(2, 2)][(2,)], 3.0*np.sin(obliquity/2)**4)
