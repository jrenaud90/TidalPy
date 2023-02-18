""" Test the `TidalPy.radial_solver.numerical.initial` for known initial functions. """

import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.radial_solver.numerical.initial import known_initial_guess_funcs

from TidalPy.radial_solver.numerical.initial.initial_solution_dynamic import (
    liquid_guess_kamata   as liquid_dynamic_compressible_kmn15,
    liquid_guess_takeuchi as liquid_dynamic_compressible_other,
    solid_guess_kamata    as solid_dynamic_compressible_kmn15,
    solid_guess_takeuchi  as solid_dynamic_compressible_other)

from TidalPy.radial_solver.numerical.initial.initial_solution_dynamic_incomp import (
    liquid_guess_kamata   as liquid_dynamic_incompressible_kmn15,
    liquid_guess_takeuchi as liquid_dynamic_incompressible_other,
    solid_guess_kamata    as solid_dynamic_incompressible_kmn15,
    solid_guess_takeuchi  as solid_dynamic_incompressible_other)

from TidalPy.radial_solver.numerical.initial.initial_solution_static import (
    liquid_guess_saito   as liquid_static_compressible_other,
    solid_guess_kamata   as solid_static_compressible_kmn15,
    solid_guess_takeuchi as solid_static_compressible_other)

from TidalPy.radial_solver.numerical.initial.initial_solution_static_incomp import (
    liquid_guess_saito   as liquid_static_incompressible_other,
    solid_guess_kamata   as solid_static_incompressible_kmn15,
    solid_guess_takeuchi as solid_static_incompressible_other)

# Kamata et al. (2015) did not provide equations for static liquid layers so use the Saito (1974) method instead.
liquid_static_compressible_kmn15   = liquid_static_compressible_other
liquid_static_incompressible_kmn15 = liquid_static_incompressible_other

# Known initial guess functions stored by: is_solid, is_static, is_incompressible, is_kamata
correct_initial_guess_funcs = {
    (True, True, True, True)    : solid_static_incompressible_kmn15,
    (True, True, True, False)   : solid_static_incompressible_other,
    (True, True, False, True)   : solid_static_compressible_kmn15,
    (True, True, False, False)  : solid_static_compressible_other,
    (True, False, True, True)   : solid_dynamic_incompressible_kmn15,
    (True, False, True, False)  : solid_dynamic_incompressible_other,
    (True, False, False, True)  : solid_dynamic_compressible_kmn15,
    (True, False, False, False) : solid_dynamic_compressible_other,
    (False, True, True, True)   : liquid_static_incompressible_kmn15,
    (False, True, True, False)  : liquid_static_incompressible_other,
    (False, True, False, True)  : liquid_static_compressible_kmn15,
    (False, True, False, False) : liquid_static_compressible_other,
    (False, False, True, True)  : liquid_dynamic_incompressible_kmn15,
    (False, False, True, False) : liquid_dynamic_incompressible_other,
    (False, False, False, True) : liquid_dynamic_compressible_kmn15,
    (False, False, False, False): liquid_dynamic_compressible_other,
    }

@pytest.mark.parametrize('is_solid', (True, False))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('is_kamata', (True, False))
def test_known_initial_funcs(is_solid, is_static, is_incompressible, is_kamata):

    # Get the expected initial function
    correct_initial_func = \
        correct_initial_guess_funcs[(is_solid, is_static, is_incompressible, is_kamata)]

    # Get the initial function according to TidalPy
    tidalpy_initial_func = \
        known_initial_guess_funcs[(is_solid, is_static, is_incompressible, is_kamata)]

    assert tidalpy_initial_func is correct_initial_func
