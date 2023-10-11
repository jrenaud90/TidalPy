""" Test the `TidalPy.radial_solver.numerical.initial` module. """

import TidalPy
TidalPy.test_mode()

def test_radial_solver_numerical_initial_package():
    """ Test radial_solver.numerical.initial import. """

    from TidalPy.radial_solver.numerical.initial import takeuchi_phi_psi, takeuchi_phi_psi_general, z_calc
    from TidalPy.radial_solver.numerical.initial import (
        solid_static_incompressible_kmn15,
        solid_static_incompressible_other,
        solid_static_compressible_kmn15,
        solid_static_compressible_other,
        solid_dynamic_incompressible_kmn15,
        solid_dynamic_incompressible_other,
        solid_dynamic_compressible_kmn15,
        solid_dynamic_compressible_other,
        liquid_static_incompressible_kmn15,
        liquid_static_incompressible_other,
        liquid_static_compressible_kmn15,
        liquid_static_compressible_other,
        liquid_dynamic_incompressible_kmn15,
        liquid_dynamic_incompressible_other,
        liquid_dynamic_compressible_kmn15,
        liquid_dynamic_compressible_other)
    from TidalPy.radial_solver.numerical.initial import known_initial_guess_funcs
    from TidalPy.radial_solver.numerical.initial import find_initial_guess

    assert True
