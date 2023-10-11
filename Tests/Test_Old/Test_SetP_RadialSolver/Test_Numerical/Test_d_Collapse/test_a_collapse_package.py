""" Test the `TidalPy.radial_solver.numerical.collapse` module. """

import TidalPy
TidalPy.test_mode()

def test_radial_solver_numerical_collapse_package():
    """ Test radial_solver.numerical.collapse import. """

    from TidalPy.radial_solver.numerical.collapse import solid_surface, static_liquid_surface, dynamic_liquid_surface
    from TidalPy.radial_solver.numerical.collapse import collapse_solutions

    assert True
