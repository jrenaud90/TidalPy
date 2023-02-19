""" Test the `TidalPy.radial_solver.numerical` module. """

import TidalPy
TidalPy.test_mode()


def test_radial_solver_numerical_package():
    """ Test radial_solver.numerical import. """

    from TidalPy.radial_solver.numerical import radial_solver, radial_solver_numba
    assert True
