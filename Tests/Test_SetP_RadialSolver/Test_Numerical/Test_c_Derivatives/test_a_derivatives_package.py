""" Test the `TidalPy.radial_solver.numerical.derivatives` module. """

import TidalPy
TidalPy.test_mode()

def test_radial_solver_numerical_derivatives_package():
    """ Test radial_solver.numerical.derivatives import. """

    from TidalPy.radial_solver.numerical.derivatives import (
        radial_derivatives_liquid_dynamic, radial_derivatives_solid_dynamic,
        radial_derivatives_liquid_static, radial_derivatives_solid_static,
        radial_derivatives_liquid_dynamic_incomp, radial_derivatives_solid_dynamic_incomp,
        radial_derivatives_liquid_static_incomp, radial_derivatives_solid_static_incomp)
    from TidalPy.radial_solver.numerical.derivatives import (
        dynamic_liquid_ode, dynamic_solid_ode,
        static_liquid_ode, static_solid_ode)
    from TidalPy.radial_solver.numerical.derivatives import known_multilayer_odes, find_ode

    assert True
