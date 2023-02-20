""" Test the `TidalPy.radial_solver.numerical.interfaces` module. """

import TidalPy
TidalPy.test_mode()

def test_radial_solver_numerical_interfaces_package():
    """ Test radial_solver.numerical.interfaces import. """

    from TidalPy.radial_solver.numerical.interfaces import (interface_LDy_LDy, interface_LSt_LSt, interface_LDy_LSt,
                                                            interface_LSt_LDy, interface_LDy_SDy, interface_LSt_SSt,
                                                            interface_LDy_SSt, interface_LSt_SDy, interface_SDy_LDy,
                                                            interface_SSt_LSt, interface_SDy_LSt, interface_SSt_LDy,
                                                            interface_SDy_SDy, interface_SSt_SSt, interface_SDy_SSt,
                                                            interface_SSt_SDy)
    from TidalPy.radial_solver.numerical.interfaces import find_interface_func
    assert True
