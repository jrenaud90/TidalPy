""" Test the `TidalPy.radial_solver.numerical.interfaces` for known interfaces. """

import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.radial_solver.numerical.interfaces import find_interface_func, known_interfaces

# Stored by (is_lower_solid, is upper_solid)
correct_interfaces = {}

from TidalPy.radial_solver.numerical.interfaces.solid_solid import both_dynamic as ss_dd
from TidalPy.radial_solver.numerical.interfaces.solid_solid import both_static as ss_ss
from TidalPy.radial_solver.numerical.interfaces.solid_solid import dynamic_static as ss_ds
from TidalPy.radial_solver.numerical.interfaces.solid_solid import static_dynamic as ss_sd

solid_solid = {
    # Lower Static, Upper Static
    (True, True) : (ss_ss, False),
    # Lower Static, Upper Dynamic
    (True, False) : (ss_sd, False),
    # Lower Dynamic, Upper Static
    (False, True) : (ss_ds, False),
    # Lower Dynamic, Upper Dynamic
    (False, False) : (ss_dd, False),
    }
correct_interfaces[(True, True)] = solid_solid

from TidalPy.radial_solver.numerical.interfaces.solid_liquid import both_dynamic as sl_dd
from TidalPy.radial_solver.numerical.interfaces.solid_liquid import both_static as sl_ss
from TidalPy.radial_solver.numerical.interfaces.solid_liquid import dynamic_static as sl_ds
from TidalPy.radial_solver.numerical.interfaces.solid_liquid import static_dynamic as sl_sd

solid_liquid = {
    # Lower Static, Upper Static
    (True, True) : (sl_ss, True),
    # Lower Static, Upper Dynamic
    (True, False) : (sl_sd, False),
    # Lower Dynamic, Upper Static
    (False, True) : (sl_ds, True),
    # Lower Dynamic, Upper Dynamic
    (False, False) : (sl_dd, False),
    }
correct_interfaces[(True, False)] = solid_liquid

from TidalPy.radial_solver.numerical.interfaces.liquid_solid import both_dynamic as ls_dd
from TidalPy.radial_solver.numerical.interfaces.liquid_solid import both_static as ls_ss
from TidalPy.radial_solver.numerical.interfaces.liquid_solid import dynamic_static as ls_ds
from TidalPy.radial_solver.numerical.interfaces.liquid_solid import static_dynamic as ls_sd

liquid_solid = {
    # Lower Static, Upper Static
    (True, True) : (ls_ss, True),
    # Lower Static, Upper Dynamic
    (True, False) : (ls_sd, True),
    # Lower Dynamic, Upper Static
    (False, True) : (ls_ds, False),
    # Lower Dynamic, Upper Dynamic
    (False, False) : (ls_dd, False),
    }
correct_interfaces[(False, True)] = liquid_solid

from TidalPy.radial_solver.numerical.interfaces.liquid_liquid import both_dynamic as ll_dd
from TidalPy.radial_solver.numerical.interfaces.liquid_liquid import both_static as ll_ss
from TidalPy.radial_solver.numerical.interfaces.liquid_liquid import dynamic_static as ll_ds
from TidalPy.radial_solver.numerical.interfaces.liquid_liquid import static_dynamic as ll_sd

liquid_liquid = {
    # Lower Static, Upper Static
    (True, True) : (ll_ss, False),
    # Lower Static, Upper Dynamic
    (True, False) : (ll_sd, True),
    # Lower Dynamic, Upper Static
    (False, True) : (ll_ds, True),
    # Lower Dynamic, Upper Dynamic
    (False, False) : (ll_dd, False),
    }
correct_interfaces[(False, False)] = liquid_liquid

@pytest.mark.parametrize('is_lower_solid', (True, False))
@pytest.mark.parametrize('is_upper_solid', (True, False))
@pytest.mark.parametrize('is_lower_static', (True, False))
@pytest.mark.parametrize('is_upper_static', (True, False))
def test_known_interfaces(is_lower_solid, is_upper_solid, is_lower_static, is_upper_static):

    tidalpy_interface, tidalpy_extra_arg = \
        known_interfaces[(is_lower_solid, is_upper_solid, is_lower_static, is_upper_static)]
    correct_interface, correct_extra_arg = \
        correct_interfaces[(is_lower_solid, is_upper_solid)][(is_lower_static, is_upper_static)]

    assert tidalpy_interface is correct_interface
    assert tidalpy_extra_arg is correct_extra_arg

@pytest.mark.parametrize('is_lower_solid', (True, False))
@pytest.mark.parametrize('is_upper_solid', (True, False))
@pytest.mark.parametrize('is_lower_static', (True, False))
@pytest.mark.parametrize('is_upper_static', (True, False))
def test_find_interface_func(is_lower_solid, is_upper_solid, is_lower_static, is_upper_static):

    # Get the expected interface function
    correct_interface, correct_extra_arg = \
        correct_interfaces[(is_lower_solid, is_upper_solid)][(is_lower_static, is_upper_static)]

    # Get the interface function and additional arguments according to TidalPy
    static_liquid_density, interface_gravity, G = 10., 20., 30.
    interface_func, extra_inputs = \
        find_interface_func(is_lower_solid, is_lower_static, is_upper_solid, is_upper_static,
                            static_liquid_density, interface_gravity, G_to_use=G)

    assert interface_func is correct_interface

    if correct_extra_arg:
        # There should be extra arguments
        assert len(extra_inputs) == 3
        assert extra_inputs[0] == interface_gravity
        assert extra_inputs[1] == static_liquid_density
        assert extra_inputs[2] == G
    else:
        # There should be no extra arguments
        assert extra_inputs == tuple()
