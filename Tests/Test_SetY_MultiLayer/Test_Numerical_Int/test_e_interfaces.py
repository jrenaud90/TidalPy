""" Tests for calculating the initial guess at the bottom of a liquid or solid layer for various types of interfaces
"""

import numpy as np
from numba.typed.typedlist import List as nbTypedList

import TidalPy
TidalPy.test_mode()

from TidalPy.constants import G
from TidalPy.tides.multilayer.numerical_int.initial_conditions import (liquid_dynamic_guess_ts72,
                                                                       liquid_static_guess_s74,
                                                                       solid_dynamic_guess_ts72,
                                                                       solid_static_guess_ts72)
from TidalPy.tides.multilayer.numerical_int.interfaces import (interface_LDy_LDy, interface_LDy_SDy, interface_LDy_SSt,
                                                               interface_LSt_LSt, interface_LSt_SDy, interface_LSt_SSt,
                                                               interface_SDy_LDy, interface_SDy_LSt, interface_SDy_SDy,
                                                               interface_SDy_SSt, interface_SSt_LDy, interface_SSt_LSt,
                                                               interface_SSt_SDy, interface_SSt_SSt)

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array = np.linspace(0., 1.e6, 11)
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.complex128)
bulk_array = 10.e10 * np.ones(10, dtype=np.float64)
radius_array_to_use = radius_array[1:]
frequency = 2. * np.pi / (86400. * 1.)

# Find lower layer guesses
# # l = 2
LDy_guess_l2 = liquid_dynamic_guess_ts72(radius_array_to_use, bulk_array, density_array, frequency, order_l=2)
LSt_guess_l2 = liquid_static_guess_s74(radius_array_to_use, order_l=2)
SSt_guess_l2 = solid_static_guess_ts72(radius_array_to_use, shear_array, bulk_array, density_array, order_l=2)
SDy_guess_l2 = solid_dynamic_guess_ts72(
    radius_array_to_use, shear_array, bulk_array, density_array, frequency,
    order_l=2
    )


def test_solid_solid_interface():
    tests = {
        interface_SSt_SSt: (SSt_guess_l2,),
        interface_SDy_SDy: (SDy_guess_l2,),
        interface_SSt_SDy: (SSt_guess_l2,),
        interface_SDy_SSt: (SDy_guess_l2,)
        }

    i = 0
    for test_func, test_input in tests.items():
        result = test_func(*test_input)
        assert type(result) in [nbTypedList, list]
        assert len(result) == 3
        for solution in range(3):
            assert result[solution].dtype == test_input[0][0].dtype

        i += 1


def test_liquid_liquid_interface():
    tests = {
        interface_LSt_LSt: (LSt_guess_l2,),
        interface_LDy_LDy: (LDy_guess_l2,),
        # TODO: these are not implemented yet
        # interface_LSt_LDy: (LSt_guess_l2,),
        # interface_LDy_LSt: (LDy_guess_l2,)
        }

    for test_func, test_input in tests.items():
        results = test_func(*test_input)

        assert type(results) in [nbTypedList, list]
        for result in results:
            assert type(result.dtype) == type(test_input[0][0].dtype)


def test_liquid_solid_interface():
    tests = {
        interface_LSt_SSt: (LSt_guess_l2, gravity_array[-1], density_array[-1]),
        interface_LDy_SDy: (LDy_guess_l2,),
        interface_LSt_SDy: (LSt_guess_l2, gravity_array[-1], density_array[-1]),
        interface_LDy_SSt: (LDy_guess_l2,)
        }

    for test_func, test_input in tests.items():
        results = test_func(*test_input)

        assert type(results) in [nbTypedList, list]
        for result in results:
            assert type(result.dtype) == type(test_input[0][0].dtype)


def test_solid_liquid_interface():
    tests = {
        interface_SSt_LSt: (SSt_guess_l2, gravity_array[-1], density_array[-1]),
        interface_SDy_LDy: (SDy_guess_l2,),
        interface_SSt_LDy: (SSt_guess_l2,),
        interface_SDy_LSt: (SDy_guess_l2, gravity_array[-1], density_array[-1])
        }

    for test_func, test_input in tests.items():
        results = test_func(*test_input)

        assert type(results) in [nbTypedList, list]
        for result in results:
            assert type(result.dtype) == type(test_input[0][0].dtype)
