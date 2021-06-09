""" Tests for calculating the initial guess at the bottom of a liquid or solid layer for various types of interfaces
"""

import numpy as np

import TidalPy
from TidalPy.constants import G
from TidalPy.tides.multilayer.numerical_int import liquid_dynamic_guess, liquid_static_guess, solid_static_guess, \
    solid_dynamic_guess
from TidalPy.tides.multilayer.numerical_int.interfaces import interface_LDy_LDy, interface_LSt_LSt, interface_SDy_SDy, \
    interface_SSt_SSt, interface_SSt_SDy, interface_SDy_SSt, interface_LDy_SDy, interface_LDy_SSt, interface_LSt_SDy, \
    interface_LSt_SSt, interface_SDy_LDy, interface_SDy_LSt, interface_SSt_LDy, interface_SSt_LSt

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()

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
LDy_guess_l2 = liquid_dynamic_guess(radius_array_to_use, bulk_array, density_array, frequency, order_l=2)
LSt_guess_l2 = liquid_static_guess(radius_array_to_use, order_l=2)
SSt_guess_l2 = solid_static_guess(radius_array_to_use, shear_array, bulk_array, density_array, order_l=2)
SDy_guess_l2 = solid_dynamic_guess(radius_array_to_use, shear_array, bulk_array, density_array, frequency, order_l=2)

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
        assert type(result) == tuple
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

    i = 0
    for test_func, test_input in tests.items():
        result = test_func(*test_input)

        if i == 0:
            # Static Result
            assert type(result) == np.ndarray
            assert result.dtype == test_input[0].dtype
        else:
            # Dynamic Result
            assert type(result) == tuple
            assert len(result) == 2
            for solution in range(2):
                assert result[solution].dtype == test_input[0][0].dtype

        i += 1

def test_liquid_solid_interface():
    tests = {
        interface_LSt_SSt: (LSt_guess_l2, gravity_array[-1], density_array[-1]),
        interface_LDy_SDy: (LDy_guess_l2,),
        interface_LSt_SDy: (LSt_guess_l2, gravity_array[-1], density_array[-1]),
        interface_LDy_SSt: (LDy_guess_l2,)
    }

    i = 0
    for test_func, test_input in tests.items():
        result = test_func(*test_input)

        if i in (0, 3):
            # Static Result
            assert type(result) == tuple
            assert len(result) == 3
            for solution in range(3):
                if i in (0, 2):
                    # Static Input
                    assert result[solution].dtype == test_input[0].dtype
                else:
                    # Dynamic Input
                    assert result[solution].dtype == test_input[0][0].dtype
        else:
            # Dynamic Result
            assert type(result) == tuple
            assert len(result) == 3
            for solution in range(3):
                if i in (0, 2):
                    # Static Input
                    assert result[solution].dtype == test_input[0].dtype
                else:
                    # Dynamic Input
                    assert result[solution].dtype == test_input[0][0].dtype

        i += 1

def test_solid_liquid_interface():
    tests = {
        interface_SSt_LSt: (SSt_guess_l2, gravity_array[-1], density_array[-1]),
        interface_SDy_LDy: (SDy_guess_l2,),
        interface_SSt_LDy: (SSt_guess_l2,),
        interface_SDy_LSt: (SDy_guess_l2, gravity_array[-1], density_array[-1])
    }

    i = 0
    for test_func, test_input in tests.items():
        result = test_func(*test_input)

        if i in (0, 3):
            # Static Result
            assert type(result) == np.ndarray
            if i in (0, 2):
                # Static Input
                assert result.dtype == test_input[0][0].dtype
            else:
                # Dynamic Input
                assert result.dtype == test_input[0][0].dtype
        else:
            # Dynamic Result
            assert type(result) == tuple
            assert len(result) == 2
            for solution in range(2):
                if i in (0, 2):
                    # Static Input
                    assert result[solution].dtype == test_input[0][0].dtype
                else:
                    # Dynamic Input
                    assert result[solution].dtype == test_input[0][0].dtype

        i += 1