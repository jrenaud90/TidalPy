import pytest
from math import isnan, isclose
from copy import copy

import numpy as np

from TidalPy.exceptions import ArgumentException
from TidalPy.rheology import Maxwell, Andrade
from TidalPy.RadialSolver.helpers import build_rs_input_homogeneous_layers, PlanetBuildData, build_rs_input_from_data

planet_radius     = 1.0e6
forcing_frequency = 1.0e-5
density_list                 = [10_000.0, 9_000.0, 7_500.0, 3_500.0, 2_000]
static_bulk_modulus_list     = [200.0e9, 150.0e9, 100.0e9, 80.0e9, 40.0e9]
static_shear_modulus_list    = [80.0e9, 60.0e9, 40.0e9, 20.0e9, 5.0e9]
bulk_viscosity_list          = [1.0e30, 1.0e28, 1.0e26, 1.0e24, 1.0e22]
shear_viscosity_list         = [1.0e26, 1.0e24, 1.0e22, 1.0e20, 1.0e18]
layer_type_list              = ['solid', 'liquid', 'solid', 'liquid', 'solid']
layer_is_static_list         = [False, True, False, True, False]
layer_is_incompressible_list = [True, False, True, False, True]
shear_rheology_model_list    = [Maxwell(), Andrade(), Maxwell(), Andrade(), Maxwell()]
bulk_rheology_model_list     = [Andrade(), Maxwell(), Andrade(), Maxwell(), Andrade()]
slices_tuple_list            = [5, 10, 11, 12, 8]

radius_fraction_lists    = [
    [1.0],
    [0.1, 1.0],
    [0.1, 0.3, 1.0],
    [0.1, 0.3, 0.45, 1.0],
    [0.1, 0.3, 0.45, 0.63, 1.0],
    ]
thickness_fraction_lists = [
    [1.0],
    [0.1, 0.9],
    [0.1, 0.2, 0.7],
    [0.1, 0.2, 0.15, 0.55],
    [0.1, 0.2, 0.15, 0.18, 0.37],
    ] 
volume_fraction_lists    = [
    [1.0],
    [0.00100000, 0.999],
    [0.00100000, 0.02600000, 0.973],
    [0.00100000, 0.02600000, 0.06412500, 0.908875],
    [0.00100000, 0.02600000, 0.06412500, 0.15892200, 0.74995300],
    ]


@pytest.mark.parametrize('use_slices_tuple', (True, False))
@pytest.mark.parametrize('num_layers', (1, 2, 3, 4, 5))
def test_build_planet_constant_layers(num_layers, use_slices_tuple):
    """ Test the `build_rs_input_homogeneous_layers` helper function. """

    slice_per_layer = 15
    if use_slices_tuple:
        slices_tuple = tuple(slices_tuple_list[:num_layers])
    else:
        slices_tuple = None
    
    common_inputs = (
        planet_radius,
        forcing_frequency, 
        tuple(density_list[:num_layers]),
        tuple(static_bulk_modulus_list[:num_layers]),
        tuple(static_shear_modulus_list[:num_layers]),
        tuple(bulk_viscosity_list[:num_layers]),
        tuple(shear_viscosity_list[:num_layers]),
        tuple(layer_type_list[:num_layers]),
        tuple(layer_is_static_list[:num_layers]),
        tuple(layer_is_incompressible_list[:num_layers]),
        tuple(shear_rheology_model_list[:num_layers]),
        tuple(bulk_rheology_model_list[:num_layers])
    )
    
    rad_frac_result = build_rs_input_homogeneous_layers(
        *common_inputs,
        radius_fraction_tuple=tuple(radius_fraction_lists[num_layers-1]),
        slices_tuple=slices_tuple,
        slice_per_layer=slice_per_layer)
    
    vol_frac_result = build_rs_input_homogeneous_layers(
        *common_inputs,
        volume_fraction_tuple=tuple(volume_fraction_lists[num_layers-1]),
        slices_tuple=slices_tuple,
        slice_per_layer=slice_per_layer)
    
    thick_frac_result = build_rs_input_homogeneous_layers(
        *common_inputs,
        thickness_fraction_tuple=tuple(thickness_fraction_lists[num_layers-1]),
        slices_tuple=slices_tuple,
        slice_per_layer=slice_per_layer)
    
    # Check that output is correct
    for res in (rad_frac_result, vol_frac_result, thick_frac_result):
        assert isinstance(res, PlanetBuildData)
        assert len(res) == 10

        # Check all fields are present
        assert res.radius_array is not None
        assert res.density_array is not None
        assert res.complex_bulk_modulus_array is not None
        assert res.complex_shear_modulus_array is not None
        assert res.frequency is not None
        assert res.planet_bulk_density is not None
        assert res.layer_types is not None
        assert res.is_static_bylayer is not None
        assert res.is_incompressible_bylayer is not None
        assert res.upper_radius_bylayer_array is not None
        
        # Unpack as a tuple
        radius_array, density_array, complex_bulk_array, complex_shear_array, freq, planet_bulk_density, \
            layer_type_tuple, layer_is_static_tuple, layer_is_incompressible_tuple, upper_radius_array = res
        
        assert type(radius_array) == np.ndarray
        assert radius_array.dtype == np.float64
        if use_slices_tuple:
            assert radius_array.size == sum(slices_tuple_list[:num_layers])
        else:
            assert radius_array.size == num_layers * slice_per_layer
        
        assert radius_array.size == density_array.size
        assert density_array.dtype == np.float64
        assert radius_array.size == complex_bulk_array.size
        assert complex_bulk_array.dtype == np.complex128
        assert radius_array.size == complex_shear_array.size
        assert complex_shear_array.dtype == np.complex128
        assert len(layer_type_tuple) == num_layers
        assert len(layer_is_static_tuple) == num_layers
        assert len(layer_is_incompressible_tuple) == num_layers
        assert type(upper_radius_array) == np.ndarray
        assert upper_radius_array.size == num_layers
        assert upper_radius_array.dtype == np.float64
        assert type(planet_bulk_density) == float

        assert layer_type_tuple == tuple(layer_type_list[:num_layers])
        assert layer_is_static_tuple == tuple(layer_is_static_list[:num_layers])
        assert layer_is_incompressible_tuple == tuple(layer_is_incompressible_list[:num_layers])

    # Check that the different input methods give the same result
    for array_i in (0, 1, 2, 3, 9):
        assert np.allclose(rad_frac_result[array_i], vol_frac_result[array_i])
        assert np.allclose(rad_frac_result[array_i], thick_frac_result[array_i])
    
    # Check floats are the same across different inputs
    for float_i in (4, 5):
        assert isclose(rad_frac_result[float_i], vol_frac_result[float_i])
        assert isclose(rad_frac_result[float_i], thick_frac_result[float_i])
    


def test_build_planet_constant_layers_exceptions():

    num_layers = 5

    common_inputs_list = [
        planet_radius,
        forcing_frequency, 
        tuple(density_list[:num_layers]),
        tuple(static_bulk_modulus_list[:num_layers]),
        tuple(static_shear_modulus_list[:num_layers]),
        tuple(bulk_viscosity_list[:num_layers]),
        tuple(shear_viscosity_list[:num_layers]),
        tuple(layer_type_list[:num_layers]),
        tuple(layer_is_static_list[:num_layers]),
        tuple(layer_is_incompressible_list[:num_layers]),
        tuple(shear_rheology_model_list[:num_layers]),
        tuple(bulk_rheology_model_list[:num_layers])
    ]

    # Check different layer lens
    with pytest.raises(ArgumentException):
        num_layers_2 = 3

        for list_change in (2, 3, 4, 5, 6, 7, 8, 9):
            inputs = copy(common_inputs_list)
            inputs[list_change] = tuple(common_inputs_list[:num_layers_2])

        result = build_rs_input_homogeneous_layers(
            *tuple(inputs),
            radius_fraction_tuple=tuple(radius_fraction_lists[num_layers-1]),
            perform_checks=True)
    
    # Check missing structure inputs
    with pytest.raises(ArgumentException):
        result = build_rs_input_homogeneous_layers(
            *tuple(common_inputs_list),
            thickness_fraction_tuple=None,
            radius_fraction_tuple=None,
            volume_fraction_tuple=None,
            perform_checks=True)
    
    # Check error is raised if incorrect
    with pytest.raises(ArgumentException):
        for list_change in (10, 11):
            inputs = copy(common_inputs_list)
            inputs[list_change] = tuple(['Maxwell', 'Andrade', 'Maxwell', 'Andrade', 'Maxwell'])
            result = build_rs_input_homogeneous_layers(
                *tuple(inputs),
                radius_fraction_tuple=tuple(radius_fraction_lists[num_layers-1]),
                perform_checks=True)
            
    # Check bad number of slices
    with pytest.raises(ArgumentException):
        result = build_rs_input_homogeneous_layers(
            *tuple(common_inputs_list),
            radius_fraction_tuple=tuple(radius_fraction_lists[num_layers-1]),
            slices_tuple=(6, 10, 3, 12, 8),
            perform_checks=True)
    with pytest.raises(ArgumentException):
        result = build_rs_input_homogeneous_layers(
            *tuple(common_inputs_list),
            radius_fraction_tuple=tuple(radius_fraction_lists[num_layers-1]),
            slice_per_layer=3,
            perform_checks=True)
    
    # Check bad structure inputs
    with pytest.raises(ArgumentException):
        # Too short
        result = build_rs_input_homogeneous_layers(
            *tuple(common_inputs_list),
            radius_fraction_tuple=(0.1, 0.3, 0.45, 0.63, 0.9),
            perform_checks=True)
    with pytest.raises(ArgumentException):
        # Too far
        result = build_rs_input_homogeneous_layers(
            *tuple(common_inputs_list),
            radius_fraction_tuple=(0.1, 0.3, 0.45, 0.63, 1.5),
            perform_checks=True)
    with pytest.raises(ArgumentException):
        # Wrong order
        result = build_rs_input_homogeneous_layers(
            *tuple(common_inputs_list),
            radius_fraction_tuple=(0.1, 0.3, 0.45, 1.0, 0.63),
            perform_checks=True)
    with pytest.raises(ValueError):
        # Too large volume
        result = build_rs_input_homogeneous_layers(
            *tuple(common_inputs_list),
            volume_fraction_tuple=(0.00100000, 0.02600000, 0.6412500, 0.15892200, 0.74995300),
            perform_checks=True)
    with pytest.raises(ValueError):
        # Too small volume
        result = build_rs_input_homogeneous_layers(
            *tuple(common_inputs_list),
            volume_fraction_tuple=(0.00100000, 0.02600000, 0.0412500, 0.15892200, 0.14995300),
            perform_checks=True)
    with pytest.raises(ValueError):
        # Too large thickness
        result = build_rs_input_homogeneous_layers(
            *tuple(common_inputs_list),
            volume_fraction_tuple=(0.1, 0.2, 0.65, 0.18, 0.37),
            perform_checks=True)
    with pytest.raises(ValueError):
        # Too small thickness
        result = build_rs_input_homogeneous_layers(
            *tuple(common_inputs_list),
            volume_fraction_tuple=(0.1, 0.2, 0.05, 0.18, 0.07),
            perform_checks=True)
    

@pytest.mark.parametrize('num_layers', (1, 2, 3, 4, 5))
@pytest.mark.parametrize('num_slices', (100, 300, 2000))
@pytest.mark.parametrize('r0', (0.0, 2.4))
def test_build_rs_input_from_data(num_layers, num_slices, r0):

    layer_upper_radius_tuple = list()
    for i in range(num_layers): 
        layer_upper_radius_tuple.append(radius_fraction_lists[num_layers - 1][i] * planet_radius)
    layer_upper_radius_tuple = tuple(layer_upper_radius_tuple)

    radius_array  = np.linspace(r0, planet_radius, num_slices)
    density_array = np.linspace(10_000.0, 4000.0, num_slices)
    static_bulk_modulus_array = np.linspace(1500.0e9, 400.0e9, num_slices)
    static_shear_modulus_array = np.linspace(500.0e9, 50.0e9, num_slices)
    bulk_viscosity_array = np.logspace(32, 26, num_slices, dtype=np.float64)
    shear_viscosity_array = np.logspace(32, 18, num_slices, dtype=np.float64)
    
    output = build_rs_input_from_data(
        forcing_frequency,
        radius_array,
        density_array,
        static_bulk_modulus_array,
        static_shear_modulus_array,
        bulk_viscosity_array,
        shear_viscosity_array,
        layer_upper_radius_tuple,
        tuple(layer_type_list[:num_layers]),
        tuple(layer_is_static_list[:num_layers]),
        tuple(layer_is_incompressible_list[:num_layers]),
        tuple(shear_rheology_model_list[:num_layers]),
        tuple(bulk_rheology_model_list[:num_layers]),
        perform_checks = True,
        warnings = False)

    # Unpack as a tuple
    radius_array_new, density_array_new, complex_bulk_array_new, complex_shear_array_new, freq, planet_bulk_density_new, \
        layer_type_tuple_new, layer_is_static_tuple_new, layer_is_incompressible_tuple_new, upper_radius_array_new = output

    assert radius_array_new.size >= radius_array.size
    assert radius_array_new.size == density_array_new.size
    assert density_array.dtype == np.float64
    assert radius_array_new.size == complex_bulk_array_new.size
    assert complex_bulk_array_new.dtype == np.complex128
    assert radius_array_new.size == complex_shear_array_new.size
    assert complex_shear_array_new.dtype == np.complex128
    assert len(layer_type_tuple_new) == num_layers
    assert len(layer_is_static_tuple_new) == num_layers
    assert len(layer_is_incompressible_tuple_new) == num_layers
    assert type(upper_radius_array_new) == np.ndarray
    assert upper_radius_array_new.size == num_layers
    assert upper_radius_array_new.dtype == np.float64
    assert type(planet_bulk_density_new) == float

    assert layer_type_tuple_new == tuple(layer_type_list[:num_layers])
    assert layer_is_static_tuple_new == tuple(layer_is_static_list[:num_layers])
    assert layer_is_incompressible_tuple_new == tuple(layer_is_incompressible_list[:num_layers])
    
    for layer_i in range(num_layers):
        if layer_i == 0:
            assert radius_array_new[0] == 0.0
        elif layer_i == (num_layers - 1):
            # One copy of upper radius at planet surface
            assert np.sum(radius_array_new == layer_upper_radius_tuple[layer_i]) == 1
        else:
            # Two copies of radius at interfaces
            assert np.sum(radius_array_new == layer_upper_radius_tuple[layer_i]) == 2

            # Check that other properties match expectations
            first_index = np.argwhere(radius_array_new == layer_upper_radius_tuple[layer_i])[0]

            # Check that radius is the same at and above
            assert radius_array_new[first_index] == layer_upper_radius_tuple[layer_i]
            assert radius_array_new[first_index + 1] == layer_upper_radius_tuple[layer_i]

            # Assert that physical properties are the same at this and below
            try:
                assert density_array_new[first_index] == density_array_new[first_index - 1]
            except:
                import pdb; pdb.set_trace()
            assert complex_bulk_array_new[first_index] == complex_bulk_array_new[first_index - 1]
            assert complex_shear_array_new[first_index] == complex_shear_array_new[first_index - 1]
            assert density_array_new[first_index] == density_array_new[first_index - 1]

            # Assert that physical properties are the same at second index and above
            assert density_array_new[first_index + 1] == density_array_new[first_index + 2]
            assert complex_bulk_array_new[first_index + 1] == complex_bulk_array_new[first_index + 2]
            assert complex_shear_array_new[first_index + 1] == complex_shear_array_new[first_index + 2]
            assert density_array_new[first_index + 1] == density_array_new[first_index + 2]


def test_build_rs_input_from_data_check_errors():

    N = 100
    radius_array = np.linspace(0., 10., N)
    density_array = np.linspace(1000., 500., N)
    static_bulk_modulus_array = np.linspace(10., 20., N)
    static_shear_modulus_array = np.linspace(10., 20., N)
    bulk_viscosity_array = np.linspace(10., 20., N)
    shear_viscosity_array = np.linspace(10., 20., N)
    layer_upper_radius_tuple = (4.0, 10.0)
    layer_types = ('solid', 'liquid')
    layer_is_statics = (True, False)
    layer_is_incompressibles = (True, False)
    shear_rheology_models = (Maxwell(), Andrade())
    bulk_rheology_models = (Andrade(), Maxwell())

    # Radius array does not end at planet radius according to upper radius list
    with pytest.raises(ArgumentException):
        output = build_rs_input_from_data(
            forcing_frequency,
            np.linspace(0., 9.5, N),
            density_array,
            static_bulk_modulus_array,
            static_shear_modulus_array,
            bulk_viscosity_array,
            shear_viscosity_array,
            layer_upper_radius_tuple,
            layer_types,
            layer_is_statics,
            layer_is_incompressibles,
            shear_rheology_models,
            bulk_rheology_models,
            perform_checks = True,
            warnings = False)
    
    # Incorrect number of layers on assumption tuples
    with pytest.raises(ArgumentException):
        output = build_rs_input_from_data(
            forcing_frequency,
            radius_array,
            density_array,
            static_bulk_modulus_array,
            static_shear_modulus_array,
            bulk_viscosity_array,
            shear_viscosity_array,
            layer_upper_radius_tuple,
            ('solid', 'liquid', 'solid'),
            layer_is_statics,
            layer_is_incompressibles,
            shear_rheology_models,
            bulk_rheology_models,
            perform_checks = True,
            warnings = False)
    
    with pytest.raises(ArgumentException):
        output = build_rs_input_from_data(
            forcing_frequency,
            radius_array,
            density_array,
            static_bulk_modulus_array,
            static_shear_modulus_array,
            bulk_viscosity_array,
            shear_viscosity_array,
            layer_upper_radius_tuple,
            layer_types,
            layer_is_statics,
            (False,),
            shear_rheology_models,
            bulk_rheology_models,
            perform_checks = True,
            warnings = False)
    
    with pytest.raises(ArgumentException):
        output = build_rs_input_from_data(
            forcing_frequency,
            radius_array,
            density_array,
            static_bulk_modulus_array,
            static_shear_modulus_array,
            bulk_viscosity_array,
            shear_viscosity_array,
            layer_upper_radius_tuple,
            layer_types,
            (True, False, True),
            layer_is_incompressibles,
            shear_rheology_models,
            bulk_rheology_models,
            perform_checks = True,
            warnings = False)
    
    # Too few slices for number of layers
    N2 = 8
    radius_array2 = np.linspace(0., 10., N2)
    density_array2 = np.linspace(1000., 500., N2)
    static_bulk_modulus_array2 = np.linspace(10., 20., N2)
    static_shear_modulus_array2 = np.linspace(10., 20., N2)
    bulk_viscosity_array2 = np.linspace(10., 20., N2)
    shear_viscosity_array2 = np.linspace(10., 20., N2)
    layer_upper_radius_tuple = (4.0, 10.0)
    layer_types = ('solid', 'liquid')
    layer_is_statics = (True, False)
    layer_is_incompressibles = (True, False)
    shear_rheology_models = (Maxwell(), Andrade())
    bulk_rheology_models = (Andrade(), Maxwell())
    with pytest.raises(ArgumentException):
        output = build_rs_input_from_data(
            forcing_frequency,
            radius_array2,
            density_array2,
            static_bulk_modulus_array2,
            static_shear_modulus_array2,
            bulk_viscosity_array2,
            shear_viscosity_array2,
            layer_upper_radius_tuple,
            layer_types,
            layer_is_statics,
            layer_is_incompressibles,
            shear_rheology_models,
            bulk_rheology_models,
            perform_checks = True,
            warnings = False)
    
    # Arrays with mismatched sizes
    with pytest.raises(ArgumentException):
        output = build_rs_input_from_data(
            forcing_frequency,
            radius_array,
            density_array2,
            static_bulk_modulus_array,
            static_shear_modulus_array,
            bulk_viscosity_array,
            shear_viscosity_array,
            layer_upper_radius_tuple,
            layer_types,
            layer_is_statics,
            layer_is_incompressibles,
            shear_rheology_models,
            bulk_rheology_models,
            perform_checks = True,
            warnings = False)

    with pytest.raises(ArgumentException):
        output = build_rs_input_from_data(
            forcing_frequency,
            radius_array,
            density_array,
            static_bulk_modulus_array2,
            static_shear_modulus_array,
            bulk_viscosity_array,
            shear_viscosity_array,
            layer_upper_radius_tuple,
            layer_types,
            layer_is_statics,
            layer_is_incompressibles,
            shear_rheology_models,
            bulk_rheology_models,
            perform_checks = True,
            warnings = False)

    with pytest.raises(ArgumentException):
        output = build_rs_input_from_data(
            forcing_frequency,
            radius_array,
            density_array,
            static_bulk_modulus_array,
            static_shear_modulus_array2,
            bulk_viscosity_array,
            shear_viscosity_array,
            layer_upper_radius_tuple,
            layer_types,
            layer_is_statics,
            layer_is_incompressibles,
            shear_rheology_models,
            bulk_rheology_models,
            perform_checks = True,
            warnings = False)
    
    with pytest.raises(ArgumentException):
        output = build_rs_input_from_data(
            forcing_frequency,
            radius_array,
            density_array,
            static_bulk_modulus_array,
            static_shear_modulus_array,
            bulk_viscosity_array2,
            shear_viscosity_array,
            layer_upper_radius_tuple,
            layer_types,
            layer_is_statics,
            layer_is_incompressibles,
            shear_rheology_models,
            bulk_rheology_models,
            perform_checks = True,
            warnings = False)
    
    with pytest.raises(ArgumentException):
        output = build_rs_input_from_data(
            forcing_frequency,
            radius_array,
            density_array,
            static_bulk_modulus_array,
            static_shear_modulus_array,
            bulk_viscosity_array,
            shear_viscosity_array2,
            layer_upper_radius_tuple,
            layer_types,
            layer_is_statics,
            layer_is_incompressibles,
            shear_rheology_models,
            bulk_rheology_models,
            perform_checks = True,
            warnings = False)
    
    # Radius array that is not in ascending order
    radius_array3 = np.linspace(0.0, 50.0, N)
    tmp = radius_array3[50]
    radius_array3[50] = radius_array3[51]
    radius_array3[51] = tmp

    with pytest.raises(ArgumentException):
        output = build_rs_input_from_data(
            forcing_frequency,
            radius_array3,
            density_array,
            static_bulk_modulus_array,
            static_shear_modulus_array,
            bulk_viscosity_array,
            shear_viscosity_array2,
            layer_upper_radius_tuple,
            layer_types,
            layer_is_statics,
            layer_is_incompressibles,
            shear_rheology_models,
            bulk_rheology_models,
            perform_checks = True,
            warnings = False)
