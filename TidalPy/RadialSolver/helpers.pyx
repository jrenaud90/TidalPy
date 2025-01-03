# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.math cimport cbrt
from libcpp cimport bool as cpp_bool
from libcpp.vector cimport vector

import numpy as np
cimport numpy as cnp
cnp.import_array()

from collections import namedtuple

from TidalPy.exceptions import ArgumentException
from TidalPy.constants cimport d_PI_DBL
from TidalPy.utilities.math.numerics cimport cf_isclose
from TidalPy.rheology.base cimport RheologyModelBase

import logging
log = logging.getLogger('TidalPy')


PlanetBuildData = namedtuple("PlanetBuildData", 
    (
        "radius_array",
        "density_array",
        "complex_bulk_modulus_array",
        "complex_shear_modulus_array",
        "frequency",
        "planet_bulk_density",
        "layer_types",
        "is_static_bylayer",
        "is_incompressible_bylayer",
        "upper_radius_bylayer_array"
    ))

def build_rs_input_homogenous_layers(
        double planet_radius,
        double forcing_frequency,
        tuple density_tuple,
        tuple static_bulk_modulus_tuple,
        tuple static_shear_modulus_tuple,
        tuple bulk_viscosity_tuple,
        tuple shear_viscosity_tuple,
        tuple layer_type_tuple,
        tuple layer_is_static_tuple,
        tuple layer_is_incompressible_tuple,
        tuple shear_rheology_model_tuple,
        tuple bulk_rheology_model_tuple,
        tuple radius_fraction_tuple = None,
        tuple thickness_fraction_tuple = None,
        tuple volume_fraction_tuple = None,
        tuple slices_tuple = None,
        size_t slice_per_layer = 10,
        cpp_bool perform_checks = True):
    
    # Optimizations
    cdef double R3 = planet_radius**3
    
    # Check that the correct number of layers were provided for each input tuple
    cdef size_t layer_i
    cdef size_t num_layers = len(density_tuple)
    if perform_checks:
        if len(static_bulk_modulus_tuple) != num_layers:
            raise ArgumentException("Unexpected length found for `static_bulk_modulus_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(static_shear_modulus_tuple) != num_layers:
            raise ArgumentException("Unexpected length found for `static_shear_modulus_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(bulk_viscosity_tuple) != num_layers:
            raise ArgumentException("Unexpected length found for `bulk_viscosity_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(shear_viscosity_tuple) != num_layers:
            raise ArgumentException("Unexpected length found for `shear_viscosity_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(layer_type_tuple) != num_layers:
            raise ArgumentException("Unexpected length found for `layer_type_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(layer_is_static_tuple) != num_layers:
            raise ArgumentException("Unexpected length found for `layer_is_static_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(layer_is_incompressible_tuple) != num_layers:
            raise ArgumentException("Unexpected length found for `layer_is_incompressible_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(shear_rheology_model_tuple) != num_layers:
            raise ArgumentException("Unexpected length found for `shear_rheology_model_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(bulk_rheology_model_tuple) != num_layers:
            raise ArgumentException("Unexpected length found for `bulk_rheology_model_tuple`. All input tuples must be the same length (equal to the number of layers).")

        if slices_tuple is not None:
            if len(slices_tuple) != num_layers:
                raise ArgumentException("Unexpected length found for `slices_tuple`. All input tuples must be the same length (equal to the number of layers).")

    # Convert to radius fractions
    cdef double thickness_frac, last_layer_frac, vol_frac, layer_radius
    cdef list thickness_fraction_list
    if thickness_fraction_tuple is None:
        if volume_fraction_tuple is None:
            if radius_fraction_tuple is None:
                raise ArgumentException("Must provide `thickness_fraction_tuple`, `radius_fraction_tuple`, or `volume_fraction_tuple`.")
            else:
                if perform_checks:
                    if len(radius_fraction_tuple) != num_layers:
                        raise ArgumentException("Unexpected length found for `radius_fraction_tuple`. All input tuples must be the same length (equal to the number of layers).")
                thickness_fraction_list = list()
                last_layer_frac = 0.0
                for layer_i in range(num_layers):
                    if perform_checks:
                        if radius_fraction_tuple[layer_i] <= last_layer_frac:
                            raise ArgumentException("`radius_fraction_tuple` entries must be ordered from smallest to largest with no repeated values.")

                    thickness_frac = (radius_fraction_tuple[layer_i] - last_layer_frac)
                    thickness_fraction_list.append(thickness_frac)
                    last_layer_frac = radius_fraction_tuple[layer_i]
                if perform_checks:
                    # Check that the last layer is at the planet surface
                    if not cf_isclose(last_layer_frac, 1.0):
                        raise ArgumentException(f"The last entry in `radius_fraction_tuple` must be equal to 1 for the top of planet (found {last_layer_frac}).")
        else:
            if perform_checks:
                if len(volume_fraction_tuple) != num_layers:
                    raise ArgumentException("Unexpected length found for `volume_fraction_tuple`. All input tuples must be the same length (equal to the number of layers).")
            thickness_fraction_list = list()
            radius_below = 0.0
            for layer_i in range(num_layers):
                vol_frac       = volume_fraction_tuple[layer_i]
                layer_radius   = cbrt(vol_frac*R3 + radius_below**3)
                thickness_frac = (layer_radius - radius_below) / planet_radius
                thickness_fraction_list.append(thickness_frac)
                radius_below = layer_radius
        thickness_fraction_tuple = tuple(thickness_fraction_list)
    elif perform_checks:
        if len(thickness_fraction_tuple) != num_layers:
                raise ArgumentException("Unexpected length found for `thickness_fraction_tuple`. All input tuples must be the same length (equal to the number of layers).")
    
    cdef double total_thick_frac = 0.0
    cdef size_t layer_slices = 0
    cdef size_t total_slices = 0
    for layer_i in range(num_layers):
        thickness_frac = thickness_fraction_tuple[layer_i]
        if thickness_frac <= 0.0:
            raise ArgumentException("Negative or zero radius fraction encountered.")
        total_thick_frac += thickness_frac
        
        if slices_tuple is not None:
            layer_slices = slices_tuple[layer_i]
        else:
            layer_slices = slice_per_layer
        if layer_slices < 5:
            raise ArgumentException(f"Layer {layer_i} has {layer_slices} slices when at least 5 are required.")
        total_slices += layer_slices

        # Check if rheology classes are the correct instances
        if perform_checks:
            if not isinstance(shear_rheology_model_tuple[layer_i], RheologyModelBase):
                raise ArgumentException(f"Layer {layer_i} shear rheology class is not an instance of `RheologyModelBase`.")
            if not isinstance(bulk_rheology_model_tuple[layer_i], RheologyModelBase):
                raise ArgumentException(f"Layer {layer_i} shear rheology class is not an instance of `RheologyModelBase`.")
    
    if not cf_isclose(total_thick_frac, 1.0):
        raise ValueError(f"Unexpected value found for total radius fraction, {total_thick_frac} (expected 1.0).")
    
    # Build required arrays
    cdef cnp.ndarray[cnp.float64_t, ndim=1] upper_radius_array     = np.full(num_layers, np.nan, dtype=np.float64, order='C')
    cdef cnp.ndarray[cnp.float64_t, ndim=1] radius_array           = np.full(total_slices, np.nan, dtype=np.float64, order='C')
    cdef cnp.ndarray[cnp.float64_t, ndim=1] density_array          = np.full(total_slices, np.nan, dtype=np.float64, order='C')
    cdef cnp.ndarray[cnp.float64_t, ndim=1] shear_viscosity_array  = np.full(total_slices, np.nan, dtype=np.float64, order='C')
    cdef cnp.ndarray[cnp.float64_t, ndim=1] bulk_viscosity_array   = np.full(total_slices, np.nan, dtype=np.float64, order='C')
    cdef cnp.ndarray[cnp.float64_t, ndim=1] shear_array            = np.full(total_slices, np.nan, dtype=np.float64, order='C')
    cdef cnp.ndarray[cnp.float64_t, ndim=1] bulk_array             = np.full(total_slices, np.nan, dtype=np.float64, order='C')
    cdef cnp.ndarray[cnp.complex128_t, ndim=1] complex_shear_array = np.full(total_slices, np.nan, dtype=np.complex128, order='C')
    cdef cnp.ndarray[cnp.complex128_t, ndim=1] complex_bulk_array  = np.full(total_slices, np.nan, dtype=np.complex128, order='C')

    cdef double* modulus_ptr                 = NULL
    cdef double* viscosity_ptr               = NULL
    cdef double complex* complex_modulus_ptr = NULL

    cdef double[::1] upper_radius_view    = upper_radius_array
    cdef double[::1] radius_view          = radius_array
    cdef double[::1] density_view         = density_array
    cdef double[::1] shear_viscosity_view = shear_viscosity_array
    cdef double[::1] bulk_viscosity_view  = bulk_viscosity_array
    cdef double[::1] shear_view           = shear_array
    cdef double[::1] bulk_view            = bulk_array
    cdef double complex[::1] complex_shear_view = complex_shear_array
    cdef double complex[::1] complex_bulk_view  = complex_bulk_array

    cdef RheologyModelBase rheo_instance

    cdef double planet_bulk_density = 0.0
    cdef double dr, layer_thickness
    cdef size_t slice_i, full_index
    cdef double layer_density, layer_static_bulk, layer_static_shear, layer_bulk_visc, layer_shear_visc

    cdef size_t first_slice_in_layer = 0
    cdef double last_layer_radius    = 0.0
    for layer_i in range(num_layers):
        if slices_tuple:
            layer_slices = slices_tuple[layer_i]
        else:
            layer_slices = slice_per_layer
        thickness_frac = thickness_fraction_tuple[layer_i]

        # Pull out other data
        layer_density      = density_tuple[layer_i]
        layer_static_bulk  = static_bulk_modulus_tuple[layer_i]
        layer_static_shear = static_shear_modulus_tuple[layer_i]
        layer_bulk_visc    = bulk_viscosity_tuple[layer_i]
        layer_shear_visc   = shear_viscosity_tuple[layer_i]

        # Build linspace radius array from bottom of layer to top, inclusive of both r-bot and r-top
        layer_thickness = thickness_frac * planet_radius
        layer_radius    = last_layer_radius + layer_thickness
        # The "- 1" in dr ensures that we build up to be inclusive of the top of the layer.
        # Earlier we ensured that layer_slices >= 5 so no worries about a divide by zero
        dr = layer_thickness / <double>(layer_slices - 1)

        # Other calculations
        upper_radius_view[layer_i] = layer_radius
        planet_bulk_density       += layer_density * (layer_radius**3 - last_layer_radius**3)

        # Populate arrays
        for slice_i in range(layer_slices):
            full_index = first_slice_in_layer + slice_i
            # Fill in radius
            radius_view[full_index] = last_layer_radius + slice_i * dr
            # Fill in other arrays (This function assumes these properties are constant throughout the layer)
            density_view[full_index]         = layer_density
            shear_viscosity_view[full_index] = layer_shear_visc
            bulk_viscosity_view[full_index]  = layer_bulk_visc
            shear_view[full_index]           = layer_static_shear
            bulk_view[full_index]            = layer_static_bulk

        # Call on rheology class instance to find complex moduli.
        # Setup pointers used by rheology class
        modulus_ptr         = &shear_view[first_slice_in_layer]
        viscosity_ptr       = &shear_viscosity_view[first_slice_in_layer]
        complex_modulus_ptr = &complex_shear_view[first_slice_in_layer]
        rheo_instance = shear_rheology_model_tuple[layer_i]
        rheo_instance._vectorize_modulus_viscosity(
            forcing_frequency,
            modulus_ptr,
            viscosity_ptr,
            complex_modulus_ptr,  # Modified variable
            <Py_ssize_t>layer_slices)
        
        # Repeat for bulk
        modulus_ptr         = &bulk_view[first_slice_in_layer]
        viscosity_ptr       = &bulk_viscosity_view[first_slice_in_layer]
        complex_modulus_ptr = &complex_bulk_view[first_slice_in_layer]
        rheo_instance = bulk_rheology_model_tuple[layer_i]
        rheo_instance._vectorize_modulus_viscosity(
            forcing_frequency,
            modulus_ptr,
            viscosity_ptr,
            complex_modulus_ptr,  # Modified variable
            <Py_ssize_t>layer_slices)

        # Prepare for next layer
        first_slice_in_layer += layer_slices
        last_layer_radius     = layer_radius
        
    # Finalize other parameters
    planet_bulk_density /= R3 
    
    # Build outputs
    rs_outputs = PlanetBuildData(
        radius_array,
        density_array,
        complex_bulk_array,
        complex_shear_array,
        forcing_frequency,
        planet_bulk_density,
        layer_type_tuple,
        layer_is_static_tuple,
        layer_is_incompressible_tuple,
        upper_radius_array
    )

    return rs_outputs


def build_rs_input_from_data(
        double forcing_frequency,
        cnp.ndarray[cnp.float64_t, ndim=1] radius_array,
        cnp.ndarray[cnp.float64_t, ndim=1] density_array,
        cnp.ndarray[cnp.float64_t, ndim=1] static_bulk_modulus_array,
        cnp.ndarray[cnp.float64_t, ndim=1] static_shear_modulus_array,
        cnp.ndarray[cnp.float64_t, ndim=1] bulk_viscosity_array,
        cnp.ndarray[cnp.float64_t, ndim=1] shear_viscosity_array,
        tuple layer_upper_radius_tuple,
        tuple layer_type_tuple,
        tuple layer_is_static_tuple,
        tuple layer_is_incompressible_tuple,
        tuple shear_rheology_model_tuple,
        tuple bulk_rheology_model_tuple,
        cpp_bool perform_checks = True,
        cpp_bool warnings = True):

    # Pull out data
    cdef size_t initial_num_slices = radius_array.size
    cdef size_t current_num_slices = initial_num_slices
    cdef double planet_radius = radius_array[initial_num_slices - 1]
    
    # Check that the correct number of layers were provided for each input tuple
    cdef size_t layer_i
    cdef size_t num_layers = len(layer_upper_radius_tuple)
    cdef double r
    if perform_checks:
        if layer_upper_radius_tuple[num_layers - 1] != planet_radius:
            raise ArgumentException("Upper radius of last layer must be equal to planet's radius (last element of `radius_array`).")
        if len(layer_type_tuple) != num_layers:
            raise ArgumentException("Unexpected length found for `layer_type_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(layer_is_static_tuple) != num_layers:
            raise ArgumentException("Unexpected length found for `layer_is_static_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(layer_is_incompressible_tuple) != num_layers:
            raise ArgumentException("Unexpected length found for `layer_is_incompressible_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(shear_rheology_model_tuple) != num_layers:
            raise ArgumentException("Unexpected length found for `shear_rheology_model_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(bulk_rheology_model_tuple) != num_layers:
            raise ArgumentException("Unexpected length found for `bulk_rheology_model_tuple`. All input tuples must be the same length (equal to the number of layers).")
        
        if len(density_array) != initial_num_slices:
            raise ArgumentException("Unexpected length found for `density_array`. All input arrays must be the same length (equal to the size of radius array).")
        if len(static_bulk_modulus_array) != initial_num_slices:
            raise ArgumentException("Unexpected length found for `static_bulk_modulus_array`. All input arrays must be the same length (equal to the size of radius array).")
        if len(static_shear_modulus_array) != initial_num_slices:
            raise ArgumentException("Unexpected length found for `static_shear_modulus_array`. All input arrays must be the same length (equal to the size of radius array).")
        if len(bulk_viscosity_array) != initial_num_slices:
            raise ArgumentException("Unexpected length found for `bulk_viscosity_array`. All input arrays must be the same length (equal to the size of radius array).")
        if len(shear_viscosity_array) != initial_num_slices:
            raise ArgumentException("Unexpected length found for `shear_viscosity_array`. All input arrays must be the same length (equal to the size of radius array).")

        for slice_i in range(initial_num_slices):
            if slice_i > 0:
                if r > radius_array[slice_i]:
                    raise ArgumentException("`radius_array` must be provided in ascending order.")
            r = radius_array[slice_i]

    # Start performing corrections on input arrays
    cdef cnp.ndarray[cnp.float64_t, ndim=1] upper_radius_array        = np.full(num_layers, np.nan, dtype=np.float64, order='C')
    cdef cnp.ndarray[cnp.float64_t, ndim=1] radius_array_use          = np.copy(radius_array)
    cdef cnp.ndarray[cnp.float64_t, ndim=1] density_array_use         = np.copy(density_array)
    cdef cnp.ndarray[cnp.float64_t, ndim=1] shear_viscosity_array_use = np.copy(shear_viscosity_array)
    cdef cnp.ndarray[cnp.float64_t, ndim=1] bulk_viscosity_array_use  = np.copy(bulk_viscosity_array)
    cdef cnp.ndarray[cnp.float64_t, ndim=1] shear_array_use           = np.copy(static_shear_modulus_array)
    cdef cnp.ndarray[cnp.float64_t, ndim=1] bulk_array_use            = np.copy(static_bulk_modulus_array)
    cdef double[::1] upper_radius_view    = upper_radius_array

    cdef vector[size_t] slices_per_layer = vector[size_t]()
    slices_per_layer.resize(num_layers)

    cdef double layer_upper_radius
    cdef double last_layer_upper_radius = 0.0
    cdef size_t first_slice_in_layer = 0
    cdef size_t layer_slices = 0
    cdef size_t insert_at_slice_i = 0
    cdef double tmp_density, tmp_shear, tmp_bulk, tmp_shear_visc, tmp_bulk_visc
    cdef double slice_volume = 0.0
    cdef double growing_mass = 0.0
    cdef double last_r3, r3
    cdef cpp_bool r_present
    for layer_i in range(num_layers):
        layer_upper_radius = layer_upper_radius_tuple[layer_i]
        upper_radius_view[layer_i] = layer_upper_radius

        # Check that this radius is provided in the radius array, twice at layer interfaces
        if layer_i == 0:
            # In bottom-most layer, check that the radius array starts at 0
            if not cf_isclose(radius_array_use[0], 0.0):
                if warnings:
                    log.warning("Radius array must start at zero, inserting r=0.0 at slice 0; other parameters will be set to their values at previous slice 0.")
                # Get other properties at previous slice 0
                tmp_density    = density_array_use[0]
                tmp_shear_visc = shear_viscosity_array_use[0]
                tmp_bulk_visc  = bulk_viscosity_array_use[0]
                tmp_shear      = shear_array_use[0]
                tmp_bulk       = bulk_array_use[0]

                # Create copies of arrays with inserted values
                radius_array_use          = np.insert(radius_array_use, 0, 0.0)
                density_array_use         = np.insert(density_array_use, 0, tmp_density)
                shear_viscosity_array_use = np.insert(shear_viscosity_array_use, 0, tmp_shear_visc)
                bulk_viscosity_array_use  = np.insert(bulk_viscosity_array_use, 0, tmp_bulk_visc)
                shear_array_use           = np.insert(shear_array_use, 0, tmp_shear)
                bulk_array_use            = np.insert(bulk_array_use, 0, tmp_bulk)
                
                # Arrays have now grown by 1
                current_num_slices += 1
        else:
            # Check that the last layer's upper radius is provided _again_ for the base of this layer.
            if not cf_isclose(radius_array_use[first_slice_in_layer], last_layer_upper_radius):
                # The base radius of this layer must be equal to the last layer's upper radius (for interfaces between layers there will be _two_ duplicate radius values)
                if warnings:
                    log.warning(f"Layer {layer_i} starts at a different radius value than the previous layer's upper radius. Interface radius values must be provided twice, once for each layer above and below the interface. This radius value will be added; other parameters will be set to their values at the slice just after this upper radius value.")
                # Get other properties at previous slice 0
                tmp_density    = density_array_use[first_slice_in_layer]
                tmp_shear_visc = shear_viscosity_array_use[first_slice_in_layer]
                tmp_bulk_visc  = bulk_viscosity_array_use[first_slice_in_layer]
                tmp_shear      = shear_array_use[first_slice_in_layer]
                tmp_bulk       = bulk_array_use[first_slice_in_layer]

                # Create copies of arrays with inserted values
                radius_array_use          = np.insert(radius_array_use, first_slice_in_layer, last_layer_upper_radius)
                density_array_use         = np.insert(density_array_use, first_slice_in_layer, tmp_density)
                shear_viscosity_array_use = np.insert(shear_viscosity_array_use, first_slice_in_layer, tmp_shear_visc)
                bulk_viscosity_array_use  = np.insert(bulk_viscosity_array_use, first_slice_in_layer, tmp_bulk_visc)
                shear_array_use           = np.insert(shear_array_use, first_slice_in_layer, tmp_shear)
                bulk_array_use            = np.insert(bulk_array_use, first_slice_in_layer, tmp_bulk)
                
                # Arrays have now grown by 1
                current_num_slices += 1
        
        # Check that the upper radius is provided in the radius array for this layer.
        r_present = False
        for slice_i in range(first_slice_in_layer, current_num_slices):
            r = radius_array_use[slice_i]
            if cf_isclose(r, layer_upper_radius):
                # This upper r is in the radius array, we are good!
                r_present = True
                # We still want to count this upper radius as part of this layer, so increment the layer slices
                layer_slices += 1
                break
            elif r > layer_upper_radius:
                # We have passed this layer's upper radius.
                break
            else:
                # Still in the layer
                layer_slices += 1

        if not r_present:
            # We need to add a radius value to the arrays at the top of this layer.
            if warnings:
                log.warning(f"Layer {layer_i} does not have its upper radius in the radius array which is required. It will be added; other parameters will be set to their values at the slice just before this upper radius value.")
            # Get other properties at slice just before where the upper radius should be added.
            insert_at_slice_i = first_slice_in_layer + layer_slices
            tmp_density       = density_array_use[insert_at_slice_i - 1]
            tmp_shear_visc    = shear_viscosity_array_use[insert_at_slice_i - 1]
            tmp_bulk_visc     = bulk_viscosity_array_use[insert_at_slice_i - 1]
            tmp_shear         = shear_array_use[insert_at_slice_i - 1]
            tmp_bulk          = bulk_array_use[insert_at_slice_i - 1]

            # Create copies of arrays with inserted values
            radius_array_use          = np.insert(radius_array_use, insert_at_slice_i, layer_upper_radius)
            density_array_use         = np.insert(density_array_use, insert_at_slice_i, tmp_density)
            shear_viscosity_array_use = np.insert(shear_viscosity_array_use, insert_at_slice_i, tmp_shear_visc)
            bulk_viscosity_array_use  = np.insert(bulk_viscosity_array_use, insert_at_slice_i, tmp_bulk_visc)
            shear_array_use           = np.insert(shear_array_use, insert_at_slice_i, tmp_shear)
            bulk_array_use            = np.insert(bulk_array_use, insert_at_slice_i, tmp_bulk)
            
            # Arrays have now grown by 1
            current_num_slices += 1
            layer_slices += 1

        # Check that there are enough slices for this layer
        if layer_slices < 5:
            raise ArgumentException(f"Layer {layer_i} has {layer_slices} slices when at least 5 are required.")
        
        # Perform other calculations for this layer
        if first_slice_in_layer == 0:
            last_r3 = 0.0
        else:
            last_r3 = radius_array_use[first_slice_in_layer - 1]**3

        for slice_i in range(first_slice_in_layer, first_slice_in_layer + layer_slices):
            r3            = radius_array_use[slice_i]**3
            tmp_density   = density_array_use[slice_i]
            slice_volume  = (4.0 / 3.0) * d_PI_DBL * (r3 - last_r3)
            growing_mass += slice_volume * tmp_density
            last_r3       = r3

        # Prepare for next layer check
        slices_per_layer[layer_i] = layer_slices
        last_layer_upper_radius   = layer_upper_radius
        first_slice_in_layer     += layer_slices
        layer_slices              = 0

        # Check if rheology classes are the correct instances
        if perform_checks:
            if not isinstance(shear_rheology_model_tuple[layer_i], RheologyModelBase):
                raise ArgumentException(f"Layer {layer_i} shear rheology class is not an instance of `RheologyModelBase`.")
            if not isinstance(bulk_rheology_model_tuple[layer_i], RheologyModelBase):
                raise ArgumentException(f"Layer {layer_i} shear rheology class is not an instance of `RheologyModelBase`.")
    
    # Get other global parameters
    cdef double planet_bulk_density = growing_mass / ((4. / 3.) * d_PI_DBL * planet_radius**3)
    
    # Build other required arrays
    cdef cnp.ndarray[cnp.complex128_t, ndim=1] complex_shear_array = np.full(current_num_slices, np.nan, dtype=np.complex128, order='C')
    cdef cnp.ndarray[cnp.complex128_t, ndim=1] complex_bulk_array  = np.full(current_num_slices, np.nan, dtype=np.complex128, order='C')

    cdef double[::1] shear_viscosity_view = shear_viscosity_array_use
    cdef double[::1] bulk_viscosity_view  = bulk_viscosity_array_use
    cdef double[::1] shear_view           = shear_array_use
    cdef double[::1] bulk_view            = bulk_array_use
    cdef double complex[::1] complex_shear_view = complex_shear_array
    cdef double complex[::1] complex_bulk_view  = complex_bulk_array

    cdef double* modulus_ptr                 = NULL
    cdef double* viscosity_ptr               = NULL
    cdef double complex* complex_modulus_ptr = NULL

    # Calculate the complex moduli for each layer
    cdef RheologyModelBase rheo_instance
    first_slice_in_layer = 0
    for layer_i in range(num_layers):
        layer_slices = slices_per_layer[layer_i]

        # Call on rheology class instance to find complex moduli.
        # Setup pointers used by rheology class
        modulus_ptr         = &shear_view[first_slice_in_layer]
        viscosity_ptr       = &shear_viscosity_view[first_slice_in_layer]
        complex_modulus_ptr = &complex_shear_view[first_slice_in_layer]
        rheo_instance       = shear_rheology_model_tuple[layer_i]
        rheo_instance._vectorize_modulus_viscosity(
            forcing_frequency,
            modulus_ptr,
            viscosity_ptr,
            complex_modulus_ptr,  # Modified variable
            <Py_ssize_t>layer_slices)
        
        # Repeat for bulk
        modulus_ptr         = &bulk_view[first_slice_in_layer]
        viscosity_ptr       = &bulk_viscosity_view[first_slice_in_layer]
        complex_modulus_ptr = &complex_bulk_view[first_slice_in_layer]
        rheo_instance       = bulk_rheology_model_tuple[layer_i]
        rheo_instance._vectorize_modulus_viscosity(
            forcing_frequency,
            modulus_ptr,
            viscosity_ptr,
            complex_modulus_ptr,  # Modified variable
            <Py_ssize_t>layer_slices)
        
        # Prepare for next layer
        first_slice_in_layer += layer_slices
    
    # Build outputs
    rs_outputs = PlanetBuildData(
        radius_array_use,
        density_array_use,
        complex_bulk_array,
        complex_shear_array,
        forcing_frequency,
        planet_bulk_density,
        layer_type_tuple,
        layer_is_static_tuple,
        layer_is_incompressible_tuple,
        upper_radius_array
    )

    return rs_outputs
