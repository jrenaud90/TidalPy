# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.math cimport cbrt
from libcpp cimport bool as cpp_bool
from libcpp.vector cimport vector

import numpy as np
cimport numpy as cnp
cnp.import_array()

from collections import namedtuple

from TidalPy.constants cimport d_PI_DBL
from TidalPy.utilities.math.numerics cimport cf_isclose
from TidalPy.rheology.base cimport RheologyModelBase

cdef double PI_4_3 = (4.0 / 3.0) * d_PI_DBL 

from TidalPy.exceptions import ArgumentException
from TidalPy.logger import get_logger
log = get_logger('TidalPy')


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

def build_rs_input_homogeneous_layers(
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
    cdef cnp.ndarray[cnp.float64_t, ndim=1] upper_radius_array     = np.empty(num_layers, dtype=np.float64, order='C')
    cdef cnp.ndarray[cnp.float64_t, ndim=1] radius_array           = np.empty(total_slices, dtype=np.float64, order='C')
    cdef cnp.ndarray[cnp.float64_t, ndim=1] density_array          = np.empty(total_slices, dtype=np.float64, order='C')
    cdef vector[double] data_array_vec  = vector[double]()
    data_array_vec.resize(total_slices * 4)
    cdef cnp.ndarray[cnp.complex128_t, ndim=1] complex_shear_array = np.empty(total_slices, dtype=np.complex128, order='C')
    cdef cnp.ndarray[cnp.complex128_t, ndim=1] complex_bulk_array  = np.empty(total_slices, dtype=np.complex128, order='C')

    cdef double complex* complex_modulus_ptr = NULL

    cdef double[::1] upper_radius_view = upper_radius_array
    cdef double[::1] radius_view       = radius_array
    cdef double[::1] density_view      = density_array
    cdef double* shear_viscosity_ptr   = &data_array_vec[0]
    cdef double* bulk_viscosity_ptr    = &data_array_vec[total_slices]
    cdef double* shear_ptr             = &data_array_vec[2 * total_slices]
    cdef double* bulk_ptr              = &data_array_vec[3 * total_slices]
    cdef double complex[::1] complex_shear_view = complex_shear_array
    cdef double complex[::1] complex_bulk_view  = complex_bulk_array

    cdef RheologyModelBase rheo_instance

    cdef double planet_bulk_density = 0.0
    cdef double dr, layer_thickness
    cdef size_t slice_i, full_index
    cdef double layer_density, layer_static_bulk, layer_static_shear, layer_bulk_visc, layer_shear_visc

    cdef size_t first_slice_in_layer = 0
    cdef double last_layer_radius    = 0.0
    cdef double last_layer_radius3   = 0.0
    cdef double layer_radius3        = 0.0

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
        layer_radius3              = layer_radius**3
        planet_bulk_density       += layer_density * (layer_radius3 - last_layer_radius3)
        last_layer_radius3         = layer_radius3

        # Populate arrays
        for slice_i in range(layer_slices):
            full_index = first_slice_in_layer + slice_i
            # Fill in radius
            if slice_i == 0:
                radius_view[full_index] = last_layer_radius
            elif slice_i == (layer_slices - 1):
                radius_view[full_index] = last_layer_radius + layer_thickness
            else:
                radius_view[full_index] = last_layer_radius + slice_i * dr
            # Fill in other arrays (This function assumes these properties are constant throughout the layer)
            density_view[full_index]         = layer_density
            shear_viscosity_ptr[full_index]  = layer_shear_visc
            bulk_viscosity_ptr[full_index]   = layer_bulk_visc
            shear_ptr[full_index]            = layer_static_shear
            bulk_ptr[full_index]             = layer_static_bulk

        # Call on rheology class instance to find complex moduli.
        # Setup pointers used by rheology class
        complex_modulus_ptr = &complex_shear_view[first_slice_in_layer]
        rheo_instance = shear_rheology_model_tuple[layer_i]
        rheo_instance._vectorize_modulus_viscosity(
            forcing_frequency,
            &shear_ptr[first_slice_in_layer],
            &shear_viscosity_ptr[first_slice_in_layer],
            complex_modulus_ptr,  # Modified variable
            <Py_ssize_t>layer_slices)
        
        # Repeat for bulk
        complex_modulus_ptr = &complex_bulk_view[first_slice_in_layer]
        rheo_instance = bulk_rheology_model_tuple[layer_i]
        rheo_instance._vectorize_modulus_viscosity(
            forcing_frequency,
            &bulk_ptr[first_slice_in_layer],
            &bulk_viscosity_ptr[first_slice_in_layer],
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
        double[::1] radius_array,
        double[::1] density_array,
        double[::1] static_bulk_modulus_array,
        double[::1] static_shear_modulus_array,
        double[::1] bulk_viscosity_array,
        double[::1] shear_viscosity_array,
        tuple layer_upper_radius_tuple,
        tuple layer_type_tuple,
        tuple layer_is_static_tuple,
        tuple layer_is_incompressible_tuple,
        tuple shear_rheology_model_tuple,
        tuple bulk_rheology_model_tuple,
        cpp_bool perform_checks = True,
        cpp_bool warnings = True):

    # Pull out data
    cdef size_t num_slices_input = radius_array.size
    cdef double planet_radius = radius_array[num_slices_input - 1]
    
    # Check that the correct number of layers were provided for each input tuple
    cdef size_t slice_i
    cdef size_t layer_i
    cdef size_t num_layers = len(layer_upper_radius_tuple)
    cdef double r
    if perform_checks:
        if layer_upper_radius_tuple[num_layers - 1] != planet_radius:
            raise ArgumentException(f"Upper radius of last layer must be equal to planet's radius (last element of `radius_array`). Expected: {layer_upper_radius_tuple[num_layers - 1]:0.3e} Actual: {planet_radius:0.3e}.")
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
        
        if density_array.size != num_slices_input:
            raise ArgumentException("Unexpected length found for `density_array`. All input arrays must be the same length (equal to the size of radius array).")
        if static_bulk_modulus_array.size != num_slices_input:
            raise ArgumentException("Unexpected length found for `static_bulk_modulus_array`. All input arrays must be the same length (equal to the size of radius array).")
        if static_shear_modulus_array.size != num_slices_input:
            raise ArgumentException("Unexpected length found for `static_shear_modulus_array`. All input arrays must be the same length (equal to the size of radius array).")
        if bulk_viscosity_array.size != num_slices_input:
            raise ArgumentException("Unexpected length found for `bulk_viscosity_array`. All input arrays must be the same length (equal to the size of radius array).")
        if shear_viscosity_array.size != num_slices_input:
            raise ArgumentException("Unexpected length found for `shear_viscosity_array`. All input arrays must be the same length (equal to the size of radius array).")

        for slice_i in range(num_slices_input):
            if slice_i > 0:
                if r > radius_array[slice_i]:
                    raise ArgumentException("`radius_array` must be provided in ascending order.")
            r = radius_array[slice_i]

    # Make a copy of the upper radius array and store it as an array
    cdef cnp.ndarray[cnp.float64_t, ndim=1] upper_radius_array = np.empty(num_layers, dtype=np.float64, order='C')
    cdef double[::1] upper_radius_view = upper_radius_array

    for layer_i in range(num_layers):
        upper_radius_view[layer_i] = layer_upper_radius_tuple[layer_i]
    
    # We need to make copies of the arrays. We do not know what the final size of the new arrays will be though because
    #  we may find that we need to add values if they are missing. Instead of performing inserts which are very 
    #  expensive as they may cause memory reallocs.
    # Let's also make 1 giant data array to again avoid a lot of allocs / reallocs
    # Max num slices = current slices + 2 for each layer assuming we have to add one to the top and bottom of each.
    cdef size_t max_possible_num_slices = num_slices_input + 2 * num_layers
    cdef cnp.ndarray[cnp.float64_t, ndim=2] data_array = np.empty((6, max_possible_num_slices), dtype=np.float64, order='C')
    cdef double[:, ::1] data_array_view = data_array


    cdef double* radius_array_use_ptr          = &data_array_view[0, 0]
    cdef double* density_array_use_ptr         = &data_array_view[1, 0]
    cdef double* shear_viscosity_array_use_ptr = &data_array_view[2, 0]
    cdef double* bulk_viscosity_array_use_ptr  = &data_array_view[3, 0]
    cdef double* shear_array_use_ptr           = &data_array_view[4, 0]
    cdef double* bulk_array_use_ptr            = &data_array_view[5, 0]

    cdef vector[size_t] slices_per_layer = vector[size_t]()
    slices_per_layer.resize(num_layers)

    cdef size_t layer_slices = 0
    cdef double radius_slicei, density_slicei, shear_slicei, bulk_slicei, shear_visc_slicei, bulk_visc_slicei
    cdef double slice_volume = 0.0
    cdef double growing_mass = 0.0
    cdef cpp_bool r_present

    # Start performing corrections on input arrays

    # Update mass to find bulk density
    cdef double r3      = 0.0
    cdef double last_r3 = 0.0

    # Step through each slice in the planet using the old arrays as a baseline.
    cdef size_t slice_i_inputs  = 0
    cdef size_t slice_i_outputs = 0
    cdef double layer_upper_radius      = 0.0
    cdef double last_layer_upper_radius = 0.0

    for layer_i in range(num_layers):
        # Prepare to parse this layer
        layer_upper_radius = upper_radius_view[layer_i]
        layer_slices       = 0

        # Get values at bottom of layer
        radius_slicei     = radius_array[slice_i_inputs]
        density_slicei    = density_array[slice_i_inputs]
        shear_slicei      = static_shear_modulus_array[slice_i_inputs]
        bulk_slicei       = static_bulk_modulus_array[slice_i_inputs]
        shear_visc_slicei = shear_viscosity_array[slice_i_inputs]
        bulk_visc_slicei  = bulk_viscosity_array[slice_i_inputs]

        # Fix array's bottom-most values
        # Check that the last layer's upper radius is provided _again_ for the base of this layer.
        # When layer_i == 0 then last_layer_upper_radius == 0.0 which we want to check for anyways
        if not cf_isclose(radius_slicei, last_layer_upper_radius):
            # The base radius of this layer must be equal to the last layer's upper radius (for interfaces between layers there will be _two_ duplicate radius values)
            if warnings:
                if layer_i == 0:
                    log.warning("Radius array must start at zero, inserting r=0.0 at slice 0; other parameters will be set to their values at previous slice 0.")
                else:
                    log.warning(f"Layer {layer_i} starts at a different radius value than the previous layer's upper radius. Interface radius values must be provided twice, once for each layer above and below the interface. This radius value will be added; other parameters will be set to their values at the slice just after this upper radius value.")
        
            # Insert a new slice that uses the same physical properties that were previously at this slice.
            radius_array_use_ptr[slice_i_outputs]          = last_layer_upper_radius
            density_array_use_ptr[slice_i_outputs]         = density_slicei
            shear_viscosity_array_use_ptr[slice_i_outputs] = shear_visc_slicei
            bulk_viscosity_array_use_ptr[slice_i_outputs]  = bulk_visc_slicei
            shear_array_use_ptr[slice_i_outputs]           = shear_slicei
            bulk_array_use_ptr[slice_i_outputs]            = bulk_slicei

            # Output arrays have now grown by 1 due to insert
            slice_i_outputs += 1
            layer_slices    += 1

            # Update physical properties
            r3            = radius_slicei**3
            slice_volume  = PI_4_3 * (r3 - last_r3)
            growing_mass += slice_volume * density_slicei
            last_r3       = r3
        
        # Now step through the rest of this layer's slices and save the values. 
        r_present = False
        for slice_i in range(slice_i_inputs, num_slices_input):
            # Pull out input values
            radius_slicei     = radius_array[slice_i]
            density_slicei    = density_array[slice_i]
            shear_slicei      = static_shear_modulus_array[slice_i]
            bulk_slicei       = static_bulk_modulus_array[slice_i]
            shear_visc_slicei = shear_viscosity_array[slice_i]
            bulk_visc_slicei  = bulk_viscosity_array[slice_i]

            if cf_isclose(radius_slicei, layer_upper_radius):
                # This upper r is in the radius array, we are good!
                r_present = True
            elif radius_slicei > layer_upper_radius:
                # We have passed this layer's upper radius.
                break
            else:
                # Still in the layer
                pass

            # Add values of this slice to our output array
            radius_array_use_ptr[slice_i_outputs]          = radius_slicei
            density_array_use_ptr[slice_i_outputs]         = density_slicei
            shear_viscosity_array_use_ptr[slice_i_outputs] = shear_visc_slicei
            bulk_viscosity_array_use_ptr[slice_i_outputs]  = bulk_visc_slicei
            shear_array_use_ptr[slice_i_outputs]           = shear_slicei
            bulk_array_use_ptr[slice_i_outputs]            = bulk_slicei

            # Output arrays have now grown by 1 due to insert
            slice_i_outputs += 1
            layer_slices    += 1
            # Update the input array slice tracker using the number of slices in the layer
            slice_i_inputs  += 1

            # Update physical properties
            r3            = radius_slicei**3
            slice_volume  = PI_4_3 * (r3 - last_r3)
            growing_mass += slice_volume * density_slicei
            last_r3       = r3

            if r_present:
                break

        # Insert any values into the top of this layer
        if not r_present:
            # We need to add a radius value to the arrays at the top of this layer.
            if warnings:
                log.warning(f"Layer {layer_i} does not have its upper radius in the radius array which is required. It will be added; other parameters will be set to their values at the slice just before this upper radius value.")
            # Get other properties at slice just before where the upper radius should be added.
            density_slicei    = density_array[slice_i - 1]
            shear_slicei      = static_shear_modulus_array[slice_i - 1]
            bulk_slicei       = static_bulk_modulus_array[slice_i - 1]
            shear_visc_slicei = shear_viscosity_array[slice_i - 1]
            bulk_visc_slicei  = bulk_viscosity_array[slice_i - 1]

            radius_array_use_ptr[slice_i_outputs]          = layer_upper_radius
            density_array_use_ptr[slice_i_outputs]         = density_slicei
            shear_viscosity_array_use_ptr[slice_i_outputs] = shear_visc_slicei
            bulk_viscosity_array_use_ptr[slice_i_outputs]  = bulk_visc_slicei
            shear_array_use_ptr[slice_i_outputs]           = shear_slicei
            bulk_array_use_ptr[slice_i_outputs]            = bulk_slicei

            # Output arrays have now grown by 1 due to insert
            slice_i_outputs += 1
            layer_slices    += 1

            # Update physical properties
            r3            = radius_slicei**3
            slice_volume  = PI_4_3 * (r3 - last_r3)
            growing_mass += slice_volume * density_slicei
            last_r3       = r3

        # Check that there are enough slices for this layer
        if layer_slices < 5:
            raise ArgumentException(f"Layer {layer_i} has {layer_slices} slices when at least 5 are required.")

        # Prepare for next layer check
        slices_per_layer[layer_i] = layer_slices
        last_layer_upper_radius   = layer_upper_radius
    
    # Get other global parameters
    cdef double planet_bulk_density = growing_mass / (PI_4_3 * planet_radius**3)
    
    # Build other required arrays
    cdef cnp.ndarray[cnp.complex128_t, ndim=1] complex_shear_array = np.empty(slice_i_outputs, dtype=np.complex128, order='C')
    cdef cnp.ndarray[cnp.complex128_t, ndim=1] complex_bulk_array  = np.empty(slice_i_outputs, dtype=np.complex128, order='C')
    cdef double complex[::1] complex_shear_view = complex_shear_array
    cdef double complex[::1] complex_bulk_view  = complex_bulk_array
    cdef double complex* complex_modulus_ptr = NULL

    # Calculate the complex moduli for each layer
    cdef RheologyModelBase rheo_instance
    cdef size_t first_slice_in_layer_output = 0
    for layer_i in range(num_layers):
        layer_slices = slices_per_layer[layer_i]

        # Check if rheology classes are the correct instances
        if perform_checks:
            if not isinstance(shear_rheology_model_tuple[layer_i], RheologyModelBase):
                raise ArgumentException(f"Layer {layer_i} shear rheology class is not an instance of `RheologyModelBase`.")
            if not isinstance(bulk_rheology_model_tuple[layer_i], RheologyModelBase):
                raise ArgumentException(f"Layer {layer_i} shear rheology class is not an instance of `RheologyModelBase`.")

        # Call on rheology class instance to find complex moduli.
        # Setup pointers used by rheology class
        complex_modulus_ptr = &complex_shear_view[first_slice_in_layer_output]
        rheo_instance       = shear_rheology_model_tuple[layer_i]
        rheo_instance._vectorize_modulus_viscosity(
            forcing_frequency,
            &shear_array_use_ptr[first_slice_in_layer_output],
            &shear_viscosity_array_use_ptr[first_slice_in_layer_output],
            complex_modulus_ptr,  # Modified variable
            <Py_ssize_t>layer_slices)
        
        # Repeat for bulk
        complex_modulus_ptr = &complex_bulk_view[first_slice_in_layer_output]
        rheo_instance       = bulk_rheology_model_tuple[layer_i]
        rheo_instance._vectorize_modulus_viscosity(
            forcing_frequency,
            &bulk_array_use_ptr[first_slice_in_layer_output],
            &bulk_viscosity_array_use_ptr[first_slice_in_layer_output],
            complex_modulus_ptr,  # Modified variable
            <Py_ssize_t>layer_slices)
        
        # Prepare for next layer
        first_slice_in_layer_output += layer_slices
    
    # Create numpy array views of the data
    cdef cnp.ndarray[cnp.float64_t, ndim=1] radius_array_use  = np.array(data_array_view[0, 0:slice_i_outputs], copy=False)
    cdef cnp.ndarray[cnp.float64_t, ndim=1] density_array_use = np.array(data_array_view[1, 0:slice_i_outputs], copy=False)
    
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
