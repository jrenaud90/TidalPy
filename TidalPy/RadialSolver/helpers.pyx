# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.math cimport cbrt
from libcpp cimport bool as cpp_bool

import numpy as np
cimport numpy as cnp
cnp.import_array()

from TidalPy.utilities.math.numerics cimport cf_isclose
from TidalPy.rheology.base cimport RheologyModelBase

def build_planet_constant_layers(
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
    
    # Convert to radius fractions
    cdef size_t layer_i
    cdef size_t num_layers = len(density_tuple)
    if perform_checks:
        if len(static_bulk_modulus_tuple) != num_layers:
            raise AttributeError("Unexpected length found for `static_bulk_modulus_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(static_shear_modulus_tuple) != num_layers:
            raise AttributeError("Unexpected length found for `static_shear_modulus_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(bulk_viscosity_tuple) != num_layers:
            raise AttributeError("Unexpected length found for `bulk_viscosity_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(shear_viscosity_tuple) != num_layers:
            raise AttributeError("Unexpected length found for `shear_viscosity_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(layer_type_tuple) != num_layers:
            raise AttributeError("Unexpected length found for `layer_type_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(layer_is_static_tuple) != num_layers:
            raise AttributeError("Unexpected length found for `layer_is_static_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(layer_is_incompressible_tuple) != num_layers:
            raise AttributeError("Unexpected length found for `layer_is_incompressible_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(shear_rheology_model_tuple) != num_layers:
            raise AttributeError("Unexpected length found for `shear_rheology_model_tuple`. All input tuples must be the same length (equal to the number of layers).")
        if len(bulk_rheology_model_tuple) != num_layers:
            raise AttributeError("Unexpected length found for `bulk_rheology_model_tuple`. All input tuples must be the same length (equal to the number of layers).")

        if slices_tuple is not None:
            if len(slices_tuple) != num_layers:
                raise AttributeError("Unexpected length found for `slices_tuple`. All input tuples must be the same length (equal to the number of layers).")

    cdef double thickness_frac, last_layer_frac, vol_frac, cross_term, layer_radius
    cdef list thickness_fraction_list
    if thickness_fraction_tuple is None:
        if volume_fraction_tuple is None:
            if radius_fraction_tuple is None:
                raise AttributeError("Must provide `thickness_fraction_tuple`, `radius_fraction_tuple`, or `volume_fraction_tuple`.")
            else:
                if perform_checks:
                    if len(radius_fraction_tuple) != num_layers:
                        raise AttributeError("Unexpected length found for `radius_fraction_tuple`. All input tuples must be the same length (equal to the number of layers).")
                thickness_fraction_list = list()
                last_layer_frac = 0.0
                for layer_i in range(num_layers):
                    if perform_checks:
                        if radius_fraction_tuple[layer_i] <= last_layer_frac:
                            raise AttributeError("`radius_fraction_tuple` entries must be ordered from smallest to largest with no repeated values.")

                    thickness_frac = (radius_fraction_tuple[layer_i] - last_layer_frac)
                    thickness_fraction_list.append(thickness_frac)
                    last_layer_frac = radius_fraction_tuple[layer_i]
                if perform_checks:
                    # Check that the last layer is at the planet surface
                    if not cf_isclose(last_layer_frac, 1.0):
                        raise AttributeError(f"The last entry in `radius_fraction_tuple` must be equal to 1 for the top of planet (found {last_layer_frac}).")
        else:
            if perform_checks:
                if len(volume_fraction_tuple) != num_layers:
                    raise AttributeError("Unexpected length found for `volume_fraction_tuple`. All input tuples must be the same length (equal to the number of layers).")
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
                raise AttributeError("Unexpected length found for `thickness_fraction_tuple`. All input tuples must be the same length (equal to the number of layers).")
    
    cdef double total_thick_frac = 0.0
    cdef size_t layer_slices = 0
    cdef size_t total_slices = 0
    for layer_i in range(num_layers):
        thickness_frac = thickness_fraction_tuple[layer_i]
        if thickness_frac <= 0.0:
            raise AttributeError("Negative or zero radius fraction encountered.")
        total_thick_frac += thickness_frac
        
        if slices_tuple is not None:
            layer_slices = slices_tuple[layer_i]
        else:
            layer_slices = slice_per_layer
        if layer_slices < 5:
            raise AttributeError(f"Layer {layer_i} has {layer_slices} slices when at least 5 are required.")
        total_slices += layer_slices

        # Check if rheology classes are the correct instances
        if perform_checks:
            if not isinstance(shear_rheology_model_tuple[layer_i], RheologyModelBase):
                raise AttributeError(f"Layer {layer_i} shear rheology class is not an instance of `RheologyModelBase`.")
            if not isinstance(bulk_rheology_model_tuple[layer_i], RheologyModelBase):
                raise AttributeError(f"Layer {layer_i} shear rheology class is not an instance of `RheologyModelBase`.")
    
    if not cf_isclose(total_thick_frac, 1.0):
        raise ValueError(f"Unexpected value found for total radius fraction, {total_thick_frac} (expected 1.0).")
    
    # Build required arrays
    cdef cnp.ndarray upper_radius_array    = np.nan * np.ones(num_layers, dtype=np.float64, order='C')
    cdef cnp.ndarray radius_array          = np.nan * np.ones(total_slices, dtype=np.float64, order='C')
    cdef cnp.ndarray density_array         = np.nan * np.ones(total_slices, dtype=np.float64, order='C')
    cdef cnp.ndarray shear_viscosity_array = np.nan * np.ones(total_slices, dtype=np.float64, order='C')
    cdef cnp.ndarray bulk_viscosity_array  = np.nan * np.ones(total_slices, dtype=np.float64, order='C')
    cdef cnp.ndarray shear_array           = np.nan * np.ones(total_slices, dtype=np.float64, order='C')
    cdef cnp.ndarray bulk_array            = np.nan * np.ones(total_slices, dtype=np.float64, order='C')
    cdef cnp.ndarray complex_shear_array   = np.nan * np.ones(total_slices, dtype=np.complex128, order='C')
    cdef cnp.ndarray complex_bulk_array    = np.nan * np.ones(total_slices, dtype=np.complex128, order='C')

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

    cdef size_t last_layer_top_slice = 0
    cdef double last_layer_radius = 0.0
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
        planet_bulk_density += layer_density * (layer_radius**3 - last_layer_radius**3)

        # Populate arrays
        for slice_i in range(layer_slices):
            full_index = last_layer_top_slice + slice_i
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
        modulus_ptr         = &shear_view[last_layer_top_slice]
        viscosity_ptr       = &shear_viscosity_view[last_layer_top_slice]
        complex_modulus_ptr = &complex_shear_view[last_layer_top_slice]
        rheo_instance = shear_rheology_model_tuple[layer_i]
        rheo_instance._vectorize_modulus_viscosity(
            forcing_frequency,
            modulus_ptr,
            viscosity_ptr,
            complex_modulus_ptr,  # Modified variable
            <Py_ssize_t>layer_slices)
        
        # Repeat for bulk
        modulus_ptr         = &bulk_view[last_layer_top_slice]
        viscosity_ptr       = &bulk_viscosity_view[last_layer_top_slice]
        complex_modulus_ptr = &complex_bulk_view[last_layer_top_slice]
        rheo_instance = bulk_rheology_model_tuple[layer_i]
        rheo_instance._vectorize_modulus_viscosity(
            forcing_frequency,
            modulus_ptr,
            viscosity_ptr,
            complex_modulus_ptr,  # Modified variable
            <Py_ssize_t>layer_slices)

        # Prepare for next layer
        last_layer_top_slice += layer_slices
        last_layer_radius = layer_radius
        
    # Finalize other parameters
    planet_bulk_density /= R3 
    
    # Build outputs
    cdef tuple rs_outputs = (
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
