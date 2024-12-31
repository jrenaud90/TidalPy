# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

cdef void cf_collapse_layer_solution(
        double complex* solution_ptr,
        double complex* constant_vector_ptr,
        double complex** storage_by_solution,
        double* layer_radius_ptr,
        double* layer_density_ptr,
        double* layer_gravity_ptr,
        double frequency_to_use,
        size_t layer_start_index,
        size_t num_layer_slices,
        size_t num_sols,
        size_t max_num_y,
        size_t num_ys,
        size_t num_output_ys,
        size_t ytype_i,
        int layer_type,
        cpp_bool layer_is_static,
        cpp_bool layer_is_incomp
        ) noexcept nogil:
    """
    cf_collapse_layer_solution

    Compute the full set of radial functions (y-values) for a given layer of a planetary model by collapsing all solutions 
    across multiple slices and solutions within the layer.

    Parameters
    ----------
    solution_ptr : double complex*
        Pointer to an array where the final computed solutions will be stored. The array must have sufficient space to 
        accommodate all solutions across all slices and y-values.
    constant_vector_ptr : double complex*
        Pointer to an array of constants used to scale the solutions from each independent solution vector.
    storage_by_solution : double complex**
        Array of pointers to arrays that store precomputed intermediate solutions for each slice in the layer.
    layer_radius_ptr : double*
        Pointer to an array containing the radii for each slice in the layer.
    layer_density_ptr : double*
        Pointer to an array containing the density for each slice in the layer.
    layer_gravity_ptr : double*
        Pointer to an array containing the gravitational acceleration for each slice in the layer.
    frequency_to_use : double
        Frequency parameter used in the calculation of \( y_3 \) for dynamic liquid layers.
    layer_start_index : size_t
        Starting index for the layer's slices in the storage arrays.
    num_layer_slices : size_t
        Number of slices within the layer.
    num_sols : size_t
        Number of independent solutions to be collapsed for each y-value.
    max_num_y : size_t
        Maximum number of y-values supported by the solution.
    num_ys : size_t
        Number of y-values defined for the current layer type.
    num_output_ys : size_t
        Number of y-values in the output array.
    ytype_i : size_t
        Index specifying the type of y-value being processed.
    layer_type : int
        Type of layer being processed:
            - 0: Solid layer
            - 1: Liquid layer
    layer_is_static : cpp_bool
        Indicates whether the layer is static.
    layer_is_incomp : cpp_bool
        Indicates whether the layer is incompressible.
        Currently, this parameter is not directly used in this function.

    Returns
    -------
    None
        The function directly modifies `solution_ptr` to store the computed solutions.

    Notes
    -----
    - **Solid Layers**:
        - All \( y \)-values are computed directly from the provided solutions and constants.
    - **Static Liquid Layers**:
        - Only \( y_5 \) is defined; all other \( y \)-values are set to `NaN`.
    - **Dynamic Liquid Layers**:
        - \( y_1, y_2, y_5, y_6 \) are computed directly.
        - \( y_3 \) is calculated after the main loop using the layer radius, gravity, and density.
    - All y-values not defined for a particular layer type are set to `NaN`.

    Raises
    ------
    None
        This function does not explicitly raise exceptions but operates under the assumption that input parameters are valid.
    """

    # Use constant vectors to find the full y from all of the solutions in this layer
    cdef size_t solution_i
    cdef size_t y_i
    cdef size_t y_rhs_i
    cdef size_t lhs_y_index
    cdef size_t slice_i
    cdef size_t slice_i_shifted
    cdef size_t slice_end
    cdef cpp_bool calculate_y3

    slice_end = layer_start_index + num_layer_slices
    calculate_y3 = False

    if (not (layer_type == 0)) and (not layer_is_static):
        # For this layer type we can also calculate y3 based on the other ys. We will do this at the end.
        calculate_y3 = True
    
    y_rhs_i = 0
    y_i = 0
    while max_num_y > y_i:
        lhs_y_index = ytype_i * max_num_y + y_i
        # Bail out early for ys that the layer type is not defined for.
        if (layer_type == 0):
            # Solid layers have the same number of ys as the max_y_num.
            y_rhs_i = y_i
        else:
            # Liquid layers have less ys than max_y_num. Set y_rhs_i to the appropriate index.
            if layer_is_static:
                # Liquid static layers only has y5 (stored at index 0).
                if y_i == 4:
                    y_rhs_i = 0
                else:
                    # All results should be NAN which the solution pointer should already be set to. Continue.
                    y_i += 1
                    continue
            else:
                # Liquid dynamic layers have y1, y2, y5, y6 (indices 0, 1, 2, 3)
                # For this layer type we can also calculate y3 based on the other ys.
                if y_i < 2:
                    y_rhs_i = y_i
                elif (y_i > 3) and (y_i < 6):
                    y_rhs_i = y_i - 2
                else:
                    # y4 == NAN for liquid dynamic layers; y3 will be calculated later.
                    y_i += 1
                    continue
        
        # Otherwise we will collapse the solutions together for each slice.
        solution_i = 0
        while num_sols > solution_i:
            # New solution; reset slice index
            slice_i = 0
            slice_i_shifted = layer_start_index
            while slice_end > slice_i_shifted:
                if solution_i == 0:
                    # Initialize values
                    solution_ptr[slice_i_shifted * num_output_ys + lhs_y_index] = \
                        (constant_vector_ptr[solution_i] *
                        storage_by_solution[solution_i][slice_i * num_ys + y_rhs_i])
                else:
                    # Add new results to old value
                    solution_ptr[slice_i_shifted * num_output_ys + lhs_y_index] += \
                        (constant_vector_ptr[solution_i] *
                        storage_by_solution[solution_i][slice_i * num_ys + y_rhs_i])

                slice_i += 1
                slice_i_shifted += 1
            solution_i += 1
        y_i += 1
    
    if calculate_y3:
        # Handle y3 for dynamic liquid layers.
        lhs_y_index = ytype_i * max_num_y
        slice_i = 0
        slice_i_shifted = layer_start_index
        while slice_end > slice_i_shifted:
            solution_ptr[slice_i_shifted * num_output_ys + (lhs_y_index + 2)] = \
                (1. / (frequency_to_use**2 * layer_radius_ptr[slice_i])) * \
                (solution_ptr[slice_i_shifted * num_output_ys + (lhs_y_index + 0)] * layer_gravity_ptr[slice_i] -
                solution_ptr[slice_i_shifted * num_output_ys + (lhs_y_index + 1)] / layer_density_ptr[slice_i] -
                solution_ptr[slice_i_shifted * num_output_ys + (lhs_y_index + 4)])

            slice_i += 1
            slice_i_shifted += 1
