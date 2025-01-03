# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from TidalPy.constants cimport d_NAN_DBL
from TidalPy.utilities.math.complex cimport cf_build_dblcmplx

cdef double complex cmplx_NAN = cf_build_dblcmplx(d_NAN_DBL, d_NAN_DBL)

cdef void cf_top_to_bottom_interface_bc(
        double complex* constant_vector_ptr,
        double complex* layer_above_constant_vector_ptr,
        double complex* uppermost_y_per_solution_ptr,
        double gravity_upper,
        double layer_above_lower_gravity,
        double density_upper,
        double layer_above_lower_density,
        int layer_type,
        int layer_above_type,
        cpp_bool layer_is_static,
        cpp_bool layer_above_is_static,
        cpp_bool layer_is_incomp,
        cpp_bool layer_above_is_incomp,
        size_t num_sols,
        size_t max_num_y
        ) noexcept nogil:
    
    # Loop variables
    cdef size_t solution_i
    
    # Interfaces are defined at the bottom of the layer in question. However, this function is calculating
    #  the transition at the top of each layer as it works its way down.
    #  So, for interface values, we actually need the ones of the layer above us.
    cdef double interface_gravity = 0.5 * (gravity_upper + layer_above_lower_gravity)
    cdef double liquid_density_at_interface = d_NAN_DBL

    cdef cpp_bool layer_is_solid = False
    cdef cpp_bool layer_above_is_solid = False
    if layer_type == 0:
        layer_is_solid = True
    if layer_above_type == 0:
        layer_above_is_solid = True

    if not layer_is_solid:
        if layer_is_static:
            liquid_density_at_interface = density_upper
        elif not layer_above_is_solid and layer_above_is_static:
            liquid_density_at_interface = layer_above_lower_density
        else:
            liquid_density_at_interface = density_upper
    elif not layer_above_is_solid:
        liquid_density_at_interface = layer_above_lower_density

    # Solve for constant vector
    cdef double complex y4_frac_1  = cmplx_NAN
    cdef double complex y4_frac_2  = cmplx_NAN
    cdef double complex gamma_1    = cmplx_NAN
    cdef double complex gamma_2    = cmplx_NAN
    cdef double complex lambda_1   = cmplx_NAN
    cdef double complex lambda_2   = cmplx_NAN
    cdef double complex lower_s1y1 = cmplx_NAN
    cdef double complex lower_s1y2 = cmplx_NAN
    cdef double complex lower_s1y5 = cmplx_NAN
    cdef double complex lower_s1y6 = cmplx_NAN
    cdef double complex lower_s2y1 = cmplx_NAN
    cdef double complex lower_s2y2 = cmplx_NAN
    cdef double complex lower_s2y5 = cmplx_NAN
    cdef double complex lower_s2y6 = cmplx_NAN

    if layer_is_solid:
        # Solid Layer
        if layer_above_is_solid:
            # Both layers are solid. Constants are the same.
            solution_i = 0
            while num_sols > solution_i:
                constant_vector_ptr[solution_i] = layer_above_constant_vector_ptr[solution_i]
                solution_i += 1
        else:
            # Create some helper functions that will be needed
            y4_frac_1 = (
                    -uppermost_y_per_solution_ptr[0 * max_num_y + 3] /
                    uppermost_y_per_solution_ptr[2 * max_num_y + 3]
                )
            y4_frac_2 = (
                    -uppermost_y_per_solution_ptr[1 * max_num_y + 3] /
                    uppermost_y_per_solution_ptr[2 * max_num_y + 3]
                )

            if layer_above_is_static:
                # Need to find 3 solid constants from 1 liquid constant
                # S74, Page 131
                constant_vector_ptr[0] = layer_above_constant_vector_ptr[0]
                # Derived by JPR based on Eq 21 (2nd line) of S74
                gamma_1 = \
                    (
                        uppermost_y_per_solution_ptr[0 * max_num_y + 1] +
                        y4_frac_1 * uppermost_y_per_solution_ptr[2 * max_num_y + 1]
                    ) - \
                    (
                        liquid_density_at_interface *
                        (
                            interface_gravity * 
                            (
                                uppermost_y_per_solution_ptr[0 * max_num_y + 0] +
                                y4_frac_1 * uppermost_y_per_solution_ptr[2 * max_num_y + 0]
                            ) -
                            (
                                uppermost_y_per_solution_ptr[0 * max_num_y + 4] +
                                y4_frac_1 * uppermost_y_per_solution_ptr[2 * max_num_y + 4]
                            )
                        )
                    )
                gamma_2 = \
                    (
                        uppermost_y_per_solution_ptr[1 * max_num_y + 1] + 
                        y4_frac_2 * uppermost_y_per_solution_ptr[2 * max_num_y + 1]
                    ) - \
                    (
                        liquid_density_at_interface *
                        (
                            interface_gravity * 
                            (
                                uppermost_y_per_solution_ptr[1 * max_num_y + 0] +
                                y4_frac_2 * uppermost_y_per_solution_ptr[2 * max_num_y + 0]
                            ) - 
                            (
                                uppermost_y_per_solution_ptr[1 * max_num_y + 4] + 
                                y4_frac_2 * uppermost_y_per_solution_ptr[2 * max_num_y + 4]
                            )
                        )
                    )

                constant_vector_ptr[1] = (-gamma_1 / gamma_2) * constant_vector_ptr[0]
                # TS72, Eq. 142 (utilizes y_4 = 0)
                constant_vector_ptr[2] = y4_frac_1 * constant_vector_ptr[0] + y4_frac_2 * constant_vector_ptr[1]

            else:
                # Need to find 3 solid constants from 2 liquid constants
                # TS72, Eq. 144
                constant_vector_ptr[0] = layer_above_constant_vector_ptr[0]
                constant_vector_ptr[1] = layer_above_constant_vector_ptr[1]
                # TS72, Eq. 142 (utilizes y_4 = 0)
                constant_vector_ptr[2] = y4_frac_1 * constant_vector_ptr[0] + y4_frac_2 * constant_vector_ptr[1]
    else:
        if layer_is_static:
            if not layer_above_is_solid:
                # Liquid layer above
                if layer_above_is_static:
                    # Both layers are static liquids. Constants are the same.
                    constant_vector_ptr[0] = layer_above_constant_vector_ptr[0]
                else:
                    # Dynamic liquid above
                    # JPR decided to follow a similar approach as Eq. 20 in S74:
                    #   Treat the lower static liquid as normal.
                    #   The upper dynamic liquid layer is treated like the solid layer in Eq. 20 except
                    #    that y_3 is undefined as is "set 3" solution mentioned in that text.
                    constant_vector_ptr[0] = layer_above_constant_vector_ptr[0]
            else:
                # Solid layer above
                # Based on S74. The constant in this layer is just equal to the constant in solution 1 of the
                #  layer above.
                constant_vector_ptr[0] = layer_above_constant_vector_ptr[0]
        else:
            if not layer_above_is_solid:
                # Liquid layer above
                if layer_above_is_static:
                    # Need to find 2 liquid (dynamic) constants from 1 liquid (static) constant
                    # S74, Page 131
                    constant_vector_ptr[0] = layer_above_constant_vector_ptr[0]
                    # Derived by JPR based on Eq 21 (2nd line) of S74
                    # Pull out ys
                    # # Solution 1
                    lower_s1y1 = uppermost_y_per_solution_ptr[0 * max_num_y + 0]
                    lower_s1y2 = uppermost_y_per_solution_ptr[0 * max_num_y + 1]
                    lower_s1y5 = uppermost_y_per_solution_ptr[0 * max_num_y + 2]
                    lower_s1y6 = uppermost_y_per_solution_ptr[0 * max_num_y + 3]
                    # # Solution 2
                    lower_s2y1 = uppermost_y_per_solution_ptr[1 * max_num_y + 0]
                    lower_s2y2 = uppermost_y_per_solution_ptr[1 * max_num_y + 1]
                    lower_s2y5 = uppermost_y_per_solution_ptr[1 * max_num_y + 2]
                    lower_s2y6 = uppermost_y_per_solution_ptr[1 * max_num_y + 3]
                    # lambda_j = (y_2j - rho * ( g * y_1j - y_5j))
                    lambda_1 = lower_s1y2 - liquid_density_at_interface * \
                            (interface_gravity * lower_s1y1 - lower_s1y5)
                    lambda_2 = lower_s2y2 - liquid_density_at_interface * \
                            (interface_gravity * lower_s2y1 - lower_s2y5)
                    constant_vector_ptr[1] = (-lambda_1 / lambda_2) * constant_vector_ptr[0]
                else:
                    # Both layers are dynamic liquids. Constants are the same.
                    solution_i = 0
                    while num_sols > solution_i:
                        constant_vector_ptr[solution_i] = layer_above_constant_vector_ptr[solution_i]
                        solution_i += 1
            else:
                # Solid layer above
                # TS72 Eqs. 148-149
                constant_vector_ptr[0] = layer_above_constant_vector_ptr[0]
                constant_vector_ptr[1] = layer_above_constant_vector_ptr[1]
