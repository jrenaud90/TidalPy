from libcpp cimport bool as bool_cpp_t

cdef extern from "love.c":
    void find_love_cf(
        double complex* complex_love_numbers_ptr,
        double complex* surface_solutions_ptr,
        double surface_gravity
        ) noexcept nogil

cdef class RadialSolverSolution():

    cdef public str message
    cdef public bool_cpp_t success

    # Result structure information
    cdef size_t num_ys
    cdef size_t num_slices
    cdef size_t total_size
    cdef size_t num_ytypes
    cdef tuple ytypes

    # Result pointers and data
    cdef double complex* full_solution_ptr
    cdef double complex[::1] full_solution_view

    # Love number information
    cdef double complex* complex_love_ptr
    cdef double complex[::1] complex_love_view


cdef RadialSolverSolution cf_radial_solver(
    const double[:] radius_array,
    const double[:] density_array,
    const double[:] gravity_array,
    const double[:] bulk_modulus_array,
    const double complex[:] complex_shear_modulus_array,
    double frequency,
    double planet_bulk_density,
    tuple is_solid_by_layer,
    tuple is_static_by_layer,
    tuple is_incompressible_by_layer,
    tuple upper_radius_by_layer,
    unsigned int degree_l = *,
    tuple solve_for = *,
    bool_cpp_t use_kamata = *,
    int integration_method = *,
    double integration_rtol = *,
    double integration_atol = *,
    bool_cpp_t scale_rtols_by_layer_type = *,
    size_t max_num_steps = *,
    size_t expected_size = *,
    size_t max_ram_MB = *,
    double max_step = *,
    bool_cpp_t limit_solution_to_radius = *,
    bool_cpp_t nondimensionalize = *,
    bool_cpp_t verbose = *,
    bool_cpp_t raise_on_fail = *
    )
