from libcpp cimport bool as bool_cpp_t

from TidalPy.radial_solver.numerical.derivatives.ode_base_x cimport RadialSolverBase


cdef class SolidDynamicCompressible(RadialSolverBase):
    pass

cdef class SolidDynamicIncompressible(RadialSolverBase):
    pass

cdef class SolidStaticCompressible(RadialSolverBase):
    pass

cdef class SolidStaticIncompressible(RadialSolverBase):
    pass

cdef class LiquidDynamicCompressible(RadialSolverBase):
    pass

cdef class LiquidDynamicIncompressible(RadialSolverBase):
    pass

cdef class LiquidStaticCompressible(RadialSolverBase):
    pass

cdef class LiquidStaticIncompressible(RadialSolverBase):
    pass

cdef RadialSolverBase build_solver(
    bool_cpp_t is_solid,
    bool_cpp_t is_static,
    bool_cpp_t is_incomp,

    # RadialSolverBase Inputs
    Py_ssize_t num_slices,
    Py_ssize_t num_ys,
    double* radius_array_ptr,
    double* density_array_ptr,
    double* gravity_array_ptr,
    double* bulk_modulus_array_ptr,
    double complex* shear_modulus_array_ptr,
    double frequency,
    unsigned int degree_l,
    double G_to_use,

    # Regular CySolver Inputs
    (double, double) t_span,
    double* y0_ptr,
    double* rtols,
    double* atols,
    unsigned char rk_method,
    double max_step,
    Py_ssize_t max_num_steps,
    Py_ssize_t expected_size,

    # Additional optional arguments for RadialSolver class
    bool_cpp_t limit_solution_to_radius
    )
