from TidalPy.RadialSolver.derivatives.base cimport RadialSolverBase


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

cdef RadialSolverBase cf_build_solver(
    int layer_type,
    bint is_static,
    bint is_incomp,

    # RadialSolverBase Inputs
    size_t num_slices,
    size_t num_ys,
    double* radius_array_ptr,
    double* density_array_ptr,
    double* gravity_array_ptr,
    double* bulk_modulus_array_ptr,
    double complex* shear_modulus_array_ptr,
    double frequency,
    unsigned char degree_l,
    double G_to_use,

    # Regular CySolver Inputs
    (double, double) t_span,
    double* y0_ptr,
    double* rtols,
    double* atols,
    unsigned char rk_method,
    double max_step,
    size_t max_num_steps,
    size_t expected_size,
    size_t max_ram_MB,

    # Additional optional arguments for RadialSolver class
    bint limit_solution_to_radius
    )
