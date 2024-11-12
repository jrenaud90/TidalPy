from CyRK cimport DiffeqFuncType, PreEvalFunc, CySolverResult

from TidalPy.Material.eos.solver cimport EOSSolutionVec


cdef struct RadialSolverDiffeqArgStruct:
    double degree_l
    double lp1
    double lm1
    double llp1
    double G
    double grav_coeff
    double frequency
    CySolverResult* eos_solution_ptr

cdef void cf_solid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* void_args_ptr,
        PreEvalFunc eos_func) noexcept nogil

cdef void cf_solid_dynamic_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* void_args_ptr,
        PreEvalFunc eos_func) noexcept nogil

cdef void cf_solid_static_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* void_args_ptr,
        PreEvalFunc eos_func) noexcept nogil

cdef void cf_solid_static_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* void_args_ptr,
        PreEvalFunc eos_func) noexcept nogil

cdef void cf_liquid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* void_args_ptr,
        PreEvalFunc eos_func) noexcept nogil

cdef void cf_liquid_dynamic_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* void_args_ptr,
        PreEvalFunc eos_func) noexcept nogil

cdef void cf_liquid_static_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* void_args_ptr,
        PreEvalFunc eos_func) noexcept nogil

cdef DiffeqFuncType cf_find_layer_diffeq(
        int layer_type,
        int layer_is_static,
        int layer_is_incomp) noexcept nogil