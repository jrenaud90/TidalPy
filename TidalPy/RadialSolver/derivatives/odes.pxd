from libcpp.memory cimport shared_ptr

from CyRK cimport DiffeqFuncType, PreEvalFunc

from TidalPy.Material.eos.eos_solution cimport EOSSolutionCC


cdef struct RadialSolverDiffeqArgStruct:
    double degree_l
    double lp1
    double lm1
    double llp1
    double G
    double grav_coeff
    double frequency
    size_t layer_index
    EOSSolutionCC* eos_solution_ptr

cdef void cf_solid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept nogil

cdef void cf_solid_dynamic_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept nogil

cdef void cf_solid_static_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept nogil

cdef void cf_solid_static_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept nogil

cdef void cf_liquid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept nogil

cdef void cf_liquid_dynamic_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept nogil

cdef void cf_liquid_static_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept nogil

cdef DiffeqFuncType cf_find_layer_diffeq(
        int layer_type,
        int layer_is_static,
        int layer_is_incomp) noexcept nogil