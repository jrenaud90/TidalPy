from CyRK cimport DiffeqFuncType, PreEvalFunc

from TidalPy.Material_x.eos.eos_solution cimport c_EOSSolution


cdef extern from "odes_.hpp" nogil:

    cdef struct c_RadialSolverArgs:
        double degree_l
        double lp1
        double lm1
        double llp1
        double G
        double grav_coeff
        double frequency
        size_t layer_index
        c_EOSSolution* eos_solution_ptr

    cdef void c_solid_dynamic_compressible(
        double* dy_ptr, double radius, double* y_ptr,
        char* args_ptr, PreEvalFunc unused) noexcept nogil

    cdef void c_solid_dynamic_incompressible(
        double* dy_ptr, double radius, double* y_ptr,
        char* args_ptr, PreEvalFunc unused) noexcept nogil

    cdef void c_solid_static_compressible(
        double* dy_ptr, double radius, double* y_ptr,
        char* args_ptr, PreEvalFunc unused) noexcept nogil

    cdef void c_solid_static_incompressible(
        double* dy_ptr, double radius, double* y_ptr,
        char* args_ptr, PreEvalFunc unused) noexcept nogil

    cdef void c_liquid_dynamic_compressible(
        double* dy_ptr, double radius, double* y_ptr,
        char* args_ptr, PreEvalFunc unused) noexcept nogil

    cdef void c_liquid_dynamic_incompressible(
        double* dy_ptr, double radius, double* y_ptr,
        char* args_ptr, PreEvalFunc unused) noexcept nogil

    cdef void c_liquid_static_incompressible(
        double* dy_ptr, double radius, double* y_ptr,
        char* args_ptr, PreEvalFunc unused) noexcept nogil

    cdef DiffeqFuncType c_find_layer_diffeq(
        int layer_type, int layer_is_static, int layer_is_incomp) noexcept nogil

    cdef size_t c_find_num_shooting_solutions(
        int layer_type, int layer_is_static, int layer_is_incomp) noexcept nogil
