from CyRK cimport PreEvalFunc, CySolverResult
from CyRK.utils.memory cimport shared_ptr, make_shared
from CyRK.utils.vector cimport vector
from libcpp cimport bool as cpp_bool

from TidalPy.Material.eos.ode cimport EOS_ODEInput

ctypedef CySolverResult* CySolveResultPtr
ctypedef vector[CySolveResultPtr] EOSSolutionVec

cdef struct GlobalEOSSolutionStorage:
    vector[double] radius
    vector[double] pressure
    vector[double] gravity
    vector[double] density
    vector[double] static_shear_modulus
    vector[double] shear_viscosity
    vector[double complex] complex_shear_modulus
    vector[double] static_bulk_modulus
    vector[double] bulk_viscosity
    vector[double complex] complex_bulk_modulus

cdef EOSSolutionVec solve_eos(
    cpp_bool* success_ptr,
    char* message_ptr,
    double* radius_array_ptr,
    size_t len_radius_array,
    double* layer_upper_radii,
    unsigned int num_layers,
    PreEvalFunc* eos_function_bylayer_ptrs,
    EOS_ODEInput** eos_input_bylayer_ptrs,
    double planet_bulk_density,
    double surface_pressure = *,
    double G_to_use = *,
    unsigned int integration_method = *,
    double rtol = *,
    double atol = *,
    double pressure_tol = *,
    unsigned int max_iters = *,
    cpp_bool verbose = *
    ) noexcept nogil
