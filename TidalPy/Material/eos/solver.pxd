from CyRK cimport CySolveOutput, PreEvalFunc
from CyRK.utils.memory cimport shared_ptr, make_shared
from CyRK.utils.vector cimport vector

from TidalPy.Material.eos.ode cimport eos_solution, EOS_ODEInput

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

ctypedef shared_ptr[GlobalEOSSolutionStorage] GlobalEOSSolution

cdef vector[CySolveOutput] solve_eos(
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
    unsigned int max_iters = *
    ) noexcept nogil
