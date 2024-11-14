from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
from libcpp cimport bool as cpp_bool

from CyRK cimport PreEvalFunc, CySolveOutput
from CyRK.utils.memory cimport shared_ptr, make_shared

from TidalPy.Material.eos.ode cimport EOS_ODEInput
from TidalPy.Material.eos.eos_solution cimport EOSSolution, EOS_DY_VALUES

ctypedef vector[CySolveOutput] CySolverResult_sptr_vec
ctypedef shared_ptr[CySolverResult_sptr_vec] CySolverResult_sptr_vec_sptr

cdef void solve_eos(
    shared_ptr[EOSSolution] eos_solution,
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
