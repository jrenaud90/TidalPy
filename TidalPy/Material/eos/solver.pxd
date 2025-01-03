from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr, unique_ptr
from libcpp cimport bool as cpp_bool

from CyRK cimport PreEvalFunc, CySolveOutput

from TidalPy.Material.eos.ode cimport EOS_ODEInput
from TidalPy.Material.eos.eos_solution cimport EOSSolutionCC, EOS_DY_VALUES

ctypedef vector[CySolveOutput] CySolverResult_sptr_vec
ctypedef shared_ptr[CySolverResult_sptr_vec] CySolverResult_sptr_vec_sptr

cdef void solve_eos(
    EOSSolutionCC* eos_solution_ptr,
    vector[PreEvalFunc] eos_function_bylayer_ptr_vec,
    vector[EOS_ODEInput] eos_input_bylayer_vec,
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
