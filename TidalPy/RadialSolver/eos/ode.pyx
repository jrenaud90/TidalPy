# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from CyRK.cy.cysolverNew cimport PreEvalFunc

from TidalPy.RadialSolver.eos.common cimport EOSOutput


cdef void eos_solution(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* input_args,
        PreEvalFunc eos_function) noexcept nogil:
    """ Solve for EOS components as a function of radius. """


    # Solve for gravity
    dy_ptr[8] = grav_coeff * eos_output.density - 2.0 * gravity * r_inverse

    # Solve for pressure
    dy_ptr[9] = -eos_output.density * gravity
