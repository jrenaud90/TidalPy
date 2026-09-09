# distutils: language = c++

cdef extern from "orbit_solver_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_OrbitState:
        double orbital_frequency
        double semi_major_axis
        double eccentricity
        double target_mass
        double host_mass

    cdef cppclass c_OrbitDerivatives:
        double da_dt
        double de_dt
        double dn_dt

    cdef cppclass c_OrbitSolver:
        c_OrbitSolver()
        double calc_da_dt(const c_OrbitState& state, double dU_dM)
        double calc_de_dt(const c_OrbitState& state, double dU_dM, double dU_dw)
        double calc_dn_dt(double orbital_frequency, double semi_major_axis, double da_dt)
        c_OrbitDerivatives calc_derivatives(const c_OrbitState& state, double dU_dM, double dU_dw)


cdef class OrbitSolver:
    cdef c_OrbitSolver _solver
