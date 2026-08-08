# distutils: language = c++

cdef extern from "spin_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_SpinConfig:
        double moment_of_inertia_factor

    cdef cppclass c_Spin:
        c_Spin()
        c_Spin(const c_SpinConfig& config)
        const c_SpinConfig& get_config()
        double calc_moment_of_inertia(
            double mass,
            double radius_outer,
            double radius_inner)
        double calc_dspin_dt(
            double host_mass,
            double dU_dO,
            double moment_of_inertia)
        double calc_synchronous_spin(double orbital_frequency)


cdef class Spin:
    cdef c_Spin _spin
