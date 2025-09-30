# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# Largely based off of the methods used in [BurnMan](https://github.com/geodynamics/burnman/blob/main/burnman/eos/birch_murnaghan.py)

from libc.math cimport pow, cbrt

cdef struct BMEOSInput:

    double K_0
    double Kprime_0
    double P_0
    double V_0


cdef double bulk_modulus(double volume, double pressure, BMEOSInput* bm_parameter_ptr) noexcept nogil:

    cdef double K_0      = bm_parameter_ptr.K_0
    cdef double Kprime_0 = bm_parameter_ptr.Kprime_0
    cdef double V_0      = bm_parameter_ptr.V_0

    cdef double x = V_0 / volume
    cdef double f = 0.5 * (pow(x, 2.0 / 3.0) - 1.0)

    return K_0 * pow(1.0 + 2.0 * f, 5.0 / 2.0) * (
        1.0
        + (3.0 * Kprime_0 - 5.0) * f
        + (27.0 / 2.0) * (Kprime_0 - 4.0) * f * f
    )


cdef double pressure(double V0_over_V, BMEOSInput* bm_parameter_ptr) noexcept nogil:

    cdef double K_0      = bm_parameter_ptr.K_0
    cdef double P_0      = bm_parameter_ptr.P_0
    cdef double Kprime_0 = bm_parameter_ptr.Kprime_0

    cdef double V0_over_V_cbrt = cbrt(V0_over_V)
    cdef double V0_over_V_cbrt2 = V0_over_V_cbrt * V0_over_V_cbrt
    cdef double V0_over_V_cbrt5 = pow(V0_over_V_cbrt, 5.0)
    cdef double V0_over_V_cbrt7 = V0_over_V_cbrt5 * V0_over_V_cbrt2

    cdef double f = (V0_over_V_cbrt7 - V0_over_V_cbrt5) * \
        (1.0 - 0.75 * (4.0 - Kprime_0) * (V0_over_V_cbrt2 - 1.0))

    return P_0 + (3.0 / 2.0) * f * K_0


def volume(pressure, BMEOSInput* bm_parameter_ptr):

    def delta_pressure(volume):
        return pressure_third_order(params["V_0"] / volume, params) - pressure

    try:
        sol = bracket(delta_pressure, params["V_0"], 1.0e-2 * params["V_0"])
    except ValueError:
        raise ValueError(
            "Cannot find a volume, perhaps you are outside of the "
            "range of validity for the equation of state?"
        )
    return opt.brentq(delta_pressure, sol[0], sol[1])
