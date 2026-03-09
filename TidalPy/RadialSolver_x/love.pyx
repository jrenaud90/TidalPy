# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as cnp
cnp.import_array()


cdef class LoveNumbers:
    
    def __cinit__(self, double complex k_, double complex h_, double complex l_):
        cdef cpp_complex[double] cpp_k = cpp_complex[double](k_.real, k_.imag)
        cdef cpp_complex[double] cpp_h = cpp_complex[double](h_.real, h_.imag)
        cdef cpp_complex[double] cpp_l = cpp_complex[double](l_.real, l_.imag)

        self._cinst = c_LoveNumbers(cpp_k, cpp_h, cpp_l)
    
    @property
    def k(self):
        return self._cinst.k
    
    @property
    def h(self):
        return self._cinst.h
    
    @property
    def l(self):
        return self._cinst.l

    @property
    def Q_k(self):
        return self._cinst.get_Q_k()
    
    @property
    def Q_h(self):
        return self._cinst.get_Q_h()
    
    @property
    def Q_l(self):
        return self._cinst.get_Q_l()
    
    @property
    def lag_k(self):
        return self._cinst.get_lag_k()
    
    @property
    def lag_h(self):
        return self._cinst.get_lag_h()
    
    @property
    def lag_l(self):
        return self._cinst.get_lag_l()


def find_love(
        double complex[::1] surface_solutions,
        double surface_gravity
        ):
    """
    Compute Love and Shida numbers from radial solution y-values at the planet surface.

    Parameters
    ----------
    surface_solutions : ndarray[complex128]
        y-values at the surface: [y1, y2, y3, y4, y5, y6].
    surface_gravity : double
        Gravitational acceleration at the surface [m s-2].

    Returns
    -------
    love : LoveNumbers
        Object containing k, h, l Love/Shida numbers with Q and lag properties.
    """
    cdef cpp_complex[double]* surface_ptr = <cpp_complex[double]*>&surface_solutions[0]

    cdef c_LoveNumbers love = c_find_love(surface_ptr, surface_gravity)

    return LoveNumbers(
        love.k,
        love.h,
        love.l,
    )
