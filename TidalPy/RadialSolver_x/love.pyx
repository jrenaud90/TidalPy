# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False


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
