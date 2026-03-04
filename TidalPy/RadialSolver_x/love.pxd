from libcpp.complex cimport complex as cpp_complex


cdef extern from "love_.hpp" nogil:

    # TODO: Move this stuff into radial solver eventually
    cdef cppclass c_LoveNumbers:
        cpp_complex[double] k
        cpp_complex[double] h
        cpp_complex[double] l
        c_LoveNumbers() except +
        c_LoveNumbers(cpp_complex[double]& k_, cpp_complex[double]& h_, cpp_complex[double]& l_) except +
        c_LoveNumbers(double k_, double h_, double l_) except +
        double get_Q_k() const
        double get_Q_h() const
        double get_Q_l() const
        double get_lag_k() const
        double get_lag_h() const
        double get_lag_l() const


cdef class LoveNumbers:
    cdef c_LoveNumbers _cinst
