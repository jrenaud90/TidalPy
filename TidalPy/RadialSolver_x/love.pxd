from libcpp.complex cimport complex as cpp_complex


cdef extern from "love_.hpp" nogil:

    cdef cppclass c_LoveNumbers:
        cpp_complex[double] k
        cpp_complex[double] h
        cpp_complex[double] l
        c_LoveNumbers() except +
        c_LoveNumbers(
            const cpp_complex[double]& k_,
            const cpp_complex[double]& h_,
            const cpp_complex[double]& l_) except +
        c_LoveNumbers(
            const double k_,
            const double h_,
            const double l_) except +
        double get_Q_k() const
        double get_Q_h() const
        double get_Q_l() const
        double get_lag_k() const
        double get_lag_h() const
        double get_lag_l() const

    cdef c_LoveNumbers c_find_love(
        cpp_complex[double]* surface_solutions_ptr,
        double surface_gravity
    ) noexcept nogil


cdef class LoveNumbers:
    cdef c_LoveNumbers _cinst
