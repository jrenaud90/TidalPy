cdef extern from "love_.cpp" nogil:

    void find_love_cf(
        double* complex_love_numbers_ptr,
        double* surface_solutions_ptr,
        double surface_gravity)
