cdef extern from "tpy_constants.hpp" nogil:
    # Math / Numerics
    double NAN_DBL
    double INF_DBL
    long double PI_LDBL
    double PI_DBL

    # Physics
    double G

    # TidalPy Specific (extremes and limits)
    double MIN_FREQUENCY
    double MAX_FREQUENCY
    double MIN_MODULUS
