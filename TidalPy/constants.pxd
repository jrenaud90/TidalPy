cdef extern from "constants_.hpp" namespace "TidalPyConstants" nogil:
    # --- Compile-time Constants (Read-Only) ---
    const double d_ppm
    const double d_ppb
    const double d_INF
    const double d_PI
    const double d_NAN
    
    const double d_DBL_MAX
    const double d_DBL_MIN
    const double d_DBL_MANT_DIGITS
    const double d_EPS
    const double d_EPS_10
    const double d_EPS_100

    const double d_MASS_SOLAR
    const double d_RADIUS_SOLAR
    const double d_LUMINOSITY_SOLAR

    const double d_MASS_TRAP1
    const double d_RADIUS_TRAP1
    const double d_LUMINOSITY_TRAP1

    const double d_MASS_EARTH
    const double d_RADIUS_EARTH

    const double d_MASS_JUPITER
    const double d_RADIUS_JUPITER

    const double d_MASS_PLUTO
    const double d_RADIUS_PLUTO

    const double d_MASS_IO
    const double d_RADIUS_IO


cdef extern from "constants_.hpp" nogil:

    cdef cppclass TidalPyConfig:
        double d_MIN_FREQUENCY
        double d_MAX_FREQUENCY
        double d_MIN_SPIN_ORBIT_DIFF
        double d_MIN_VISCOSITY
        double d_MIN_MODULUS
        double d_MIN_THICKNESS
        double d_G
        double d_AU
        double d_SBC
        double d_R
        double d_K_BOLTZMAN
        double d_TEST_CONST
    
    cdef TidalPyConfig* tidalpy_config_ptr

# Expose the API function
cdef TidalPyConfig* get_shared_config_address()