from libcpp cimport bool as cpp_bool

cdef extern from "constants_.hpp" namespace "TidalPyConstants" nogil:
    # --- Compile-time Constants (Read-Only) ---
    const double d_ppm
    const double d_ppb
    const double d_INF
    const double d_PI
    const double d_NAN
    const double d_SECONDS_PER_MYR
    
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
        double   d_MIN_FREQUENCY
        double   d_MAX_FREQUENCY
        double   d_MIN_SPIN_ORBIT_DIFF
        double   d_MIN_VISCOSITY
        double   d_MIN_MODULUS
        double   d_MIN_SOLID_RIGIDITY
        double   d_MIN_THICKNESS
        double   d_NUMERICAL_FLOOR
        double   d_LAYER_CONTINUITY_RTOL
        double   d_MAX_START_RADIUS_FRAC
        double   d_MIN_SURFACE_RCOND
        double   d_FREQUENCY_MATCH_RTOL
        double   d_MIN_NUSSELT
        double   d_MAX_EOS_MASS_RATIO
        double   d_EOS_INVERT_RTOL
        int      d_EOS_INVERT_MAX_ITERS
        int      d_TIDES_3D_LATITUDE_NODES
        int      d_TIDES_3D_LONGITUDE_NODES
        int      d_TIDES_3D_RADIAL_SLICES
        int      d_EOS_SOLVER_METHOD
        double   d_EOS_SOLVER_RTOL
        double   d_EOS_SOLVER_ATOL
        double   d_EOS_SOLVER_PRESSURE_TOL
        int      d_EOS_SOLVER_MAX_ITERS
        int      d_EOS_SOLVER_SLICES_PER_LAYER
        cpp_bool d_EOS_SOLVER_NONDIMENSIONALIZE
        cpp_bool d_EOS_SOLVER_SOLVE_TEMPERATURE
        int      d_RADIAL_SOLVER_METHOD
        double   d_RADIAL_SOLVER_RTOL
        double   d_RADIAL_SOLVER_ATOL
        cpp_bool d_RADIAL_SOLVER_USE_KAMATA
        double   d_RADIAL_SOLVER_START_RADIUS_TOL
        cpp_bool d_RADIAL_SOLVER_SCALE_RTOLS
        int      d_RADIAL_SOLVER_MAX_NUM_STEPS
        int      d_RADIAL_SOLVER_EXPECTED_SIZE
        int      d_RADIAL_SOLVER_MAX_RAM_MB
        cpp_bool d_RADIAL_SOLVER_NONDIMENSIONALIZE
        double d_G
        double d_AU
        double d_SBC
        double d_R
        double d_K_BOLTZMANN
        double d_TEST_CONST
    
    cdef TidalPyConfig* tidalpy_config_ptr

    void set_tidalpy_config_ptr(TidalPyConfig* ptr)

cdef TidalPyConfig* get_shared_config_address()
