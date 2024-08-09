cdef extern from "tpy_constants.hpp" nogil:
    # Numerical constants
    cdef const double NAN_DBL
    cdef const double INF_DBL
    cdef const long double PI_LDBL
    cdef const double PI_DBL

    # Astrophysical Constants
    cdef const double au # Au m-1
    cdef const double G  # m^3 kg^-1 s^-2

    # Sun
    cdef const double mass_solar        # [kg]
    cdef const double radius_solar      # [m]
    cdef const double luminosity_solar  # [Watts]

    # Earth
    cdef const double mass_earth      # [kg]
    cdef const double radius_earth    # [m]

    # Jupiter
    cdef const double mass_jupiter    # [kg]
    cdef const double radius_jupiter  # [m]

    # Pluto
    cdef const double mass_pluto    # [kg]
    cdef const double radius_pluto  # [m]

    # Io
    cdef const double mass_io    # [kg]
    cdef const double radius_io  # [m]

    # Other Physics
    cdef const double molar_gas_constant         # J mol^-1 K^-1
    cdef const double stefan_boltzmann_constant  # W m^-2 K^-4

    # Mathematics
    cdef const double ppm
    cdef const double ppb

    # SI prefixes
    cdef const double quetta
    cdef const double ronna
    cdef const double yotta
    cdef const double zetta
    cdef const double exa
    cdef const double peta
    cdef const double tera
    cdef const double giga
    cdef const double mega
    cdef const double kilo
    cdef const double hecto
    cdef const double deka
    cdef const double deci
    cdef const double centi
    cdef const double milli
    cdef const double micro
    cdef const double nano
    cdef const double pico
    cdef const double femto
    cdef const double atto
    cdef const double zepto
    cdef const double yocto
    cdef const double ronto
    cdef const double quecto

    # Extremes
    # Forcing Frequency Extremes
    # Assume that any forcing period larger than a Gyr leads to a zero in frequency.
    # Converting to frequency is roughly 1.0e-17 rads s-1
    cdef const double MIN_FREQUENCY
    # Assume max frequency is for a forcing period of 1 micro-second
    cdef const double MAX_FREQUENCY

    #/ Shear/Bulk Modulus Extremes
    cdef const double MIN_MODULUS


# Alias Names
cdef double Au  = au
cdef double sbc = stefan_boltzmann_constant
cdef double SBC = sbc
cdef double newtons_constant = G

cdef double M_sol   = mass_solar
cdef double M_earth = mass_earth
cdef double M_jup   = mass_jupiter
cdef double M_pluto = mass_pluto

cdef double R_sol   = radius_solar
cdef double R_earth = radius_earth
cdef double R_jup   = radius_jupiter
cdef double R_pluto = radius_pluto

cdef double L_sol   = luminosity_solar
