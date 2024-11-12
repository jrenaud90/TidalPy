# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# Extremes
# Forcing Frequency Extremes
# Assume that any forcing period larger than a Myr leads to a zero in frequency.
# Converting to frequency is roughly 1.0e-14 rads s-1
cdef double MIN_FREQUENCY = 1.0e-14
# Assume max frequency is for a forcing period of 1 micro-second
cdef double MAX_FREQUENCY = 1.0e8

# Shear/Bulk Modulus Extremes
cdef double MIN_MODULUS = 1.0e-3

# Mathematics
cdef double ppm = 1.e-6
cdef double ppb = 1.e-9

from libc.float cimport DBL_MAX, DBL_MIN, DBL_MANT_DIG
from libc.math cimport M_PI
cdef double PI_DBL = M_PI

from libcpp.limits cimport numeric_limits
cdef double INF_DBL = numeric_limits[double].infinity()
cdef double EPS_DBL = numeric_limits[double].epsilon()
cdef double NAN_DBL = numeric_limits[double].quiet_NaN()

# Astrophysical Constants
from scipy.constants import G as G_
from scipy.constants import au as au_
from scipy.constants import Stefan_Boltzmann
from scipy.constants import R as R_

# Sun
cdef double mass_solar       = 1.988435e30  # [kg]
cdef double radius_solar     = 6.957e8  # [m]
cdef double luminosity_solar = 3.848e26  # [Watts]

# TRAPPIST-1
cdef double mass_trap1       = 0.0898 * mass_solar  # [kg]
cdef double radius_trap1     = 0.1192 * radius_solar  # [m]
cdef double luminosity_trap1 = 0.000553 * luminosity_solar # [Watts]

# Earth
cdef double mass_earth   = 5.9721986e24  # [kg]
cdef double radius_earth = 6.371008e6  # [m]

# Jupiter
cdef double mass_jupiter   = 1.89813e27  # [kg]
cdef double radius_jupiter = 6.9911e7  # [m]

# Pluto
cdef double mass_pluto   = 1.309e22  # [kg]
cdef double radius_pluto = 1.1899e6  # [m]

# Io
cdef double mass_io   = 8.9298e22  # [kg]
cdef double radius_io = 1.82149e6  # [m]

# Alias Names
cdef double G = G_
cdef double R = R_
cdef double Au = au_
cdef double sbc = Stefan_Boltzmann
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

cdef double L_sol = luminosity_solar
