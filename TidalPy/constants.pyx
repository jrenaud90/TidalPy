# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# In order for both cython and python to see these variables we have to declare two. the "d_" prefix are pure C
# These are then read into python objects which share the same name except without the "d_" prefix

# Extremes
# Forcing Frequency Extremes
# Assume that any forcing period larger than a Myr leads to a zero in frequency.
# Converting to frequency is roughly 1.0e-14 rads s-1
import TidalPy

cdef double d_MIN_FREQUENCY = 1.0e-14
if TidalPy.config:
    d_MIN_FREQUENCY = TidalPy.config['tides']['modes']['minimum_frequency']
MIN_FREQUENCY = d_MIN_FREQUENCY

# Assume max frequency is for a forcing period of 1 micro-second
cdef double d_MAX_FREQUENCY = 1.0e8
if TidalPy.config:
    d_MAX_FREQUENCY = TidalPy.config['tides']['modes']['maximum_frequency']
MAX_FREQUENCY = d_MAX_FREQUENCY

# Shear/Bulk Modulus Extremes
cdef double d_MIN_MODULUS = 1.0e-3
MIN_MODULUS = d_MIN_MODULUS

# Mathematics
cdef double d_ppm = 1.e-6
ppm = d_ppm
cdef double d_ppb = 1.e-9
ppb = d_ppb

from libc.float cimport DBL_MAX as _DBL_MAX
from libc.float cimport DBL_MIN as _DBL_MIN
from libc.float cimport DBL_MANT_DIG as _DBL_MANT_DIG
from libc.math cimport M_PI as _M_PI
cdef double d_DBL_MAX      = _DBL_MAX
cdef double d_DBL_MIN      = _DBL_MIN
cdef double d_DBL_MANT_DIG = _DBL_MANT_DIG
cdef double d_PI_DBL       = _M_PI
DBL_MAX      = d_DBL_MAX
DBL_MIN      = d_DBL_MIN
DBL_MANT_DIG = d_DBL_MANT_DIG
PI_DBL       = d_PI_DBL

from libcpp.limits cimport numeric_limits
cdef double d_INF_DBL = numeric_limits[double].infinity()
INF_DBL = d_INF_DBL
cdef double d_EPS_DBL = numeric_limits[double].epsilon()
EPS_DBL = d_EPS_DBL
cdef double d_EPS_DBL_10  = 10 * d_EPS_DBL
EPS_DBL_10  = d_EPS_DBL_10
cdef double d_EPS_DBL_100 = 100 * d_EPS_DBL
EPS_DBL_100 = d_EPS_DBL_100
cdef double d_NAN_DBL = numeric_limits[double].quiet_NaN()
NAN_DBL = d_NAN_DBL

# Astrophysical Constants
from scipy.constants import G as G_
from scipy.constants import au as au_
from scipy.constants import Stefan_Boltzmann
from scipy.constants import R as R_

# Sun
cdef double d_mass_solar       = 1.988435e30  # [kg]
mass_solar = d_mass_solar
cdef double d_radius_solar     = 6.957e8  # [m]
radius_solar = d_radius_solar
cdef double d_luminosity_solar = 3.848e26  # [Watts]
luminosity_solar = d_luminosity_solar

# TRAPPIST-1
cdef double d_mass_trap1       = 0.0898 * d_mass_solar  # [kg]
mass_trap1 = d_mass_trap1
cdef double d_radius_trap1     = 0.1192 * d_radius_solar  # [m]
radius_trap1 = d_radius_trap1
cdef double d_luminosity_trap1 = 0.000553 * d_luminosity_solar # [Watts]
luminosity_trap1 = d_luminosity_trap1

# Earth
cdef double d_mass_earth   = 5.9721986e24  # [kg]
mass_earth = d_mass_earth
cdef double d_radius_earth = 6.371008e6  # [m]
radius_earth = d_radius_earth

# Jupiter
cdef double d_mass_jupiter   = 1.89813e27  # [kg]
mass_jupiter = d_mass_jupiter
cdef double d_radius_jupiter = 6.9911e7  # [m]
radius_jupiter = d_radius_jupiter

# Pluto
cdef double d_mass_pluto   = 1.309e22  # [kg]
mass_pluto = d_mass_pluto
cdef double d_radius_pluto = 1.1899e6  # [m]
radius_pluto = d_radius_pluto

# Io
cdef double d_mass_io   = 8.9298e22  # [kg]
mass_io = d_mass_io
cdef double d_radius_io = 1.82149e6  # [m]
radius_io = d_radius_io

# Alias Names
cdef double d_G = G_
G = d_G
cdef double d_R = R_
R = d_R
cdef double d_Au = au_
Au = d_Au
au = d_Au
cdef double d_sbc = Stefan_Boltzmann
sbc = d_sbc
cdef double d_SBC = sbc
SBC = d_SBC
cdef double d_newtons_constant = G
newtons_constant = d_newtons_constant

cdef double d_M_sol   = mass_solar
M_sol = d_M_sol
cdef double d_M_earth = mass_earth
M_earth = d_M_earth
cdef double d_M_jup   = mass_jupiter
M_jup = d_M_jup
cdef double d_M_pluto = mass_pluto
M_pluto = d_M_pluto

cdef double d_R_sol   = radius_solar
R_sol = d_R_sol
cdef double d_R_earth = radius_earth
R_earth = d_R_earth
cdef double d_R_jup   = radius_jupiter
R_jup = d_R_jup
cdef double d_R_pluto = radius_pluto
R_pluto = d_R_pluto

cdef double d_L_sol = luminosity_solar
L_sol = d_L_sol
