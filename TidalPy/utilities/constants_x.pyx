# distutils: language = c
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# Astrophysical Constants
cdef double G = 6.67430e-11

# Extremes
# Forcing Frequency Extremes
# Assume that any forcing period larger than a Gyr leads to a zero in frequency.
# Converting to frequency is roughly 1.0e-17 rads s-1
cdef double MIN_FREQUENCY = 1.0e-17
# Assume max frequency is for a forcing period of 1 micro-second
cdef double MAX_FREQUENCY = 1.0e8

# Shear/Bulk Modulus Extremes
cdef double MIN_MODULUS = 1.0e-3