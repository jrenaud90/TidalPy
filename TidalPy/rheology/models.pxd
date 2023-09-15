from TidalPy.rheology.common cimport RheologyModelBase

########################################################################################################################
######################################## New Rheological Models Go Below Here! #########################################
########################################################################################################################

cdef class ElasticRheology(RheologyModelBase):
    pass

cdef class NewtonRheology(RheologyModelBase):
    pass

cdef class MaxwellRheology(RheologyModelBase):
    pass

cdef class VoigtRheology(RheologyModelBase):
    cdef double voigt_modulus_scale
    cdef double voigt_viscosity_scale

cdef class BurgersRheology(RheologyModelBase):
    cdef double voigt_modulus_scale
    cdef double voigt_viscosity_scale

cdef class AndradeRheology(RheologyModelBase):
    cdef double alpha
    cdef double alpha_factorial
    cdef double zeta
    cdef double complex sine_term

cdef class SundbergCooperRheology(RheologyModelBase):
    cdef double voigt_modulus_scale
    cdef double voigt_viscosity_scale
    cdef double alpha
    cdef double alpha_factorial
    cdef double zeta
    cdef double complex sine_term
