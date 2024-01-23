from TidalPy.rheology.base cimport RheologyModelBase

########################################################################################################################
######################################## New Rheological Models Go Below Here! #########################################
########################################################################################################################

cdef class Elastic(RheologyModelBase):
    pass

cdef class Newton(RheologyModelBase):
    pass

cdef class Maxwell(RheologyModelBase):
    pass

cdef class Voigt(RheologyModelBase):
    cdef double voigt_modulus_scale
    cdef double voigt_viscosity_scale

cdef class Burgers(RheologyModelBase):
    cdef double voigt_modulus_scale
    cdef double voigt_viscosity_scale

cdef class Andrade(RheologyModelBase):
    cdef double alpha
    cdef double alpha_factorial
    cdef double zeta
    cdef double complex sine_term

cdef class SundbergCooper(RheologyModelBase):
    cdef double voigt_modulus_scale
    cdef double voigt_viscosity_scale
    cdef double alpha
    cdef double alpha_factorial
    cdef double zeta
    cdef double complex sine_term
