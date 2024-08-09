from CyRK.utils.memory cimport shared_ptr, make_shared
from CyRK.utils.vector cimport vector

from TidalPy.Material.eos.ode cimport eos_solution, EOS_ODEInput

cdef struct GlobalEOSSolutionStorage:
    vector[double] radius
    vector[double] pressure
    vector[double] gravity
    vector[double] density
    vector[double] static_shear_modulus
    vector[double] shear_viscosity
    vector[double complex] complex_shear_modulus
    vector[double] static_bulk_modulus
    vector[double] bulk_viscosity
    vector[double complex] complex_bulk_modulus

ctypedef shared_ptr[GlobalEOSSolutionStorage] GlobalEOSSolution
