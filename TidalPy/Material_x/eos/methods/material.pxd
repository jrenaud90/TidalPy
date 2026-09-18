# distutils: language = c++
"""Cython declarations for the material-EOS pre-evaluation used by the world EOS solve."""

from TidalPy.Material_x.eos.material_eos cimport c_MaterialEOSBase


cdef extern from "material_.hpp" nogil:

    cdef struct c_MaterialEOSInput:
        c_MaterialEOSBase* eos_model_ptr
        double             temperature
        double             length_scale
        double             pascal_scale
        double             density_scale

    cdef void c_preeval_material_eos(
            char* preeval_output,
            double radius,
            double* radial_solutions,
            char* preeval_input
            ) noexcept nogil
