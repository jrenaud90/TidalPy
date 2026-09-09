# distutils: language = c++
"""
material.pxd
Cython declarations for the material-EOS pre-evaluation used by the world-level
EOS solve. Exposes the per-layer input struct and the CyRK-compatible pre-eval
function that dispatches to a layer's c_MaterialEOSBase model.
"""

from TidalPy.Material_x.eos.material_eos cimport c_MaterialEOSBase


cdef extern from "material_.hpp" nogil:

    cdef struct c_MaterialEOSInput:
        c_MaterialEOSBase* eos_model_ptr
        double             temperature_k

    cdef void c_preeval_material_eos(
            # Values that will be updated by the function
            char* preeval_output,
            # Input that is used by the pre-eval
            double radius,
            double* radial_solutions,
            char* preeval_input
            ) noexcept nogil
