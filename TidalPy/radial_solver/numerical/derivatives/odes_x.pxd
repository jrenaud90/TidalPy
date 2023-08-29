from TidalPy.radial_solver.numerical.derivatives.ode_base_x cimport RadialSolverBase


cdef class SolidDynamicCompressible(RadialSolverBase):
    pass

cdef class SolidDynamicIncompressible(RadialSolverBase):
    pass

cdef class SolidStaticCompressible(RadialSolverBase):
    pass

cdef class SolidStaticIncompressible(RadialSolverBase):
    pass

cdef class LiquidDynamicCompressible(RadialSolverBase):
    pass

cdef class LiquidDynamicIncompressible(RadialSolverBase):
    pass

cdef class LiquidStaticCompressible(RadialSolverBase):
    pass

cdef class LiquidStaticIncompressible(RadialSolverBase):
    pass
