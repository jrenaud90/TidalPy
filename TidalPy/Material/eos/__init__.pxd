# Common functions and structures
from TidalPy.Material.eos.ode cimport EOS_ODEInput, eos_diffeq
from TidalPy.Material.eos.solver cimport EOSSolutionVec, solve_eos

# Specific EOS functions and structures
# Interpolate
from TidalPy.Material.eos.methods.interpolate cimport InterpolateEOSInput, preeval_interpolate

