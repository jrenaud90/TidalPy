# Common functions and structures
from TidalPy.Material.eos.common cimport EOSOutput
from TidalPy.Material.eos.ode cimport EOS_ODEInput, eos_diffeq
from TidalPy.Material.eos.solver cimport EOSSolutionVec, GlobalEOSSolutionStorage, solve_eos

# Specific EOS functions and structures
# Interpolate
from TidalPy.Material.eos.methods.interpolate cimport InterpolateEOSInput, preeval_interpolate

