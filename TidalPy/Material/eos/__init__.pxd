# Common functions and structures
from TidalPy.Material.eos.eos_solution cimport EOSSolutionCC
from TidalPy.Material.eos.ode cimport EOS_ODEInput, eos_diffeq
from TidalPy.Material.eos.solver cimport solve_eos

# Specific EOS functions and structures
# Interpolate
from TidalPy.Material.eos.methods.interpolate cimport InterpolateEOSInput, preeval_interpolate
