# Common functions and structures
from TidalPy.Material_x.eos.eos_solution cimport c_EOSSolution
from TidalPy.Material_x.eos.ode cimport c_EOS_ODEInput, c_eos_diffeq, C_EOS_Y_VALUES, C_EOS_EXTRA_VALUES, C_EOS_DY_VALUES
from TidalPy.Material_x.eos.solver cimport c_solve_eos

# Specific EOS functions and structures
# Interpolate
from TidalPy.Material_x.eos.methods.interpolate cimport c_InterpolateEOSInput, c_preeval_interpolate
