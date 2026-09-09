# distutils: language = c++
"""
potential_3d.pxd
Cython declarations for the dynamic Kaula 3D tidal-potential engine (potential_3d_.hpp).

The engine returns each active mode's degree, signed forcing frequency, and the COMPLEX potential
angular-factor amplitude (the mode's e^{i omega t} pulled out) used by the secular 3D heating.
"""

from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex


cdef extern from "potential_3d_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_PotentialPointC:
        cpp_complex[double] U
        cpp_complex[double] dU_dtheta
        cpp_complex[double] dU_dphi
        cpp_complex[double] d2U_dtheta2
        cpp_complex[double] d2U_dphi2
        cpp_complex[double] d2U_dtheta_dphi

    cdef cppclass c_TidalPotential3DMode:
        int degree_l
        double mode_frequency
        c_PotentialPointC potential

    vector[c_TidalPotential3DMode] c_tidal_potential_3d_modes(
        double planet_radius,
        double semi_major_axis,
        double orbital_frequency,
        double spin_frequency,
        double obliquity,
        double eccentricity,
        double host_mass,
        double G_to_use,
        int min_degree_l,
        int max_degree_l,
        int obliquity_truncation,
        int eccentricity_truncation,
        double colatitude,
        double longitude,
        int* error_code) except +
