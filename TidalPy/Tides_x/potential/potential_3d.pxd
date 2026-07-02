# distutils: language = c++
"""
potential_3d.pxd
Cython declarations for the dynamic Kaula 3D tidal-potential engine (potential_3d_.hpp).
"""

from libcpp.vector cimport vector


cdef extern from "potential_point_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_PotentialPoint:
        double U
        double dU_dtheta
        double dU_dphi
        double d2U_dtheta2
        double d2U_dphi2
        double d2U_dtheta_dphi


cdef extern from "potential_3d_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_TidalPotential3DMode:
        int degree_l
        double mode_frequency
        c_PotentialPoint potential

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
        double time,
        int* error_code) except +
