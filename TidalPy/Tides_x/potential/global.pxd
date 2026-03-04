

from libc.stdint cimport int16_t

from TidalPy.utilities.lookups cimport c_IntMap, c_Key2, c_Key4
from TidalPy.Tides_x.potential.potential_common cimport c_ModeMap, c_UniqueFreqIndexMap, c_UniqueFreqMap


cdef extern from "global_.hpp" nogil:
    struct c_GlobalPotentialResultAtMode:
        double dU_dM
        double dU_dw
        double dU_dO
        double E_dot
        c_GlobalPotentialResultAtMode(
            double dU_dM_,
            double dU_dw_,
            double dU_dO_,
            double E_dot_)


    struct c_GlobalPotentialStorage:
        c_ModeMap mode_map
        c_UniqueFreqIndexMap unique_freq_index_map
        c_UniqueFreqMap unique_freq_map
        c_IntMap[c_Key4, c_GlobalPotentialResultAtMode] potential_map
        int error_code
        int working_on_l

    c_GlobalPotentialStorage c_global_potential(
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
        int eccentricity_truncation
    )
