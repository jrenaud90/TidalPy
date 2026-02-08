# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from TidalPy.utilities.lookups import IntMap2, IntMap4

def global_potential(
        double planet_radius,
        double semi_major_axis,
        double orbital_frequency,
        double spin_frequency,
        double obliquity,
        double eccentricity,
        double host_mass,
        double G_to_use,
        int min_degree_l=2,
        int max_degree_l=2,
        object obliquity_truncation='gen',
        int eccentricity_truncation=6
    ):

    # Clean up non-C inputs
    cdef int i_obliquity_truncation = 0
    if isinstance(obliquity_truncation, str):
        if obliquity_truncation.lower() in ('gen', 'general'):
            i_obliquity_truncation = 10
        elif obliquity_truncation.lower() in ('off',):
            i_obliquity_truncation = 0
        else:
            try:
                i_obliquity_truncation = int(obliquity_truncation)
            except ValueError:
                raise ValueError("Unexpected obliquity truncation encountered.")
    elif isinstance(obliquity_truncation, int):
        i_obliquity_truncation = obliquity_truncation
    if i_obliquity_truncation not in (0, 2, 4, 10):
        raise NotImplementedError("Unsupported obliquity truncation encountered.")

    # Run C++ code
    cdef c_GlobalPotentialResult c_result = c_global_potential(
        planet_radius,
        semi_major_axis,
        orbital_frequency,
        spin_frequency,
        obliquity,
        eccentricity,
        host_mass,
        G_to_use,
        min_degree_l,
        max_degree_l,
        i_obliquity_truncation,
        eccentricity_truncation
    )

        struct PotentialResultAtMode:
        double dU_dM
        double dU_dw
        double dU_dO
        double E_dot
        PotentialResultAtMode(
            double dU_dM_,
            double dU_dw_,
            double dU_dO_,
            double E_dot_)


    struct c_GlobalPotentialResult:
        ModeMap mode_map
        cUniqueFreqIndexMap unique_freq_index_map
        c_UniqueFreqMap unique_freq_map
        c_IntMap[c_Key4, PotentialResultAtMode] potential_map
        int error_code
        int working_on_l
