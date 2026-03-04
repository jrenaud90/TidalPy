# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from TidalPy.constants cimport tidalpy_config_ptr, get_shared_config_address
from TidalPy.Tides_x.potential.potential_common cimport ModeMap, UniqueFrequencyMap
from TidalPy.Tides_x.potential.potential_common import ModeMap, UniqueFrequencyMap

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

    # Ensure the global config pointer is initialized before calling C++ code.
    global tidalpy_config_ptr
    if tidalpy_config_ptr == NULL:
        tidalpy_config_ptr = get_shared_config_address()

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
    cdef c_GlobalPotentialStorage c_result = c_global_potential(
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

    # Check for errors
    if c_result.error_code != 0:
        if c_result.error_code == -20:
            raise NotImplementedError(
                f"Global potential error code -20: Could not find l,m coefficient "
                f"(working on degree l={c_result.working_on_l}). Perhaps unsupported degree l.")
        elif c_result.error_code == -1:
            raise NotImplementedError(
                f"Global potential error code -1: Unsupported obliquity truncation "
                f"(working on degree l={c_result.working_on_l}).")
        elif c_result.error_code == -2:
            raise NotImplementedError(
                f"Global potential error code -2: Unsupported degree l={c_result.working_on_l} "
                f"for obliquity functions.")
        else:
            raise RuntimeError(
                f"Unknown global potential error code: {c_result.error_code} "
                f"(working on degree l={c_result.working_on_l}).")

    # Convert C++ results to Python-accessible objects.
    # Mode map
    cdef ModeMap mode_map = ModeMap()
    mode_map._cinst = c_result.mode_map

    # Unique frequency index map
    cdef UniqueFrequencyMap unique_freq_index_map = UniqueFrequencyMap()
    unique_freq_index_map._cinst = c_result.unique_freq_index_map

    # Unique frequency map (vector of c_FrequencyStorage -> list of (frequency, num_instances))
    cdef list unique_freq_list = []
    cdef size_t i
    for i in range(c_result.unique_freq_map.size()):
        unique_freq_list.append(
            (c_result.unique_freq_map[i].frequency, c_result.unique_freq_map[i].num_instances)
        )

    # Potential map (c_IntMap<c_Key4, c_GlobalPotentialResultAtMode> -> dict)
    cdef dict potential_dict = dict()
    cdef c_Key4 pkey
    for i in range(c_result.potential_map.size()):
        pkey = c_result.potential_map.data[i].first
        potential_dict[(pkey.a, pkey.b, pkey.c, pkey.d)] = (
            c_result.potential_map.data[i].second.dU_dM,
            c_result.potential_map.data[i].second.dU_dw,
            c_result.potential_map.data[i].second.dU_dO,
            c_result.potential_map.data[i].second.E_dot
        )

    return mode_map, unique_freq_index_map, unique_freq_list, potential_dict
