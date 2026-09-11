# distutils: language = c++
"""Cython declarations for the native `_x` radial-solver input builders (build_inputs_.hpp)."""

from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex
from libcpp.vector cimport vector

from TidalPy.rheology_x.rheology cimport c_RheologyBase


cdef extern from "build_inputs_.hpp" namespace "tidalpy" nogil:

    cdef size_t C_RS_MIN_SLICES_PER_LAYER

    cdef cppclass c_RadialSolverInputs:
        vector[double] radius_m
        vector[double] density_kg_m3
        vector[cpp_complex[double]] complex_bulk_modulus_pa
        vector[cpp_complex[double]] complex_shear_modulus_pa
        vector[double] upper_radius_bylayer_m
        vector[size_t] slices_bylayer
        double forcing_frequency_rad_s
        double planet_bulk_density_kg_m3

    void c_thickness_from_radius_fractions(
        const vector[double]& radius_fraction_bylayer,
        vector[double]& out_thickness_fraction_bylayer) except +

    void c_thickness_from_volume_fractions(
        double planet_radius_m,
        const vector[double]& volume_fraction_bylayer,
        vector[double]& out_thickness_fraction_bylayer) except +

    void c_build_rs_input_homogeneous_layers(
        double planet_radius_m,
        double forcing_frequency_rad_s,
        const vector[double]& density_bylayer,
        const vector[double]& static_bulk_modulus_bylayer,
        const vector[double]& static_shear_modulus_bylayer,
        const vector[double]& bulk_viscosity_bylayer,
        const vector[double]& shear_viscosity_bylayer,
        const vector[double]& thickness_fraction_bylayer,
        const vector[size_t]& slices_bylayer,
        const vector[const c_RheologyBase*]& shear_rheology_bylayer,
        const vector[const c_RheologyBase*]& bulk_rheology_bylayer,
        c_RadialSolverInputs& out) except +

    void c_build_rs_input_from_data(
        double forcing_frequency_rad_s,
        const vector[double]& radius_m,
        const vector[double]& density_kg_m3,
        const vector[double]& static_bulk_modulus_pa,
        const vector[double]& static_shear_modulus_pa,
        const vector[double]& bulk_viscosity_pas,
        const vector[double]& shear_viscosity_pas,
        const vector[double]& layer_upper_radius_bylayer_m,
        const vector[const c_RheologyBase*]& shear_rheology_bylayer,
        const vector[const c_RheologyBase*]& bulk_rheology_bylayer,
        cpp_bool warnings,
        c_RadialSolverInputs& out) except +


# Resolve a rheology argument (one model, one model name, or a per-layer sequence of either) into
# per-layer C++ pointers. Returns the Python objects that own those pointers; keep the returned list
# alive for as long as the pointers are used.
cdef list cy_resolve_rheology_bylayer(
    object models, size_t num_layers, str argument_name, vector[const c_RheologyBase*]& out_ptrs)
