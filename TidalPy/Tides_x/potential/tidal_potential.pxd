# distutils: language = c++
"""
tidal_potential.pxd
Cython declarations for TidalPy's tidal-potential truncation model hierarchy.

Exports the C++ models, the config struct, the per-mode output structs, and the enum factory so other
extensions (the 3D dissipation kernel / the world) can cimport and build or evaluate a tidal potential.
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp cimport bool as cpp_bool

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase


cdef extern from "potential_point_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_PotentialPoint:
        double U
        double dU_dtheta
        double dU_dphi
        double d2U_dtheta2
        double d2U_dphi2
        double d2U_dtheta_dphi


cdef extern from "tidal_potential_base_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_TidalPotentialState:
        double orbital_frequency
        double spin_frequency
        double eccentricity
        double obliquity
        double host_mass
        double semi_major_axis

    cdef cppclass c_TidalPotentialModeSet:
        int num_modes
        double mode_frequency[9]
        c_PotentialPoint potential[9]

    cdef cppclass c_TidalPotentialBase(c_PhysicsBase):
        c_TidalPotentialModeSet calc_modes(
            const c_TidalPotentialState& state,
            double radius,
            double colatitude,
            double longitude,
            double time) const
        int num_modes() const


cdef extern from "tidal_potential_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_TidalPotentialConfig:
        cpp_bool use_static

    cdef cppclass c_SyncLowEPotential(c_TidalPotentialBase):
        c_SyncLowEPotential() except +
        c_SyncLowEPotential(const c_TidalPotentialConfig& cfg) except +

    cdef cppclass c_NSRModesPotential(c_TidalPotentialBase):
        c_NSRModesPotential() except +
        c_NSRModesPotential(const c_TidalPotentialConfig& cfg) except +
        cpp_bool get_use_static() const

    cdef enum class c_TidalPotentialModel:
        SyncLowE
        NSRModes

    c_TidalPotentialModel c_tidal_potential_model_from_name(const string& model_name) except +
    unique_ptr[c_TidalPotentialBase] c_find_tidal_potential(
        c_TidalPotentialModel model, const c_TidalPotentialConfig& cfg) except +


cdef class TidalPotentialBase(PhysicsBase):
    cdef unique_ptr[c_TidalPotentialBase] _potential_ptr   # owns the most-derived C++ model object
    cpdef dict get_config_dict(self)


cdef class SyncLowEPotential(TidalPotentialBase):
    cdef c_SyncLowEPotential* _sync_ptr                    # non-owning; ownership via _potential_ptr


cdef class NSRModesPotential(TidalPotentialBase):
    cdef c_NSRModesPotential* _nsr_ptr                     # non-owning; ownership via _potential_ptr
    cpdef dict get_config_dict(self)
