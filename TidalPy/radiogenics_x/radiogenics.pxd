# distutils: language = c++
"""
radiogenics.pxd
Cython declarations for TidalPy's radiogenics model hierarchy.

Exports the three C++ radiogenics models, the combined config struct, and the
Python wrapper classes so other extensions can cimport and build/attach
radiogenics models at C speed.

Usage::

    from TidalPy.radiogenics_x.radiogenics cimport (
        RadiogenicsBase, IsotopeRadiogenics,
        c_RadiogenicsBase, c_IsotopeRadiogenics, c_RadiogenicsConfig)
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "radiogenics_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_RadiogenicsBase(c_PhysicsBase):
        double calc_heating(double time_s, double mass_kg) const
        void calc_heating_vectorize_time(
            const vector[double]& time_s,
            double mass_kg,
            vector[double]& out_heating) except +
        void calc_heating_vectorize_mass(
            double time_s,
            const vector[double]& mass_kg,
            vector[double]& out_heating) except +
        void calc_heating_vectorize_all(
            const vector[double]& time_s,
            const vector[double]& mass_kg,
            vector[double]& out_heating) except +


cdef extern from "radiogenics_.hpp" namespace "tidalpy" nogil:

    # A single radioactive isotope (value type; no base class).
    cdef cppclass c_Isotope:
        c_Isotope() except +
        c_Isotope(string name, double hpr_w_kg, double half_life_s,
                  double mass_frac, double concentration) except +
        string name
        double heat_production_w_kg
        double half_life_s
        double mass_frac
        double concentration
        double decay_constant() const
        double specific_heating(double time_s, double ref_time_s) const

    # A named, literature-sourced set of isotopes plus its reference time.
    cdef cppclass c_IsotopeDataset:
        vector[c_Isotope] isotopes
        double ref_time_s

    cdef cppclass c_RadiogenicsConfig:
        vector[c_Isotope] isotopes
        double fixed_heat_production_w_kg
        double average_half_life_s
        double ref_time_s

    # Built-in isotope dataset catalog (raises ValueError on unknown name).
    c_IsotopeDataset c_get_isotope_dataset(const string& name) except +
    vector[string] c_isotope_dataset_names()

    cdef cppclass c_OffRadiogenics(c_RadiogenicsBase):
        c_OffRadiogenics() except +
        c_OffRadiogenics(const c_RadiogenicsConfig& cfg) except +

    cdef cppclass c_IsotopeRadiogenics(c_RadiogenicsBase):
        c_IsotopeRadiogenics() except +
        c_IsotopeRadiogenics(const c_RadiogenicsConfig& cfg) except +
        const vector[c_Isotope]& get_isotopes() const
        double get_ref_time()      const
        size_t get_num_isotopes()  const

    cdef cppclass c_FixedRadiogenics(c_RadiogenicsBase):
        c_FixedRadiogenics() except +
        c_FixedRadiogenics(const c_RadiogenicsConfig& cfg) except +
        double get_fixed_heat_production() const
        double get_average_half_life()     const
        double get_ref_time()              const

    # Enum naming each radiogenics model.
    cdef enum class c_RadiogenicsModel:
        Off
        Isotope
        Fixed

    # Map a name/alias to the enum (raises ValueError on unknown name).
    c_RadiogenicsModel c_radiogenics_model_from_name(const string& model_name) except +

    # Build the model named by the enum; returns an owning unique_ptr.
    unique_ptr[c_RadiogenicsBase] c_find_radiogenics(
        c_RadiogenicsModel model, const c_RadiogenicsConfig& cfg) except +


# =====================================================================================================================
# Cython wrapper class declarations
# =====================================================================================================================
cdef class RadiogenicsBase(PhysicsBase):
    cdef unique_ptr[c_RadiogenicsBase] _radiogenics_ptr   # owns the most-derived C++ model object
    cpdef dict get_config_dict(self)


cdef class OffRadiogenics(RadiogenicsBase):
    pass


cdef class IsotopeRadiogenics(RadiogenicsBase):
    cdef c_IsotopeRadiogenics* _isotope_ptr   # non-owning; ownership via RadiogenicsBase._radiogenics_ptr
    cpdef dict get_config_dict(self)


cdef class FixedRadiogenics(RadiogenicsBase):
    cdef c_FixedRadiogenics* _fixed_ptr   # non-owning; ownership via RadiogenicsBase._radiogenics_ptr
    cpdef dict get_config_dict(self)
