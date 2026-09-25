# distutils: language = c++

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase


cdef extern from "radiogenics_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_RadiogenicsBase(c_PhysicsBase):
        double calc_heating(double time, double mass) const
        void calc_heating_vectorize(
            const vector[double]& time,
            const vector[double]& mass,
            vector[double]& out_heating) except +


cdef extern from "radiogenics_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_Isotope:
        c_Isotope() except +
        c_Isotope(
            string name,
            double hpr,
            double half_life,
            double mass_frac,
            double concentration) except +
        string name
        double heat_production
        double half_life
        double mass_frac
        double concentration
        double decay_constant() const
        double specific_heating(double time, double ref_time) const

    cdef cppclass c_IsotopeDataset:
        vector[c_Isotope] isotopes
        double ref_time

    cdef cppclass c_RadiogenicsConfig:
        vector[c_Isotope] isotopes
        double fixed_heat_production
        double average_half_life
        double ref_time

    # Raises ValueError on an unknown dataset name.
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

    cdef enum class c_RadiogenicsModel:
        Off
        Isotope
        Fixed

    # Raises ValueError on an unknown name.
    c_RadiogenicsModel c_radiogenics_model_from_name(const string& model_name) except +

    unique_ptr[c_RadiogenicsBase] c_find_radiogenics(
        c_RadiogenicsModel model, const c_RadiogenicsConfig& cfg) except +


cdef class RadiogenicsBase(PhysicsBase):
    cdef unique_ptr[c_RadiogenicsBase] _radiogenics_ptr   # owns the most-derived C++ model object
    cdef void _adopt(self, unique_ptr[c_RadiogenicsBase]& model) noexcept


cdef class OffRadiogenics(RadiogenicsBase):
    pass


cdef class IsotopeRadiogenics(RadiogenicsBase):
    pass


cdef class FixedRadiogenics(RadiogenicsBase):
    pass
