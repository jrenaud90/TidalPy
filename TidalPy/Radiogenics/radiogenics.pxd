# distutils: language = c++

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities.classes.classes cimport PhysicsBase, c_PhysicsBase, c_ParamMap


cdef extern from "radiogenics_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_RadiogenicsBase(c_PhysicsBase):
        double calc_heating(double time, double mass) const
        double get_ref_time() const
        void calc_heating_vectorize(
            const vector[double]& time,
            const vector[double]& mass,
            vector[double]& out_heating) except +


cdef extern from "radiogenics_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_Isotope:
        string name
        double heat_production
        double half_life
        double mass_frac
        double concentration

    cdef cppclass c_IsotopeDataset:
        vector[c_Isotope] isotopes
        double ref_time

    # Raises ValueError on an unknown dataset name.
    c_IsotopeDataset c_get_isotope_dataset(const string& name) except +
    vector[string] c_isotope_dataset_names()

    cdef cppclass c_IsotopeRadiogenics(c_RadiogenicsBase):
        size_t get_num_isotopes() const
        const vector[string]& get_isotope_names() const

    unique_ptr[c_RadiogenicsBase] c_find_radiogenics(const string& model_name, const c_ParamMap& params) except +
    unique_ptr[c_RadiogenicsBase] c_make_isotope_radiogenics(
        const c_ParamMap& params, const vector[string]& isotope_names) except +
    string c_radiogenics_canonical_name(const string& model_name) except +
    vector[string] c_radiogenics_model_names() except +


cdef class RadiogenicsBase(PhysicsBase):
    cdef c_RadiogenicsBase* _radiogenics(self) except NULL


cdef class OffRadiogenics(RadiogenicsBase):
    pass


cdef class IsotopeRadiogenics(RadiogenicsBase):
    pass


cdef class FixedRadiogenics(RadiogenicsBase):
    pass
