# distutils: language = c++

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase
from TidalPy.viscosity_x.viscosity cimport c_ViscosityBase
from TidalPy.partial_melt_x.partial_melt cimport c_PartialMeltBase


cdef extern from "material_eos_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_MaterialState:
        double density
        double melt_fraction
        double shear_modulus
        double bulk_modulus
        double shear_viscosity
        double bulk_viscosity

    cdef cppclass c_MaterialEOSBase(c_PhysicsBase):
        double calc_density(double pressure, double temperature, double radius) const
        double calc_bulk_modulus(double pressure, double temperature, double radius) const
        void calc_material_state(
            double pressure, double temperature, cpp_bool thermal_density, double radius, c_MaterialState& out) const
        double get_thermal_expansion() const
        double get_reference_temperature() const
        double get_tabulated_shear_modulus(double radius) const
        double get_tabulated_bulk_modulus(double radius) const
        double get_tabulated_shear_viscosity(double radius) const
        double get_tabulated_bulk_viscosity(double radius) const
        # Static constants and the shear law
        double get_shear_modulus_static() const
        double get_bulk_modulus_static() const
        double get_shear_viscosity_static() const
        double get_bulk_viscosity_static() const
        double get_shear_modulus_pressure_derivative() const
        double get_shear_modulus_temperature_derivative() const
        double get_shear_modulus_reference_temperature() const
        double get_thermal_conductivity() const
        double get_heat_capacity() const
        double calc_thermal_diffusivity(double density) const
        void set_shear_modulus_static(double value)
        void set_bulk_modulus_static(double value)
        void set_shear_viscosity_static(double value)
        void set_bulk_viscosity_static(double value)
        # Owned viscosity and partial-melt models
        void set_shear_viscosity(unique_ptr[c_ViscosityBase] model)
        void set_bulk_viscosity(unique_ptr[c_ViscosityBase] model)
        void set_partial_melt(unique_ptr[c_PartialMeltBase] model)
        c_ViscosityBase*   get_shear_viscosity_model() const
        c_ViscosityBase*   get_bulk_viscosity_model() const
        c_PartialMeltBase* get_partial_melt_model() const

    cdef cppclass c_MaterialEOSConfig:
        double reference_density
        double reference_bulk_modulus
        double bulk_modulus_derivative
        double thermal_expansion
        double reference_temperature
        double invert_rtol
        int    invert_max_iters
        double shear_modulus_static
        double bulk_modulus_static
        double shear_viscosity_static
        double bulk_viscosity_static
        double shear_modulus_pressure_derivative
        double shear_modulus_temperature_derivative
        double shear_modulus_reference_temperature
        double thermal_conductivity
        double heat_capacity
        vector[double] radius
        vector[double] density
        vector[double] shear_modulus
        vector[double] bulk_modulus
        vector[double] shear_viscosity
        vector[double] bulk_viscosity

    cdef cppclass c_ConstantDensityEOS(c_MaterialEOSBase):
        c_ConstantDensityEOS() except +
        c_ConstantDensityEOS(const c_MaterialEOSConfig& cfg) except +
        double get_reference_density() const

    cdef cppclass c_BirchMurnaghanEOS(c_MaterialEOSBase):
        c_BirchMurnaghanEOS() except +
        c_BirchMurnaghanEOS(const c_MaterialEOSConfig& cfg) except +
        double get_reference_density()       const
        double get_reference_bulk_modulus()  const
        double get_bulk_modulus_derivative() const
        double get_invert_rtol()             const
        int    get_invert_max_iters()        const

    cdef cppclass c_VinetEOS(c_MaterialEOSBase):
        c_VinetEOS() except +
        c_VinetEOS(const c_MaterialEOSConfig& cfg) except +
        double get_reference_density()       const
        double get_reference_bulk_modulus()  const
        double get_bulk_modulus_derivative() const
        double get_invert_rtol()             const
        int    get_invert_max_iters()        const

    cdef cppclass c_InterpolatedEOS(c_MaterialEOSBase):
        c_InterpolatedEOS() except +
        c_InterpolatedEOS(const c_MaterialEOSConfig& cfg) except +
        size_t get_num_points() const
        cpp_bool has_shear_modulus() const
        cpp_bool has_bulk_modulus() const
        cpp_bool has_shear_viscosity() const
        cpp_bool has_bulk_viscosity() const

    double eos_bm_pressure(double eta, double K0, double K0_prime)
    double eos_vinet_pressure(double eta, double K0, double K0_prime)

    cdef enum class c_MaterialEOSModel:
        Constant
        BirchMurnaghan
        Vinet
        Interpolated

    c_MaterialEOSModel c_material_eos_model_from_name(const string& model_name) except +
    unique_ptr[c_MaterialEOSBase] c_find_material_eos(
        c_MaterialEOSModel model, const c_MaterialEOSConfig& cfg) except +


# The EOS model's own config entries plus a sub-table for each attached viscosity or partial-melt model.
# Shared with the layer wrappers, which emit it as their ``material`` table.
cdef dict cy_material_config(const c_MaterialEOSBase* eos_ptr)


cdef class MaterialEOSBase(PhysicsBase):
    cdef unique_ptr[c_MaterialEOSBase] _eos_ptr   # owns the most-derived C++ model; the typed pointers below do not
    cdef c_MaterialEOSBase* _model(self) except NULL


cdef class ConstantDensityEOS(MaterialEOSBase):
    cdef c_ConstantDensityEOS* _constant_ptr


cdef class BirchMurnaghanEOS(MaterialEOSBase):
    cdef c_BirchMurnaghanEOS* _bm_ptr


cdef class VinetEOS(MaterialEOSBase):
    cdef c_VinetEOS* _vinet_ptr


cdef class InterpolatedEOS(MaterialEOSBase):
    cdef c_InterpolatedEOS* _interp_ptr
