#pragma once

#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <memory>
#include <string>

#include "cysolution.hpp"  // Part of the CyRK python package. Header should be included in setup using CyRK.get_include()

#include "nondimensional_.hpp" // Part of the TidalPy.utilities module
#include "constants_.hpp"      // Part of the TidalPy

#include "ode_.hpp" // For C_EOS_Y_VALUES, C_EOS_EXTRA_VALUES, C_EOS_DY_VALUES
#include "../../utilities/arrays/interp_.hpp"  // cf_binary_search_with_guess, cf_interp, cf_interp_complex


/// Helper: Check if two doubles are approximately equal.
inline bool c_eos_isclose(double a, double b)
{
    // TODO: Make the numerics.pyx utility a C header that we can import here instead of reimplementing it.
    const double rtol = 1.0e-9;
    const double atol = 0.0;

    if (std::isnan(a))
    {
        return false;
    }

    if (std::isnan(b))
    {
        return false;
    }

    if (a == b)
    {
        return true;
    }

    const double lhs = std::fabs(a - b);
    const double rhs = std::fmax(rtol * std::fmax(std::fabs(a), std::fabs(b)), atol);

    return lhs <= rhs;
}


/// C++ class storing the equation of state integration results for a layered planet.
///
/// Stores CyRK integration results for each layer, and provides methods for interpolating
/// the full planet's gravity, pressure, mass, moment of inertia, density, and complex moduli.
class c_EOSSolution
{

// Attributes
protected:

public:
    int iterations             = -1;
    int error_code             = -100;
    int nondim_status          = 0;
    int solution_nondim_status = 0;
    bool success          = false;
    bool max_iters_hit    = false;
    bool radius_array_set  = false;
    bool other_vecs_set    = false;
    bool p_use_array_interp = false;  // set by inject_from_world_eos; call() uses array interpolation

    // Optional dense structure source (non-owning). When set together with p_use_array_interp, the solved structure
    // dependent variables (gravity/pressure/mass/moi) are read from this external dense EOS solution rather than the
    // stored arrays, preserving the dense (polynomial) accuracy of the integrated EOS structure. Density and the
    // complex moduli remain array-interpolated (density is algebraic; the complex moduli come from the rheology, not
    // the EOS solve). Used by the world Love-number solve, which reuses the world's already-solved (dense) EOS while
    // supplying frequency-dependent complex moduli. The source is dimensional (SI); the scales convert to this
    // solution's (non-dim) units: source_radius = this_radius * p_structure_length_scale; this_value = source_value /
    // p_structure_*_scale.
    const c_EOSSolution* p_structure_dense_source = nullptr;
    double p_structure_length_scale  = 1.0;
    double p_structure_gravity_scale = 1.0;
    double p_structure_pascal_scale  = 1.0;
    double p_structure_mass_scale    = 1.0;
    double p_structure_moi_scale     = 1.0;

    std::string message         = "No Message Set.";
    size_t current_layers_saved = 0;
    size_t num_layers           = 0;
    size_t radius_array_size    = 0;
    size_t num_cyolver_calls    = 0;

    double pressure_error   = TidalPyConstants::d_NAN;
    double surface_gravity  = TidalPyConstants::d_NAN;
    double surface_pressure = TidalPyConstants::d_NAN;
    double central_pressure = TidalPyConstants::d_NAN;
    double radius           = TidalPyConstants::d_NAN;
    double mass             = TidalPyConstants::d_NAN;
    double moi              = TidalPyConstants::d_NAN;

    double redim_length_scale  = TidalPyConstants::d_NAN;
    double redim_gravity_scale = TidalPyConstants::d_NAN;
    double redim_mass_scale    = TidalPyConstants::d_NAN;
    double redim_density_scale = TidalPyConstants::d_NAN;
    double redim_moi_scale     = TidalPyConstants::d_NAN;
    double redim_pascal_scale  = TidalPyConstants::d_NAN;

    // Store results from CyRK's cysolve_ivp.
    std::vector<double> upper_radius_bylayer_vec  = std::vector<double>();
    std::vector<size_t> steps_taken_vec           = std::vector<size_t>();
    std::vector<std::unique_ptr<CySolverResult>> cysolver_results_uptr_bylayer_vec = std::vector<std::unique_ptr<CySolverResult>>();

    // Copy of user-provided radius array
    std::vector<double> radius_array_vec = std::vector<double>();

    // Store interpolated arrays based on provided radius array
    std::vector<double> gravity_array_vec  = std::vector<double>();
    std::vector<double> pressure_array_vec = std::vector<double>();
    std::vector<double> mass_array_vec     = std::vector<double>();
    std::vector<double> moi_array_vec      = std::vector<double>();
    std::vector<double> density_array_vec  = std::vector<double>();
    std::vector<std::complex<double>> complex_shear_array_vec = std::vector<std::complex<double>>();
    std::vector<std::complex<double>> complex_bulk_array_vec  = std::vector<std::complex<double>>();
    // Static (real) shear/bulk viscosity [Pa s] vs radius (EOS-model extra outputs).
    std::vector<double> shear_viscosity_array_vec = std::vector<double>();
    std::vector<double> bulk_viscosity_array_vec  = std::vector<double>();

// Methods
protected:

public:

    virtual ~c_EOSSolution()
    {
        // Reset each unique pointer in the cysolver vector
        for (size_t i = 0; i < this->cysolver_results_uptr_bylayer_vec.size(); i++)
        {
            this->cysolver_results_uptr_bylayer_vec[i]->dense_vec.clear();
            // Reset this class reference to the unique pointer.
            this->cysolver_results_uptr_bylayer_vec[i].reset();
        }
        // Clear all vectors
        this->cysolver_results_uptr_bylayer_vec.clear();
        this->upper_radius_bylayer_vec.clear();
        this->radius_array_vec.clear();
        this->gravity_array_vec.clear();
        this->pressure_array_vec.clear();
        this->mass_array_vec.clear();
        this->moi_array_vec.clear();
        this->density_array_vec.clear();
        this->complex_shear_array_vec.clear();
        this->complex_bulk_array_vec.clear();
        this->shear_viscosity_array_vec.clear();
        this->bulk_viscosity_array_vec.clear();
    }

    c_EOSSolution()
    {
    }

    c_EOSSolution(
            double* upper_radius_bylayer_ptr,
            size_t num_layers_,
            double* radius_array_ptr,
            size_t radius_array_size_
        ) :
            error_code(0),
            current_layers_saved(0),
            num_layers(num_layers_)
    {
        this->cysolver_results_uptr_bylayer_vec.reserve(num_layers_);
        this->upper_radius_bylayer_vec.resize(num_layers_);

        // Store the upper radius of each layer to make it easier to call interpolators later
        std::memcpy(this->upper_radius_bylayer_vec.data(), upper_radius_bylayer_ptr, this->num_layers * sizeof(double));

        // Use provided radius array to setup various storage vectors
        this->change_radius_array(radius_array_ptr, radius_array_size_);
    }


    /// Save a CyRK solver result for one layer.
    void save_cyresult(std::unique_ptr<CySolverResult> new_cysolver_result_uptr)
    {
        this->cysolver_results_uptr_bylayer_vec.push_back(std::move(new_cysolver_result_uptr));
        this->current_layers_saved++;
    }


    /// Record the number of integration steps taken for a layer.
    void save_steps_taken(size_t steps_taken)
    {
        this->steps_taken_vec.push_back(steps_taken);
        this->num_cyolver_calls++;
    }


    /// Fallback: interpolate all EOS outputs directly from the stored array vectors.
    /// Used when p_use_array_interp is true (set by inject_from_world_eos).
    /// Arrays must already be in the units the caller expects; no nondim scaling is applied.
    void _call_interp_arrays(const double radius_val, double* y_interp_ptr) const noexcept
    {
        const size_t n  = this->radius_array_size;
        // cf_interp/cf_binary_search_with_guess take non-const double* but only read the data.
        double* r = const_cast<double*>(this->radius_array_vec.data());
        double  r_val = radius_val;   // mutable copy for cf_interp's desired_x_ptr arg

        // Initial index guess from normalised position in the radius range.
        const double r_left  = r[0];
        const double r_right = r[n > 0 ? n - 1 : 0];
        size_t j = 0;
        if (r_right > r_left)
        {
            const double frac = (radius_val - r_left) / (r_right - r_left);
            const double clamped = frac < 0.0 ? 0.0 : (frac > 1.0 ? 1.0 : frac);
            j = static_cast<size_t>(static_cast<double>(n) * clamped);
            if (j >= n) j = n - 1;
        }
        int b_code = 0;
        j = cf_binary_search_with_guess(radius_val, r, n, j, &b_code);

        // Interpolate each quantity in the order that matches the CySolverResult layout:
        //   [0] gravity, [1] pressure, [2] mass, [3] moi, [4] density,
        //   [5] shear_real, [6] shear_imag, [7] bulk_real, [8] bulk_imag
        cf_interp(&r_val, r, const_cast<double*>(this->gravity_array_vec.data()),  n, &j, &y_interp_ptr[0]);
        cf_interp(&r_val, r, const_cast<double*>(this->pressure_array_vec.data()), n, &j, &y_interp_ptr[1]);
        cf_interp(&r_val, r, const_cast<double*>(this->mass_array_vec.data()),     n, &j, &y_interp_ptr[2]);
        cf_interp(&r_val, r, const_cast<double*>(this->moi_array_vec.data()),      n, &j, &y_interp_ptr[3]);
        cf_interp(&r_val, r, const_cast<double*>(this->density_array_vec.data()),  n, &j, &y_interp_ptr[4]);

        double shear_result[2] = {0.0, 0.0};
        cf_interp_complex(
            radius_val, r,
            const_cast<double*>(reinterpret_cast<const double*>(this->complex_shear_array_vec.data())),
            n, &j, shear_result);
        y_interp_ptr[5] = shear_result[0];
        y_interp_ptr[6] = shear_result[1];

        double bulk_result[2] = {0.0, 0.0};
        cf_interp_complex(
            radius_val, r,
            const_cast<double*>(reinterpret_cast<const double*>(this->complex_bulk_array_vec.data())),
            n, &j, bulk_result);
        y_interp_ptr[7] = bulk_result[0];
        y_interp_ptr[8] = bulk_result[1];
    }


    /// Interpolate at a single radius using the CySolverResult from a specific layer.
    void call(
        const size_t layer_index,
        const double radius_val,
        double* y_interp_ptr) const
    {
        if (this->p_use_array_interp) [[unlikely]]
        {
            // EOS was injected from a pre-solved world; use simple array interpolation.
            // Arrays are already in the units inject_from_world_eos provided.
            this->_call_interp_arrays(radius_val, y_interp_ptr);

            // Override the SOLVED structure dependent variables (gravity/pressure/mass/moi) with the dense-output
            // values from the source EOS solution, preserving the polynomial accuracy of the integrated structure.
            // Density (y[4]) and the complex moduli (y[5..8]) keep their array-interpolated values.
            if (this->p_structure_dense_source) [[unlikely]]
            {
                double src_out[C_EOS_DY_VALUES];
                const double src_radius = radius_val * this->p_structure_length_scale;
                this->p_structure_dense_source->call(layer_index, src_radius, src_out);
                y_interp_ptr[0] = src_out[0] / this->p_structure_gravity_scale;   // gravity
                y_interp_ptr[1] = src_out[1] / this->p_structure_pascal_scale;    // pressure
                y_interp_ptr[2] = src_out[2] / this->p_structure_mass_scale;      // mass
                y_interp_ptr[3] = src_out[3] / this->p_structure_moi_scale;       // moment of inertia
            }
            return;
        }

        if (layer_index < this->current_layers_saved) [[likely]]
        {
            this->cysolver_results_uptr_bylayer_vec[layer_index]->call(radius_val, y_interp_ptr);

            if (this->nondim_status == 1)
            {
                // Acceleration due to Gravity
                y_interp_ptr[0] *= this->redim_gravity_scale;

                // Pressure
                y_interp_ptr[1] *= this->redim_pascal_scale;

                // Total mass
                y_interp_ptr[2] *= this->redim_mass_scale;

                // Moment of inertia
                y_interp_ptr[3] *= this->redim_moi_scale;

                // Density
                y_interp_ptr[4] *= this->redim_density_scale;

                // Shear modulus (real and complex)
                y_interp_ptr[5] *= this->redim_pascal_scale;
                y_interp_ptr[6] *= this->redim_pascal_scale;

                // Bulk modulus (real and complex)
                y_interp_ptr[7] *= this->redim_pascal_scale;
                y_interp_ptr[8] *= this->redim_pascal_scale;
            }
            else if (this->nondim_status == -1)
            {
                // Acceleration due to Gravity
                y_interp_ptr[0] /= this->redim_gravity_scale;

                // Pressure
                y_interp_ptr[1] /= this->redim_pascal_scale;

                // Total mass
                y_interp_ptr[2] /= this->redim_mass_scale;

                // Moment of inertia
                y_interp_ptr[3] /= this->redim_moi_scale;

                // Density
                y_interp_ptr[4] /= this->redim_density_scale;

                // Shear modulus (real and complex)
                y_interp_ptr[5] /= this->redim_pascal_scale;
                y_interp_ptr[6] /= this->redim_pascal_scale;

                // Bulk modulus (real and complex)
                y_interp_ptr[7] /= this->redim_pascal_scale;
                y_interp_ptr[8] /= this->redim_pascal_scale;
            }
        }
        else
        {
            throw std::out_of_range("Layer index out of range.");
        }
    }


    /// Prepare storage vectors for a new or changed radius array.
    void change_radius_array(
        double* new_radius_ptr,
        size_t new_radius_size)
    {
        this->radius_array_size = new_radius_size;
        if (this->radius_array_set)
        {
            // Radius array was already set before. Reset the storage array vectors.
            this->radius_array_vec.clear();
            this->gravity_array_vec.clear();
            this->pressure_array_vec.clear();
            this->mass_array_vec.clear();
            this->moi_array_vec.clear();
            this->density_array_vec.clear();
            this->complex_shear_array_vec.clear();
            this->complex_bulk_array_vec.clear();
            this->shear_viscosity_array_vec.clear();
            this->bulk_viscosity_array_vec.clear();
            for (size_t i = 0; i < this->cysolver_results_uptr_bylayer_vec.size(); i++)
            {
                this->cysolver_results_uptr_bylayer_vec[i]->dense_vec.clear();
                this->cysolver_results_uptr_bylayer_vec[i].reset();
            }
            this->cysolver_results_uptr_bylayer_vec.clear();
            this->current_layers_saved = 0;

            // Indicate that all vectors are cleared.
            this->other_vecs_set = false;
        }
        this->radius_array_set = true;

        // Reserve capacity in vectors for new radius array size
        this->gravity_array_vec.reserve(this->radius_array_size);
        this->pressure_array_vec.reserve(this->radius_array_size);
        this->mass_array_vec.reserve(this->radius_array_size);
        this->moi_array_vec.reserve(this->radius_array_size);
        this->density_array_vec.reserve(this->radius_array_size);
        this->complex_shear_array_vec.reserve(this->radius_array_size);
        this->complex_bulk_array_vec.reserve(this->radius_array_size);
        this->shear_viscosity_array_vec.reserve(this->radius_array_size);
        this->bulk_viscosity_array_vec.reserve(this->radius_array_size);

        // Copy over the radius array values
        this->radius_array_vec.resize(this->radius_array_size);
        std::memcpy(this->radius_array_vec.data(), new_radius_ptr, new_radius_size * sizeof(double));

        // Get constants
        this->radius = this->radius_array_vec.back();
    }

    /// Run full planet interpolation through each layer using the stored radius array.
    void interpolate_full_planet()
    {
        this->solution_nondim_status = this->nondim_status;

        if (this->current_layers_saved == 0)
        {
            throw std::runtime_error("No layers have been saved. Can not perform interpolation.");
        }

        size_t current_layer_index        = 0;
        double current_layer_upper_radius = this->upper_radius_bylayer_vec[0];

        double y_interp_arr[C_EOS_DY_VALUES];
        double* y_interp_ptr = &y_interp_arr[0];

        bool ready_for_next_layer = false;

        size_t interface_check = 0;

        for (size_t radius_i = 0; radius_i < this->radius_array_size; radius_i++)
        {
            const double radius_val = this->radius_array_vec[radius_i];

            if (c_eos_isclose(radius_val, current_layer_upper_radius))
            {
                // At the layer's radius. We want to capture it once at interfaces
                // (there will be two of the same radii for interface layers)
                if (interface_check == 1)
                {
                    // No longer in current layer
                    ready_for_next_layer = true;
                }
                interface_check++;
            }
            else if (radius_val > current_layer_upper_radius)
            {
                // No longer in current layer
                ready_for_next_layer = true;
            }

            if (ready_for_next_layer)
            {
                current_layer_index++;
                if (current_layer_index > (this->num_layers - 1))
                {
                    // Outside of the planet
                    break;
                }
                else
                {
                    interface_check = 0;
                    current_layer_upper_radius = this->upper_radius_bylayer_vec[current_layer_index];
                    ready_for_next_layer       = false;
                }
            }

            // Call interpolate using temp array as holder.
            this->cysolver_results_uptr_bylayer_vec[current_layer_index]->call(radius_val, y_interp_ptr);

            // Store results
            this->gravity_array_vec.push_back(y_interp_ptr[0]);
            this->pressure_array_vec.push_back(y_interp_ptr[1]);
            this->mass_array_vec.push_back(y_interp_ptr[2]);
            this->moi_array_vec.push_back(y_interp_ptr[3]);
            this->density_array_vec.push_back(y_interp_ptr[4]);

            // CyRK only deals with doubles so we need to convert the 2 doubles into 1 complex for each 
            //  modulus array.
            this->complex_shear_array_vec.push_back(std::complex<double>(y_interp_ptr[5], y_interp_ptr[6]));
            this->complex_bulk_array_vec.push_back(std::complex<double>(y_interp_ptr[7], y_interp_ptr[8]));

            // Static viscosities (real) carried as EOS extra outputs 9 and 10.
            this->shear_viscosity_array_vec.push_back(y_interp_ptr[9]);
            this->bulk_viscosity_array_vec.push_back(y_interp_ptr[10]);

            // Record central pressure
            if (current_layer_index == 0 && radius_i == 0)
            {
                this->central_pressure = y_interp_ptr[1];
            }
        }

        // At the end we should be at the planet's surface. Use the last interpolated values to set surface constants
        this->surface_gravity  = y_interp_ptr[0];
        this->surface_pressure = y_interp_ptr[1];
        this->mass             = y_interp_ptr[2];
        this->moi              = y_interp_ptr[3];

        // Finished
        this->other_vecs_set = true;
    }

    /// Populate all structure arrays directly from pre-solved world EOS data.
    ///
    /// Bypasses CySolverResult — use when a LayeredWorld has already run solve_eos
    /// and the result needs to be handed into a c_RadialSolutionStorage without
    /// re-integrating the ODE.  All arrays are assumed to be in SI (MKS) units.
    void inject_from_world_eos(
        const double* radius_ptr,
        const double* gravity_ptr,
        const double* pressure_ptr,
        const double* mass_ptr,
        const double* moi_ptr,
        const double* density_ptr,
        const std::complex<double>* complex_shear_ptr,
        const std::complex<double>* complex_bulk_ptr,
        size_t n)
    {
        if (n == 0)
        {
            throw std::invalid_argument("inject_from_world_eos: array length n must be > 0.");
        }

        this->radius_array_size = n;

        this->radius_array_vec.assign(radius_ptr,        radius_ptr        + n);
        this->gravity_array_vec.assign(gravity_ptr,      gravity_ptr       + n);
        this->pressure_array_vec.assign(pressure_ptr,    pressure_ptr      + n);
        this->mass_array_vec.assign(mass_ptr,            mass_ptr          + n);
        this->moi_array_vec.assign(moi_ptr,              moi_ptr           + n);
        this->density_array_vec.assign(density_ptr,      density_ptr       + n);
        this->complex_shear_array_vec.assign(complex_shear_ptr, complex_shear_ptr + n);
        this->complex_bulk_array_vec.assign(complex_bulk_ptr,   complex_bulk_ptr  + n);

        this->radius           = radius_ptr[n - 1];
        this->surface_gravity  = gravity_ptr[n - 1];
        this->surface_pressure = pressure_ptr[n - 1];
        this->mass             = mass_ptr[n - 1];
        this->moi              = moi_ptr[n - 1];
        this->central_pressure = pressure_ptr[0];

        this->nondim_status          = 0;
        this->solution_nondim_status = 0;
        this->success                = true;
        this->error_code             = 0;
        this->radius_array_set       = true;
        this->other_vecs_set         = true;
        this->p_use_array_interp     = true;
    }


    /// Handle dimensionalization/redimensionalization of solution data.
    void dimensionalize_data(
        c_NonDimensionalScales* nondim_scales,
        bool redimensionalize)
    {
        // Save scalers
        this->redim_length_scale  = nondim_scales->length_conversion;
        this->redim_gravity_scale = nondim_scales->length_conversion / nondim_scales->second2_conversion;
        this->redim_mass_scale    = nondim_scales->mass_conversion;
        this->redim_density_scale = nondim_scales->density_conversion;
        this->redim_moi_scale     = nondim_scales->mass_conversion * nondim_scales->length_conversion * nondim_scales->length_conversion;
        this->redim_pascal_scale  = nondim_scales->pascal_conversion;

        // Figure out how to set flags used during eos solution calls
        if (this->solution_nondim_status == 0)
        {
            // EOS was not explicitly non-dim'd or re-dim'd when solution was found.
            if (redimensionalize)
            {
                // User is now "redimensionalizing." Thus we assume that the solution was non-dim'd
                this->nondim_status = 1;
            }
            else
            {
                // User is now "dimensionalizing." Thus we assume that the solution was dim'd
                this->nondim_status = -1;
            }
        }
        else
        {
            // TODO: Deal with this case! For now push the problem to the user when they try to call...
            throw std::runtime_error("Unsupported dimensionalization encountered.");
        }

        // Update other constants
        if (redimensionalize)
        {
            this->pressure_error *= this->redim_pascal_scale;
        }
        else
        {
            this->pressure_error /= this->redim_pascal_scale;
        }

        // Update layer data
        for (size_t layer_i = 0; layer_i < this->num_layers; layer_i++)
        {
            if (redimensionalize)
            {
                this->upper_radius_bylayer_vec[layer_i] *= this->redim_length_scale;
            }
            else
            {
                this->upper_radius_bylayer_vec[layer_i] /= this->redim_length_scale;
            }
        }

        if (this->other_vecs_set)
        {
            // Work through arrays
            for (size_t slice_i = 0; slice_i < this->radius_array_size; slice_i++)
            {
                if (redimensionalize)
                {
                    this->radius_array_vec[slice_i]        *= this->redim_length_scale;
                    this->gravity_array_vec[slice_i]       *= this->redim_gravity_scale;
                    this->pressure_array_vec[slice_i]      *= this->redim_pascal_scale;
                    this->mass_array_vec[slice_i]          *= this->redim_mass_scale;
                    this->moi_array_vec[slice_i]           *= this->redim_moi_scale;
                    this->density_array_vec[slice_i]       *= this->redim_density_scale;
                    this->complex_shear_array_vec[slice_i] *= this->redim_pascal_scale;
                    this->complex_bulk_array_vec[slice_i]  *= this->redim_pascal_scale;
                }
                else
                {
                    this->radius_array_vec[slice_i]        /= this->redim_length_scale;
                    this->gravity_array_vec[slice_i]       /= this->redim_gravity_scale;
                    this->pressure_array_vec[slice_i]      /= this->redim_pascal_scale;
                    this->mass_array_vec[slice_i]          /= this->redim_mass_scale;
                    this->moi_array_vec[slice_i]           /= this->redim_moi_scale;
                    this->density_array_vec[slice_i]       /= this->redim_density_scale;
                    this->complex_shear_array_vec[slice_i] /= this->redim_pascal_scale;
                    this->complex_bulk_array_vec[slice_i]  /= this->redim_pascal_scale;
                }
            }

            // Update global constants
            this->radius           = this->radius_array_vec[this->radius_array_size - 1];
            this->surface_gravity  = this->gravity_array_vec[this->radius_array_size - 1];
            this->surface_pressure = this->pressure_array_vec[this->radius_array_size - 1];
            this->mass             = this->mass_array_vec[this->radius_array_size - 1];
            this->moi              = this->moi_array_vec[this->radius_array_size - 1];
            this->central_pressure = this->pressure_array_vec[0];
        }
    }
};
