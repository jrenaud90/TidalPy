#pragma once

#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <memory>
#include <string>

#include "cysolution.hpp"  // CyRK: CySolverResult

#include "../../Utilities_x/dimensions/nondimensional_.hpp" // c_NonDimensionalScales
#include "constants_.hpp"

#include "ode_.hpp" // C_EOS_Y_VALUES, C_EOS_EXTRA_VALUES, C_EOS_DY_VALUES
#include "../../utilities/arrays/interp_.hpp"  // c_binary_search_with_guess, c_interp, c_interp_complex
#include "../../Utilities_x/math_x/numerics_.hpp"  // c_isclose


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

    // Optional non-owning dense structure source (the world Love solve). With p_use_array_interp, gravity, pressure,
    // mass, and moi are read from it at dense accuracy while density and the complex moduli stay array-interpolated.
    // The source is SI: source_radius = this_radius * p_structure_length_scale, this_value = source_value / scale.
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

    // Radial density derivative per slice (array units per unit radius). When present the density lookup is cubic
    // Hermite (slope continuous across slices, so the radial integrator sees no kinks). Filled by
    // inject_from_world_eos; empty otherwise.
    std::vector<double> density_slope_array_vec = std::vector<double>();

    // Per-layer EOS functions and arguments saved by c_solve_eos. The retained integrators carry only the four
    // structure variables; density, moduli, and viscosities are evaluated on demand from the interpolated state.
    std::vector<PreEvalFunc>    eos_function_bylayer_vec = std::vector<PreEvalFunc>();
    std::vector<c_EOS_ODEInput> eos_input_bylayer_vec    = std::vector<c_EOS_ODEInput>();

    // Per-layer partition of the radius array (see update_slice_partition). An interface radius is the last slice
    // of the lower layer and the first slice of the upper one, so a lookup for a layer must stay inside its slices.
    std::vector<size_t> first_slice_bylayer_vec = std::vector<size_t>();
    std::vector<size_t> num_slices_bylayer_vec  = std::vector<size_t>();

// Methods
protected:

public:

    virtual ~c_EOSSolution()
    {
        for (size_t i = 0; i < this->cysolver_results_uptr_bylayer_vec.size(); i++)
        {
            this->cysolver_results_uptr_bylayer_vec[i]->dense_vec.clear();
            this->cysolver_results_uptr_bylayer_vec[i].reset();
        }
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
        std::memcpy(this->upper_radius_bylayer_vec.data(), upper_radius_bylayer_ptr, this->num_layers * sizeof(double));
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


    /// Keep the per-layer EOS evaluation functions and their arguments so the density, moduli, and viscosities can
    /// be evaluated at any radius after the solve. The copies ask the functions for every output.
    void save_eos_functions(
        const std::vector<PreEvalFunc>& eos_function_bylayer,
        const std::vector<c_EOS_ODEInput>& eos_input_bylayer)
    {
        this->eos_function_bylayer_vec = eos_function_bylayer;
        this->eos_input_bylayer_vec    = eos_input_bylayer;
        for (c_EOS_ODEInput& input : this->eos_input_bylayer_vec)
        {
            input.update_bulk  = true;
            input.update_shear = true;
            input.final_solve  = true;
        }
    }


    /// Partition the radius array by layer: a layer's slices run from its first slice through the first copy of
    /// its upper radius, the second copy of an interface radius starting the next layer. Both copies share a radius
    /// but carry their own layer's density and moduli, so every array lookup must stay inside one layer's slices.
    void update_slice_partition() noexcept
    {
        const size_t n = this->radius_array_vec.size();
        this->first_slice_bylayer_vec.assign(this->num_layers, 0);
        this->num_slices_bylayer_vec.assign(this->num_layers, 0);
        if (n == 0 || this->upper_radius_bylayer_vec.size() < this->num_layers)
        {
            return;
        }
        size_t next_first = 0;
        for (size_t layer_i = 0; layer_i < this->num_layers; ++layer_i)
        {
            const double layer_upper = this->upper_radius_bylayer_vec[layer_i];
            size_t count          = 0;
            size_t interface_hits = 0;
            for (size_t slice_i = next_first; slice_i < n; ++slice_i)
            {
                const double radius_check = this->radius_array_vec[slice_i];
                if (c_isclose(radius_check, layer_upper, 1.0e-9, 0.0))
                {
                    if (++interface_hits > 1)
                    {
                        break;
                    }
                }
                else if (radius_check > layer_upper)
                {
                    break;
                }
                ++count;
            }
            this->first_slice_bylayer_vec[layer_i] = next_first;
            this->num_slices_bylayer_vec[layer_i]  = count;
            next_first += count;
        }
    }


    /// Linear interpolation of gravity and density inside one layer's slices, in array units. Used where the dense
    /// output is unavailable (the dynamic-liquid y3 reconstruction).
    void interp_structure_in_layer(
        const size_t layer_index,
        const double radius_val,
        double* gravity_ptr,
        double* density_ptr) const noexcept
    {
        size_t first = 0;
        size_t n     = this->radius_array_vec.size();
        if (layer_index < this->num_slices_bylayer_vec.size() && this->num_slices_bylayer_vec[layer_index] > 1)
        {
            first = this->first_slice_bylayer_vec[layer_index];
            n     = this->num_slices_bylayer_vec[layer_index];
        }
        if (n == 0)
        {
            *gravity_ptr = TidalPyConstants::d_NAN;
            *density_ptr = TidalPyConstants::d_NAN;
            return;
        }
        double* radius_data_ptr = const_cast<double*>(this->radius_array_vec.data()) + first;
        double  radius_query    = radius_val;
        size_t  j               = this->p_seed_index(radius_val, radius_data_ptr, n);
        int     b_code          = 0;
        j = c_binary_search_with_guess(radius_val, radius_data_ptr, n, j, &b_code);
        c_interp(&radius_query, radius_data_ptr, const_cast<double*>(this->gravity_array_vec.data()) + first, n, &j,
                 gravity_ptr);
        c_interp(&radius_query, radius_data_ptr, const_cast<double*>(this->density_array_vec.data()) + first, n, &j,
                 density_ptr);
    }


    /// Radius in solve units for an SI radius: the retained integrators live in the units the solve ran in, so a
    /// re-dimensionalized solution converts an SI query back before evaluating them.
    double convert_radius_si_to_solve(const double radius_si) const noexcept
    {
        return (this->nondim_status == 1) ? radius_si / this->redim_length_scale : radius_si;
    }


protected:
    /// Index guess for a binary search from the query's fractional position in the array.
    static size_t p_seed_index(const double radius_val, const double* radius_data_ptr, const size_t n) noexcept
    {
        const double r_left  = radius_data_ptr[0];
        const double r_right = radius_data_ptr[n - 1];
        size_t j = 0;
        if (r_right > r_left)
        {
            const double frac    = (radius_val - r_left) / (r_right - r_left);
            const double clamped = frac < 0.0 ? 0.0 : (frac > 1.0 ? 1.0 : frac);
            j = static_cast<size_t>(static_cast<double>(n) * clamped);
            if (j >= n)
            {
                j = n - 1;
            }
        }
        return j;
    }


    /// Evaluate one layer's retained integrator at a radius in solve units, with no rescaling: the four structure
    /// variables from the dense output, then the density, moduli, and viscosities from the layer's EOS function at
    /// that state. Writes C_EOS_DY_VALUES doubles; the extra outputs are NaN when no EOS function was saved.
    void p_evaluate_solver(const size_t layer_index, const double radius_val, double* y_interp_ptr) const
    {
        this->cysolver_results_uptr_bylayer_vec[layer_index]->call(radius_val, y_interp_ptr);
        if (layer_index < this->eos_function_bylayer_vec.size()
            && this->eos_function_bylayer_vec[layer_index] != nullptr)
        {
            c_EOSOutput eos_output;
            // The EOS functions take their arguments through non-const pointers but leave this solution unchanged.
            char* input_ptr = const_cast<char*>(
                reinterpret_cast<const char*>(&this->eos_input_bylayer_vec[layer_index]));
            this->eos_function_bylayer_vec[layer_index](
                reinterpret_cast<char*>(&eos_output), radius_val, y_interp_ptr, input_ptr);
            y_interp_ptr[4]  = eos_output.density;
            y_interp_ptr[5]  = eos_output.shear_modulus.real();
            y_interp_ptr[6]  = eos_output.shear_modulus.imag();
            y_interp_ptr[7]  = eos_output.bulk_modulus.real();
            y_interp_ptr[8]  = eos_output.bulk_modulus.imag();
            y_interp_ptr[9]  = eos_output.shear_viscosity;
            y_interp_ptr[10] = eos_output.bulk_viscosity;
        }
        else
        {
            for (size_t value_i = C_EOS_Y_VALUES; value_i < C_EOS_DY_VALUES; ++value_i)
            {
                y_interp_ptr[value_i] = TidalPyConstants::d_NAN;
            }
        }
    }


    /// Apply the solution's dimensional state to the first `count` outputs of an evaluation (structure variables,
    /// density, and moduli). The viscosities (indices 9 and 10) are always SI and are left alone.
    void p_rescale_outputs(double* y_interp_ptr, const size_t count) const noexcept
    {
        if (this->nondim_status == 0)
        {
            return;
        }
        double scales[9] = {
            this->redim_gravity_scale, this->redim_pascal_scale, this->redim_mass_scale, this->redim_moi_scale,
            this->redim_density_scale, this->redim_pascal_scale, this->redim_pascal_scale, this->redim_pascal_scale,
            this->redim_pascal_scale};
        const size_t limit = (count < 9) ? count : 9;
        if (this->nondim_status == 1)
        {
            for (size_t value_i = 0; value_i < limit; ++value_i)
            {
                y_interp_ptr[value_i] *= scales[value_i];
            }
        }
        else
        {
            for (size_t value_i = 0; value_i < limit; ++value_i)
            {
                y_interp_ptr[value_i] /= scales[value_i];
            }
        }
    }

public:


    /// Interpolate every EOS output from the stored arrays, searching only the calling layer's slices (see
    /// update_slice_partition). Used when p_use_array_interp is true (set by inject_from_world_eos). The arrays
    /// stay in the units they were injected in; no scaling is applied.
    void _call_interp_arrays(const size_t layer_index, const double radius_val, double* y_interp_ptr) const noexcept
    {
        size_t first = 0;
        size_t n     = this->radius_array_size;
        if (layer_index < this->num_slices_bylayer_vec.size() && this->num_slices_bylayer_vec[layer_index] > 1)
        {
            first = this->first_slice_bylayer_vec[layer_index];
            n     = this->num_slices_bylayer_vec[layer_index];
        }
        if (n == 0)
        {
            // NaN rather than reading uninitialized memory.
            for (size_t value_i = 0; value_i < C_EOS_DY_VALUES; ++value_i)
            {
                y_interp_ptr[value_i] = TidalPyConstants::d_NAN;
            }
            return;
        }
        // c_interp/c_binary_search_with_guess take non-const double* but only read the data.
        double* radius_data_ptr = const_cast<double*>(this->radius_array_vec.data()) + first;
        double  radius_query = radius_val;   // mutable copy for c_interp's desired_x_ptr arg

        size_t j      = this->p_seed_index(radius_val, radius_data_ptr, n);
        int    b_code = 0;
        j = c_binary_search_with_guess(radius_val, radius_data_ptr, n, j, &b_code);

        // Density: cubic Hermite between the bracketing slices when the slopes are stored, else linear below.
        const bool hermite_density = (this->density_slope_array_vec.size() == this->radius_array_size) && (n >= 2);
        if (hermite_density)
        {
            size_t left = (b_code == -1) ? 0 : j;
            if (left >= n - 1)
            {
                left = n - 2;
            }
            const double* density_ptr = this->density_array_vec.data() + first;
            const double* slope_ptr   = this->density_slope_array_vec.data() + first;
            const double  r_left      = radius_data_ptr[left];
            const double  span        = radius_data_ptr[left + 1] - r_left;
            double t = (span > 0.0) ? (radius_val - r_left) / span : 0.0;
            t = (t < 0.0) ? 0.0 : ((t > 1.0) ? 1.0 : t);
            const double t2 = t * t;
            const double t3 = t2 * t;
            y_interp_ptr[4] =
                (2.0 * t3 - 3.0 * t2 + 1.0) * density_ptr[left]
                + (t3 - 2.0 * t2 + t) * span * slope_ptr[left]
                + (-2.0 * t3 + 3.0 * t2) * density_ptr[left + 1]
                + (t3 - t2) * span * slope_ptr[left + 1];
        }

        // Output layout matches the CySolverResult: [0] gravity, [1] pressure, [2] mass, [3] moi, [4] density,
        // [5] shear_real, [6] shear_imag, [7] bulk_real, [8] bulk_imag, [9] shear viscosity, [10] bulk viscosity.
        c_interp(
            &radius_query,
            radius_data_ptr,
            const_cast<double*>(this->gravity_array_vec.data()) + first,
            n,
            &j,
            &y_interp_ptr[0]);
        c_interp(
            &radius_query,
            radius_data_ptr,
            const_cast<double*>(this->pressure_array_vec.data()) + first,
            n,
            &j,
            &y_interp_ptr[1]);
        c_interp(
            &radius_query,
            radius_data_ptr,
            const_cast<double*>(this->mass_array_vec.data()) + first,
            n,
            &j,
            &y_interp_ptr[2]);
        c_interp(
            &radius_query,
            radius_data_ptr,
            const_cast<double*>(this->moi_array_vec.data()) + first,
            n,
            &j,
            &y_interp_ptr[3]);
        if (!hermite_density)
        {
            c_interp(
                &radius_query,
                radius_data_ptr,
                const_cast<double*>(this->density_array_vec.data()) + first,
                n,
                &j,
                &y_interp_ptr[4]);
        }

        double shear_result[2] = {0.0, 0.0};
        c_interp_complex(
            radius_val,
            radius_data_ptr,
            const_cast<double*>(reinterpret_cast<const double*>(this->complex_shear_array_vec.data() + first)),
            n,
            &j,
            shear_result);
        y_interp_ptr[5] = shear_result[0];
        y_interp_ptr[6] = shear_result[1];

        double bulk_result[2] = {0.0, 0.0};
        c_interp_complex(
            radius_val,
            radius_data_ptr,
            const_cast<double*>(reinterpret_cast<const double*>(this->complex_bulk_array_vec.data() + first)),
            n,
            &j,
            bulk_result);
        y_interp_ptr[7] = bulk_result[0];
        y_interp_ptr[8] = bulk_result[1];

        // Viscosities are NaN when not stored.
        if (this->shear_viscosity_array_vec.size() == this->radius_array_size
            && this->bulk_viscosity_array_vec.size() == this->radius_array_size)
        {
            c_interp(
                &radius_query,
                radius_data_ptr,
                const_cast<double*>(this->shear_viscosity_array_vec.data()) + first,
                n,
                &j,
                &y_interp_ptr[9]);
            c_interp(
                &radius_query,
                radius_data_ptr,
                const_cast<double*>(this->bulk_viscosity_array_vec.data()) + first,
                n,
                &j,
                &y_interp_ptr[10]);
        }
        else
        {
            y_interp_ptr[9]  = TidalPyConstants::d_NAN;
            y_interp_ptr[10] = TidalPyConstants::d_NAN;
        }
    }


    /// Evaluate every EOS output at a single radius in solve units for a specific layer, writing C_EOS_DY_VALUES
    /// doubles: gravity, pressure, mass, moment of inertia, density, shear modulus (real, imaginary), bulk modulus
    /// (real, imaginary), shear viscosity, bulk viscosity. A re-dimensionalized solution returns SI values for all
    /// but the viscosities, which are SI in every state.
    void call(
        const size_t layer_index,
        const double radius_val,
        double* y_interp_ptr) const
    {
        if (this->p_use_array_interp) [[unlikely]]
        {
            this->_call_interp_arrays(layer_index, radius_val, y_interp_ptr);

            // The structure variables come from the dense source at its accuracy; density and moduli stay
            // array-interpolated.
            if (this->p_structure_dense_source) [[unlikely]]
            {
                double src_out[C_EOS_Y_VALUES];
                const double src_radius = radius_val * this->p_structure_length_scale;
                this->p_structure_dense_source->call_y_si(layer_index, src_radius, src_out);
                y_interp_ptr[0] = src_out[0] / this->p_structure_gravity_scale;   // gravity
                y_interp_ptr[1] = src_out[1] / this->p_structure_pascal_scale;    // pressure
                y_interp_ptr[2] = src_out[2] / this->p_structure_mass_scale;      // mass
                y_interp_ptr[3] = src_out[3] / this->p_structure_moi_scale;       // moment of inertia
            }
            return;
        }

        if (layer_index >= this->current_layers_saved) [[unlikely]]
        {
            throw std::out_of_range("Layer index out of range.");
        }
        this->p_evaluate_solver(layer_index, radius_val, y_interp_ptr);
        this->p_rescale_outputs(y_interp_ptr, C_EOS_DY_VALUES);
    }


    /// The four structure variables (gravity, pressure, mass, moment of inertia) at a radius in solve units,
    /// without evaluating the layer's EOS function. Writes C_EOS_Y_VALUES doubles.
    void call_y(
        const size_t layer_index,
        const double radius_val,
        double* y_interp_ptr) const
    {
        if (this->p_use_array_interp) [[unlikely]]
        {
            double full_out[C_EOS_DY_VALUES];
            this->call(layer_index, radius_val, full_out);
            for (size_t value_i = 0; value_i < C_EOS_Y_VALUES; ++value_i)
            {
                y_interp_ptr[value_i] = full_out[value_i];
            }
            return;
        }
        if (layer_index >= this->current_layers_saved) [[unlikely]]
        {
            throw std::out_of_range("Layer index out of range.");
        }
        this->cysolver_results_uptr_bylayer_vec[layer_index]->call(radius_val, y_interp_ptr);
        this->p_rescale_outputs(y_interp_ptr, C_EOS_Y_VALUES);
    }


    /// `call` for an SI radius [m]: the radius is converted into solve units first.
    void call_si(const size_t layer_index, const double radius_si, double* y_interp_ptr) const
    {
        this->call(layer_index, this->convert_radius_si_to_solve(radius_si), y_interp_ptr);
    }


    /// `call_y` for an SI radius [m]: the radius is converted into solve units first.
    void call_y_si(const size_t layer_index, const double radius_si, double* y_interp_ptr) const
    {
        this->call_y(layer_index, this->convert_radius_si_to_solve(radius_si), y_interp_ptr);
    }


    /// Prepare storage vectors for a new or changed radius array.
    void change_radius_array(
        double* new_radius_ptr,
        size_t new_radius_size)
    {
        this->radius_array_size = new_radius_size;
        if (this->radius_array_set)
        {
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
            this->other_vecs_set = false;
        }
        this->radius_array_set = true;

        this->gravity_array_vec.reserve(this->radius_array_size);
        this->pressure_array_vec.reserve(this->radius_array_size);
        this->mass_array_vec.reserve(this->radius_array_size);
        this->moi_array_vec.reserve(this->radius_array_size);
        this->density_array_vec.reserve(this->radius_array_size);
        this->complex_shear_array_vec.reserve(this->radius_array_size);
        this->complex_bulk_array_vec.reserve(this->radius_array_size);
        this->shear_viscosity_array_vec.reserve(this->radius_array_size);
        this->bulk_viscosity_array_vec.reserve(this->radius_array_size);

        this->radius_array_vec.resize(this->radius_array_size);
        std::memcpy(this->radius_array_vec.data(), new_radius_ptr, new_radius_size * sizeof(double));

        this->radius = this->radius_array_vec.back();
        this->update_slice_partition();
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

            if (c_isclose(radius_val, current_layer_upper_radius))
            {
                // An interface radius appears twice; the first copy belongs to the lower layer.
                if (interface_check == 1)
                {
                    ready_for_next_layer = true;
                }
                interface_check++;
            }
            else if (radius_val > current_layer_upper_radius)
            {
                ready_for_next_layer = true;
            }

            if (ready_for_next_layer)
            {
                current_layer_index++;
                if (current_layer_index > (this->num_layers - 1))
                {
                    break;
                }
                else
                {
                    interface_check = 0;
                    current_layer_upper_radius = this->upper_radius_bylayer_vec[current_layer_index];
                    ready_for_next_layer       = false;
                }
            }

            this->p_evaluate_solver(current_layer_index, radius_val, y_interp_ptr);

            this->gravity_array_vec.push_back(y_interp_ptr[0]);
            this->pressure_array_vec.push_back(y_interp_ptr[1]);
            this->mass_array_vec.push_back(y_interp_ptr[2]);
            this->moi_array_vec.push_back(y_interp_ptr[3]);
            this->density_array_vec.push_back(y_interp_ptr[4]);

            // CyRK carries the moduli as real and imaginary double pairs.
            this->complex_shear_array_vec.push_back(std::complex<double>(y_interp_ptr[5], y_interp_ptr[6]));
            this->complex_bulk_array_vec.push_back(std::complex<double>(y_interp_ptr[7], y_interp_ptr[8]));

            this->shear_viscosity_array_vec.push_back(y_interp_ptr[9]);
            this->bulk_viscosity_array_vec.push_back(y_interp_ptr[10]);

            if (current_layer_index == 0 && radius_i == 0)
            {
                this->central_pressure = y_interp_ptr[1];
            }
        }

        // The last interpolated values are the surface values.
        this->surface_gravity  = y_interp_ptr[0];
        this->surface_pressure = y_interp_ptr[1];
        this->mass             = y_interp_ptr[2];
        this->moi              = y_interp_ptr[3];

        this->other_vecs_set = true;
    }

    /// Populate the structure arrays from a world's solved EOS without re-integrating, in whatever units the
    /// arrays are given in. A density_slope_ptr (radial density derivative per slice) switches the density lookup
    /// to cubic Hermite; null keeps it linear.
    void inject_from_world_eos(
        const double* radius_ptr,
        const double* gravity_ptr,
        const double* pressure_ptr,
        const double* mass_ptr,
        const double* moi_ptr,
        const double* density_ptr,
        const std::complex<double>* complex_shear_ptr,
        const std::complex<double>* complex_bulk_ptr,
        size_t n,
        const double* density_slope_ptr = nullptr)
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
        if (density_slope_ptr != nullptr)
        {
            this->density_slope_array_vec.assign(density_slope_ptr, density_slope_ptr + n);
        }
        else
        {
            this->density_slope_array_vec.clear();
        }

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
        this->update_slice_partition();
    }


    /// Scale the structure variables, density, and moduli into or out of non-dimensional units. The viscosity
    /// arrays are SI in every state and are left alone.
    void dimensionalize_data(
        c_NonDimensionalScales* nondim_scales,
        bool redimensionalize)
    {
        this->redim_length_scale  = nondim_scales->length_conversion;
        this->redim_gravity_scale = nondim_scales->length_conversion / nondim_scales->second2_conversion;
        this->redim_mass_scale    = nondim_scales->mass_conversion;
        this->redim_density_scale = nondim_scales->density_conversion;
        this->redim_moi_scale     = nondim_scales->mass_conversion * nondim_scales->length_conversion * nondim_scales->length_conversion;
        this->redim_pascal_scale  = nondim_scales->pascal_conversion;

        // A solution that was neither non-dimensionalized nor re-dimensionalized at solve time is assumed to be in
        // the state the requested direction implies.
        if (this->solution_nondim_status == 0)
        {
            if (redimensionalize)
            {
                this->nondim_status = 1;
            }
            else
            {
                this->nondim_status = -1;
            }
        }
        else
        {
            // TODO: Handle repeated same-direction dimensionalization requests instead of raising.
            throw std::runtime_error("Unsupported dimensionalization encountered.");
        }

        if (redimensionalize)
        {
            this->pressure_error *= this->redim_pascal_scale;
        }
        else
        {
            this->pressure_error /= this->redim_pascal_scale;
        }

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

            this->radius           = this->radius_array_vec[this->radius_array_size - 1];
            this->surface_gravity  = this->gravity_array_vec[this->radius_array_size - 1];
            this->surface_pressure = this->pressure_array_vec[this->radius_array_size - 1];
            this->mass             = this->mass_array_vec[this->radius_array_size - 1];
            this->moi              = this->moi_array_vec[this->radius_array_size - 1];
            this->central_pressure = this->pressure_array_vec[0];
        }
    }
};
