#pragma once

#include <algorithm>
#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <memory>
#include <string>

#include "cysolution.hpp"  // CyRK: CySolverResult

#include "../../Utilities_x/dimensions/nondimensional_.hpp"
#include "constants_.hpp"

#include "ode_.hpp"
#include "../../utilities/arrays/interp_.hpp"
#include "../../Utilities_x/math_x/numerics_.hpp"
#include "../../Utilities_x/arrays/layer_partition_.hpp"  // c_partition_radius_by_layer


/// Dense-output read of a finished CyRK integration that never hands back unwritten memory.
///
/// CyRK's `call` returns an error without writing when the query lies outside the integrated domain, which a
/// radius computed by arithmetic can miss by a rounding step at a layer end. A query within the layer
/// continuity tolerance of an end is evaluated at that end; anything further out, a non-finite query, or a
/// failed call leaves every output NaN.
///
/// Parameters
/// ----------
/// result_ptr : Integration with dense output; null gives NaN.
/// time_value : Query in the integration's own independent-variable units.
/// y_out_ptr : Receives `num_y` values.
/// num_y : Number of values CyRK writes for this integration.
///
/// Returns
/// -------
/// True when `y_out_ptr` holds a valid evaluation.
inline bool c_call_dense_checked(
        CySolverResult* result_ptr,
        double time_value,
        double* y_out_ptr,
        const size_t num_y) noexcept
{
    for (size_t y_i = 0; y_i < num_y; ++y_i) { y_out_ptr[y_i] = TidalPyConstants::d_NAN; }
    if (!result_ptr || !result_ptr->config_uptr || !std::isfinite(time_value)) { return false; }

    const double domain_low  = std::min(result_ptr->config_uptr->t_start, result_ptr->config_uptr->t_end);
    const double domain_high = std::max(result_ptr->config_uptr->t_start, result_ptr->config_uptr->t_end);
    const double end_rtol    = tidalpy_config_ptr ? tidalpy_config_ptr->d_LAYER_CONTINUITY_RTOL : 0.0;
    const double slack       = end_rtol * std::max(std::abs(domain_low), std::abs(domain_high));
    if (time_value < domain_low)
    {
        if (time_value < domain_low - slack) { return false; }
        time_value = domain_low;
    }
    else if (time_value > domain_high)
    {
        if (time_value > domain_high + slack) { return false; }
        time_value = domain_high;
    }

    if (result_ptr->call(time_value, y_out_ptr) != CyrkErrorCodes::NO_ERROR)
    {
        for (size_t y_i = 0; y_i < num_y; ++y_i) { y_out_ptr[y_i] = TidalPyConstants::d_NAN; }
        return false;
    }
    return true;
}


/// Equation-of-state integration results for a layered planet.
///
/// Holds the CyRK results of each layer and interpolates the planet's gravity, pressure, mass, moment of
/// inertia, density, and unrelaxed moduli at any radius. The viscoelastic response at a forcing frequency
/// belongs to a rheology rather than to the equation of state and is reached through call_material.
class c_EOSSolution
{

// Attributes
protected:

public:
    int iterations             = -1;
    int error_code             = -100;
    int nondim_status          = 0;
    int solution_nondim_status = 0;
    bool success               = false;
    bool max_iters_hit         = false;
    bool radius_array_set      = false;
    bool other_vecs_set        = false;

    // Between this solution's units and the provider's SI: provider_radius = this_radius * length scale,
    // this_value = provider_value / scale. All one in a solve that ran dimensional.
    double p_structure_length_scale  = 1.0;
    double p_structure_gravity_scale = 1.0;
    double p_structure_pascal_scale  = 1.0;
    double p_structure_density_scale = 1.0;

    // Optional state provider, installed by the world Love solve. One call at an SI radius for one layer
    // fills the whole evaluation layout plus the complex moduli at the solve's forcing frequency. It stays
    // installed after the solve, which is what lets the solution answer at any radius afterwards.
    using MaterialEval = std::function<void(
        size_t layer_index,
        double radius_si,
        double* state_out,
        std::complex<double>& shear_out,
        std::complex<double>& bulk_out)>;
    MaterialEval p_material_eval;

    // Whatever the retained integrators' diffeq arguments point into (the per-layer material inputs and models,
    // the heat sources), co-owned so that every later dense call reads what the solve used, however long the
    // solution is kept and whatever its owner does next.
    std::shared_ptr<const void> input_keepalive;

    std::string message         = "No Message Set.";
    size_t current_layers_saved = 0;
    size_t num_layers           = 0;
    size_t radius_array_size    = 0;
    size_t num_cysolver_calls    = 0;

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
    // One retained integrator per radial segment, ascending. A layer is one segment unless its temperature
    // profile has a kink, so the mapping below finds a layer's segments.
    std::vector<std::unique_ptr<CySolverResult>> cysolver_results_uptr_vec =
        std::vector<std::unique_ptr<CySolverResult>>();
    std::vector<double> segment_upper_radius_vec   = std::vector<double>();
    std::vector<size_t> first_segment_bylayer_vec  = std::vector<size_t>();
    std::vector<size_t> num_segments_bylayer_vec   = std::vector<size_t>();
    // Used when the solve did not carry temperature.
    std::vector<double> segment_temperature_vec    = std::vector<double>();
    // C_EOS_THERMAL_Y_VALUES when the solve integrated temperature and heat flow.
    size_t num_y_solved = C_EOS_Y_VALUES;

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
    // Static (real) shear/bulk viscosity [Pa s] vs radius.
    std::vector<double> shear_viscosity_array_vec = std::vector<double>();
    std::vector<double> bulk_viscosity_array_vec  = std::vector<double>();
    // Temperature [K] and heat flow [W] vs radius.
    std::vector<double> temperature_array_vec     = std::vector<double>();
    std::vector<double> heat_flow_array_vec       = std::vector<double>();

    // Saved by c_solve_eos. The retained integrators carry only the four structure variables; density,
    // moduli, and viscosities are evaluated on demand from the interpolated state.
    std::vector<PreEvalFunc>    eos_function_bylayer_vec = std::vector<PreEvalFunc>();
    std::vector<c_EOS_ODEInput> eos_input_bylayer_vec    = std::vector<c_EOS_ODEInput>();

    // An interface radius is the last slice of the lower layer and the first slice of the upper one, so a
    // lookup for a layer must stay inside its own slices.
    std::vector<size_t> first_slice_bylayer_vec = std::vector<size_t>();
    std::vector<size_t> num_slices_bylayer_vec  = std::vector<size_t>();

// Methods
protected:

public:

    virtual ~c_EOSSolution()
    {
        this->clear_segments();
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
        this->cysolver_results_uptr_vec.reserve(num_layers_);
        this->upper_radius_bylayer_vec.resize(num_layers_);
        std::memcpy(this->upper_radius_bylayer_vec.data(), upper_radius_bylayer_ptr, this->num_layers * sizeof(double));
        this->set_segments_from_layers();
        this->change_radius_array(radius_array_ptr, radius_array_size_);
    }


    /// Save a CyRK solver result for one segment, in ascending radius.
    void save_cyresult(std::unique_ptr<CySolverResult> new_cysolver_result_uptr)
    {
        this->cysolver_results_uptr_vec.push_back(std::move(new_cysolver_result_uptr));
        this->current_layers_saved++;
    }

    /// Release the retained integrators and their dense output.
    void clear_segments() noexcept
    {
        for (size_t i = 0; i < this->cysolver_results_uptr_vec.size(); i++)
        {
            if (this->cysolver_results_uptr_vec[i])
            {
                this->cysolver_results_uptr_vec[i]->dense_vec.clear();
                this->cysolver_results_uptr_vec[i].reset();
            }
        }
        this->cysolver_results_uptr_vec.clear();
    }

    /// The layout of a solve that does not carry temperature.
    void set_segments_from_layers()
    {
        this->segment_upper_radius_vec  = this->upper_radius_bylayer_vec;
        this->first_segment_bylayer_vec.resize(this->num_layers);
        this->num_segments_bylayer_vec.assign(this->num_layers, 1);
        this->segment_temperature_vec.assign(this->num_layers, TidalPyConstants::d_NAN);
        for (size_t layer_i = 0; layer_i < this->num_layers; ++layer_i)
        {
            this->first_segment_bylayer_vec[layer_i] = layer_i;
        }
    }

    void set_segments(const std::vector<c_EOSSegment>& segment_vec)
    {
        const size_t num_segments = segment_vec.size();
        this->segment_upper_radius_vec.resize(num_segments);
        this->segment_temperature_vec.resize(num_segments);
        this->first_segment_bylayer_vec.assign(this->num_layers, 0);
        this->num_segments_bylayer_vec.assign(this->num_layers, 0);
        for (size_t segment_i = 0; segment_i < num_segments; ++segment_i)
        {
            const c_EOSSegment& segment = segment_vec[segment_i];
            this->segment_upper_radius_vec[segment_i] = segment.upper_radius;
            this->segment_temperature_vec[segment_i]  = segment.start_temperature;
            const size_t layer_i = segment.layer_index;
            if (layer_i >= this->num_layers) { continue; }
            if (this->num_segments_bylayer_vec[layer_i] == 0)
            {
                this->first_segment_bylayer_vec[layer_i] = segment_i;
            }
            this->num_segments_bylayer_vec[layer_i]++;
        }
    }

    /// The last segment at or above the radius; the layer's first as a fallback.
    size_t segment_index(const size_t layer_index, const double radius_val) const noexcept
    {
        if (layer_index >= this->num_segments_bylayer_vec.size()) { return layer_index; }
        const size_t first = this->first_segment_bylayer_vec[layer_index];
        const size_t count = this->num_segments_bylayer_vec[layer_index];
        if (count == 0) { return first; }
        for (size_t offset = 0; offset < count - 1; ++offset)
        {
            if (radius_val <= this->segment_upper_radius_vec[first + offset]) { return first + offset; }
        }
        return first + count - 1;
    }


    void save_steps_taken(size_t steps_taken)
    {
        this->steps_taken_vec.push_back(steps_taken);
        this->num_cysolver_calls++;
    }


    /// Kept so the density, moduli, and viscosities can be evaluated at any radius after the solve. The
    /// copies ask the functions for every output.
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
        }
    }


    /// A layer's slices run from its first through the first copy of its upper radius; the second copy starts
    /// the next layer. Both copies share a radius but carry their own layer's density and moduli, so every
    /// array lookup must stay inside one layer's slices.
    void update_slice_partition() noexcept
    {
        if (this->upper_radius_bylayer_vec.size() < this->num_layers)
        {
            this->first_slice_bylayer_vec.assign(this->num_layers, 0);
            this->num_slices_bylayer_vec.assign(this->num_layers, 0);
            return;
        }
        tidalpy::c_partition_radius_by_layer(
            this->radius_array_vec.data(),
            this->radius_array_vec.size(),
            this->upper_radius_bylayer_vec.data(),
            this->num_layers,
            this->first_slice_bylayer_vec,
            this->num_slices_bylayer_vec);
    }




    /// The retained integrators live in the units the solve ran in, so a re-dimensionalized solution
    /// converts an SI query back before evaluating them.
    double convert_radius_si_to_solve(const double radius_si) const noexcept
    {
        return (this->nondim_status == 1) ? radius_si / this->redim_length_scale : radius_si;
    }


protected:
    /// No rescaling: the four structure variables from the dense output, then the density, moduli, and
    /// viscosities from the layer's EOS function at that state. The extra outputs are NaN when no EOS
    /// function was saved.
    void p_evaluate_solver(
        const size_t layer_index,
        const double radius_val,
        double* y_interp_ptr,
        c_EOSMaterialState* material_out = nullptr) const
    {
        // The integrator writes num_y_solved values into a buffer of its own, and the structure variables
        // are copied out: the evaluation layout uses slots 4 and 5 for the density and the shear modulus.
        const size_t segment_i = this->segment_index(layer_index, radius_val);
        double state_arr[C_EOS_THERMAL_Y_VALUES];
        const bool state_found = (segment_i < this->cysolver_results_uptr_vec.size())
            && c_call_dense_checked(
                this->cysolver_results_uptr_vec[segment_i].get(), radius_val, &state_arr[0], C_EOS_THERMAL_Y_VALUES);
        if (!state_found)
        {
            // Outside the layer's solved domain: every output is unknown.
            for (size_t value_i = 0; value_i < C_EOS_DY_VALUES; ++value_i)
            {
                y_interp_ptr[value_i] = TidalPyConstants::d_NAN;
            }
            if (material_out)
            {
                material_out->gravity       = TidalPyConstants::d_NAN;
                material_out->density       = TidalPyConstants::d_NAN;
                material_out->shear_modulus = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
                material_out->bulk_modulus  = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
            }
            return;
        }
        for (size_t value_i = 0; value_i < C_EOS_Y_VALUES; ++value_i)
        {
            y_interp_ptr[value_i] = state_arr[value_i];
        }
        if (material_out)
        {
            material_out->gravity = y_interp_ptr[0];
        }
        if (this->num_y_solved >= C_EOS_THERMAL_Y_VALUES)
        {
            y_interp_ptr[C_EOS_TEMPERATURE_INDEX] = state_arr[4];
            y_interp_ptr[C_EOS_HEAT_FLOW_INDEX]   = state_arr[5];
        }
        else
        {
            // A solve without temperature reports each segment's uniform value and no heat flow.
            y_interp_ptr[C_EOS_TEMPERATURE_INDEX] = (segment_i < this->segment_temperature_vec.size())
                ? this->segment_temperature_vec[segment_i] : TidalPyConstants::d_NAN;
            y_interp_ptr[C_EOS_HEAT_FLOW_INDEX]   = 0.0;
        }
        if (layer_index < this->eos_function_bylayer_vec.size()
            && this->eos_function_bylayer_vec[layer_index] != nullptr)
        {
            c_EOSOutput eos_output;
            // The EOS functions take non-const pointers but leave this solution unchanged.
            char* input_ptr = const_cast<char*>(
                reinterpret_cast<const char*>(&this->eos_input_bylayer_vec[layer_index]));
            // The EOS function reads the state layout, where a thermal solve keeps its temperature at index
            // 4; in the evaluation layout that slot is the density this call is about to fill.
            this->eos_function_bylayer_vec[layer_index](
                reinterpret_cast<char*>(&eos_output), radius_val, &state_arr[0], input_ptr);
            y_interp_ptr[C_EOS_DENSITY_INDEX]         = eos_output.density;
            y_interp_ptr[C_EOS_SHEAR_MODULUS_INDEX]   = eos_output.shear_modulus.real();
            y_interp_ptr[C_EOS_BULK_MODULUS_INDEX]    = eos_output.bulk_modulus.real();
            y_interp_ptr[C_EOS_SHEAR_VISCOSITY_INDEX] = eos_output.shear_viscosity;
            y_interp_ptr[C_EOS_BULK_VISCOSITY_INDEX]  = eos_output.bulk_viscosity;
            y_interp_ptr[C_EOS_MELT_FRACTION_INDEX]   = eos_output.melt_fraction;
            if (material_out)
            {
                material_out->density       = eos_output.density;
                material_out->shear_modulus = eos_output.shear_modulus;
                material_out->bulk_modulus  = eos_output.bulk_modulus;
            }
        }
        else
        {
            // The temperature and heat flow set above stand; only the material outputs are unknown.
            for (size_t value_i = C_EOS_Y_VALUES; value_i < C_EOS_TEMPERATURE_INDEX; ++value_i)
            {
                y_interp_ptr[value_i] = TidalPyConstants::d_NAN;
            }
            y_interp_ptr[C_EOS_MELT_FRACTION_INDEX] = TidalPyConstants::d_NAN;
        }
    }


    /// Applies to the structure variables, density, and the two static moduli. The viscosities and
    /// everything past them are SI in every state and are left alone.
    void p_rescale_outputs(double* y_interp_ptr, const size_t count) const noexcept
    {
        if (this->nondim_status == 0)
        {
            return;
        }
        const size_t scaled_values = C_EOS_SHEAR_VISCOSITY_INDEX;  // slots 0 through the last modulus
        double scales[C_EOS_SHEAR_VISCOSITY_INDEX] = {
            this->redim_gravity_scale, this->redim_pascal_scale, this->redim_mass_scale, this->redim_moi_scale,
            this->redim_density_scale, this->redim_pascal_scale, this->redim_pascal_scale};
        const size_t limit = (count < scaled_values) ? count : scaled_values;
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




    /// Every EOS output at one radius in solve units, in the evaluation layout of eos_layout_.hpp. A
    /// re-dimensionalized solution returns SI for all but the viscosities, which are SI in every state.
    /// Frequency independent throughout; a viscoelastic response comes from call_material.
    void call_nondim(
        const size_t layer_index,
        const double radius_val,
        double* y_interp_ptr) const
    {
        // Provider mode: this solution stores no grid at all, and the provider answers in SI at the exact
        // radius asked for. Only frequency-independent values belong in this layout; the provider's
        // complex moduli reach the solver through call_material.
        if (this->p_material_eval)
        {
            std::complex<double> shear_unused;
            std::complex<double> bulk_unused;
            this->p_material_eval(
                layer_index, radius_val * this->p_structure_length_scale, y_interp_ptr, shear_unused, bulk_unused);
            return;
        }

        if (layer_index >= this->current_layers_saved) [[unlikely]]
        {
            throw std::out_of_range("Layer index out of range.");
        }
        this->p_evaluate_solver(layer_index, radius_val, y_interp_ptr);
        this->p_rescale_outputs(y_interp_ptr, C_EOS_DY_VALUES);
    }


    /// The radial solver's read at an integration radius: gravity, density, and the complex moduli.
    void call_material(
        const size_t layer_index,
        const double radius_val,
        c_EOSMaterialState& out) const
    {
        out = c_EOSMaterialState();

        // One provider call reads the world's solved EOS and applies the rheology at the exact radius asked
        // for; its SI answers are scaled into this solution's units.
        if (this->p_material_eval)
        {
            double state[C_EOS_DY_VALUES];
            this->p_material_eval(
                layer_index, radius_val * this->p_structure_length_scale, state, out.shear_modulus, out.bulk_modulus);
            out.gravity        = state[0] / this->p_structure_gravity_scale;
            out.density        = state[C_EOS_DENSITY_INDEX] / this->p_structure_density_scale;
            out.shear_modulus /= this->p_structure_pascal_scale;
            out.bulk_modulus  /= this->p_structure_pascal_scale;
            return;
        }

        // Otherwise: this solution's own retained integrators plus the layer's EOS model.
        if (layer_index >= this->current_layers_saved) [[unlikely]]
        {
            throw std::out_of_range("Layer index out of range.");
        }
        double scratch[C_EOS_DY_VALUES];
        this->p_evaluate_solver(layer_index, radius_val, scratch, &out);
        // The same rescale the evaluation layout applies, on the four values this carries.
        if (this->nondim_status == 1)
        {
            out.gravity       *= this->redim_gravity_scale;
            out.density       *= this->redim_density_scale;
            out.shear_modulus *= this->redim_pascal_scale;
            out.bulk_modulus  *= this->redim_pascal_scale;
        }
        else if (this->nondim_status != 0)
        {
            out.gravity       /= this->redim_gravity_scale;
            out.density       /= this->redim_density_scale;
            out.shear_modulus /= this->redim_pascal_scale;
            out.bulk_modulus  /= this->redim_pascal_scale;
        }
    }


    /// `call_material` for an SI radius [m].
    void call_material_si(const size_t layer_index, const double radius_si, c_EOSMaterialState& out) const
    {
        this->call_material(layer_index, this->convert_radius_si_to_solve(radius_si), out);
    }


    /// The four structure variables alone, without evaluating the layer's EOS function.
    void call_y(
        const size_t layer_index,
        const double radius_val,
        double* y_interp_ptr) const
    {
        if (layer_index >= this->current_layers_saved) [[unlikely]]
        {
            throw std::out_of_range("Layer index out of range.");
        }
        double state_arr[C_EOS_THERMAL_Y_VALUES];
        const size_t segment_i = this->segment_index(layer_index, radius_val);
        if (segment_i < this->cysolver_results_uptr_vec.size())
        {
            c_call_dense_checked(
                this->cysolver_results_uptr_vec[segment_i].get(), radius_val, &state_arr[0], C_EOS_THERMAL_Y_VALUES);
        }
        else
        {
            for (size_t value_i = 0; value_i < C_EOS_THERMAL_Y_VALUES; ++value_i)
            {
                state_arr[value_i] = TidalPyConstants::d_NAN;
            }
        }
        for (size_t value_i = 0; value_i < C_EOS_Y_VALUES; ++value_i)
        {
            y_interp_ptr[value_i] = state_arr[value_i];
        }
        this->p_rescale_outputs(y_interp_ptr, C_EOS_Y_VALUES);
    }


    /// `call_nondim` for an SI radius [m].
    void call_si(const size_t layer_index, const double radius_si, double* y_interp_ptr) const
    {
        this->call_nondim(layer_index, this->convert_radius_si_to_solve(radius_si), y_interp_ptr);
    }


    /// `call_y` for an SI radius [m].
    void call_y_si(const size_t layer_index, const double radius_si, double* y_interp_ptr) const
    {
        this->call_y(layer_index, this->convert_radius_si_to_solve(radius_si), y_interp_ptr);
    }


    /// The innermost layer whose upper radius reaches it. An interface radius belongs to the lower layer,
    /// matching the convention the rest of this solution uses.
    size_t layer_at_radius_si(const double radius_si) const noexcept
    {
        const double radius_solve = this->convert_radius_si_to_solve(radius_si);
        for (size_t layer_i = 0; layer_i < this->num_layers; ++layer_i)
        {
            if (radius_solve <= this->upper_radius_bylayer_vec[layer_i])
            {
                return layer_i;
            }
        }
        return (this->num_layers > 0) ? (this->num_layers - 1) : 0;
    }




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
            this->temperature_array_vec.clear();
            this->heat_flow_array_vec.clear();
            this->clear_segments();
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

    /// Interpolate the whole planet onto the stored radius array, layer by layer.
    void interpolate_full_planet()
    {
        this->solution_nondim_status = this->nondim_status;

        if (this->current_layers_saved == 0)
        {
            throw std::runtime_error("No layers have been saved. Can not perform interpolation.");
        }

        if (this->radius_array_size == 0)
        {
            throw std::runtime_error("No radius slices have been set. Can not perform interpolation.");
        }

        size_t current_layer_index        = 0;
        double current_layer_upper_radius = this->upper_radius_bylayer_vec[0];

        // Zero initialized because the surface values are read out after the loop: a loop that breaks on
        // its first pass would otherwise leave the planet's mass and moi holding stack garbage.
        double y_interp_arr[C_EOS_DY_VALUES] = {0.0};
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

            // The evaluation layout carries only the unrelaxed moduli, so the complex ones come separately.
            c_EOSMaterialState material_state;
            this->p_evaluate_solver(current_layer_index, radius_val, y_interp_ptr, &material_state);

            this->gravity_array_vec.push_back(y_interp_ptr[0]);
            this->pressure_array_vec.push_back(y_interp_ptr[1]);
            this->mass_array_vec.push_back(y_interp_ptr[2]);
            this->moi_array_vec.push_back(y_interp_ptr[3]);
            this->density_array_vec.push_back(y_interp_ptr[C_EOS_DENSITY_INDEX]);

            this->complex_shear_array_vec.push_back(material_state.shear_modulus);
            this->complex_bulk_array_vec.push_back(material_state.bulk_modulus);

            this->shear_viscosity_array_vec.push_back(y_interp_ptr[C_EOS_SHEAR_VISCOSITY_INDEX]);
            this->bulk_viscosity_array_vec.push_back(y_interp_ptr[C_EOS_BULK_VISCOSITY_INDEX]);
            this->temperature_array_vec.push_back(y_interp_ptr[C_EOS_TEMPERATURE_INDEX]);
            this->heat_flow_array_vec.push_back(y_interp_ptr[C_EOS_HEAT_FLOW_INDEX]);

            if (current_layer_index == 0 && radius_i == 0)
            {
                this->central_pressure = y_interp_ptr[1];
            }
        }

        this->surface_gravity  = y_interp_ptr[0];
        this->surface_pressure = y_interp_ptr[1];
        this->mass             = y_interp_ptr[2];
        this->moi              = y_interp_ptr[3];

        this->other_vecs_set = true;
    }



    /// The viscosity arrays are SI in every state and are left alone.
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

        // A solution neither non-dimensionalized nor re-dimensionalized at solve time is taken to be in the
        // state the requested direction implies.
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

        // A solution may carry its scalars without the sampled arrays (a radial solver's stand-in for a world EOS),
        // so the loop is bounded by what every array actually holds.
        const size_t num_slices = std::min({
            this->radius_array_size, this->radius_array_vec.size(), this->gravity_array_vec.size(),
            this->pressure_array_vec.size(), this->mass_array_vec.size(), this->moi_array_vec.size(),
            this->density_array_vec.size(), this->complex_shear_array_vec.size(), this->complex_bulk_array_vec.size()});
        if (this->other_vecs_set && (num_slices > 0))
        {
            for (size_t slice_i = 0; slice_i < num_slices; slice_i++)
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

            this->radius           = this->radius_array_vec[num_slices - 1];
            this->surface_gravity  = this->gravity_array_vec[num_slices - 1];
            this->surface_pressure = this->pressure_array_vec[num_slices - 1];
            this->mass             = this->mass_array_vec[num_slices - 1];
            this->moi              = this->moi_array_vec[num_slices - 1];
            this->central_pressure = this->pressure_array_vec[0];
        }
        else
        {
            // No arrays to read the scalars back from: convert them in place.
            const auto convert = [redimensionalize](double& value, double scale) {
                value = redimensionalize ? value * scale : value / scale;
            };
            convert(this->radius,           this->redim_length_scale);
            convert(this->surface_gravity,  this->redim_gravity_scale);
            convert(this->surface_pressure, this->redim_pascal_scale);
            convert(this->central_pressure, this->redim_pascal_scale);
            convert(this->mass,             this->redim_mass_scale);
            convert(this->moi,              this->redim_moi_scale);
        }
    }
};
