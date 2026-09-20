// rs_solution_.hpp: radial solver solution storage.
#pragma once

#include <cstring>
#include <cmath>
#include <array>
#include <complex>
#include <vector>
#include <memory>
#include <string>
#include <functional>

#include "love_.hpp"
#include "rs_constants_.hpp"
#include "../Material_x/eos/eos_solution_.hpp"   // also provides CyRK's CySolverResult (complete type)
#include "../Material_x/eos/methods/interpolate_.hpp"  // c_InterpolateEOSInput (persisted standalone EOS args)
#include "../../constants_.hpp"
#include "../Utilities_x/dimensions/nondimensional_.hpp"
#include "../utilities/arrays/interp_.hpp"        // c_interp / c_binary_search_with_guess (shared array-interp)


// Error Codes:
// -1 : Equation of State storage (c_EOSSolution) could not be initialized.
// -2 : (set by python wrapper) Unknown / Unsupported boundary condition provided.
// -5 : There was a problem with the inputs to radial solver
//
// -1X : Error in shooting method
// -10 : Error in finding starting conditions
// -11 : Numerical integration failed
// -12 : Error using ZGESV solver with boundary condition
//
// -2X : Error in propagation matrix method
// -20 : Unknown core starting conditions
// -21 : Error using ZGESV solver with boundary condition

class c_RadialSolutionStorage
{
public:
    bool success        = false;
    int error_code      = -100;
    int degree_l        = 0;
    std::string message = "No Message Set.";
    size_t num_ytypes   = 0;
    size_t num_slices   = 0;
    size_t num_layers   = 0;
    size_t total_size   = 0;

    // Equation of state solution
    std::unique_ptr<c_EOSSolution> eos_solution_uptr = nullptr;

    // Radial solution results (stores double-pairs for complex values)
    std::vector<double> full_solution_vec = std::vector<double>();

    // Love number attributes (stores double-pairs for complex values)
    std::vector<c_LoveNumbers> complex_love_vec = std::vector<c_LoveNumbers>();

    // Surface y-solution (SI), laid out [ytype * C_MAX_NUM_Y + y_index], cached by find_love once per solve.
    std::vector<std::complex<double>> p_surface_y_si = std::vector<std::complex<double>>();

    // The surface boundary conditions this solve produced blocks for, in order (tidal = 1, free = 0,
    // loading = 2).
    std::vector<int> p_bc_models = std::vector<int>();

    // Diagnostic data
    std::vector<size_t> shooting_method_steps_taken_vec = std::vector<size_t>();

    // Worst-case error amplification of the surface boundary condition solve across ytypes (shooting method
    // only; stays 0 for the matrix method). See c_estimate_surface_amplification in boundaries_.hpp.
    double surface_amplification = 0.0;

    // ================================================================================================================
    // Interpolant-based (shooting) solution
    // ================================================================================================================
    // When true, get_radial_solution evaluates the per-(layer, solution) dense CyRK interpolants at any radius and
    // collapses them with the constants below; the matrix method leaves it false and fills full_solution_vec, which
    // get_radial_solution then interpolates linearly.
    bool p_uses_interpolants = false;

    // The propagation matrix's own radius grid (solve units), laid down inside c_matrix_propagate and kept only
    // because full_solution_vec is interpolated against it. Empty after a shooting solve, which grids nothing.
    std::vector<double> p_matrix_radius_solve = std::vector<double>();

    // Rebuilds the complex moduli this solve used, from a layer's static modulus and viscosity. The world
    // installs it holding shared ownership of the layers' rheologies and the solved frequency, so it keeps
    // working after those layers are gone; the rheology models are pure in their three arguments, which is what
    // makes that safe. Type-erased because the rheology classes live above this header, the same reason
    // c_EOSSolution::MaterialEval is. Empty when the solve was handed its moduli rather than deriving them (the
    // supplied-moduli path), and then the complex getters report NaN while the static ones still answer.
    using ComplexModuliEval = std::function<void(
        size_t layer_index,
        double static_shear, double shear_viscosity,
        double static_bulk,  double bulk_viscosity,
        std::complex<double>& shear_out, std::complex<double>& bulk_out)>;
    ComplexModuliEval p_complex_moduli_eval;
    double p_love_frequency_si = TidalPyConstants::d_NAN;

    /// The complex shear and bulk moduli [Pa] at an SI radius, rebuilt the way the solve built them: the layer's
    /// rheology applied to the static modulus and viscosity the solved EOS reports there, at the solved
    /// frequency. Both NaN when this solution carries no rheology.
    void get_complex_moduli_si(
        const double radius_si,
        std::complex<double>& shear_out,
        std::complex<double>& bulk_out) const
    {
        const std::complex<double> cNAN(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
        shear_out = cNAN;
        bulk_out  = cNAN;

        size_t layer_i = 0;
        double solve_r = 0.0;
        if (!this->p_locate_eos(radius_si, layer_i, solve_r)) { return; }

        if (this->p_complex_moduli_eval)
        {
            double state[C_EOS_DY_VALUES];
            this->eos_solution_uptr->call(layer_i, solve_r, &state[0]);
            this->p_complex_moduli_eval(
                layer_i,
                state[C_EOS_SHEAR_MODULUS_INDEX], state[C_EOS_SHEAR_VISCOSITY_INDEX],
                state[C_EOS_BULK_MODULUS_INDEX],  state[C_EOS_BULK_VISCOSITY_INDEX],
                shear_out, bulk_out);
            return;
        }

        // No rheology to apply, because this solve was provided arrays of moduli. So instead we will use
        // linear interpolation on those input arrays to return the moduli if requested.
        c_EOSMaterialState material_state;
        this->eos_solution_uptr->call_material(layer_i, solve_r, material_state);
        shear_out = material_state.shear_modulus;
        bulk_out  = material_state.bulk_modulus;
    }

    // Dense CyRK results [layer][solution]; owns the force-retained integrators.
    std::vector<std::vector<std::unique_ptr<CySolverResult>>> p_interp_by_layer_sol;

    // Collapse constants [ytype][layer][solution], at most 3 solutions per layer.
    std::vector<std::vector<std::array<std::complex<double>, 3>>> p_constants_by_ytype_layer;

    // Per-layer metadata for the collapse, in solve units. char rather than bool for a stable data().
    std::vector<int>    p_layer_types        = std::vector<int>();
    std::vector<char>   p_layer_is_static    = std::vector<char>();
    std::vector<char>   p_layer_is_incomp    = std::vector<char>();
    std::vector<size_t> p_num_sols_by_layer  = std::vector<size_t>();
    std::vector<double> p_upper_radii_solve  = std::vector<double>();   // layer upper radii (solve units)
    size_t p_start_layer_i          = 0;
    double p_starting_radius_solve  = 0.0;   // radii below this return NaN
    double p_frequency_solve        = 0.0;   // forcing frequency (solve units) for y3 reconstruction

    // Dimensional context: radius_solve = r_si / p_length_conv, and the scales re-dimensionalize the solve-unit y to
    // SI. Identity when the solve ran in SI.
    double p_length_conv  = 1.0;
    double p_disp_scale   = 1.0;   // y1, y3
    double p_stress_scale = 1.0;   // y2, y4
    double p_pot_scale    = 1.0;   // y6   (y5 is unitless)
    // Whether the EOS arrays are still non-dim; when re-dimensionalized (export mode) the gravity and density read
    // during the dynamic-liquid y3 reconstruction are divided back into solve units.
    bool   p_eos_is_nondim = false;
    double p_grav_conv     = 1.0;
    double p_dens_conv     = 1.0;

    // Persisted EOS inputs for the standalone shooting path, in EOS solve units. The dense EOS re-evaluation
    // references them through c_InterpolateEOSInput, so they must outlive c_radial_solver. Unused by the world path.
    std::vector<double> p_eos_in_radius_nd  = std::vector<double>();
    std::vector<double> p_eos_in_density_nd = std::vector<double>();
    std::vector<std::complex<double>> p_eos_in_bulk_nd  = std::vector<std::complex<double>>();
    std::vector<std::complex<double>> p_eos_in_shear_nd = std::vector<std::complex<double>>();
    std::vector<c_InterpolateEOSInput> p_eos_interp_inputs = std::vector<c_InterpolateEOSInput>();

    c_RadialSolutionStorage() = default;

    c_RadialSolutionStorage(
        size_t num_ytypes,
        double* upper_radius_bylayer_ptr,
        size_t num_layers,
        double* radius_array_ptr,
        size_t size_radius_array,
        int degree_l) :
            success(false),
            error_code(0),
            degree_l(degree_l),
            num_ytypes(num_ytypes),
            num_slices(size_radius_array),
            num_layers(num_layers)
    {
        this->eos_solution_uptr = std::make_unique<c_EOSSolution>(
            upper_radius_bylayer_ptr,
            num_layers,
            radius_array_ptr,
            this->num_slices
            );

        // Up to 3 solutions per layer.
        this->shooting_method_steps_taken_vec.resize(3 * this->num_layers);
        for (size_t layer_i = 0; layer_i < this->num_layers; ++layer_i)
        {
            this->shooting_method_steps_taken_vec[3 * layer_i]     = 0;
            this->shooting_method_steps_taken_vec[3 * layer_i + 1] = 0;
            this->shooting_method_steps_taken_vec[3 * layer_i + 2] = 0;
        }

        if (this->eos_solution_uptr.get())
        {
            this->change_radius_array(
                radius_array_ptr,
                size_radius_array,
                false  // not an array change
                );

            this->message = "Radial solution storage initialized successfully.";
        }
        else
        {
            this->error_code = -1;
            this->message = "c_RadialSolutionStorage:: Could not initialize equation of state storage.";
        }
    }

    virtual ~c_RadialSolutionStorage()
    {
        this->eos_solution_uptr.reset();
    }

    c_EOSSolution* get_eos_solution_ptr()
    {
        return this->eos_solution_uptr.get();
    }

    void change_radius_array(
        double* new_radius_array_ptr,
        size_t new_size_radius_array,
        bool array_changed)
    {
        if (this->error_code == 0)
        {
            if (array_changed)
            {
                if (this->eos_solution_uptr.get())
                {
                    this->eos_solution_uptr->change_radius_array(new_radius_array_ptr, new_size_radius_array);
                }

                this->message = "Radius array changed. Radial solution reset.";
                this->success = false;
            }

            this->num_slices = new_size_radius_array;
            this->total_size = static_cast<size_t>(C_MAX_NUM_Y_REAL) * this->num_slices * this->num_ytypes;

            this->full_solution_vec.resize(this->total_size);
            this->complex_love_vec.resize(this->num_ytypes);
        }
    }

    void find_love()
    {
        if (!(this->success && this->eos_solution_uptr->success && this->error_code == 0)) [[unlikely]]
            return;

        // c_find_love is scale-invariant (the displacement and gravity scales cancel), so the solve-unit surface y and
        // surface gravity give the same k, h, l as the SI pair. The SI surface y is cached here for both methods.
        std::complex<double> surface_solutions[C_MAX_NUM_Y];
        this->p_surface_y_si.assign(
            this->num_ytypes * C_MAX_NUM_Y, std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN)
        );

        if (this->p_uses_interpolants)
        {
            const double surface_r_solve =
                this->p_upper_radii_solve.empty() ? 0.0 : this->p_upper_radii_solve.back();
            for (size_t ytype_i = 0; ytype_i < this->num_ytypes; ++ytype_i)
            {
                this->eval_solveunits(surface_r_solve, ytype_i, surface_solutions);
                this->complex_love_vec[ytype_i] =
                    c_find_love(surface_solutions, this->eos_solution_uptr->surface_gravity);
                this->cache_surface_y(ytype_i, surface_solutions);
            }
            return;
        }

        // Matrix path: the surface slice of the grid (real/imag pairs).
        const size_t top_slice_i   = this->num_slices - 1;
        const size_t num_output_ys = C_MAX_NUM_Y_REAL * this->num_ytypes;
        for (size_t ytype_i = 0; ytype_i < this->num_ytypes; ++ytype_i)
        {
            for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i)
            {
                const size_t lhs_y_index = ytype_i * C_MAX_NUM_Y_REAL + y_i * 2;
                double real_part = this->full_solution_vec[top_slice_i * num_output_ys + lhs_y_index];
                double imag_part = this->full_solution_vec[top_slice_i * num_output_ys + lhs_y_index + 1];
                surface_solutions[y_i] = std::complex<double>(real_part, imag_part);
            }
            this->complex_love_vec[ytype_i] =
                c_find_love(surface_solutions, this->eos_solution_uptr->surface_gravity);
            this->cache_surface_y(ytype_i, surface_solutions);
        }
    }

    // Re-dimensionalize a solve-unit y1..y6 vector to SI in place (y5 is unitless).
    void apply_redimensionalization(std::complex<double>* y6) const noexcept
    {
        y6[0] *= this->p_disp_scale;    // y1
        y6[2] *= this->p_disp_scale;    // y3
        y6[1] *= this->p_stress_scale;  // y2
        y6[3] *= this->p_stress_scale;  // y4
        y6[5] *= this->p_pot_scale;     // y6
    }

    // Cache the SI surface y for one ytype from a solve-unit surface y vector.
    void cache_surface_y(size_t ytype_i, const std::complex<double>* surface_solve_units) noexcept
    {
        std::complex<double> y_si[C_MAX_NUM_Y];
        for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i) y_si[y_i] = surface_solve_units[y_i];
        this->apply_redimensionalization(y_si);
        for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i)
            this->p_surface_y_si[ytype_i * C_MAX_NUM_Y + y_i] = y_si[y_i];
    }

    void dimensionalize_data(
        c_NonDimensionalScales* nondim_scales,
        bool redimensionalize,
        bool include_eos = true)
    {
        // include_eos = false leaves the EOS arrays non-dim for reuse across frequency solves (the world path).
        double* full_solution_ptr       = this->full_solution_vec.data();
        c_EOSSolution* eos_solution_ptr = this->get_eos_solution_ptr();
        if (include_eos)
            eos_solution_ptr->dimensionalize_data(nondim_scales, redimensionalize);

        const double displacement_scale = (nondim_scales->second2_conversion / nondim_scales->length_conversion);
        const double stress_scale       = (nondim_scales->mass_conversion / nondim_scales->length3_conversion);
        const double potential_scale    = (1.0 / nondim_scales->length_conversion);

        // The shooting grid is unfilled (re-dimensionalized on the fly by get_radial_solution) and must not be scaled.
        if (this->success && !this->p_uses_interpolants)
        {
            for (size_t solver_i = 0; solver_i < this->num_ytypes; ++solver_i)
            {
                const size_t bc_stride = solver_i * C_MAX_NUM_Y_REAL;
                for (size_t slice_i = 0; slice_i < this->num_slices; ++slice_i)
                {
                    const size_t slice_stride = bc_stride + slice_i * C_MAX_NUM_Y_REAL * this->num_ytypes;
                    // Real and imaginary pairs; y5 (8, 9) is unitless.
                    // y1
                    full_solution_ptr[slice_stride + 0] *= displacement_scale;
                    full_solution_ptr[slice_stride + 1] *= displacement_scale;
                    // y3
                    full_solution_ptr[slice_stride + 4] *= displacement_scale;
                    full_solution_ptr[slice_stride + 5] *= displacement_scale;
                    // y2
                    full_solution_ptr[slice_stride + 2] *= stress_scale;
                    full_solution_ptr[slice_stride + 3] *= stress_scale;
                    // y4
                    full_solution_ptr[slice_stride + 6] *= stress_scale;
                    full_solution_ptr[slice_stride + 7] *= stress_scale;
                    // y6
                    full_solution_ptr[slice_stride + 10] *= potential_scale;
                    full_solution_ptr[slice_stride + 11] *= potential_scale;
                }
            }
        }
    }

    // ================================================================================================================
    // Interpolant-based radial-solution evaluation (the dense calling system)
    // ================================================================================================================

    // Clear the stored shooting interpolants, constants, and metadata before a re-solve.
    void reset_interpolant_storage() noexcept
    {
        this->p_uses_interpolants = false;
        this->p_interp_by_layer_sol.clear();
        this->p_constants_by_ytype_layer.clear();
        this->p_layer_types.clear();
        this->p_layer_is_static.clear();
        this->p_layer_is_incomp.clear();
        this->p_num_sols_by_layer.clear();
        this->p_upper_radii_solve.clear();
        this->p_start_layer_i         = 0;
        this->p_starting_radius_solve = 0.0;
        this->p_frequency_solve       = 0.0;
    }

    // Record the SI-to-solve-unit radius map, the y re-dimensionalization scales, and the EOS unit state.
    void set_dimensional_context(
            double length_conv,
            double disp_scale,
            double stress_scale,
            double pot_scale,
            bool eos_is_nondim,
            double grav_conv,
            double dens_conv) noexcept
    {
        this->p_length_conv   = length_conv;
        this->p_disp_scale    = disp_scale;
        this->p_stress_scale  = stress_scale;
        this->p_pot_scale     = pot_scale;
        this->p_eos_is_nondim = eos_is_nondim;
        this->p_grav_conv     = grav_conv;
        this->p_dens_conv     = dens_conv;
    }

    // Collapsed y1..y6 (solve units) at a solve-unit radius, shooting path only. Returns false and NaN-fills out6
    // if unsolved, below the starting radius, or out of range.
    bool eval_solveunits(double radius_solve, size_t ytype_i, std::complex<double>* out6) const
    {
        const std::complex<double> cNAN(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
        for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i) out6[y_i] = cNAN;

        if (!this->p_uses_interpolants || !this->success) return false;
        if (ytype_i >= this->num_ytypes) return false;
        if (!(radius_solve >= this->p_starting_radius_solve)) return false;   // also rejects NaN

        // Locate the layer. Upper radii ascend; an interface radius belongs to the lower of the two layers, so the
        // first layer whose upper radius is >= the query wins. A tiny relative slack absorbs the exact-surface case.
        size_t target_layer_i = this->num_layers;
        for (size_t layer_i = this->p_start_layer_i; layer_i < this->num_layers; ++layer_i)
        {
            const double upper = this->p_upper_radii_solve[layer_i];
            if (radius_solve <= upper * (1.0 + 1.0e-12) + 1.0e-300) { target_layer_i = layer_i; break; }
        }
        if (target_layer_i >= this->num_layers) target_layer_i = this->num_layers - 1;   // clamp slight surface overshoot

        const size_t num_sols = this->p_num_sols_by_layer[target_layer_i];
        const int  layer_type = this->p_layer_types[target_layer_i];
        const bool is_static  = this->p_layer_is_static[target_layer_i] != 0;
        if (num_sols == 0 || num_sols > 3) return false;

        // CyRK writes two reals per complex y.
        const size_t num_ys = 2 * num_sols;
        std::complex<double> ysol[3][C_MAX_NUM_Y];
        double real_out[2 * C_MAX_NUM_Y];
        for (size_t sol_i = 0; sol_i < num_sols; ++sol_i)
        {
            // CySolverResult::call is non-const; unique_ptr::get() yields a non-const pointer from this const method.
            CySolverResult* interp = this->p_interp_by_layer_sol[target_layer_i][sol_i].get();
            if (!interp) return false;
            interp->call(radius_solve, real_out);
            for (size_t y_i = 0; y_i < num_ys; ++y_i)
                ysol[sol_i][y_i] = std::complex<double>(real_out[2 * y_i], real_out[2 * y_i + 1]);
        }

        const std::array<std::complex<double>, 3>& constants =
            this->p_constants_by_ytype_layer[ytype_i][target_layer_i];

        // Collapse: out6[y] = sum_sol const[sol] * ysol[sol][mapped_y]. The y-index mapping mirrors
        // c_collapse_layer_solution (liquid layers store fewer ys); undefined ys are left NaN.
        const bool calculate_y3 = (layer_type != 0) && (!is_static);   // dynamic liquid reconstructs y3 below
        for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i)
        {
            size_t y_rhs_i;
            if (layer_type == 0)
            {
                y_rhs_i = y_i;          // solid: all 6 ys
            }
            else if (is_static)
            {
                if (y_i == 4)
                {
                    y_rhs_i = 0;        // static liquid: only y5 (stored at index 0)   
                }
                else continue;
            }
            else
            {
                if (y_i < 2)
                {
                    y_rhs_i = y_i;      // dynamic liquid: y1, y2 (0, 1)
                }
                else if (y_i > 3 && y_i < 6)
                {
                    y_rhs_i = y_i - 2;  // y5, y6 (2, 3)
                }
                else continue;          // y3 reconstructed, y4 undefined
            }
            std::complex<double> acc(0.0, 0.0);
            for (size_t sol_i = 0; sol_i < num_sols; ++sol_i)
                acc += constants[sol_i] * ysol[sol_i][y_rhs_i];
            out6[y_i] = acc;
        }

        if (calculate_y3)
        {
            // y3 = (1/(w^2 r)) (y1 g - y2/rho - y5) in solve units; gravity and density come from this layer's own
            // slices so an interface radius takes this layer's density rather than the neighbor's.
            const double eos_r = this->p_eos_is_nondim ? radius_solve : radius_solve * this->p_length_conv;
            double g_solve = 0.0, rho_solve = 0.0;
            this->eos_solution_uptr->interp_structure_in_layer(target_layer_i, eos_r, &g_solve, &rho_solve);
            if (!this->p_eos_is_nondim) { g_solve /= this->p_grav_conv; rho_solve /= this->p_dens_conv; }
            const double w = this->p_frequency_solve;
            out6[2] = (1.0 / (w * w * radius_solve)) * (out6[0] * g_solve - out6[1] / rho_solve - out6[4]);
        }
        return true;
    }

    // Collapsed y1..y6 (SI) at an SI radius for one ytype: the shooting path evaluates the dense interpolants, the
    // matrix path linearly interpolates its grid. Returns false and NaN-fills out6 on failure or out of range.
    bool get_radial_solution(double radius_si, size_t ytype_i, std::complex<double>* out6) const
    {
        const std::complex<double> cNAN(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);

        if (this->p_uses_interpolants)
        {
            const double radius_solve = radius_si / this->p_length_conv;
            if (!this->eval_solveunits(radius_solve, ytype_i, out6))
            {
                for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i) out6[y_i] = cNAN;
                return false;
            }
            this->apply_redimensionalization(out6);
            return true;
        }

        // Matrix path: full_solution_vec is already SI.
        for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i) out6[y_i] = cNAN;
        if (!this->success || ytype_i >= this->num_ytypes || this->num_slices < 2) return false;

        // The matrix method's grid
        const std::vector<double>& rad = this->p_matrix_radius_solve;
        const double eos_r = this->p_eos_is_nondim ? (radius_si / this->p_length_conv) : radius_si;
        const size_t n = this->num_slices;
        if (n == 0 || rad.size() < n) return false;
        if (eos_r < rad[0] || eos_r > rad[n - 1]) return false;

        // Locate bracketing slices [j, j+1].
        size_t j = 0;
        while (j + 1 < n && rad[j + 1] < eos_r) ++j;
        if (j + 1 >= n) j = n - 2;
        const double r0 = rad[j], r1 = rad[j + 1];
        const double frac = (r1 > r0) ? (eos_r - r0) / (r1 - r0) : 0.0;

        const size_t num_output_ys = C_MAX_NUM_Y_REAL * this->num_ytypes;   // doubles per slice
        const double* fv = this->full_solution_vec.data();
        for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i)
        {
            const size_t off0 = j * num_output_ys + ytype_i * C_MAX_NUM_Y_REAL + y_i * 2;
            const size_t off1 = off0 + num_output_ys;
            const std::complex<double> v0(fv[off0], fv[off0 + 1]);
            const std::complex<double> v1(fv[off1], fv[off1 + 1]);
            out6[y_i] = v0 + (v1 - v0) * frac;
        }
        return true;
    }

    // y1..y6 at n SI radii for one ytype; out holds n * C_MAX_NUM_Y complex values laid out [radius][y].
    void get_radial_solution_array(
        const double* radii_si, size_t n, size_t ytype_i, std::complex<double>* out) const
    {
        for (size_t i = 0; i < n; ++i)
            this->get_radial_solution(radii_si[i], ytype_i, &out[i * C_MAX_NUM_Y]);
    }

    // SI surface y1..y6 for one ytype from the find_love cache; false and NaN-filled if unavailable.
    bool get_surface_y(size_t ytype_i, std::complex<double>* out6) const
    {
        const std::complex<double> cNAN(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
        if (!this->success || ytype_i >= this->num_ytypes
            || this->p_surface_y_si.size() < (ytype_i + 1) * C_MAX_NUM_Y)
        {
            for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i) out6[y_i] = cNAN;
            return false;
        }
        for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i)
            out6[y_i] = this->p_surface_y_si[ytype_i * C_MAX_NUM_Y + y_i];
        return true;
    }

    // The layer holding an SI radius and that radius in the interpolant's solve-unit domain. False when there is
    // nothing to read. Shared by the EOS readers below so they always agree on which layer a radius belongs to.
    bool p_locate_eos(double radius_si, size_t& layer_out, double& solve_radius_out) const
    {
        if (!this->eos_solution_uptr || !this->success) return false;
        solve_radius_out = radius_si / this->p_length_conv;
        layer_out = (this->num_layers == 0) ? 0 : this->num_layers - 1;
        for (size_t layer_i = 0; layer_i < this->num_layers; ++layer_i)
        {
            const double upper = (layer_i < this->p_upper_radii_solve.size())
                ? this->p_upper_radii_solve[layer_i]
                : TidalPyConstants::d_INF;
            if (solve_radius_out <= upper * (1.0 + 1.0e-12) + 1.0e-300)
            {
                layer_out = layer_i;
                break;
            }
        }
        return true;
    }

    // Dense EOS evaluation at an SI radius: the radius is converted into the interpolant's solve-unit domain, the
    // layer located there, and eos->call re-dimensionalizes the outputs. out holds C_EOS_DY_VALUES doubles in the
    // evaluation layout of eos_layout_.hpp, which is frequency-independent: [0] gravity [1] pressure [2] mass
    // [3] moi [4] density [5] shear modulus [6] bulk modulus [7,8] viscosities [9] temperature [10] heat flow
    // [11] melt fraction.
    bool get_eos_si(double radius_si, double* out) const
    {
        size_t target_layer_i = 0;
        double solve_r        = 0.0;
        if (!this->p_locate_eos(radius_si, target_layer_i, solve_r)) return false;
        this->eos_solution_uptr->call(target_layer_i, solve_r, out);
        return true;
    }

    // The complex moduli the EOS itself carries at an SI radius, which is a real answer only on a path that was
    // handed its moduli as arrays (see get_complex_moduli_si).
    bool get_eos_material_si(double radius_si, c_EOSMaterialState& out) const
    {
        size_t target_layer_i = 0;
        double solve_r        = 0.0;
        if (!this->p_locate_eos(radius_si, target_layer_i, solve_r)) return false;
        this->eos_solution_uptr->call_material(target_layer_i, solve_r, out);
        return true;
    }

    // Fill the SI grid from get_radial_solution over the EOS radius grid, for the array-returning standalone API.
    void sample_onto_grid()
    {
        if (!this->success) return;
        const c_EOSSolution* eos = this->eos_solution_uptr.get();
        const std::vector<double>& rad = eos->radius_array_vec;
        const size_t n = this->num_slices;
        const double length_conv = this->p_eos_is_nondim ? this->p_length_conv : 1.0;  // EOS-units radius -> SI
        const size_t num_output_ys = C_MAX_NUM_Y_REAL * this->num_ytypes;
        std::complex<double> out6[C_MAX_NUM_Y];
        for (size_t slice_i = 0; slice_i < n; ++slice_i)
        {
            const double radius_si = rad[slice_i] * length_conv;
            for (size_t ytype_i = 0; ytype_i < this->num_ytypes; ++ytype_i)
            {
                this->get_radial_solution(radius_si, ytype_i, out6);
                for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i)
                {
                    const size_t base = slice_i * num_output_ys + ytype_i * C_MAX_NUM_Y_REAL + y_i * 2;
                    this->full_solution_vec[base]     = out6[y_i].real();
                    this->full_solution_vec[base + 1] = out6[y_i].imag();
                }
            }
        }
    }
};
