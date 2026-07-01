// rs_solution_.hpp - Radial solver solution storage class
// Ported from TidalPy/RadialSolver/rs_solution_.hpp + rs_solution_.cpp
#pragma once

#include <cstring>
#include <cmath>
#include <array>
#include <complex>
#include <vector>
#include <memory>
#include <string>

#include "love_.hpp"
#include "rs_constants_.hpp"
#include "../Material_x/eos/eos_solution_.hpp"   // also provides CyRK's CySolverResult (complete type)
#include "../Material_x/eos/methods/interpolate_.hpp"  // c_InterpolateEOSInput (persisted standalone EOS args)
#include "../../constants_.hpp"
#include "../utilities/dimensions/nondimensional_.hpp"
#include "../utilities/arrays/interp_.hpp"        // cf_interp / cf_binary_search_with_guess (shared array-interp)


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

    // Cached surface y-solution (SI), laid out [ytype * C_MAX_NUM_Y + y_index]. Computed once per solve (both the
    // shooting and matrix methods) inside find_love, because the surface y - and the Love numbers derived from it -
    // is always needed. Interior y-values stay on demand via get_radial_solution.
    std::vector<std::complex<double>> p_surface_y_si = std::vector<std::complex<double>>();

    // Diagnostic data
    std::vector<size_t> shooting_method_steps_taken_vec = std::vector<size_t>();

    // ================================================================================================================
    // Interpolant-based (shooting) solution
    //
    // When p_uses_interpolants is true, the y-solution is served on demand from per-(layer, solution) dense CyRK
    // interpolants plus the collapse constants below, instead of the gridded full_solution_vec. get_radial_solution
    // evaluates the interpolants at any radius (no linear interpolation between grid slices), then collapses and
    // re-dimensionalizes. The matrix method leaves this false and keeps filling full_solution_vec; its
    // get_radial_solution linearly interpolates that grid.
    // ================================================================================================================
    bool p_uses_interpolants = false;

    // Dense CyRK results: [layer][solution]. Owns the force-retained integrators (dense output captured).
    std::vector<std::vector<std::unique_ptr<CySolverResult>>> p_interp_by_layer_sol;

    // Collapse constants: [ytype][layer][solution] (<= 3 solutions per layer). Found from the surface boundary
    // condition and propagated down through the interfaces (exactly as the gridded collapse does), but stored so any
    // radius can be collapsed on demand.
    std::vector<std::vector<std::array<std::complex<double>, 3>>> p_constants_by_ytype_layer;

    // Per-layer metadata needed to collapse at a radius. Radii/frequency are in solve units (whatever the shooting
    // solver ran in: non-dim in the world path, SI when nondimensionalize is off). char (not bool) for a stable data().
    std::vector<int>    p_layer_types        = std::vector<int>();
    std::vector<char>   p_layer_is_static    = std::vector<char>();
    std::vector<char>   p_layer_is_incomp    = std::vector<char>();
    std::vector<size_t> p_num_sols_by_layer  = std::vector<size_t>();
    std::vector<double> p_upper_radii_solve  = std::vector<double>();   // layer upper radii (solve units)
    size_t p_start_layer_i          = 0;
    double p_starting_radius_solve  = 0.0;   // radii below this return NaN
    double p_frequency_solve        = 0.0;   // forcing frequency (solve units) for y3 reconstruction

    // Dimensional context: maps an external SI radius into solve units (radius_solve = r_si / p_length_conv) and
    // re-dimensionalizes the solve-unit y-solution to SI. All default to identity (used when the solve ran in SI).
    double p_length_conv  = 1.0;
    double p_disp_scale   = 1.0;   // y1, y3
    double p_stress_scale = 1.0;   // y2, y4
    double p_pot_scale    = 1.0;   // y6   (y5 is unitless)
    // EOS unit state. The EOS is queried for gravity/density during y3 reconstruction of dynamic-liquid layers; when
    // the EOS has been re-dimensionalized (export mode) its outputs are divided back into solve units.
    bool   p_eos_is_nondim = false;
    double p_grav_conv     = 1.0;
    double p_dens_conv     = 1.0;

    // Persisted EOS inputs for the standalone (c_radial_solver) shooting path. The EOS cysolver's dense extra-output
    // re-invoke (density + complex moduli) reads its stored c_InterpolateEOSInput, which only references, never
    // copies, these arrays. They must therefore (a) outlive c_radial_solver (the original per-call locals dangled,
    // corrupting any post-solve eos->call), and (b) stay in the EOS solve units (non-dim when the solve ran non-dim),
    // so c_EOSSolution::call's nondim_status redim restores SI. get_eos_si then evaluates the dense interpolant at a
    // solve-unit radius and gets correct SI back. The world path never fills these (it injects array-interp EOS data).
    std::vector<double> p_eos_in_radius_nd  = std::vector<double>();
    std::vector<double> p_eos_in_density_nd = std::vector<double>();
    std::vector<std::complex<double>> p_eos_in_bulk_nd  = std::vector<std::complex<double>>();
    std::vector<std::complex<double>> p_eos_in_shear_nd = std::vector<std::complex<double>>();
    std::vector<c_InterpolateEOSInput> p_eos_interp_inputs = std::vector<c_InterpolateEOSInput>();

    // Default constructor
    c_RadialSolutionStorage() = default;

    // Main constructor
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
        // Create equation of state class instance
        this->eos_solution_uptr = std::make_unique<c_EOSSolution>(
            upper_radius_bylayer_ptr,
            num_layers,
            radius_array_ptr,
            this->num_slices
            );

        // Setup diagnostic array. Size = max possible number of solutions per layer (3) * num layers
        this->shooting_method_steps_taken_vec.resize(3 * this->num_layers);
        for (size_t layer_i = 0; layer_i < this->num_layers; ++layer_i)
        {
            this->shooting_method_steps_taken_vec[3 * layer_i]     = 0;
            this->shooting_method_steps_taken_vec[3 * layer_i + 1] = 0;
            this->shooting_method_steps_taken_vec[3 * layer_i + 2] = 0;
        }

        // Setup radius array based vectors
        if (this->eos_solution_uptr.get())
        {
            this->change_radius_array(
                radius_array_ptr,
                size_radius_array,
                false  // We are in initialization, this is not an array change.
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

            // Love number structure for each ytype
            this->complex_love_vec.resize(this->num_ytypes);
        }
    }

    void find_love()
    {
        if (!(this->success && this->eos_solution_uptr->success && this->error_code == 0)) [[unlikely]]
            return;

        // c_find_love is scale-invariant: it uses (surface y, surface gravity) and the displacement/gravity scales
        // cancel exactly, so the solve-unit (non-dim, pre-redimensionalization) pair gives the same k, h, l as the SI
        // pair. find_love runs before any re-dimensionalization, so surface_gravity is in solve units here.
        //
        // The surface y is always needed (it yields the Love numbers), so it is cached here once per solve - in SI -
        // for both methods. p_surface_y_si feeds get_surface_y / get_love_surface_y without recomputation; interior
        // y-values remain on demand via get_radial_solution.
        std::complex<double> surface_solutions[C_MAX_NUM_Y];
        this->p_surface_y_si.assign(
            this->num_ytypes * C_MAX_NUM_Y, std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN)
        );

        if (this->p_uses_interpolants)
        {
            // Shooting (interpolant) path: evaluate the collapsed surface y from the dense interpolants.
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

        // Matrix path: read the surface slice of the gridded full_solution_vec (stored as real/imag pairs).
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

    // Re-dimensionalize a solve-unit y1..y6 vector to SI in place (y5 is unitless). Shared by get_radial_solution and
    // the surface-y cache so both apply identical scaling.
    void apply_redimensionalization(std::complex<double>* y6) const noexcept
    {
        y6[0] *= this->p_disp_scale;    // y1
        y6[2] *= this->p_disp_scale;    // y3
        y6[1] *= this->p_stress_scale;  // y2
        y6[3] *= this->p_stress_scale;  // y4
        y6[5] *= this->p_pot_scale;     // y6
    }

    // Store the SI surface y for one ytype into the cache from a solve-unit surface y vector.
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
        // Perform dimensionalization on the EOS solution first.
        // include_eos = false redimensionalizes only the y-solution, leaving the EOS
        // arrays untouched (used by the cached world radial solver, which keeps its
        // EOS arrays permanently non-dim and reuses them across frequency solves).
        double* full_solution_ptr       = this->full_solution_vec.data();
        c_EOSSolution* eos_solution_ptr = this->get_eos_solution_ptr();
        if (include_eos)
            eos_solution_ptr->dimensionalize_data(nondim_scales, redimensionalize);

        const double displacement_scale = (nondim_scales->second2_conversion / nondim_scales->length_conversion);
        const double stress_scale       = (nondim_scales->mass_conversion / nondim_scales->length3_conversion);
        const double potential_scale    = (1.0 / nondim_scales->length_conversion);

        // The shooting method serves its y-solution from the dense interpolants and re-dimensionalizes on the fly in
        // get_radial_solution, so the gridded full_solution_vec is not filled and must not be scaled here. (The EOS
        // arrays above are still re-dimensionalized when include_eos is set.)
        if (this->success && !this->p_uses_interpolants)
        {
            for (size_t solver_i = 0; solver_i < this->num_ytypes; ++solver_i)
            {
                const size_t bc_stride = solver_i * C_MAX_NUM_Y_REAL;
                for (size_t slice_i = 0; slice_i < this->num_slices; ++slice_i)
                {
                    const size_t slice_stride = bc_stride + slice_i * C_MAX_NUM_Y_REAL * this->num_ytypes;
                    // y1 (real and imag)
                    full_solution_ptr[slice_stride + 0] *= displacement_scale;
                    full_solution_ptr[slice_stride + 1] *= displacement_scale;

                    // y3 (real and imag)
                    full_solution_ptr[slice_stride + 4] *= displacement_scale;
                    full_solution_ptr[slice_stride + 5] *= displacement_scale;

                    // y2 (real and imag)
                    full_solution_ptr[slice_stride + 2] *= stress_scale;
                    full_solution_ptr[slice_stride + 3] *= stress_scale;

                    // y4 (real and imag)
                    full_solution_ptr[slice_stride + 6] *= stress_scale;
                    full_solution_ptr[slice_stride + 7] *= stress_scale;

                    // y5 is unitless - no conversion needed

                    // y6 (real and imag)
                    full_solution_ptr[slice_stride + 10] *= potential_scale;
                    full_solution_ptr[slice_stride + 11] *= potential_scale;
                }
            }
        }
    }

    // ================================================================================================================
    // Interpolant-based radial-solution evaluation (the dense calling system)
    // ================================================================================================================

    // Clear any previously stored shooting interpolants / constants / metadata. Called at the start of a shooting
    // solve so a re-solve (e.g. a new forcing frequency) does not accumulate stale interpolators.
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

    // Record how to map an external SI radius into solve units and re-dimensionalize the solve-unit y-solution. Set by
    // the solve wrapper once the non-dimensionalization scales are known. eos_is_nondim reflects whether the EOS
    // arrays are currently non-dim (fast world path) or have been re-dimensionalized to SI (export path).
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

    // Evaluate the collapsed y1..y6 at a radius given in solve units, writing solve-unit complex values into out6
    // (length C_MAX_NUM_Y). Returns false (and NaN-fills) if unsolved, below the starting radius, or out of range.
    // Shooting (interpolant) path only.
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
        for (size_t layer_i = this->p_start_layer_i; layer_i < this->num_layers; ++L)
        {
            const double upper = this->p_upper_radii_solve[L];
            if (radius_solve <= upper * (1.0 + 1.0e-12) + 1.0e-300) { target_layer_i = layer_i; break; }
        }
        if (target_layer_i >= this->num_layers) target_layer_i = this->num_layers - 1;   // clamp slight surface overshoot

        const size_t num_sols = this->p_num_sols_by_layer[target_layer_i];
        const int  layer_type = this->p_layer_types[target_layer_i];
        const bool is_static  = this->p_layer_is_static[target_layer_i] != 0;
        if (num_sols == 0 || num_sols > 3) return false;

        // Evaluate each independent solution's dense interpolant at this radius (CyRK writes 2x real per complex y).
        const size_t num_ys = 2 * num_sols;            // y-values defined for this layer type
        std::complex<double> ysol[3][C_MAX_NUM_Y];
        double real_out[2 * C_MAX_NUM_Y];              // up to num_ys complex -> 2*num_ys doubles
        for (size_t sol_i = 0; sol_i < num_sols; ++sol_i)
        {
            // unique_ptr::get() yields a non-const CySolverResult* even from this const method; call() is non-const
            // (mirrors how c_EOSSolution::call invokes its stored interpolants).
            CySolverResult* interp = this->p_interp_by_layer_sol[target_layer_i][sol_i].get();
            if (!interp) return false;
            interp->call(radius_solve, real_out);
            for (size_t y_i = 0; y_i < num_ys; ++y_i)
                ysol[sol_i][y_i] = std::complex<double>(real_out[2 * y_i], real_out[2 * y_i + 1]);
        }

        const std::array<std::complex<double>, 3>& constants =
            this->p_constants_by_ytype_layer[ytype_i][layer_i];

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
            // y3 = (1/(w^2 r)) (y1 g - y2/rho - y5), all in solve units. Gravity/density are read from the EOS's
            // structural arrays via the shared array-interp utility (cf_interp + cf_binary_search_with_guess) - the
            // same arrays the gridded collapse used (layer_gravity_ptr/layer_density_ptr), so the dense y3 reproduces
            // the grid y3. (eos->call would route through the cysolver dense output, whose unit handling differs across
            // solve paths.) Mirrors c_EOSSolution::_call_interp_arrays.
            const double eos_r = this->p_eos_is_nondim ? radius_solve : radius_solve * this->p_length_conv;
            c_EOSSolution* eos = this->eos_solution_uptr.get();
            double*      r_arr = eos->radius_array_vec.data();
            const size_t n     = eos->radius_array_vec.size();
            size_t j_guess = 0;
            if (n > 1 && r_arr[n - 1] > r_arr[0]) {
                double frac = (eos_r - r_arr[0]) / (r_arr[n - 1] - r_arr[0]);
                frac = frac < 0.0 ? 0.0 : (frac > 1.0 ? 1.0 : frac);
                j_guess = static_cast<size_t>(static_cast<double>(n) * frac);
                if (j_guess >= n) j_guess = n - 1;
            }
            int b_code = 0;
            size_t j   = cf_binary_search_with_guess(eos_r, r_arr, n, j_guess, &b_code);
            double desired = eos_r, g_solve = 0.0, rho_solve = 0.0;
            cf_interp(&desired, r_arr, eos->gravity_array_vec.data(), n, &j, &g_solve);
            cf_interp(&desired, r_arr, eos->density_array_vec.data(), n, &j, &rho_solve);
            if (!this->p_eos_is_nondim) { g_solve /= this->p_grav_conv; rho_solve /= this->p_dens_conv; }
            const double w = this->p_frequency_solve;
            out6[2] = (1.0 / (w * w * radius_solve)) * (out6[0] * g_solve - out6[1] / rho_solve - out6[4]);
        }
        return true;
    }

    // Evaluate the collapsed y1..y6 at an external SI radius for one ytype, writing SI complex values into out6.
    // Dispatches on the solve method: shooting evaluates the dense interpolants (accurate at any radius); the matrix
    // method linearly interpolates its constructed grid. Returns false (and NaN-fills) on failure / out of range.
    bool get_radial_solution(double radius_si, size_t ytype_i, std::complex<double>* out6) const
    {
        const std::complex<double> cNAN(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
        
        // Shooting method path
        if (this->p_uses_interpolants)
        {
            const double radius_solve = radius_si / this->p_length_conv;
            if (!this->eval_solveunits(radius_solve, ytype_i, out6))
            {
                for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i) out6[y_i] = cNAN;
                return false;
            }
            // Re-dimensionalize solve-unit ys to SI (y5 is unitless).
            this->apply_redimensionalization(out6);
            return true;
        }

        // Matrix path: full_solution_vec is already SI. Linearly interpolate it on the EOS radius grid.
        for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i) out6[y_i] = cNAN;
        if (!this->success || ytype_i >= this->num_ytypes || this->num_slices < 2) return false;

        const c_EOSSolution* eos = this->eos_solution_uptr.get();
        const std::vector<double>& rad = eos->radius_array_vec;     // EOS units (non-dim or SI per p_eos_is_nondim)
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

    // Vectorized convenience: evaluate y1..y6 at n SI radii for one ytype. out must hold n * C_MAX_NUM_Y complex
    // values, laid out [radius][y]. Each radius is evaluated independently via get_radial_solution.
    void get_radial_solution_array(
        const double* radii_si, size_t n, size_t ytype_i, std::complex<double>* out) const
    {
        for (size_t i = 0; i < n; ++i)
            this->get_radial_solution(radii_si[i], ytype_i, &out[i * C_MAX_NUM_Y]);
    }

    // SI y1..y6 at the planet surface for one ytype, served from the cache that find_love fills once per solve
    // (method-agnostic). Returns false (and NaN-fills) if unsolved or the cache is unavailable.
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

    // SI-aware dense EOS evaluation at an arbitrary radius. Routes through the solution's own dense EOS
    // interpolant (not a re-interpolation of the gridded arrays): converts the external SI radius into the
    // interpolant's solve-unit (non-dim) domain, locates the layer in solve units, then lets eos->call
    // redimensionalize the outputs back to SI (the cysolver interpolant lives in the non-dim radius domain,
    // so an SI radius would extrapolate). out must hold C_EOS_DY_VALUES doubles:
    //   [0] gravity [1] pressure [2] mass [3] moi [4] density [5,6] shear re/im [7,8] bulk re/im.
    bool get_eos_si(double radius_si, double* out) const
    {
        if (!this->eos_solution_uptr || !this->success) return false;
        const double solve_r = radius_si / this->p_length_conv;
        size_t layer_i = (this->num_layers == 0) ? 0 : this->num_layers - 1;
        for (size_t layer_i = 0; layer_i < this->num_layers; ++layer_i)
        {
            const double upper = (layer_i < this->p_upper_radii_solve.size())
                ? this->p_upper_radii_solve[layer_i]
                : TidalPyConstants::d_INF;
            if (solve_r <= upper * (1.0 + 1.0e-12) + 1.0e-300)
            {
                layer_i = layer_i;
                break;
            }
        }
        this->eos_solution_uptr->call(layer_i, solve_r, out);
        return true;
    }

    // Fill the gridded full_solution_vec (SI) by sampling get_radial_solution over the EOS radius grid. Used by the
    // standalone solver to keep the legacy array-returning API working after a shooting solve; the world path does
    // not call this (it queries get_radial_solution at the radii it needs).
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
