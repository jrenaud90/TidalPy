// rs_solution_.hpp: radial solver solution storage.
#pragma once

#include <algorithm>
#include <cstdio>
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
#include "layer_kind_.hpp"   // c_layer_layout
#include "matrix_types/solid_matrix_.hpp"      // c_fundamental_matrix, to continue a propagation-matrix solution
#include "../Material/eos/eos_solution_.hpp"   // also provides CyRK's CySolverResult (complete type)
#include "../constants_.hpp"
#include "../Utilities/dimensions/nondimensional_.hpp"


// Error Codes:
// -1 : Equation of State storage (c_EOSSolution) could not be initialized.
// -2 : (set by python wrapper) Unknown / Unsupported boundary condition provided.
// -5 : There was a problem with the inputs to radial solver (including a starting radius outside the planet, above
//      [numerical] max_start_radius_fraction of it, or with no radial slice above it in its layer)
//
// -1X : Error in shooting method
// -10 : Error in finding starting conditions
// -11 : Numerical integration failed
// -12 : The surface boundary condition solve returned non-finite constants
// -13 : The surface boundary condition system is singular to working precision (surface_rcond below
//       [numerical] minimum_surface_rcond), or singular by construction: degree 1 with every integrated layer
//       static, where a rigid translation meets every surface condition (both methods)
// -14 : Unknown surface boundary condition model, or an unsupported number of them (both methods)
// -15 : A compressible layer that reads the bulk modulus has a non-positive one
//
// -2X : Error in propagation matrix method
// -20 : Unknown core starting conditions
// -21 : The surface boundary condition solve returned non-finite constants
// -22 : A core starting condition other than the regular solution (core_model 1 to 4) with a manual starting radius
// -23 : The world's structure is not the single solid, static, incompressible layer the method requires

class c_RadialSolutionStorage
{
public:
    bool success        = false;
    int error_code      = -100;
    int degree_l        = 0;
    std::string message = "No Message Set.";
    size_t num_ytypes   = 0;

    // Marks the solve failed with an error code (the list above) and a message, printed as well when verbose. Returns
    // the code, so a solver can `return storage->fail(...)`.
    int fail(int code, const std::string& text, bool verbose) noexcept
    {
        this->success    = false;
        this->error_code = code;
        this->message    = text;
        if (verbose) { std::printf("%s", this->message.c_str()); }
        return code;
    }

    size_t num_slices   = 0;
    size_t num_layers   = 0;
    size_t total_size   = 0;

    // Equation of state solution
    std::unique_ptr<c_EOSSolution> eos_solution_uptr = nullptr;

    // Radial solution results (stores double-pairs for complex values)
    std::vector<double> full_solution_vec = std::vector<double>();

    // Love number attributes (stores double-pairs for complex values)
    std::vector<c_LoveNumbers> complex_love_vec = std::vector<c_LoveNumbers>();

    // Surface y-solution (SI), laid out [ytype * C_MAX_NUM_Y + y_index]; find_love caches it once per solve.
    std::vector<std::complex<double>> p_surface_y_si = std::vector<std::complex<double>>();

    // The surface boundary conditions this solve produced blocks for, in order.
    std::vector<int> p_bc_models = std::vector<int>();

    // Diagnostic data
    std::vector<size_t> shooting_method_steps_taken_vec = std::vector<size_t>();

    // Worst-case error amplification of the surface boundary-condition solve across ytypes; shooting method
    // only, and 0 for the matrix method. See c_estimate_surface_amplification in boundaries_.hpp.
    double surface_amplification = 0.0;

    // Equilibrated reciprocal condition number of the surface boundary-condition system, the rank measure the
    // amplification cannot give; shooting method only, NaN for the matrix method or before the collapse. See
    // c_estimate_surface_rcond in boundaries_.hpp.
    double surface_rcond = TidalPyConstants::d_NAN;

    // When true, get_radial_solution evaluates the per-(layer, solution) dense CyRK interpolants at any
    // radius and collapses them with the constants below. The matrix method leaves it false and fills
    // full_solution_vec on its grid, and get_radial_solution continues its propagation to the radius asked for.
    bool p_uses_interpolants = false;

    // A shooting solve run for its Love numbers alone (love_only) keeps no dense interpolants, only each surface-layer
    // solution's y at the surface, [solution * C_MAX_NUM_Y + y] (solve units), which find_love collapses. Its radial
    // functions anywhere below the surface are unavailable.
    bool p_love_only = false;
    std::vector<std::complex<double>> p_surface_top_y = std::vector<std::complex<double>>();

    // The propagation matrix's own radius grid (solve units). Empty after a shooting solve, which grids nothing.
    std::vector<double> p_matrix_radius_solve = std::vector<double>();

    // What the propagation matrix needs to continue its solution to any radius, so get_radial_solution evaluates
    // the propagation itself rather than interpolating between grid radii (solve units, SVC16 y order and signs).
    // Slice i is the shell (r_{i-1}, r_i] of slice i's material, and inside it y(r) = Y_i(r) B_i c, with Y_i the
    // fundamental matrix of that material, B_i = Y_i(r_{i-1})^-1 P_{i-1} the 6 x 3 block laid out
    // [slice * 18 + row * 3 + column], and c the surface constants of a ytype, [ytype * 3 + column]. Below the seed
    // radius, one slice under the first propagated one, a regular seed continues as y(r) = Y_seed(r)[:, 0:3] c.
    std::vector<std::complex<double>> p_matrix_shell_coeffs = std::vector<std::complex<double>>();
    std::vector<std::complex<double>> p_matrix_constants    = std::vector<std::complex<double>>();
    std::vector<double>               p_matrix_density      = std::vector<double>();
    std::vector<double>               p_matrix_gravity      = std::vector<double>();
    std::vector<std::complex<double>> p_matrix_shear        = std::vector<std::complex<double>>();
    size_t p_matrix_first_slice   = 0;       // first propagated slice; at least 2
    bool   p_matrix_regular_core  = false;   // the seed is the regular solution, so it continues to r = 0
    double p_matrix_G             = 0.0;     // gravitational constant (solve units)
    int    p_matrix_degree_l      = 0;

    // A static liquid surface layer: the surface y2 of each ytype (its boundary condition) and the liquid's
    // density and gravity there, which give the surface y1 = (y2 / rho + y5) / g (solve units). Empty otherwise.
    std::vector<std::complex<double>> p_static_surface_y2 = std::vector<std::complex<double>>();
    double p_static_surface_density = TidalPyConstants::d_NAN;
    double p_static_surface_gravity = TidalPyConstants::d_NAN;

    // Forcing frequency [rad s-1] of the last solve; NaN before one.
    double p_love_frequency_si = TidalPyConstants::d_NAN;

    /// The complex moduli [Pa] at an SI radius, from the same provider call the integrator made: a rheology
    /// solve reports its rheology at the solved frequency, a supplied-moduli solve the profile it was given.
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

        // call_material answers in solve units, which is what the integrator wants; this is a readout, so
        // convert back to SI. The scale is one in a solve that ran dimensional.
        c_EOSMaterialState material_state;
        this->eos_solution_uptr->call_material(layer_i, solve_r, material_state);
        const double pascal_scale = this->eos_solution_uptr->p_structure_pascal_scale;
        shear_out = material_state.shear_modulus * pascal_scale;
        bulk_out  = material_state.bulk_modulus * pascal_scale;
    }

    // Dense CyRK results [layer][solution]; owns the force-retained integrators.
    std::vector<std::vector<std::unique_ptr<CySolverResult>>> p_interp_by_layer_sol;

    // Collapse constants [ytype][layer][solution], at most 3 solutions per layer.
    std::vector<std::vector<std::array<std::complex<double>, 3>>> p_constants_by_ytype_layer;

    // Per-layer metadata for the collapse. char rather than bool for a stable data().
    std::vector<int>    p_layer_types        = std::vector<int>();
    std::vector<char>   p_layer_is_static    = std::vector<char>();
    std::vector<size_t> p_num_sols_by_layer  = std::vector<size_t>();
    std::vector<double> p_upper_radii_solve  = std::vector<double>();   // layer upper radii (solve units)
    std::vector<double> p_sample_radius_si   = std::vector<double>();   // radii the result grid is sampled on [m]
    size_t p_start_layer_i          = 0;
    double p_starting_radius_solve  = 0.0;   // radii below this return NaN
    double p_frequency_solve        = 0.0;   // forcing frequency (solve units) for y3 reconstruction

    // radius_solve = r_si / p_length_conv, and the scales re-dimensionalize a solve-unit y to SI. Identity
    // when the solve ran in SI. Set once per solve by set_dimensional_context.
    double p_length_conv  = 1.0;
    double p_disp_scale   = 1.0;   // y1, y3
    double p_stress_scale = 1.0;   // y2, y4
    double p_pot_scale    = 1.0;   // y6   (y5 is unitless)
    // Whether the EOS arrays are still non-dim.
    bool   p_eos_is_nondim = false;

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

        // Up to 3 solutions per layer, zero until a shooting solve counts them.
        this->shooting_method_steps_taken_vec.resize(3 * this->num_layers);

        if (this->eos_solution_uptr.get())
        {
            this->total_size = static_cast<size_t>(C_MAX_NUM_Y_REAL) * this->num_slices * this->num_ytypes;
            this->full_solution_vec.resize(this->total_size);
            this->complex_love_vec.resize(this->num_ytypes);

            this->message = "Radial solution storage initialized successfully.";
        }
        else
        {
            this->error_code = -1;
            this->message = "c_RadialSolutionStorage:: Could not initialize equation of state storage.";
        }
    }

    ~c_RadialSolutionStorage() = default;

    c_EOSSolution* get_eos_solution_ptr()
    {
        return this->eos_solution_uptr.get();
    }

    void find_love()
    {
        if (!(this->success && this->eos_solution_uptr->success && this->error_code == 0)) [[unlikely]]
            return;

        // c_find_love is scale invariant, the displacement and gravity scales cancelling, so the solve-unit
        // surface y and gravity give the same k, h, l as the SI pair. The SI surface y is cached here.
        std::complex<double> surface_solutions[C_MAX_NUM_Y];
        this->p_surface_y_si.assign(
            this->num_ytypes * C_MAX_NUM_Y, std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN)
        );

        if (this->p_uses_interpolants)
        {
            const double surface_r_solve =
                this->p_upper_radii_solve.empty() ? 0.0 : this->p_upper_radii_solve.back();
            c_RadialBasis basis;
            const bool basis_found = this->p_evaluate_basis(surface_r_solve, false, basis);
            for (size_t ytype_i = 0; ytype_i < this->num_ytypes; ++ytype_i)
            {
                this->p_collapse_basis(basis, basis_found, ytype_i, surface_r_solve, surface_solutions);
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

    // Re-dimensionalize a solve-unit y1..y6 vector to SI in place; y5 is unitless.
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
        // include_eos = false leaves the EOS arrays non-dim for reuse across frequency solves. The solution
        // itself takes the scales set_dimensional_context stored for this solve.
        double* full_solution_ptr       = this->full_solution_vec.data();
        c_EOSSolution* eos_solution_ptr = this->get_eos_solution_ptr();
        if (include_eos)
            eos_solution_ptr->dimensionalize_data(nondim_scales, redimensionalize);

        const double displacement_scale = this->p_disp_scale;
        const double stress_scale       = this->p_stress_scale;
        const double potential_scale    = this->p_pot_scale;

        // The shooting grid is unfilled, re-dimensionalized on the fly, so it must not be scaled.
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


    void reset_interpolant_storage() noexcept
    {
        this->p_uses_interpolants = false;
        this->p_love_only         = false;
        this->p_surface_top_y.clear();
        this->p_interp_by_layer_sol.clear();
        this->p_constants_by_ytype_layer.clear();
        this->p_layer_types.clear();
        this->p_layer_is_static.clear();
        this->p_num_sols_by_layer.clear();
        this->p_upper_radii_solve.clear();
        this->p_start_layer_i         = 0;
        this->p_starting_radius_solve = 0.0;
        this->p_frequency_solve       = 0.0;
        this->p_static_surface_y2.clear();
        this->p_matrix_shell_coeffs.clear();
        this->p_matrix_constants.clear();
        this->p_matrix_density.clear();
        this->p_matrix_gravity.clear();
        this->p_matrix_shear.clear();
        this->p_matrix_first_slice  = 0;
        this->p_matrix_regular_core = false;
    }

    // The scales that re-dimensionalize this solve's y to SI, from the non-dimensional scales it ran in; null for a
    // solve that ran in SI. A non-dim solve leaves the EOS arrays non-dim until an export converts them.
    void set_dimensional_context(const c_NonDimensionalScales* nondim_scales) noexcept
    {
        if (nondim_scales)
        {
            this->p_length_conv   = nondim_scales->length_conversion;
            this->p_disp_scale    = nondim_scales->second2_conversion / nondim_scales->length_conversion;
            this->p_stress_scale  = nondim_scales->mass_conversion / nondim_scales->length3_conversion;
            this->p_pot_scale     = 1.0 / nondim_scales->length_conversion;
            this->p_eos_is_nondim = true;
        }
        else
        {
            this->p_length_conv   = 1.0;
            this->p_disp_scale    = 1.0;
            this->p_stress_scale  = 1.0;
            this->p_pot_scale     = 1.0;
            this->p_eos_is_nondim = false;
        }
    }

    // The independent solutions of the layer that owns a radius, evaluated there once so every ytype can be
    // collapsed from them.
    struct c_RadialBasis
    {
        size_t layer_i    = 0;
        size_t num_sols   = 0;
        int    layer_type = 0;
        bool   is_static  = false;
        // [solution][y], in the layer's own y storage order; CyRK writes two reals per complex y.
        std::complex<double> ysol[3][C_MAX_NUM_Y];
        // Material at the radius, read only for a dynamic liquid's y3 reconstruction.
        double gravity = TidalPyConstants::d_NAN;
        double density = TidalPyConstants::d_NAN;
    };

    // Evaluate the basis at a solve-unit radius; false when unsolved, below the starting radius, or out of range.
    // A radius on a layer interface belongs to the lower layer, unless `upper_at_interface`, which takes the layer
    // above: the second copy of an interface radius in a layered radius array. A radius above the surface has no
    // solution.
    bool p_evaluate_basis(double radius_solve, bool upper_at_interface, c_RadialBasis& basis) const
    {
        if (!this->p_uses_interpolants || !this->success) return false;
        if (!(radius_solve >= this->p_starting_radius_solve)) return false;   // also rejects NaN

        // Upper radii ascend and an interface radius belongs to the lower layer, so the first layer whose
        // upper radius reaches the query wins (the first whose upper radius is above it, for an upper copy). A
        // little slack absorbs the exact-surface case.
        const double interface_rtol = 1.0e-12;
        size_t target_layer_i = this->num_layers;
        for (size_t layer_i = this->p_start_layer_i; layer_i < this->num_layers; ++layer_i)
        {
            const double upper = this->p_upper_radii_solve[layer_i];
            const bool in_layer = upper_at_interface
                ? (radius_solve < upper * (1.0 - interface_rtol))
                : (radius_solve <= upper * (1.0 + interface_rtol) + 1.0e-300);
            if (in_layer) { target_layer_i = layer_i; break; }
        }
        if (target_layer_i >= this->num_layers)
        {
            // The surface's own upper copy has no layer above it: it is the top layer's. Anything higher is outside.
            const double surface = this->p_upper_radii_solve.empty() ? 0.0 : this->p_upper_radii_solve.back();
            if (!upper_at_interface || !(radius_solve <= surface * (1.0 + interface_rtol) + 1.0e-300)) return false;
            target_layer_i = this->num_layers - 1;
        }

        basis.layer_i    = target_layer_i;
        basis.num_sols   = this->p_num_sols_by_layer[target_layer_i];
        basis.layer_type = this->p_layer_types[target_layer_i];
        basis.is_static  = this->p_layer_is_static[target_layer_i] != 0;
        if (basis.num_sols == 0 || basis.num_sols > 3) return false;

        const size_t num_ys = 2 * basis.num_sols;
        if (this->p_love_only)
        {
            // Only the surface values were kept.
            const double surface = this->p_upper_radii_solve.back();
            if ((target_layer_i + 1 != this->num_layers)
                || !(std::fabs(radius_solve - surface) <= interface_rtol * surface + 1.0e-300)
                || (this->p_surface_top_y.size() < basis.num_sols * C_MAX_NUM_Y)) return false;
            for (size_t sol_i = 0; sol_i < basis.num_sols; ++sol_i)
                for (size_t y_i = 0; y_i < num_ys; ++y_i)
                    basis.ysol[sol_i][y_i] = this->p_surface_top_y[sol_i * C_MAX_NUM_Y + y_i];
        }
        double real_out[2 * C_MAX_NUM_Y] = {};
        for (size_t sol_i = 0; (sol_i < basis.num_sols) && !this->p_love_only; ++sol_i)
        {
            // CySolverResult::call is non-const, so get() is used to escape this method's constness.
            CySolverResult* interp = this->p_interp_by_layer_sol[target_layer_i][sol_i].get();
            if (!c_call_dense_checked(interp, radius_solve, real_out, 2 * num_ys)) return false;
            for (size_t y_i = 0; y_i < num_ys; ++y_i)
                basis.ysol[sol_i][y_i] = std::complex<double>(real_out[2 * y_i], real_out[2 * y_i + 1]);
        }

        if ((basis.layer_type != 0) && (!basis.is_static))
        {
            // Gravity and density are asked of the layer that owns this radius, so an interface takes this layer's
            // density, not the neighbor's. This goes through call_material rather than the stored slice arrays
            // because a world-attached solve grids nothing: those arrays are empty there, and reading them gave an
            // uninitialized density that sent y3 to infinity through a whole dynamic liquid layer.
            c_EOSMaterialState material_state;
            this->eos_solution_uptr->call_material(target_layer_i, radius_solve, material_state);
            basis.gravity = material_state.gravity;
            basis.density = material_state.density;
        }
        return true;
    }

    // Collapse an evaluated basis into y1..y6 (solve units) for one ytype; NaN-filled when the basis was not found or
    // the ytype is out of range. Liquid layers store fewer ys, so undefined ys stay NaN.
    bool p_collapse_basis(
            const c_RadialBasis& basis,
            bool basis_found,
            size_t ytype_i,
            double radius_solve,
            std::complex<double>* out6) const
    {
        const std::complex<double> cNAN(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
        for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i) out6[y_i] = cNAN;
        if (!basis_found || (ytype_i >= this->num_ytypes)) return false;

        const std::array<std::complex<double>, 3>& constants =
            this->p_constants_by_ytype_layer[ytype_i][basis.layer_i];

        // out6[y] = sum_sol const[sol] * ysol[sol][slot of y], for each y the layer kind stores (c_layer_layout); a
        // dynamic liquid's y3 is reconstructed below, and the ys a kind does not store stay NaN.
        const bool calculate_y3 = (basis.layer_type != 0) && (!basis.is_static);   // dynamic liquid, below
        const c_LayerKindLayout& layout = c_layer_layout(basis.layer_type, basis.is_static);
        for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i)
        {
            const size_t y_rhs_i = layout.slot_of_full_y(y_i);
            if (y_rhs_i == C_Y_NOT_STORED) continue;
            std::complex<double> acc(0.0, 0.0);
            for (size_t sol_i = 0; sol_i < basis.num_sols; ++sol_i)
                acc += constants[sol_i] * basis.ysol[sol_i][y_rhs_i];
            out6[y_i] = acc;
        }

        if (calculate_y3)
        {
            // y3 = (1/(w^2 r)) (y1 g - y2/rho - y5) in solve units.
            const double w = this->p_frequency_solve;
            out6[2] = (1.0 / (w * w * radius_solve))
                * (out6[0] * basis.gravity - out6[1] / basis.density - out6[4]);
        }
        else if ((basis.layer_type != 0) && basis.is_static && (basis.layer_i + 1 == this->num_layers) &&
                 (ytype_i < this->p_static_surface_y2.size()) && !this->p_upper_radii_solve.empty())
        {
            // The free surface of a static liquid top layer: y2 is the boundary condition and y2 = rho (g y1 - y5).
            const double interface_rtol = 1.0e-12;
            if (radius_solve >= this->p_upper_radii_solve.back() * (1.0 - interface_rtol))
            {
                out6[1] = this->p_static_surface_y2[ytype_i];
                out6[0] = (out6[1] / this->p_static_surface_density + out6[4]) / this->p_static_surface_gravity;
            }
        }
        return true;
    }

    // Collapsed y1..y6 in solve units, shooting path only. False and NaN-filled when unsolved, below the
    // starting radius, or out of range. get_radial_solution is the SI form. `upper_at_interface` is as in
    // p_evaluate_basis.
    bool get_radial_solution_nondim(
            double radius_solve,
            size_t ytype_i,
            std::complex<double>* out6,
            bool upper_at_interface = false) const
    {
        c_RadialBasis basis;
        const bool basis_found =
            (ytype_i < this->num_ytypes) && this->p_evaluate_basis(radius_solve, upper_at_interface, basis);
        return this->p_collapse_basis(basis, basis_found, ytype_i, radius_solve, out6);
    }

    // Collapsed y1..y6 (SI) for one ytype: the shooting path evaluates the dense interpolants, the matrix
    // path continues its propagation inside the slice (p_matrix_y_solve). False and NaN-filled on failure or out of
    // range.
    bool get_radial_solution(
            double radius_si,
            size_t ytype_i,
            std::complex<double>* out6,
            bool upper_at_interface = false) const
    {
        const std::complex<double> cNAN(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);

        if (this->p_uses_interpolants)
        {
            const double radius_solve = radius_si / this->p_length_conv;
            if (!this->get_radial_solution_nondim(radius_solve, ytype_i, out6, upper_at_interface))
            {
                for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i) out6[y_i] = cNAN;
                return false;
            }
            this->apply_redimensionalization(out6);
            return true;
        }

        // Matrix path: the propagation continued to the radius, which the grid stays in solve units for.
        if (!this->success || !this->p_matrix_y_solve(radius_si / this->p_length_conv, ytype_i, out6))
        {
            for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i) out6[y_i] = cNAN;
            return false;
        }
        this->apply_redimensionalization(out6);
        return true;
    }

    // The propagation-matrix solution y1..y6 (solve units, TS72 convention) at a solve-unit radius, from the
    // propagation itself: inside slice i, y(r) = Y_i(r) B_i c (see p_matrix_shell_coeffs), which reproduces the grid
    // at both ends of the slice, and below the seed radius the regular solution of the seed's material. Gravity
    // inside a slice is linear between its end radii, and inside the seed radius it is that of a uniform sphere,
    // g r / r_seed, as the regular seed assumes. False and NaN-filled where there is no solution: outside the body,
    // before a solve, or below the seed radius when the core seed is not the regular solution.
    bool p_matrix_y_solve(double radius_solve, size_t ytype_i, std::complex<double>* out6) const
    {
        const std::complex<double> cNAN(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
        for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i) out6[y_i] = cNAN;

        const std::vector<double>& rad = this->p_matrix_radius_solve;
        const size_t n = rad.size();
        if ((n < 2) || (this->p_matrix_first_slice < 1) || (this->p_matrix_first_slice >= n) ||
            (ytype_i >= this->num_ytypes) || (this->p_matrix_constants.size() < 3 * (ytype_i + 1)) ||
            (this->p_matrix_shell_coeffs.size() < 18 * n) || (this->p_matrix_density.size() < n) ||
            (this->p_matrix_gravity.size() < n) || (this->p_matrix_shear.size() < n))
        {
            return false;
        }
        const size_t seed_i = this->p_matrix_first_slice - 1;
        // Also rejects NaN; a rounding step above the surface is the surface.
        const double surface = rad[n - 1];
        if (!(radius_solve >= 0.0) || !(radius_solve <= surface * (1.0 + 1.0e-12) + 1.0e-300)) return false;
        double radius_here = std::fmin(radius_solve, surface);

        const std::complex<double>* constants = &this->p_matrix_constants[3 * ytype_i];
        size_t material_i     = seed_i;
        double gravity_here   = TidalPyConstants::d_NAN;
        std::complex<double> column_weights[6];
        if (radius_here <= rad[seed_i])
        {
            if (!this->p_matrix_regular_core) return false;
            // y = Y_seed(r)[:, 0:3] c; the other columns are singular at the center and weigh nothing.
            gravity_here = (rad[seed_i] > 0.0) ? this->p_matrix_gravity[seed_i] * radius_here / rad[seed_i] : 0.0;
            for (size_t column_i = 0; column_i < 6; ++column_i)
                column_weights[column_i] = (column_i < 3) ? constants[column_i] : std::complex<double>(0.0, 0.0);
        }
        else
        {
            // The slice whose shell (r_{i-1}, r_i] holds the radius.
            size_t slice_i = seed_i + 1;
            while ((slice_i + 1 < n) && (rad[slice_i] < radius_here)) ++slice_i;
            material_i = slice_i;
            const double radius_lower = rad[slice_i - 1];
            const double radius_upper = rad[slice_i];
            const double fraction =
                (radius_upper > radius_lower) ? (radius_here - radius_lower) / (radius_upper - radius_lower) : 1.0;
            gravity_here = this->p_matrix_gravity[slice_i - 1] +
                fraction * (this->p_matrix_gravity[slice_i] - this->p_matrix_gravity[slice_i - 1]);
            // B_i c
            const std::complex<double>* shell_coeffs = &this->p_matrix_shell_coeffs[18 * slice_i];
            for (size_t row_i = 0; row_i < 6; ++row_i)
            {
                std::complex<double> weight(0.0, 0.0);
                for (size_t column_i = 0; column_i < 3; ++column_i)
                    weight += shell_coeffs[row_i * 3 + column_i] * constants[column_i];
                column_weights[row_i] = weight;
            }
        }

        // The fundamental matrix of the shell's material at this radius.
        double density_here = this->p_matrix_density[material_i];
        std::complex<double> shear_here = this->p_matrix_shear[material_i];
        std::complex<double> fundamental[36];
        c_fundamental_matrix(
            0,
            1,
            &radius_here,
            &density_here,
            &gravity_here,
            &shear_here,
            fundamental,
            nullptr,
            nullptr,
            this->p_matrix_degree_l,
            this->p_matrix_G);

        std::complex<double> y_svc[6];
        for (size_t row_i = 0; row_i < 6; ++row_i)
        {
            std::complex<double> value(0.0, 0.0);
            for (size_t column_i = 0; column_i < 6; ++column_i)
            {
                // Skipped, not multiplied, so an infinite irregular column at the center cannot give 0 * inf.
                if (column_weights[column_i] == std::complex<double>(0.0, 0.0)) continue;
                value += fundamental[row_i * 6 + column_i] * column_weights[column_i];
            }
            y_svc[row_i] = value;
        }

        // SVC16 to TS72 convention (B13 Eq. 7): swap y2 and y3, negate y5 and y6.
        out6[0] = y_svc[0];
        out6[1] = y_svc[2];
        out6[2] = y_svc[1];
        out6[3] = y_svc[3];
        out6[4] = -y_svc[4];
        out6[5] = -y_svc[5];
        return true;
    }

    // y1..y6 at n SI radii for one ytype; out holds n * C_MAX_NUM_Y complex values laid out [radius][y].
    void get_radial_solution_array(
        const double* radii_si, size_t n, size_t ytype_i, std::complex<double>* out) const
    {
        for (size_t i = 0; i < n; ++i)
            this->get_radial_solution(radii_si[i], ytype_i, &out[i * C_MAX_NUM_Y]);
    }

    // SI surface y1..y6 for one ytype from the find_love cache.
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

    // The layer holding an SI radius, and that radius in the interpolant's solve-unit domain. Shared by the
    // EOS readers below so they always agree on which layer a radius belongs to.
    bool p_locate_eos(double radius_si, size_t& layer_out, double& solve_radius_out) const
    {
        if (!this->eos_solution_uptr || !this->success) return false;
        solve_radius_out = radius_si / this->p_length_conv;
        if (!(solve_radius_out >= 0.0)) return false;   // also rejects NaN
        // Upper radii ascend, so the first layer whose top reaches the query holds it. Past the last top,
        // give or take a rounding step, is outside the body.
        for (size_t layer_i = 0; layer_i < this->num_layers; ++layer_i)
        {
            const double upper = (layer_i < this->p_upper_radii_solve.size())
                ? this->p_upper_radii_solve[layer_i]
                : TidalPyConstants::d_INF;
            if (solve_radius_out <= upper * (1.0 + 1.0e-12) + 1.0e-300)
            {
                layer_out = layer_i;
                return true;
            }
        }
        return false;
    }

    // Dense EOS evaluation at an SI radius, in the frequency-independent evaluation layout of
    // eos_layout_.hpp: [0] gravity [1] pressure [2] mass [3] moi [4] density [5] shear modulus
    // [6] bulk modulus [7,8] viscosities [9] temperature [10] heat flow [11] melt fraction.
    bool get_eos_si(double radius_si, double* out) const
    {
        size_t target_layer_i = 0;
        double solve_r        = 0.0;
        if (!this->p_locate_eos(radius_si, target_layer_i, solve_r)) return false;
        this->eos_solution_uptr->call_nondim(target_layer_i, solve_r, out);
        return true;
    }

    // Fill the SI grid the array-returning API reports (`result`) at the given radii [m], ascending, and keep the
    // radii so a plot can pair the two. A radius listed twice (an interface, as a layered radius array gives it)
    // takes the lower layer at its first copy and the upper layer at its second. Shooting path only: the matrix
    // path's grid is its solution.
    void sample_onto_radii(const double* radius_si, size_t n)
    {
        if (!this->success || !this->p_uses_interpolants || this->p_love_only) return;
        this->num_slices = n;
        this->total_size = static_cast<size_t>(C_MAX_NUM_Y_REAL) * n * this->num_ytypes;
        this->full_solution_vec.assign(this->total_size, TidalPyConstants::d_NAN);
        this->p_sample_radius_si.assign(radius_si, radius_si + n);
        const size_t num_output_ys = C_MAX_NUM_Y_REAL * this->num_ytypes;
        std::complex<double> out6[C_MAX_NUM_Y];
        c_RadialBasis basis;
        for (size_t slice_i = 0; slice_i < n; ++slice_i)
        {
            const bool upper_copy = (slice_i > 0) && (radius_si[slice_i] == radius_si[slice_i - 1]);
            const double radius_solve = radius_si[slice_i] / this->p_length_conv;
            const bool basis_found = this->p_evaluate_basis(radius_solve, upper_copy, basis);
            for (size_t ytype_i = 0; ytype_i < this->num_ytypes; ++ytype_i)
            {
                if (this->p_collapse_basis(basis, basis_found, ytype_i, radius_solve, out6))
                    this->apply_redimensionalization(out6);
                for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i)
                {
                    const size_t base = slice_i * num_output_ys + ytype_i * C_MAX_NUM_Y_REAL + y_i * 2;
                    this->full_solution_vec[base]     = out6[y_i].real();
                    this->full_solution_vec[base + 1] = out6[y_i].imag();
                }
            }
        }
    }

    // The same over the solution's own EOS radius grid, for a solution with no caller-supplied radii (a world's
    // released solution).
    void sample_onto_grid()
    {
        if (!this->success || !this->p_uses_interpolants || this->p_love_only) return;
        const c_EOSSolution* eos = this->eos_solution_uptr.get();
        const std::vector<double>& rad = eos->radius_array_vec;
        const size_t n = std::min(this->num_slices, rad.size());
        const double length_conv = this->p_eos_is_nondim ? this->p_length_conv : 1.0;  // EOS-units radius -> SI
        std::vector<double> radius_si(n);
        for (size_t slice_i = 0; slice_i < n; ++slice_i) { radius_si[slice_i] = rad[slice_i] * length_conv; }
        this->sample_onto_radii(radius_si.data(), n);
    }

    // The radii [m] the `result` grid was sampled on; empty when it holds the matrix path's own grid.
    const std::vector<double>& get_sample_radii_si() const noexcept { return this->p_sample_radius_si; }

    // Whether the solve kept only what its Love numbers need (see p_love_only): no radial functions below the surface.
    bool get_love_only() const noexcept { return this->p_love_only; }
};
