// world_radial_solver_.hpp: cached whole-planet Love-number solver owned by c_LayeredWorld.
//
// build_cache stores the frequency-independent setup once per EOS solve (non-dimensional scales and structure
// arrays, the c_ShootingInputs / c_MatrixInputs structs, the reused c_RadialSolutionStorage); solve then only
// non-dimensionalizes the complex moduli the world filled for the forcing frequency, runs the shooting or
// propagation-matrix method, and extracts the Love numbers. The helper holds no pointer back to the world, which
// keeps this header free of structures_x includes.
#pragma once

#include <cmath>
#include <complex>
#include <memory>
#include <string>
#include <vector>

#include "constants_.hpp"                                   // TidalPyConstants
#include "../Utilities_x/math_x/numerics_.hpp"                  // c_isclose
#include "../Utilities_x/arrays/layer_partition_.hpp"           // c_partition_radius_by_layer
#include "../Utilities_x/dimensions/nondimensional_.hpp"    // c_NonDimensionalScales
#include "rs_constants_.hpp"
#include "rs_solution_.hpp"
#include "love_.hpp"
#include "shooting_.hpp"
#include "matrix_.hpp"


// =====================================================================================================================
// Per-solver input structs, built once by build_cache; each solve updates only the per-call knobs
// =====================================================================================================================

struct c_ShootingInputs {
    // layer_types: 0 = solid, 1 = liquid. bool[] because std::vector<bool> is bit-packed and the solver wants bool*.
    std::vector<int>        layer_types;
    std::unique_ptr<bool[]> is_static;
    std::unique_ptr<bool[]> is_incompressible;
    size_t                  num_layers = 0;

    // Per-layer slice partitioning over the (non-dim) radius grid.
    std::vector<size_t> first_slice_index_by_layer;
    std::vector<size_t> num_slices_by_layer;

    // Surface boundary conditions to solve for, in order (tidal = 1, free = 0, loading = 2). One solve produces a
    // block of radial functions per entry. The independent solutions do not depend on the boundary condition, so
    // n conditions cost one integration and n surface solves, not n integrations.
    std::vector<int> bc_models = {1};

    // Non-dim planet scalars.
    double planet_bulk_density = 0.0;
    double G                   = 0.0;
    int    degree_l            = 2;

    // Shooting-method knobs (per-call, set from the runtime config).
    bool      use_kamata          = false;
    double    starting_radius     = 0.0;          // non-dim
    double    start_radius_tol    = 1.0e-4;
    ODEMethod integration_method  = ODEMethod::DOP853;
    double    integration_rtol    = 1.0e-5;
    double    integration_atol    = 1.0e-7;
    bool      scale_rtols         = true;
    size_t    max_num_steps       = 500000;
    size_t    expected_size       = 500;
    size_t    max_ram_MB          = 500;
    double    max_step            = 0.0;
    bool      warnings            = true;
};

// Propagation-matrix-method inputs (only valid for a single solid, static, incompressible layer).
struct c_MatrixInputs {
    size_t num_layers          = 1;
    // The method lays its own grid down inside c_matrix_propagate; this is the only thing left to say about it.
    size_t slices_per_layer    = 0;
    std::vector<int> bc_models = {1};
    double planet_bulk_density = 0.0;
    double G                   = 0.0;
    int    degree_l            = 2;
    double starting_radius     = 0.0;   // non-dim
    double start_radius_tol    = 1.0e-4;
    int    core_model          = 0;
};

// Struct-based wrappers around the positional solvers.
inline int c_shooting_solve(
        c_RadialSolutionStorage* storage, c_ShootingInputs& in, double frequency, bool verbose) noexcept {
    return c_shooting_solver(
            storage,
            frequency,
            in.planet_bulk_density,
            in.layer_types.data(),
            in.is_static.get(),
            in.is_incompressible.get(),
            in.first_slice_index_by_layer,
            in.num_slices_by_layer,
            in.bc_models.size(),
            in.bc_models.data(),
            in.G,
            in.degree_l,
            in.use_kamata,
            in.starting_radius,
            in.start_radius_tol,
            in.integration_method,
            in.integration_rtol,
            in.integration_atol,
            in.scale_rtols,
            in.max_num_steps,
            in.expected_size,
            in.max_ram_MB,
            in.max_step,
            verbose,
            in.warnings);
}

inline int c_matrix_solve(
        c_RadialSolutionStorage* storage, c_MatrixInputs& in, double frequency, bool verbose) noexcept {
    return c_matrix_propagate(
        storage,
        frequency,
        in.planet_bulk_density,
        in.slices_per_layer,
        in.bc_models.size(),
        in.bc_models.data(),
        in.G,
        in.degree_l,
        in.starting_radius,
        in.start_radius_tol,
        in.core_model,
        verbose
    );
}


// =====================================================================================================================
// Runtime (per-solve) configuration
// =====================================================================================================================
struct c_LoveSolveRuntimeConfig {
    double    frequency       = 1.0e-5;             // [rad/s]; the only physically per-call quantity
    std::vector<int> bc_models = {1};               // tidal = 1, free = 0, loading = 2; one output block each
    bool      use_prop_matrix = false;              // false = shooting method, true = propagation matrix
    int       core_model      = 0;                  // propagation-matrix core starting condition (0-4)
    bool      use_kamata      = false;
    double    starting_radius = 0.0;                // [m]; 0 -> auto
    double    start_radius_tol = 1.0e-4;
    ODEMethod integration_method = ODEMethod::DOP853;
    double    rtol               = 1.0e-5;
    double    atol               = 1.0e-7;
    bool      scale_rtols        = true;
    size_t    max_num_steps      = 500000;
    size_t    expected_size      = 500;
    size_t    max_ram_MB         = 500;
    double    max_step           = 0.0;
    bool      verbose            = false;
    bool      warnings           = true;              // enables the surface conditioning diagnostic
    bool      redim_eos_arrays   = false;             // export mode: also redimensionalize the EOS arrays
};


// =====================================================================================================================
// c_WorldRadialSolver
// =====================================================================================================================
class c_WorldRadialSolver {
public:
    c_WorldRadialSolver() = default;
    ~c_WorldRadialSolver() = default;

    // Cache signature check (layer count, slice count, degree, nondim flag, output-block count) so the world can
    // skip a rebuild. The block count is in the signature because it sizes the storage, so asking for a different
    // set of boundary conditions has to rebuild.
    bool cache_matches(
        size_t n_layers,
        size_t total_slices,
        int degree_l,
        bool nondimensionalize,
        size_t num_ytypes) const noexcept
    {
        return this->p_cache_valid
            && this->p_n_layers     == n_layers
            && this->p_total_slices == total_slices
            && this->p_degree_l     == degree_l
            && this->p_nondim       == nondimensionalize
            && this->p_num_ytypes   == num_ytypes;
    }

    // The layer flags are user-mutable without an EOS re-solve, so a cache hit must also confirm them.
    bool layer_flags_match(
        const int* layer_types,
        const bool* is_static,
        const bool* is_incompressible,
        size_t n_layers) const noexcept
    {
        const c_ShootingInputs& shoot = this->p_shooting_inputs;
        if (shoot.layer_types.size() != n_layers) { return false; }
        if (!shoot.is_static || !shoot.is_incompressible) { return false; }
        for (size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
            if (shoot.layer_types[layer_i]       != layer_types[layer_i])       { return false; }
            if (shoot.is_static[layer_i]         != is_static[layer_i])         { return false; }
            if (shoot.is_incompressible[layer_i] != is_incompressible[layer_i]) { return false; }
        }
        return true;
    }

    void invalidate() noexcept { this->p_cache_valid = false; }

    c_RadialSolutionStorage* get_storage() const noexcept { return this->p_storage.get(); }

    bool get_solved() const noexcept { return this->p_solved; }

    std::unique_ptr<c_RadialSolutionStorage> release_storage() noexcept {
        this->p_cache_valid = false;
        return std::move(this->p_storage);
    }

    // Install the per-layer state provider that both methods read gravity, density, and the complex moduli from, at
    // the exact radius asked for. Set per Love solve, because the callable carries that solve's forcing frequency,
    // and left in place afterwards so the solution can still be read. See c_EOSSolution::MaterialEval.
    void set_material_eval(c_EOSSolution::MaterialEval eval) {
        if (this->p_storage) {
            this->p_storage->get_eos_solution_ptr()->p_material_eval = std::move(eval);
        }
    }
    size_t total_slices() const noexcept { return this->p_total_slices; }

    // Frequency-independent setup from SI inputs. Returns false (with an error on the storage) if the layer slice
    // partition is invalid.
    bool build_cache(
        const std::vector<double>& radius_si,
        const std::vector<double>& density_si,
        const std::vector<double>& gravity_si,
        const std::vector<double>& pressure_si,
        const std::vector<double>& mass_si,
        const std::vector<double>& moi_si,
        const std::vector<double>& upper_radii_si,
        const int* layer_types,
        const bool* is_static,
        const bool* is_incompressible,
        size_t n_layers,
        double planet_radius,
        double bulk_density,
        int degree_l,
        bool nondimensionalize,
        const c_EOSSolution* world_eos_ptr = nullptr,
        size_t num_ytypes = 1)
    {
        const size_t total_slices = radius_si.size();
        if (num_ytypes == 0) { num_ytypes = 1; }

        this->p_n_layers           = n_layers;
        this->p_total_slices       = total_slices;
        this->p_degree_l           = degree_l;
        this->p_nondim             = nondimensionalize;
        this->p_num_ytypes         = num_ytypes;
        this->p_surface_gravity_si = gravity_si[total_slices - 1];

        // unique_ptr because c_NonDimensionalScales is not assignable.
        this->p_non_dim_uptr = std::make_unique<c_NonDimensionalScales>(planet_radius, bulk_density);

        const double length_conv  = this->p_non_dim_uptr->length_conversion;
        const double mass_conv     = this->p_non_dim_uptr->mass_conversion;
        const double sec2_conv     = this->p_non_dim_uptr->second2_conversion;
        const double pascal_conv   = this->p_non_dim_uptr->pascal_conversion;
        const double density_conv  = this->p_non_dim_uptr->density_conversion;
        const double gravity_conv  = length_conv / sec2_conv;
        const double moi_conv       = mass_conv * length_conv * length_conv;
        const double G_si           = c_get_G();

        this->p_upper_radii_nd.assign(upper_radii_si.begin(), upper_radii_si.end());
        double G_nd   = G_si;
        double rho_nd = bulk_density;
        if (nondimensionalize) {
            for (size_t layer_i = 0; layer_i < n_layers; ++layer_i)
                this->p_upper_radii_nd[layer_i] /= length_conv;
            G_nd   = G_si / (this->p_non_dim_uptr->length3_conversion / (mass_conv * sec2_conv));
            rho_nd = bulk_density / density_conv;
        }

        // The non-dim radius grid is built here and handed to the storage, which keeps it only to size its output
        // sampling. No other structure array is copied: everything the solve reads comes from the world's solved
        // EOS and the layer's models, through the dense source and provider installed below.
        this->p_radius_nd = radius_si;
        if (nondimensionalize) {
            for (size_t slice_i = 0; slice_i < total_slices; ++slice_i) {
                this->p_radius_nd[slice_i] /= length_conv;
            }
        }

        // (Re)build the reusable solution storage with the non-dim radius grid.
        this->p_storage = std::make_unique<c_RadialSolutionStorage>(
            num_ytypes,
            this->p_upper_radii_nd.data(),
            n_layers,
            this->p_radius_nd.data(),
            total_slices,
            degree_l);

        // Per-layer slice partitioning over the non-dim grid (interface radii appear in two layers).
        std::vector<size_t> first_slice_idx;
        std::vector<size_t> num_slices;
        tidalpy::c_partition_radius_by_layer(
            this->p_radius_nd.data(),
            total_slices,
            this->p_upper_radii_nd.data(),
            n_layers,
            first_slice_idx,
            num_slices);
        for (size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
            if (num_slices[layer_i] < 5) {
                this->p_storage->error_code = -5;
                this->p_storage->message    = "TidalPy: at least 5 slices per layer required";
                this->p_storage->success    = false;
                this->p_cache_valid         = false;
                return false;
            }
        }

        // Populate the shooting-method inputs (structural fields; per-call knobs set in solve()).
        c_ShootingInputs& shoot = this->p_shooting_inputs;
        shoot.num_layers          = n_layers;
        shoot.layer_types.assign(layer_types, layer_types + n_layers);
        shoot.is_static           = std::make_unique<bool[]>(n_layers);
        shoot.is_incompressible   = std::make_unique<bool[]>(n_layers);
        for (size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
            shoot.is_static[layer_i]         = is_static[layer_i];
            shoot.is_incompressible[layer_i] = is_incompressible[layer_i];
        }
        shoot.first_slice_index_by_layer = first_slice_idx;
        shoot.num_slices_by_layer        = num_slices;
        shoot.planet_bulk_density        = rho_nd;
        shoot.G                          = G_nd;
        shoot.degree_l                   = degree_l;

        // Populate the propagation-matrix inputs (structural fields). It builds its own grid, so all it needs is how
        // fine that grid should be; the count matches what the caller asked for so the output sampling is unchanged.
        c_MatrixInputs& mat = this->p_matrix_inputs;
        mat.num_layers          = n_layers;
        mat.slices_per_layer    = (n_layers > 0) ? total_slices / n_layers : 0;
        mat.planet_bulk_density = rho_nd;
        mat.G                   = G_nd;
        mat.degree_l            = degree_l;

        // The solve reads no stored profile, but the solution still reports the planet's scalars, and find_love
        // runs only on a solution marked solved. These are the same values inject_from_world_eos used to take off
        // the ends of the arrays it copied, in the units the methods work in.
        c_EOSSolution* storage_eos = this->p_storage->get_eos_solution_ptr();
        const auto to_nd = [nondimensionalize](double value, double conv) {
            return nondimensionalize ? value / conv : value;
        };
        storage_eos->radius                 = to_nd(radius_si[total_slices - 1],   length_conv);
        storage_eos->surface_gravity        = to_nd(gravity_si[total_slices - 1],  gravity_conv);
        storage_eos->surface_pressure       = to_nd(pressure_si[total_slices - 1], pascal_conv);
        storage_eos->central_pressure       = to_nd(pressure_si[0],                pascal_conv);
        storage_eos->mass                   = to_nd(mass_si[total_slices - 1],     mass_conv);
        storage_eos->moi                    = to_nd(moi_si[total_slices - 1],      moi_conv);
        storage_eos->nondim_status          = 0;
        storage_eos->solution_nondim_status = 0;
        storage_eos->success                = true;
        storage_eos->error_code             = 0;
        storage_eos->radius_array_set       = true;
        storage_eos->other_vecs_set         = true;

        // This storage's EOS stands in for the world's, so it also reports the world's solve diagnostics. Without
        // this they stay at their "never solved" defaults, which an exported solution would report as its own.
        if (world_eos_ptr) {
            storage_eos->iterations         = world_eos_ptr->iterations;
            storage_eos->pressure_error     = world_eos_ptr->pressure_error;
            storage_eos->max_iters_hit      = world_eos_ptr->max_iters_hit;
            storage_eos->message            = world_eos_ptr->message;
            storage_eos->steps_taken_vec    = world_eos_ptr->steps_taken_vec;
            storage_eos->num_cysolver_calls = world_eos_ptr->num_cysolver_calls;
        }

        // The state provider answers in SI at an SI radius; the scales convert the non-dim shooting radius up and
        // its answers back down.
        storage_eos->p_structure_length_scale  = nondimensionalize ? length_conv  : 1.0;
        storage_eos->p_structure_gravity_scale = nondimensionalize ? gravity_conv : 1.0;
        storage_eos->p_structure_pascal_scale  = nondimensionalize ? pascal_conv  : 1.0;
        storage_eos->p_structure_density_scale = nondimensionalize ? density_conv : 1.0;

        this->p_cache_valid = true;
        return true;
    }

    // Frequency-dependent step; expects the SI complex-moduli scratch to be filled for rt.frequency.
    void solve(const c_LoveSolveRuntimeConfig& rt) {
        c_RadialSolutionStorage* storage = this->p_storage.get();
        storage->success    = false;
        storage->error_code = 0;
        storage->p_bc_models = rt.bc_models;
        storage->p_love_frequency_si = rt.frequency;
        this->p_solved      = false;

        // Propagation matrix is only valid for a single solid, static, incompressible layer.
        if (rt.use_prop_matrix) {
            const bool ok = (this->p_n_layers == 1)
                && (this->p_shooting_inputs.layer_types[0] == 0)
                && this->p_shooting_inputs.is_static[0]
                && this->p_shooting_inputs.is_incompressible[0];
            if (!ok) {
                storage->error_code = -20;
                storage->message    = "TidalPy: the propagation-matrix method requires a single solid, static, "
                                      "incompressible layer.";
                storage->success    = false;
                this->p_solved      = false;
                return;
            }
        }

        const double freq_nd = this->p_nondim
            ? rt.frequency * this->p_non_dim_uptr->second_conversion
            : rt.frequency;
        const double start_r = this->p_nondim
            ? rt.starting_radius / this->p_non_dim_uptr->length_conversion
            : rt.starting_radius;

        if (rt.use_prop_matrix) {
            c_MatrixInputs& mat = this->p_matrix_inputs;
            mat.bc_models        = rt.bc_models;
            mat.starting_radius  = start_r;
            mat.start_radius_tol = rt.start_radius_tol;
            mat.core_model       = rt.core_model;
            c_matrix_solve(storage, mat, freq_nd, rt.verbose);
        } else {
            c_ShootingInputs& shoot = this->p_shooting_inputs;
            shoot.bc_models          = rt.bc_models;
            shoot.use_kamata         = rt.use_kamata;
            shoot.starting_radius    = start_r;
            shoot.start_radius_tol   = rt.start_radius_tol;
            shoot.integration_method = rt.integration_method;
            shoot.integration_rtol   = rt.rtol;
            shoot.integration_atol   = rt.atol;
            shoot.scale_rtols        = rt.scale_rtols;
            shoot.max_num_steps      = rt.max_num_steps;
            shoot.expected_size      = rt.expected_size;
            shoot.max_ram_MB         = rt.max_ram_MB;
            shoot.max_step           = rt.max_step;
            shoot.warnings           = rt.warnings;
            c_shooting_solve(storage, shoot, freq_nd, rt.verbose);
        }

        // Dimensional context for get_radial_solution; the EOS arrays are still non-dim here.
        if (storage->success) {
            if (this->p_nondim) {
                const double length_conv  = this->p_non_dim_uptr->length_conversion;
                const double sec2_conv    = this->p_non_dim_uptr->second2_conversion;
                const double disp_scale   = sec2_conv / length_conv;
                const double stress_scale =
                    this->p_non_dim_uptr->mass_conversion / this->p_non_dim_uptr->length3_conversion;
                const double pot_scale    = 1.0 / length_conv;
                const double grav_conv    = length_conv / sec2_conv;
                const double dens_conv    = this->p_non_dim_uptr->density_conversion;
                storage->set_dimensional_context(
                    length_conv,
                    disp_scale,
                    stress_scale,
                    pot_scale,
                    /*eos_is_nondim=*/true,
                    grav_conv,
                    dens_conv);
            } else {
                storage->set_dimensional_context(
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    /*eos_is_nondim=*/false,
                    1.0,
                    1.0);
            }
        }

        // Love numbers from the non-dim solution: k = y5 - 1, h = y1 g, l = y3 g, and the displacement and gravity
        // scales cancel. The storage's surface_gravity stays non-dim for the next solve.
        if (storage->success)
            storage->find_love();

        // Export mode fills the SI grid for the array-returning standalone API while the EOS is still non-dim.
        if (rt.redim_eos_arrays && storage->success && storage->p_uses_interpolants)
            storage->sample_onto_grid();

        // Re-dimensionalize the matrix y-grid and, in export mode only, the EOS arrays; the fast path keeps them
        // non-dim so the cache survives for the next frequency.
        if (this->p_nondim && storage->success) {
            storage->dimensionalize_data(
                this->p_non_dim_uptr.get(),
                true,
                /*include_eos=*/rt.redim_eos_arrays
            );

            if (rt.redim_eos_arrays) {
                storage->get_eos_solution_ptr()->surface_gravity = this->p_surface_gravity_si;
                storage->p_eos_is_nondim = false;
            }
        }
        // Export mode consumes the cache (the EOS arrays are now SI).
        if (rt.redim_eos_arrays)
            this->p_cache_valid = false;

        this->p_solved = storage->success;
    }

    // Cached state (frequency-independent unless noted).
    bool   p_cache_valid  = false;
    size_t p_n_layers     = 0;
    size_t p_total_slices = 0;
    int    p_degree_l     = 0;
    bool   p_nondim       = true;
    bool   p_solved       = false;
    size_t p_num_ytypes   = 1;
    double p_surface_gravity_si = TidalPyConstants::d_NAN;

    std::unique_ptr<c_NonDimensionalScales> p_non_dim_uptr;

    c_ShootingInputs p_shooting_inputs;
    c_MatrixInputs   p_matrix_inputs;

    std::vector<double> p_upper_radii_nd;
    // The non-dim radius grid, kept only to size and place the storage's output sampling.
    std::vector<double> p_radius_nd;

    std::unique_ptr<c_RadialSolutionStorage> p_storage;
};
