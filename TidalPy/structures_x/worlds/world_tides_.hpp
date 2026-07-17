#pragma once
/*
 * world_tides_.hpp — out-of-line definition of c_LayeredWorld::calc_tides.
 *
 * c_LayeredWorld extends the common analytic tide path (c_BaseWorld::calc_tides,
 * world_tides_base_.hpp) with two layered-world-only capabilities:
 *   - the rheology model: -Im[k_l(omega)] from the world radial solver run at each unique
 *     tidal frequency (the EOS must be solved first); and
 *   - per-layer heating distribution by each layer's tidal_scale_method.
 *
 * It runs the global-potential engine, collapses (analytic or rheology), stores the result,
 * then distributes the heat to the layers. The tide-model holder / config / result state is
 * inherited from c_BaseWorld.
 *
 * This header pulls in the heavy global-potential tables; force-include it in the layered
 * (and gas-giant) world extension only.
 */

#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "layered_.hpp"
#include "../../Tides_x/classes/tide_collapse_.hpp"      // c_global_potential, c_collapse_global_tides
#include "../../Tides_x/classes/tide_.hpp"               // c_RheologyTide (3D orchestration target)
#include "../../Tides_x/potential/potential_3d_.hpp"     // c_tidal_potential_3d_modes (dynamic Kaula engine)
#include "../../Tides_x/multilayer/kernel_.hpp"          // strain/stress/heating kernel + c_StrainRadialCoeffs

namespace tidalpy {

inline void c_LayeredWorld::calc_tides(const c_TideSolveConfig& state) {
    if (!this->p_tide) {
        throw std::runtime_error(
            "TidalPy: no tide model attached to the world — call set_tide_model() first");
    }

    const double planet_radius = this->get_radius();
    const double G_to_use = (tidalpy_config_ptr != nullptr) ? tidalpy_config_ptr->d_G : 6.674015e-11;
    const c_TideConfig& tcfg = this->p_tide_config;

    // Model-independent per-mode terms + the unique-frequency maps.
    c_GlobalPotentialStorage potential = c_global_potential(
        planet_radius,
        state.semi_major_axis,
        state.orbital_frequency,
        state.spin_frequency,
        state.obliquity,
        state.eccentricity,
        state.host_mass,
        G_to_use,
        tcfg.min_degree_l,
        tcfg.max_degree_l,
        tcfg.obliquity_truncation,
        tcfg.eccentricity_truncation
    );

    if (potential.error_code != 0) {
        this->p_tides_solved           = false;
        this->p_tide_result            = c_GlobalTideResult();
        this->p_tide_result.error_code = potential.error_code;
        throw std::runtime_error("TidalPy: global potential failed during calc_tides");
    }

    // Collapse the per-mode potential terms with the tide model's dissipation multiplier.
    this->p_tide_solver_love.clear();
    if (this->p_tide->needs_radial_solve()) {
        // Rheology model: the per-mode -Im[k_l(omega)] comes from the world radial solver,
        // which needs the EOS solved first.
        if (!this->p_eos_solved || !this->p_eos_solution) {
            this->p_tides_solved = false;
            throw std::runtime_error(
                "TidalPy: the rheology tide model needs the EOS solved first — call "
                "solve_eos() before calc_tides()");
        }

        // Solve the Love numbers once per unique (degree_l, frequency) pair, caching by that
        // pair so modes that share a degree and frequency reuse one radial solve. Then record
        // each active mode's Love numbers keyed by its (l, m, p, q).
        c_IntMap<c_Key2, tidalpy::c_LoveNumbers> love_by_l_freq;
        c_LoveSolveConfig love_cfg;
        for (const auto& mode_entry : potential.potential_map) {
            const c_Key4& lmpq_key = mode_entry.first;
            const int degree_l     = static_cast<int>(lmpq_key.a);

            bool found = false;
            const std::size_t freq_index = potential.unique_freq_index_map.get(found, lmpq_key);
            if (!found) {
                // Inactive (zero-frequency) mode; contributes nothing.
                continue;
            }
            const double frequency = potential.unique_freq_map[freq_index].frequency;

            c_Key2 lf_key(static_cast<int16_t>(degree_l), static_cast<int16_t>(freq_index));
            bool cached = false;
            love_by_l_freq.get(cached, lf_key);
            if (!cached) {
                love_cfg.degree_l        = degree_l;
                love_cfg.frequency_rad_s = frequency;
                this->solve_love_numbers(love_cfg);
                if (!this->get_love_success()) {
                    this->p_tides_solved = false;
                    throw std::runtime_error(
                        "TidalPy: radial-solver Love-number solve failed during calc_tides: "
                        + this->get_love_message());
                }
                tidalpy::c_LoveNumbers solved_love;
                solved_love.k = this->get_love_number_k(0);
                solved_love.h = this->get_love_number_h(0);
                solved_love.l = this->get_love_number_l(0);
                love_by_l_freq.set(lf_key, solved_love);
            }

            bool have = false;
            tidalpy::c_LoveNumbers mode_love = love_by_l_freq.get(have, lf_key);
            this->p_tide_solver_love.set(lmpq_key, mode_love);
        }

        this->p_tide_result = c_collapse_global_tides(potential, *this->p_tide, &this->p_tide_solver_love);
    } else {
        // Analytic models: no radial-solver Love numbers needed.
        this->p_tide_result = c_collapse_global_tides(potential, *this->p_tide, nullptr);
    }
    this->p_tides_solved = true;

    // Distribute heating to the layers by their effective tidal scale (0 for non-tidal
    // layers), and store the result on each layer so layer.get_tidal_heating() reports it.
    const double planet_volume =
        (4.0 / 3.0) * TidalPyConstants::d_PI * planet_radius * planet_radius * planet_radius;
    const std::size_t n_layers = this->p_layers.size();
    this->p_layer_tidal_heating.assign(n_layers, 0.0);
    for (std::size_t i = 0; i < n_layers; ++i) {
        c_BaseLayer* layer = this->p_layers[i].get();
        const double scale = this->effective_tidal_scale(layer, planet_volume, state);
        const double heat  = this->p_tide_result.tidal_heating * scale;
        this->p_layer_tidal_heating[i] = heat;
        layer->set_tidal_heating(heat);
    }
}

// Effective per-layer tidal-heating scale for the layer's tidal_scale_method (0 for a
// non-tidal layer).
//   user_provided   : the layer's tidal_scale field.
//   volume_fraction : layer volume / planet volume.
//   tidal_timescale : a log-Gaussian bell curve in the layer's Maxwell time tau = eta/mu
//                     (from its static shear modulus + viscosity) about the tidal forcing
//                     period 2*pi/|orbital_frequency|; width [decades] from the tide config.
//                     Returns 0 for a geometry-only layer or when mu/eta/forcing are unusable.
inline double c_LayeredWorld::effective_tidal_scale(
        const c_BaseLayer* layer, double planet_volume, const c_TideSolveConfig& state) const {
    if (!layer->get_is_tidal()) {
        return 0.0;
    }
    switch (layer->get_tidal_scale_method()) {
        case c_TidalScaleMethod::user_provided:
            return layer->get_tidal_scale();
        case c_TidalScaleMethod::volume_fraction:
            return (planet_volume > TidalPyConstants::d_EPS)
                 ? layer->get_volume() / planet_volume : 0.0;
        case c_TidalScaleMethod::tidal_timescale: {
            const auto* phys = dynamic_cast<const c_PhysicsLayer*>(layer);
            if (phys == nullptr) {
                return 0.0;   // geometry-only layer has no Maxwell time
            }
            const double shear_modulus  = phys->get_shear_modulus_static();
            const double shear_viscosity = phys->get_shear_viscosity_static();
            const double orbital_freq   = std::abs(state.orbital_frequency);
            if (shear_modulus  <= TidalPyConstants::d_EPS
             || shear_viscosity <= TidalPyConstants::d_EPS
             || orbital_freq    <= TidalPyConstants::d_EPS) {
                return 0.0;
            }
            const double maxwell_time   = shear_viscosity / shear_modulus;          // [s]
            const double forcing_period = 2.0 * TidalPyConstants::d_PI / orbital_freq;  // [s]
            double width = this->p_tide_config.tidal_timescale_width_decades;
            if (width <= TidalPyConstants::d_EPS) {
                width = 1.0;
            }
            const double z = std::log10(maxwell_time / forcing_period) / width;
            return std::exp(-0.5 * z * z);
        }
        default:
            return 0.0;
    }
}


// =====================================================================================================================
// On-demand 3D tidal heating
// =====================================================================================================================

// The 3D orchestration lives on the rheology tide model (the only TideBase with a depth-resolved
// solution); it calls the world's members directly (no callbacks). Defined here, in the world extension,
// where c_LayeredWorld + the kernel/potential headers are complete and CyRK lives, so every radial solve
// + dense call stays in its owning extension.
inline double c_RheologyTide::calc_3d_tidal_heating(
        c_LayeredWorld& world,
        const c_TideSolveConfig& state,
        double radius,
        double colatitude) const {
    const c_TideConfig& tide_cfg = world.get_tide_config();
    const double surface_radius = world.get_radius();
    const double G_to_use =
        (tidalpy_config_ptr != nullptr) ? tidalpy_config_ptr->d_G : 6.67430e-11;

    int engine_error = 0;
    const std::vector<c_TidalPotential3DMode> modes = c_tidal_potential_3d_modes(
        surface_radius, state.semi_major_axis, state.orbital_frequency, state.spin_frequency,
        state.obliquity, state.eccentricity, state.host_mass, G_to_use,
        tide_cfg.min_degree_l, tide_cfg.max_degree_l,
        tide_cfg.obliquity_truncation, tide_cfg.eccentricity_truncation,
        colatitude, 0.0, &engine_error);
    if (engine_error != 0) {
        throw std::runtime_error(
            "TidalPy: tidal potential engine failed during secular 3D tidal heating (error "
            + std::to_string(engine_error) + "); check degree/truncation levels");
    }

    bool is_solid = true;
    bool is_incompressible = false;
    const auto* physics_layer = dynamic_cast<const c_PhysicsLayer*>(world.find_layer_for_radius(radius));
    if (physics_layer != nullptr) {
        is_solid = physics_layer->get_is_solid();
        is_incompressible = physics_layer->get_is_incompressible();
    }

    const double min_freq =
        (tidalpy_config_ptr != nullptr) ? tidalpy_config_ptr->d_MIN_SPIN_ORBIT_DIFF : 1.0e-9;

    double heating = 0.0;
    c_LoveSolveConfig love_cfg;
    for (const c_TidalPotential3DMode& mode : modes) {
        const double frequency = std::abs(mode.mode_frequency);
        if (frequency <= min_freq) {
            continue;
        }

        love_cfg.degree_l = mode.degree_l;
        love_cfg.frequency_rad_s = frequency;
        world.solve_love_numbers(love_cfg);
        if (!world.get_love_success()) {
            throw std::runtime_error(
                "TidalPy: radial solve failed during secular 3D tidal heating: " + world.get_love_message());
        }

        const ::c_RadialSolutionStorage* storage = world.get_love_storage();
        std::complex<double> y_at_r[C_MAX_NUM_Y];
        if (storage == nullptr || !storage->get_radial_solution(radius, 0, y_at_r)) {
            return TidalPyConstants::d_NAN;
        }
        if (!std::isfinite(y_at_r[0].real()) || !std::isfinite(y_at_r[1].real())
         || !std::isfinite(y_at_r[2].real()) || !std::isfinite(y_at_r[3].real())) {
            return TidalPyConstants::d_NAN;
        }

        const std::complex<double> shear = world.calc_complex_shear_modulus(radius, frequency);
        const std::complex<double> bulk  = world.calc_complex_bulk_modulus(radius, frequency);

        const tides::c_StrainRadialCoeffs coeffs = tides::c_compute_strain_radial_coeffs(
            y_at_r[0], y_at_r[1], y_at_r[2], y_at_r[3], shear, bulk, radius,
            static_cast<double>(mode.degree_l), is_solid, is_incompressible);
        if (!coeffs.valid) {
            continue;   // liquid / center: no shear dissipation contribution
        }
        tides::c_Tensor6 strain;
        tides::c_Tensor6 stress;
        tides::c_compute_strain_stress(coeffs, mode.potential, colatitude, strain, stress);

        heating += 0.5 * frequency * tides::c_volumetric_heating_signed(stress, strain);
    }
    return heating;
}

// World delegation: validate preconditions, map the solve state into the potential's state struct, and
// hand off to the rheology tide model's 3D orchestration.
inline double c_LayeredWorld::get_3d_tidal_heating(
        const c_TideSolveConfig& state,
        double radius,
        double colatitude) {
    if (!this->p_tide) {
        throw std::runtime_error(
            "TidalPy: no tide model attached to the world — call set_tide_model() first");
    }
    auto* rheology = dynamic_cast<c_RheologyTide*>(this->p_tide.get());
    if (rheology == nullptr) {
        throw std::runtime_error(
            "TidalPy: 3D tidal heating requires the rheology tide model (the analytic cpl/ctl/ctl_q "
            "models have no depth-resolved radial solution)");
    }
    if (!this->p_eos_solved || !this->p_eos_solution) {
        throw std::runtime_error(
            "TidalPy: 3D tidal heating needs the EOS solved first — call solve_eos()");
    }
    return rheology->calc_3d_tidal_heating(*this, state, radius, colatitude);
}

// World delegation for the batch path: same preconditions as the scalar get_3d_tidal_heating.
inline void c_LayeredWorld::get_3d_tidal_heating_array(
        const c_TideSolveConfig& state,
        const double* radii,
        const double* colatitudes,
        size_t num_points,
        double* out_heating) {
    if (!this->p_tide) {
        throw std::runtime_error(
            "TidalPy: no tide model attached to the world — call set_tide_model() first");
    }
    auto* rheology = dynamic_cast<c_RheologyTide*>(this->p_tide.get());
    if (rheology == nullptr) {
        throw std::runtime_error(
            "TidalPy: 3D tidal heating requires the rheology tide model (the analytic cpl/ctl/ctl_q "
            "models have no depth-resolved radial solution)");
    }
    if (!this->p_eos_solved || !this->p_eos_solution) {
        throw std::runtime_error(
            "TidalPy: 3D tidal heating needs the EOS solved first — call solve_eos()");
    }
    rheology->calc_3d_tidal_heating_batch(*this, state, radii, colatitudes, num_points, out_heating);
}

// Vectorized batch form of the secular 3D heating. Builds the position-independent mode list once and
// amortizes the world radial (Love-number) solve across all query points: the radial solve depends on
// (l, frequency) only, so consecutive modes that share (degree_l, |frequency|) reuse the same solve.
// out_heating[i] accumulates each active mode's (freq/2) Im(sigma_c : conj(eps_c)) at point i; a point
// whose radius has no depth-resolved solution (liquid / center / below the solver start) is set to NaN.
inline void c_RheologyTide::calc_3d_tidal_heating_batch(
        c_LayeredWorld& world,
        const c_TideSolveConfig& state,
        const double* radii,
        const double* colatitudes,
        size_t num_points,
        double* out_heating) const {
    const c_TideConfig& tide_cfg = world.get_tide_config();
    const double surface_radius = world.get_radius();
    const double G_to_use =
        (tidalpy_config_ptr != nullptr) ? tidalpy_config_ptr->d_G : 6.67430e-11;

    int engine_error = 0;
    const std::vector<c_TidalPotential3DModeCoeff> modes = c_tidal_potential_3d_mode_coeffs(
        surface_radius, state.semi_major_axis, state.orbital_frequency, state.spin_frequency,
        state.obliquity, state.eccentricity, state.host_mass, G_to_use,
        tide_cfg.min_degree_l, tide_cfg.max_degree_l,
        tide_cfg.obliquity_truncation, tide_cfg.eccentricity_truncation, &engine_error);
    if (engine_error != 0) {
        throw std::runtime_error(
            "TidalPy: tidal potential engine failed during secular 3D tidal heating (error "
            + std::to_string(engine_error) + "); check degree/truncation levels");
    }

    for (size_t i = 0; i < num_points; ++i) {
        out_heating[i] = 0.0;
    }
    // A point whose radius has no depth-resolved solution is invalid for every mode; mark and skip it.
    std::vector<unsigned char> point_invalid(num_points, 0);

    const double min_freq =
        (tidalpy_config_ptr != nullptr) ? tidalpy_config_ptr->d_MIN_SPIN_ORBIT_DIFF : 1.0e-9;

    c_LoveSolveConfig love_cfg;
    int last_degree_l = -1;
    double last_frequency = -1.0;
    for (const c_TidalPotential3DModeCoeff& mode : modes) {
        const double frequency = std::abs(mode.mode_frequency);
        if (frequency <= min_freq) {
            continue;
        }

        // Solve the radial response only when (degree_l, frequency) changes from the previous mode
        // (the engine emits modes grouped by degree, so shared frequencies cluster).
        if (mode.degree_l != last_degree_l || frequency != last_frequency) {
            love_cfg.degree_l = mode.degree_l;
            love_cfg.frequency_rad_s = frequency;
            world.solve_love_numbers(love_cfg);
            if (!world.get_love_success()) {
                throw std::runtime_error(
                    "TidalPy: radial solve failed during secular 3D tidal heating: "
                    + world.get_love_message());
            }
            last_degree_l = mode.degree_l;
            last_frequency = frequency;
        }

        const ::c_RadialSolutionStorage* storage = world.get_love_storage();
        if (storage == nullptr) {
            throw std::runtime_error("TidalPy: missing radial solution during secular 3D tidal heating");
        }
        const double degree_l_d = static_cast<double>(mode.degree_l);

        for (size_t i = 0; i < num_points; ++i) {
            if (point_invalid[i]) {
                continue;
            }
            const double radius = radii[i];
            const double colatitude = colatitudes[i];

            std::complex<double> y_at_r[C_MAX_NUM_Y];
            if (!storage->get_radial_solution(radius, 0, y_at_r)
                || !std::isfinite(y_at_r[0].real()) || !std::isfinite(y_at_r[1].real())
                || !std::isfinite(y_at_r[2].real()) || !std::isfinite(y_at_r[3].real())) {
                point_invalid[i] = 1;
                continue;
            }

            bool is_solid = true;
            bool is_incompressible = false;
            const auto* physics_layer =
                dynamic_cast<const c_PhysicsLayer*>(world.find_layer_for_radius(radius));
            if (physics_layer != nullptr) {
                is_solid = physics_layer->get_is_solid();
                is_incompressible = physics_layer->get_is_incompressible();
            }

            const std::complex<double> shear = world.calc_complex_shear_modulus(radius, frequency);
            const std::complex<double> bulk  = world.calc_complex_bulk_modulus(radius, frequency);

            const tides::c_StrainRadialCoeffs coeffs = tides::c_compute_strain_radial_coeffs(
                y_at_r[0], y_at_r[1], y_at_r[2], y_at_r[3], shear, bulk, radius,
                degree_l_d, is_solid, is_incompressible);
            if (!coeffs.valid) {
                continue;   // liquid / center: no shear dissipation contribution at this radius
            }

            const c_PotentialPointC potential = c_eval_potential_point_3d(mode, colatitude, 0.0);
            tides::c_Tensor6 strain;
            tides::c_Tensor6 stress;
            tides::c_compute_strain_stress(coeffs, potential, colatitude, strain, stress);

            out_heating[i] += 0.5 * frequency * tides::c_volumetric_heating_signed(stress, strain);
        }
    }

    for (size_t i = 0; i < num_points; ++i) {
        if (point_invalid[i]) {
            out_heating[i] = TidalPyConstants::d_NAN;
        }
    }
}

// =====================================================================================================================
// Collapsed (summed / averaged) 3D tidal heating
// =====================================================================================================================

// Gauss-Legendre nodes/weights on [-1, 1] (Newton-Raphson on the Legendre polynomial). Used for the
// colatitude integral via the substitution x = cos(theta): int_0^pi f(theta) sin(theta) dtheta =
// sum_i w[i] f(acos(x[i])). Spectral accuracy for the smooth secular density (the 1/sin theta factors
// cancel in the physical density), so a modest order integrates the latitude sum to ~machine precision.
inline void c_gauss_legendre_nodes(
        int num_nodes,
        std::vector<double>& nodes,
        std::vector<double>& weights) {
    nodes.resize(num_nodes);
    weights.resize(num_nodes);
    const double pi = TidalPyConstants::d_PI;
    const int half = (num_nodes + 1) / 2;
    for (int i = 0; i < half; ++i) {
        double x = std::cos(pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(num_nodes) + 0.5));
        double dp = 1.0;
        for (int iter = 0; iter < 100; ++iter) {
            double p0 = 1.0;   // P_0
            double p1 = x;     // P_1
            for (int k = 2; k <= num_nodes; ++k) {
                const double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / static_cast<double>(k);
                p0 = p1;
                p1 = p2;
            }
            dp = static_cast<double>(num_nodes) * (x * p1 - p0) / (x * x - 1.0);
            const double dx = -p1 / dp;
            x += dx;
            if (std::abs(dx) <= 1.0e-15) { break; }
        }
        nodes[i]                = -x;
        nodes[num_nodes - 1 - i] = x;
        const double w = 2.0 / ((1.0 - x * x) * dp * dp);
        weights[i]                = w;
        weights[num_nodes - 1 - i] = w;
    }
}

// Produce the 3D tidal heating as a full grid over (radius, colatitude, longitude[, time]) or reduced
// (integrated) along any spatial dimension. orbit_averaged=true gives the secular density h_bar(r,theta)
// (longitude/time independent); orbit_averaged=false gives the instantaneous power sigma_ij(t) eps_dot_ij(t)
// at each user time. See c_Heating3DCollapseConfig / c_Heating3DCollapsed for the flags + output conventions.
inline c_Heating3DCollapsed c_RheologyTide::calc_3d_tidal_heating_collapsed(
        c_LayeredWorld& world,
        const c_TideSolveConfig& state,
        const double* radii,
        size_t num_radii,
        const double* colatitudes,
        size_t num_colatitudes,
        const double* longitudes,
        size_t num_longitudes,
        const double* times,
        size_t num_times,
        const c_Heating3DCollapseConfig& cfg) const {
    c_Heating3DCollapsed result;
    const double two_pi = 2.0 * TidalPyConstants::d_PI;
    const double G_to_use = (tidalpy_config_ptr != nullptr) ? tidalpy_config_ptr->d_G : 6.67430e-11;
    const double min_freq = (tidalpy_config_ptr != nullptr) ? tidalpy_config_ptr->d_MIN_SPIN_ORBIT_DIFF : 1.0e-9;
    const double surface_radius = world.get_radius();
    const c_TideConfig& tide_cfg = world.get_tide_config();
    const bool instantaneous = !cfg.orbit_averaged;
    const bool any_summed = cfg.latitude_summed || cfg.longitude_summed || cfg.radial_summed;

    // Position-independent mode list (built once, reused across every grid point).
    int engine_error = 0;
    const std::vector<c_TidalPotential3DModeCoeff> modes = c_tidal_potential_3d_mode_coeffs(
        surface_radius, state.semi_major_axis, state.orbital_frequency, state.spin_frequency,
        state.obliquity, state.eccentricity, state.host_mass, G_to_use,
        tide_cfg.min_degree_l, tide_cfg.max_degree_l,
        tide_cfg.obliquity_truncation, tide_cfg.eccentricity_truncation, &engine_error);
    if (engine_error != 0) {
        throw std::runtime_error(
            "TidalPy: tidal potential engine failed during 3D tidal heating (error "
            + std::to_string(engine_error) + "); check degree/truncation levels");
    }

    // Axis Grids
    // Radius: user array unless summed (then a per-layer trapezoid grid; r_wsum carries the r^2 Jacobian).
    std::vector<double> r_grid, r_wsum;
    std::vector<size_t> r_layer;
    const size_t num_layers = world.get_num_layers();
    if (cfg.radial_summed) {
        const int slices = (cfg.radial_slices > 1) ? cfg.radial_slices : 80;
        for (size_t layer_i = 0; layer_i < num_layers; ++layer_i) {
            const c_BaseLayer* layer = world.get_layer(layer_i);
            const double r_inner = layer->get_radius_inner();
            const double r_outer = layer->get_radius_outer();
            const double dr = (r_outer - r_inner) / static_cast<double>(slices);
            for (int s = 0; s <= slices; ++s) {
                const double rr = r_inner + dr * static_cast<double>(s);
                r_grid.push_back(rr);
                r_wsum.push_back(((s == 0 || s == slices) ? 0.5 * dr : dr) * rr * rr);  // trapezoid x r^2
                r_layer.push_back(layer_i);
            }
        }
    } else {
        r_grid.assign(radii, radii + num_radii);
    }
    // Colatitude: user array unless summed (Gauss-Legendre in cos theta; the weight absorbs sin theta).
    std::vector<double> th_grid, th_wsum;
    if (cfg.latitude_summed) {
        const int num_nodes = (cfg.latitude_nodes > 1) ? cfg.latitude_nodes : 48;
        std::vector<double> gl_x;
        c_gauss_legendre_nodes(num_nodes, gl_x, th_wsum);
        th_grid.resize(num_nodes);
        for (int i = 0; i < num_nodes; ++i) { th_grid[i] = std::acos(gl_x[i]); }
    } else {
        th_grid.assign(colatitudes, colatitudes + num_colatitudes);
    }
    // Longitude: user array unless summed. Secular -> analytic 2*pi (single point). Instantaneous ->
    // periodic trapezoid over [0, 2*pi) (uniform weight, no endpoint halving since the field is periodic).
    std::vector<double> ph_grid, ph_wsum;
    if (cfg.longitude_summed) {
        if (instantaneous) {
            const int num_nodes = (cfg.longitude_nodes > 1) ? cfg.longitude_nodes : 64;
            const double dph = two_pi / static_cast<double>(num_nodes);
            for (int k = 0; k < num_nodes; ++k) { ph_grid.push_back(k * dph); ph_wsum.push_back(dph); }
        } else {
            ph_grid.push_back(0.0);
            ph_wsum.push_back(two_pi);
        }
    } else {
        ph_grid.assign(longitudes, longitudes + num_longitudes);
    }
    // Time: user array (instantaneous only); a single dummy sample when orbit-averaged.
    std::vector<double> t_grid;
    if (instantaneous) { t_grid.assign(times, times + num_times); }
    else { t_grid.push_back(0.0); }

    const size_t nr = r_grid.size();
    const size_t nth = th_grid.size();
    const size_t nph = ph_grid.size();
    const size_t nt = t_grid.size();

    const bool surv_r = !cfg.radial_summed;
    const bool surv_th = !cfg.latitude_summed;
    const bool surv_ph = !cfg.longitude_summed;
    const bool surv_t = instantaneous;

    // Reported axes (surviving axes carry their grid; summed axes report empty).
    if (surv_r)  { result.radii = r_grid; }
    if (surv_th) { result.colatitudes = th_grid; }
    if (surv_ph) { result.longitudes = ph_grid; }
    if (surv_t)  { result.times = t_grid; }
    result.all_spatial_summed = cfg.radial_summed && cfg.latitude_summed && cfg.longitude_summed;
    result.n_layers = num_layers;
    result.n_times = surv_t ? nt : 1;
    if (surv_r)  { result.shape.push_back(nr); }
    if (surv_th) { result.shape.push_back(nth); }
    if (surv_ph) { result.shape.push_back(nph); }
    if (surv_t)  { result.shape.push_back(nt); }
    if (nr == 0 || nth == 0 || nph == 0 || nt == 0) {
        return result;
    }

    size_t out_size = 1;
    for (size_t s : result.shape) { out_size *= s; }
    result.values.assign(out_size, 0.0);
    if (result.all_spatial_summed) {
        result.layer_totals.assign(num_layers * nt, 0.0);
    }

    // Row-major flat index over the surviving axes, in the fixed order [radius, colatitude, longitude, time].
    auto surv_index = [&](size_t ir, size_t ith, size_t iph, size_t it) -> size_t {
        size_t idx = 0;
        if (surv_r)  { idx = idx * nr + ir; }
        if (surv_th) { idx = idx * nth + ith; }
        if (surv_ph) { idx = idx * nph + iph; }
        if (surv_t)  { idx = idx * nt + it; }
        return idx;
    };
    // Combined per-point weight: summed axes -> integration weight (with Jacobian); surviving spatial axes
    // -> their Jacobian when any axis is summed, else 1 (raw density). Longitude/time Jacobian is 1.
    auto combined_weight = [&](size_t ir, size_t ith, size_t iph) -> double {
        const double wr = cfg.radial_summed ? r_wsum[ir]
                                            : (any_summed ? r_grid[ir] * r_grid[ir] : 1.0);
        const double wth = cfg.latitude_summed ? th_wsum[ith]
                                               : (any_summed ? std::sin(th_grid[ith]) : 1.0);
        const double wph = cfg.longitude_summed ? ph_wsum[iph] : 1.0;
        return wr * wth * wph;
    };

    // per-(unique frequency group) radial coefficients
    // The radial solve depends on (degree l, |frequency|) only, so group the active modes and solve once
    // per group, evaluating the strain radial coefficients at every radius (reused across theta/phi/time).
    struct FreqGroup { int degree_l; double frequency; };
    std::vector<FreqGroup> groups;
    std::vector<int> mode_group(modes.size(), -1);
    for (size_t m = 0; m < modes.size(); ++m) {
        const double frequency = std::abs(modes[m].mode_frequency);
        if (frequency <= min_freq) { continue; }
        int found = -1;
        for (size_t g = 0; g < groups.size(); ++g) {
            if (groups[g].degree_l == modes[m].degree_l && groups[g].frequency == frequency) {
                found = static_cast<int>(g);
                break;
            }
        }
        if (found < 0) {
            found = static_cast<int>(groups.size());
            groups.push_back(FreqGroup{modes[m].degree_l, frequency});
        }
        mode_group[m] = found;
    }

    std::vector<std::vector<tides::c_StrainRadialCoeffs>> group_coeffs(
        groups.size(), std::vector<tides::c_StrainRadialCoeffs>(nr));
    std::vector<unsigned char> radius_solve_failed(nr, 0);  // no depth-resolved solution at this radius
    c_LoveSolveConfig love_cfg;
    for (size_t g = 0; g < groups.size(); ++g) {
        love_cfg.degree_l = groups[g].degree_l;
        love_cfg.frequency_rad_s = groups[g].frequency;
        world.solve_love_numbers(love_cfg);
        if (!world.get_love_success()) {
            throw std::runtime_error(
                "TidalPy: radial solve failed during 3D tidal heating: " + world.get_love_message());
        }
        const ::c_RadialSolutionStorage* storage = world.get_love_storage();
        if (storage == nullptr) {
            throw std::runtime_error("TidalPy: missing radial solution during 3D tidal heating");
        }
        const double degree_l_d = static_cast<double>(groups[g].degree_l);
        for (size_t ir = 0; ir < nr; ++ir) {
            const double radius = r_grid[ir];
            std::complex<double> y_at_r[C_MAX_NUM_Y];
            if (!storage->get_radial_solution(radius, 0, y_at_r)
                || !std::isfinite(y_at_r[0].real()) || !std::isfinite(y_at_r[1].real())
                || !std::isfinite(y_at_r[2].real()) || !std::isfinite(y_at_r[3].real())) {
                radius_solve_failed[ir] = 1;
                group_coeffs[g][ir].valid = false;
                continue;
            }
            bool is_solid = true;
            bool is_incompressible = false;
            const auto* physics_layer =
                dynamic_cast<const c_PhysicsLayer*>(world.find_layer_for_radius(radius));
            if (physics_layer != nullptr) {
                is_solid = physics_layer->get_is_solid();
                is_incompressible = physics_layer->get_is_incompressible();
            }
            const std::complex<double> shear = world.calc_complex_shear_modulus(radius, groups[g].frequency);
            const std::complex<double> bulk  = world.calc_complex_bulk_modulus(radius, groups[g].frequency);
            group_coeffs[g][ir] = tides::c_compute_strain_radial_coeffs(
                y_at_r[0], y_at_r[1], y_at_r[2], y_at_r[3], shear, bulk, radius,
                degree_l_d, is_solid, is_incompressible);
        }
    }

    // Evaluate and Reduce
    const double nan_v = TidalPyConstants::d_NAN;
    if (!instantaneous) {
        // Secular density h_bar(r, theta): sum over modes of (freq/2) Im(sigma_c : conj(eps_c)).
        // Longitude/time independent, so evaluate once per (r, theta) and broadcast over phi.
        for (size_t ir = 0; ir < nr; ++ir) {
            for (size_t ith = 0; ith < nth; ++ith) {
                double density = 0.0;
                if (!radius_solve_failed[ir]) {
                    for (size_t m = 0; m < modes.size(); ++m) {
                        const int g = mode_group[m];
                        if (g < 0 || !group_coeffs[g][ir].valid) { continue; }
                        const c_PotentialPointC potential =
                            c_eval_potential_point_3d(modes[m], th_grid[ith], 0.0);
                        tides::c_Tensor6 strain, stress;
                        tides::c_compute_strain_stress(group_coeffs[g][ir], potential, th_grid[ith],
                                                       strain, stress);
                        density += 0.5 * groups[g].frequency * tides::c_volumetric_heating_signed(stress, strain);
                    }
                }
                for (size_t iph = 0; iph < nph; ++iph) {
                    if (!any_summed) {
                        result.values[surv_index(ir, ith, iph, 0)] =
                            radius_solve_failed[ir] ? nan_v : density;
                    } else if (!radius_solve_failed[ir]) {
                        const double contrib = density * combined_weight(ir, ith, iph);
                        result.values[surv_index(ir, ith, iph, 0)] += contrib;
                        if (result.all_spatial_summed) {
                            result.layer_totals[r_layer[ir] * nt + 0] += contrib;
                        }
                    }
                }
            }
        }
    } else {
        // Instantaneous power sigma_ij(t) eps_dot_ij(t). Build each mode's complex stress/strain amplitude
        // at (r, theta, phi), then evolve in time. A mode with signed omega < 0 uses the conjugated
        // potential phasor at +|omega| (reproduces the true real field, cross terms included).
        std::vector<tides::c_Tensor6> mode_stress, mode_strain;
        std::vector<double> mode_freq;
        for (size_t ir = 0; ir < nr; ++ir) {
            for (size_t ith = 0; ith < nth; ++ith) {
                for (size_t iph = 0; iph < nph; ++iph) {
                    if (radius_solve_failed[ir]) {
                        if (!any_summed) {
                            for (size_t it = 0; it < nt; ++it) {
                                result.values[surv_index(ir, ith, iph, it)] = nan_v;
                            }
                        }
                        continue;
                    }
                    mode_stress.clear();
                    mode_strain.clear();
                    mode_freq.clear();
                    for (size_t m = 0; m < modes.size(); ++m) {
                        const int g = mode_group[m];
                        if (g < 0 || !group_coeffs[g][ir].valid) { continue; }
                        c_PotentialPointC potential =
                            c_eval_potential_point_3d(modes[m], th_grid[ith], ph_grid[iph]);
                        if (modes[m].mode_frequency < 0.0) {
                            potential.U = std::conj(potential.U);
                            potential.dU_dtheta = std::conj(potential.dU_dtheta);
                            potential.dU_dphi = std::conj(potential.dU_dphi);
                            potential.d2U_dtheta2 = std::conj(potential.d2U_dtheta2);
                            potential.d2U_dphi2 = std::conj(potential.d2U_dphi2);
                            potential.d2U_dtheta_dphi = std::conj(potential.d2U_dtheta_dphi);
                        }
                        tides::c_Tensor6 strain, stress;
                        tides::c_compute_strain_stress(group_coeffs[g][ir], potential, th_grid[ith],
                                                       strain, stress);
                        mode_stress.push_back(stress);
                        mode_strain.push_back(strain);
                        mode_freq.push_back(groups[g].frequency);
                    }
                    for (size_t it = 0; it < nt; ++it) {
                        const double time = t_grid[it];
                        double power = 0.0;
                        for (size_t k = 0; k < 6; ++k) {
                            double sigma = 0.0;
                            double eps_dot = 0.0;
                            for (size_t mm = 0; mm < mode_freq.size(); ++mm) {
                                const double omega = mode_freq[mm];
                                const double cos_wt = std::cos(omega * time);
                                const double sin_wt = std::sin(omega * time);
                                const std::complex<double>& sc = mode_stress[mm].c[k];
                                const std::complex<double>& ec = mode_strain[mm].c[k];
                                sigma   += sc.real() * cos_wt - sc.imag() * sin_wt;
                                eps_dot += -omega * (ec.real() * sin_wt + ec.imag() * cos_wt);
                            }
                            power += ((k < 3) ? 1.0 : 2.0) * sigma * eps_dot;
                        }
                        if (!any_summed) {
                            result.values[surv_index(ir, ith, iph, it)] = power;
                        } else {
                            const double contrib = power * combined_weight(ir, ith, iph);
                            result.values[surv_index(ir, ith, iph, it)] += contrib;
                            if (result.all_spatial_summed) {
                                result.layer_totals[r_layer[ir] * nt + it] += contrib;
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}

// World delegation for the collapse path: same preconditions as the scalar get_3d_tidal_heating.
inline c_Heating3DCollapsed c_LayeredWorld::calc_3d_tides(
        const c_TideSolveConfig& state,
        const double* radii,
        size_t num_radii,
        const double* colatitudes,
        size_t num_colatitudes,
        const double* longitudes,
        size_t num_longitudes,
        const double* times,
        size_t num_times,
        const c_Heating3DCollapseConfig& cfg) {
    if (!this->p_tide) {
        throw std::runtime_error(
            "TidalPy: no tide model attached to the world — call set_tide_model() first");
    }
    auto* rheology = dynamic_cast<c_RheologyTide*>(this->p_tide.get());
    if (rheology == nullptr) {
        throw std::runtime_error(
            "TidalPy: 3D tidal heating requires the rheology tide model (the analytic cpl/ctl/ctl_q "
            "models have no depth-resolved radial solution)");
    }
    if (!this->p_eos_solved || !this->p_eos_solution) {
        throw std::runtime_error(
            "TidalPy: 3D tidal heating needs the EOS solved first — call solve_eos()");
    }
    return rheology->calc_3d_tidal_heating_collapsed(
        *this, state, radii, num_radii, colatitudes, num_colatitudes,
        longitudes, num_longitudes, times, num_times, cfg);
}

} // namespace tidalpy
