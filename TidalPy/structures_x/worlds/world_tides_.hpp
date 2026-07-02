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

} // namespace tidalpy
