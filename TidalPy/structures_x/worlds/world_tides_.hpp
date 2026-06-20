#pragma once
/*
 * world_tides_.hpp — out-of-line definition of c_LayeredWorld::calc_tides.
 *
 * Runs the global-potential engine (eccentricity/obliquity functions + tidal potential of
 * Renaud et al. 2021) for the world's stored tide config + the supplied orbital/spin state,
 * then collapses the per-mode terms with the attached tide model into the world's total
 * tidal heating + the three orbital potential derivatives, and distributes per-layer heating
 * by each layer's tidal_scale.
 *
 * This header pulls in the heavy global-potential tables (eccentricity/obliquity). Include it
 * in ONE extension only (the world/layered extension) so those tables don't compile into
 * every translation unit that includes layered_.hpp.
 *
 * The rheology model's per-mode radial-solver Love numbers are wired in a follow-up; for now
 * calc_tides handles the analytic models (cpl/ctl/ctl_q) and raises for rheology.
 */

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "layered_.hpp"
#include "../../Tides_x/classes/tide_collapse_.hpp"   // c_global_potential, c_collapse_global_tides

namespace tidalpy {

inline void c_LayeredWorld::calc_tides(const c_TideSolveConfig& state) {
    if (!this->p_tide) {
        throw std::runtime_error(
            "TidalPy: no tide model attached to the world — call set_tide_model() first");
    }

    // The rheology model needs per-mode radial-solver Love numbers (the world radial solver
    // run at each unique tidal frequency). That coupling lands in a follow-up step; for now
    // the analytic models are fully supported.
    if (this->p_tide->needs_radial_solve()) {
        throw std::runtime_error(
            "TidalPy: the rheology tide model's radial-solver coupling is not yet wired into "
            "calc_tides; use an analytic tide model (cpl/ctl/ctl_q), or call solve_love_numbers "
            "directly, for now");
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

    // Analytic collapse (no radial-solver Love numbers needed).
    this->p_tide_result  = c_collapse_global_tides(potential, *this->p_tide, nullptr);
    this->p_tides_solved = true;

    // Distribute heating to the layers by their tidal scale (0 for non-tidal layers).
    const std::size_t n_layers = this->p_layers.size();
    this->p_layer_tidal_heating.assign(n_layers, 0.0);
    for (std::size_t i = 0; i < n_layers; ++i) {
        const c_BaseLayer* layer = this->p_layers[i].get();
        const double scale = layer->get_is_tidal() ? layer->get_tidal_scale() : 0.0;
        this->p_layer_tidal_heating[i] = this->p_tide_result.tidal_heating * scale;
    }
}

} // namespace tidalpy
