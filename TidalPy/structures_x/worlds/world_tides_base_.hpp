#pragma once
/*
 * world_tides_base_.hpp — out-of-line definition of c_BaseWorld::calc_tides (the analytic
 * global tidal-dissipation path, common to all world types).
 *
 * Runs the global-potential engine (eccentricity/obliquity functions + tidal potential of
 * Renaud et al. 2021) for the world's stored tide config + the supplied orbital/spin state,
 * then collapses the per-mode terms with the attached analytic tide model (cpl/ctl/ctl_q)
 * into the world's total tidal heating + the three orbital potential derivatives.
 *
 * This is the path a layerless world (e.g. a star) uses. The rheology model needs the radial
 * solver and is handled by c_LayeredWorld::calc_tides (world_tides_.hpp), which hides this one.
 *
 * This header pulls in the heavy global-potential tables; force-include it in the base-world
 * extension only (other extensions reach c_BaseWorld::calc_tides through inheritance).
 */

#include <stdexcept>

#include "base_.hpp"
#include "../../Tides_x/classes/tide_collapse_.hpp"   // c_global_potential, c_collapse_global_tides

namespace tidalpy {

inline void c_BaseWorld::calc_tides(const c_TideSolveConfig& state) {
    if (!this->p_tide) {
        throw std::runtime_error(
            "TidalPy: no tide model attached to the world — call set_tide_model() first");
    }

    // The rheology model needs the world radial solver, which only a layered world has.
    if (this->p_tide->needs_radial_solve()) {
        throw std::runtime_error(
            "TidalPy: the rheology tide model is only supported on a layered world (it needs "
            "the radial solver); use an analytic model (cpl/ctl/ctl_q) on this world type");
    }

    const double planet_radius = this->get_radius();
    const double G_to_use = (tidalpy_config_ptr != nullptr) ? tidalpy_config_ptr->d_G : 6.674015e-11;
    const c_TideConfig& tcfg = this->p_tide_config;

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

    this->p_tide_solver_love.clear();
    this->p_tide_result  = c_collapse_global_tides(potential, *this->p_tide, nullptr);
    this->p_tides_solved = true;
}

} // namespace tidalpy
