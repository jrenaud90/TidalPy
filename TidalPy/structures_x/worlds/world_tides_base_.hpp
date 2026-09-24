#pragma once
/*
 * world_tides_base_.hpp: out-of-line definition of c_BaseWorld::calc_tides, the analytic global tidal
 * dissipation path shared by every world type.
 *
 * Runs the global-potential engine (eccentricity and obliquity functions, tidal potential of Renaud et al.
 * 2021) for the world's tide config and the supplied orbital and spin state, then collapses the per-mode terms
 * with the attached analytic tide model (cpl, ctl, ctl_q) into the total tidal heating and the three orbital
 * potential derivatives. The rheology model needs the radial solver, so c_LayeredWorld::calc_tides
 * (world_tides_.hpp) hides this one; it runs the engine through the same c_world_global_potential.
 *
 * This header pulls in the heavy global-potential tables; force-include it in the base-world extension only
 * (other extensions reach c_BaseWorld::calc_tides through inheritance, and world_tides_.hpp includes it).
 */

#include <stdexcept>

#include "base_.hpp"
#include "../../Tides_x/classes/tide_collapse_.hpp"   // c_global_potential, c_collapse_global_tides

namespace tidalpy {

// The global potential of a world's tide config at an orbital state. On failure the world's tide state, passed in as
// tides_solved and tide_result, is marked unsolved with the engine's error code, and std::runtime_error is thrown.
inline c_GlobalPotentialStorage c_world_global_potential(
        const c_BaseWorld& world,
        const c_TideSolveConfig& state,
        bool& tides_solved,
        c_GlobalTideResult& tide_result) {
    const double planet_radius = world.get_radius();
    const double G_to_use = c_get_G();
    const c_TideConfig& tcfg = world.get_tide_config();

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
        tcfg.eccentricity_truncation,
        tcfg.eccentricity_exact_tolerance
    );

    if (potential.error_code != 0) {
        tides_solved           = false;
        tide_result            = c_GlobalTideResult();
        tide_result.error_code = potential.error_code;
        throw std::runtime_error("TidalPy: global potential failed during calc_tides");
    }
    return potential;
}

inline void c_BaseWorld::calc_tides(const c_TideSolveConfig& state) {
    if (!this->p_tide) {
        throw std::runtime_error(
            "TidalPy: no tide model attached to the world: call set_tide_model() first");
    }

    // The rheology model needs the world radial solver, which only a layered world has.
    if (this->p_tide->needs_radial_solve()) {
        throw std::runtime_error(
            "TidalPy: the rheology tide model is only supported on a layered world (it needs "
            "the radial solver); use an analytic model (cpl/ctl/ctl_q) on this world type");
    }

    this->p_check_tide_state(state);
    const c_GlobalPotentialStorage potential =
        c_world_global_potential(*this, state, this->p_tides_solved, this->p_tide_result);

    this->p_tide_solver_love.clear();
    this->p_tide_result  = c_collapse_global_tides(potential, *this->p_tide, nullptr);
    this->p_tides_solved = true;
}

} // namespace tidalpy
