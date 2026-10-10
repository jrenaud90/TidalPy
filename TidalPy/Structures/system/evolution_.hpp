#pragma once
/*
 * evolution_.hpp - the coupled orbital, spin, and thermal evolution of two worlds, one the tidal host of the other
 * (usually each other's), integrated with CyRK in C++ (System.evolve).
 *
 * The state is the shared orbit (a / a0 - 1, e - e_ref), each dissipating body's spin, and, with evolve_thermal, the
 * temperature of every layer of each body that has layers (over its initial value). The semi-major axis and the
 * eccentricity are carried as their changes, so their relative tolerances apply to the changes: the orbit changes by
 * fractions of a percent over Gyr, and CyRK's implicit methods tie their Newton tolerance to the smallest relative
 * tolerance, which a tolerance on a or e itself would have to make tiny. The eccentricity's reference e_ref is its
 * initial value; a segment restarts when e falls to a tenth of it and takes that value as the new reference, so a
 * circularizing orbit keeps its accuracy relative to e as e decays by decades (down to ten times its absolute
 * tolerance). Every right-hand side solves each body's tides at
 * the state, so the spin follows the tidal torque alone; a body with no tide model is rigid (its spin frequency stays
 * as it is) and carries no spin state.
 *
 * A spin is carried as its offset delta = s - j / m from the nearest spin-orbit commensurability j / m (s the spin
 * frequency over the mean motion, m at most the body's largest tidal degree, since a mode of order m is resonant where
 * m s is an integer). The offset reaches the mode frequencies exactly (c_mode_frequency), so a spin held within
 * rounding of a lock (Charon's synchronous offset is about 1e-16) keeps a smooth torque, and its error weight is
 * relative to its distance from the commensurability. With the low-frequency continuation of the tides
 * (c_BaseWorld::calc_continuation_frequency) the torque is a smooth odd function through each commensurability, so a
 * lock is an ordinary stable equilibrium of a stiff variable: the implicit integrator holds it with long steps, and
 * capture, passage, and release follow from the equations with no special handling. Three terminal events, each a
 * function of the offset alone, restart the integration:
 *   - the offset reaches 0.6 of the gap to the neighboring commensurability, past their midpoint: the reference
 *     moves to that neighbor (the margin keeps a spin that starts at the midpoint from rebasing at every step);
 *   - the offset crosses zero (armed only when the segment starts at least d_EVOLVE_ARM_FACTOR delta_lin away): the
 *     integrator restarts at the commensurability, so a lock narrower than its step (Charon's inner peak sits within
 *     1e-9 of synchrony) is not stepped over;
 *   - a disarmed offset moves past twice the arming distance: the crossing event is armed again (the margin keeps a
 *     spin that sits at the arming distance from restarting at every step).
 * delta_lin = omega_c / n is the offset at which the slowest resonant mode reaches the body's continuation frequency
 * omega_c (c_BaseWorld::calc_continuation_frequency).
 *
 * The Jacobian is formed by differences with steps scaled to each variable (a spin column moves only that body's tides,
 * a temperature column that body's EOS and tides), since an integrator's own differences scale with the absolute
 * tolerance, far below the spin offsets' scale. A spin column is a central difference, since a lock sits on a sharply
 * curved stretch of the torque. A body's tides are cached on the exact inputs they depend on (the orbit, its own
 * offset, and its EOS state), so a column reuses the other body's.
 *
 * Assumptions:
 *   - The worlds have no permanent (triaxial) figure and no rotational flattening; the spins follow the tidal torques.
 *   - Layer boundaries are fixed; a two-phase layer melts and freezes inside its own radii (its zones move).
 *   - A thermal body's structure is solved without its tidal heat (that heat enters only its temperature rates).
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cysolve.hpp"     // CyRK: baseline_cysolve_ivp_noreturn, CySolverResult
#include "c_events.hpp"    // CyRK: Event

#include "system_.hpp"     // c_System, c_TidalDissipation, c_WorldEvolution
#include "constants_.hpp"  // TidalPyConstants, tidalpy_config_ptr
#include "../../Tides/potential/potential_common_.hpp"   // c_min_frequency
#include "../../Utilities/classes/model_names_.hpp"   // c_to_lower
#include "../../Utilities/logging/logger_.hpp"

namespace tidalpy {

// The driver integrates in Myr, so the rates it compares are of order one.
inline constexpr double d_EVOLVE_TIME_UNIT = TidalPyConstants::d_SECONDS_PER_MYR;   // [s]
// A segment's crossing event is armed when its offset starts at least this many delta_lin from the commensurability.
inline constexpr double d_EVOLVE_ARM_FACTOR = 10.0;
// Difference steps of the Jacobian: a spin offset by this fraction of max(|delta|, delta_lin), a / a0 and e
// by this fraction of their size (e of at least d_EVOLVE_JACOBIAN_E_FLOOR), a layer temperature by this fraction of
// it. The temperature step stays well above the noise the rates carry from the EOS and Love solves (a third of the
// change a 1e-3 K step makes in a held spin's torque), while the viscosities change by under a percent across it.
inline constexpr double d_EVOLVE_JACOBIAN_SPIN_STEP = 1.0e-3;
inline constexpr double d_EVOLVE_JACOBIAN_ORBIT_STEP = 1.0e-7;
inline constexpr double d_EVOLVE_JACOBIAN_E_FLOOR = 1.0e-6;
inline constexpr double d_EVOLVE_JACOBIAN_TEMPERATURE_STEP = 1.0e-4;
// The spin offset's absolute tolerance over delta_lin.
inline constexpr double d_EVOLVE_SPIN_ATOL_FACTOR = 0.1;
// A run's first step over the span. A segment starts at no more than d_EVOLVE_MAX_FIRST_STEP_FRACTION of what is left
// of the span, and one after a rebase starts at d_EVOLVE_RESTART_STEP_FACTOR times the last step of the one before (a
// rebase changes no rate), one after a spin crossed its commensurability (where it may be captured and relax within a
// few steps of the old size) at d_EVOLVE_CROSSING_STEP_FACTOR times it, and one after a failed evaluation
// d_EVOLVE_RETRY_STEP_FACTOR times shorter than that one did, d_EVOLVE_MAX_FAILED_SEGMENTS times in a row at most.
inline constexpr double d_EVOLVE_FIRST_STEP_FRACTION = 1.0e-10;
inline constexpr double d_EVOLVE_MAX_FIRST_STEP_FRACTION = 0.25;
inline constexpr double d_EVOLVE_RESTART_STEP_FACTOR = 0.5;
inline constexpr double d_EVOLVE_CROSSING_STEP_FACTOR = 1.0e-3;
inline constexpr double d_EVOLVE_RETRY_STEP_FACTOR = 1.0e-3;
inline constexpr std::size_t C_EVOLVE_MAX_FAILED_SEGMENTS = 5;
// Segments in a row that end where they began before a run is stopped as stalled.
inline constexpr std::size_t C_EVOLVE_MAX_STALLED_SEGMENTS = 10;
// The offset, as a fraction of the gap to the neighboring commensurability, at which a spin's reference moves.
inline constexpr double d_EVOLVE_REBASE_SPACINGS = 0.6;
// The fraction of a segment's reference eccentricity at which a decaying e restarts the integration.
inline constexpr double d_EVOLVE_ECCENTRICITY_REBASE = 0.1;
inline constexpr double d_EVOLVE_MAX_ECCENTRICITY = 1.0 - 1.0e-12;
// Absolute tolerance on a / a0 - 1, and on each scaled layer temperature over its relative tolerance.
inline constexpr double d_EVOLVE_SEMI_MAJOR_AXIS_ATOL = 1.0e-12;
inline constexpr double d_EVOLVE_THERMAL_ATOL_FACTOR = 1.0e-3;
inline constexpr std::size_t C_EVOLVE_EXPECTED_STEPS = 256;
inline constexpr std::size_t C_EVOLVE_MAX_RAM_MB     = 2000;
// The least wall time [s] between two progress reports (c_PairEvolveSettings::progress).
inline constexpr double d_EVOLVE_PROGRESS_INTERVAL = 0.2;

// The spin-orbit commensurability j / m (1 <= m <= max_order) nearest a spin ratio, in lowest terms (a tie goes to the
// smaller m, so a fraction that reduces is never chosen over its reduced form).
inline std::pair<int, int> c_nearest_commensurability(double spin_ratio, int max_order) {
    std::pair<int, int> nearest{static_cast<int>(std::lround(spin_ratio)), 1};
    double distance = std::abs(spin_ratio - nearest.first);
    for (int order = 2; order <= max_order; ++order) {
        const int multiple = static_cast<int>(std::lround(order * spin_ratio));
        const double order_distance = std::abs(spin_ratio - static_cast<double>(multiple) / order);
        if (order_distance < distance) {
            nearest  = {multiple, order};
            distance = order_distance;
        }
    }
    return nearest;
}

// The distance from commensurability `numerator` / `denominator` to its neighbor above (side > 0) or below.
inline double c_commensurability_gap(int numerator, int denominator, int max_order, int side) {
    double gap = std::numeric_limits<double>::infinity();
    for (long long order = 1; order <= max_order; ++order) {
        // The multiple of 1 / order next beyond the reference on that side, exactly.
        const long long scaled = order * numerator;
        long long multiple = scaled / denominator;   // Truncated toward zero
        if (side > 0) {
            while (multiple * denominator <= scaled) { ++multiple; }
        } else {
            while (multiple * denominator >= scaled) { --multiple; }
        }
        const long long gap_numerator = std::abs(multiple * denominator - scaled);
        gap = std::min(gap, static_cast<double>(gap_numerator) / static_cast<double>(order * denominator));
    }
    return gap;
}

// Restores a world's own radial-solver settings when a run ends, however it ends.
struct c_RestoreRadialOverrides {
    c_BaseWorld*            world_ptr = nullptr;
    c_RadialSolverOverrides saved;
    ~c_RestoreRadialOverrides() {
        if (this->world_ptr == nullptr) { return; }
        try { this->world_ptr->set_radial_solver_overrides(this->saved); } catch (...) {}
    }
};

class c_EvolveTimeout : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// =====================================================================================================================
// Settings and results
// =====================================================================================================================
// The settings of a pair evolution ([evolution] in the TidalPy configuration; System.evolve documents each).
struct c_PairEvolveSettings {
    bool        evolve_thermal       = true;
    std::string method               = "LSODA";
    double      semi_major_axis_rtol = 1.0e-5;   // on the change a / a0 - 1
    double      eccentricity_rtol    = 1.0e-4;
    double      eccentricity_atol    = 1.0e-8;
    double      spin_rtol            = 1.0e-3;
    double      thermal_rtol         = 1.0e-5;
    double      radial_rtol          = 1.0e-10;   // each body's Love solves during the run
    double      radial_atol          = 1.0e-10;
    double      max_wall_time        = TidalPyConstants::d_NAN;   // [s]; not finite: no cap
    // An optional progress report, called with the time reached [s] as the integration advances, at most every
    // d_EVOLVE_PROGRESS_INTERVAL of wall time, and once at the end. It must not throw; returning false stops the run
    // as max_wall_time does, with what was integrated.
    bool  (*progress)(void* context, double time) = nullptr;
    void* progress_context = nullptr;
};

// One body at the start of a segment: its reference commensurability j / m and whether its crossing event is armed.
struct c_BodySegment {
    double reference_ratio = TidalPyConstants::d_NAN;
    bool   rigid           = false;
    bool   armed           = false;
};

// One integration segment: its start time [s], each body's reference, and how it ended when it did not reach the end
// (empty when it did).
struct c_PairSegment {
    double                     time = TidalPyConstants::d_NAN;
    std::vector<c_BodySegment> bodies;
    std::string                ended;
};

// One body's evolution at every stored step. `temperature` is step-major (steps x num_layers) and empty when the body's
// thermal state was not evolved.
struct c_BodyEvolutionRecord {
    std::size_t                world_index     = 0;
    bool                       rigid           = false;
    bool                       thermal         = false;
    std::size_t                num_layers      = 0;
    std::vector<double>        spin_ratio;       // spin frequency over the mean motion
    std::vector<double>        spin_frequency;   // [rad s-1]
    std::vector<double>        tidal_heating;    // [W]
    std::vector<double>        da_dt;            // this body's contribution [m s-1]
    std::vector<double>        de_dt;            // this body's contribution [s-1]
    std::vector<double>        dspin_dt;         // [rad s-2]
    std::vector<double>        temperature;      // [K]
    std::size_t                num_tide_solves = 0;
    std::size_t                num_eos_solves  = 0;
};

// The pair's evolution at every stored step: the shared time [s], orbit, summed orbital rates, each body's record (the
// world first, then its partner), the segments, and the outcome.
struct c_PairEvolutionRecord {
    std::vector<double>                time;
    std::vector<double>                semi_major_axis;   // [m]
    std::vector<double>                eccentricity;
    std::vector<double>                da_dt;             // [m s-1]
    std::vector<double>                de_dt;             // [s-1]
    std::vector<double>                dn_dt;             // [rad s-2]
    std::vector<c_BodyEvolutionRecord> bodies;
    std::vector<c_PairSegment>         segments;
    std::size_t                        num_rhs_calls = 0;
    std::size_t                        num_jacobians = 0;
    bool                               success = false;
    std::string                        message;
    double                             elapsed = 0.0;     // wall clock [s]
};

// =====================================================================================================================
// The driver
// =====================================================================================================================
// Integrates the shared orbit, the spins, and the layer temperatures of a pair of worlds over a span in CyRK segments
// (see the file header). Construct, then run.
class c_PairEvolver {
public:
    // The pair of `world_index` and its tidal host. Throws std::invalid_argument for a world with no tidal host or no
    // usable orbit, a dissipating world without a finite spin of at least zero, with evolve_thermal a world with layers
    // whose temperatures are not finite and positive, or settings out of range.
    c_PairEvolver(c_System* system_ptr, std::size_t world_index, const c_PairEvolveSettings& settings)
            : p_system_ptr(system_ptr), p_settings(settings) {
        if (system_ptr == nullptr) {
            throw std::invalid_argument("TidalPy: System.evolve needs a system.");
        }
        system_ptr->check_index(world_index);
        if (!system_ptr->has_tidal_host(world_index)) {
            throw std::invalid_argument("TidalPy: System.evolve needs a world with a tidal host.");
        }
        this->p_check_settings();
        this->p_method = p_evolve_method(settings.method);
        this->p_orbit_index = world_index;
        this->p_orbital_frequency0 = system_ptr->calc_orbital_frequency(world_index);
        if (!(std::isfinite(this->p_orbital_frequency0) && (this->p_orbital_frequency0 > 0.0))) {
            throw std::invalid_argument(
                "TidalPy: System.evolve needs a usable orbit (a positive semi-major axis) about the tidal host.");
        }
        this->p_semi_major_axis0  = system_ptr->get_semi_major_axis(world_index);
        this->p_orbital_frequency = this->p_orbital_frequency0;

        const std::size_t indices[2] = {
            world_index, static_cast<std::size_t>(system_ptr->get_tidal_host_index(world_index))};
        std::size_t slot = 2;
        for (std::size_t b = 0; b < 2; ++b) {
            c_EvolveBody& body = this->p_bodies[b];
            body.index      = indices[b];
            body.world_ptr  = system_ptr->p_worlds[indices[b]].get();
            body.rigid      = !body.world_ptr->get_tide_model_set();
            body.num_layers = body.world_ptr->get_num_layers();
            body.thermal    = settings.evolve_thermal && (body.num_layers > 0);
            body.spin_frequency0 = body.world_ptr->get_spin_frequency();
            if (!body.rigid && !(std::isfinite(body.spin_frequency0) && (body.spin_frequency0 >= 0.0))) {
                throw std::invalid_argument(
                    "TidalPy: System.evolve needs each dissipating world spinning prograde (a finite spin_frequency of "
                    "at least zero); '" + body.world_ptr->get_name() + "' does not.");
            }
            body.max_order = std::max(body.world_ptr->get_tide_config().max_degree_l, 1);
            const double spin_ratio = body.spin_frequency0 / this->p_orbital_frequency0;
            this->p_set_reference(body, c_nearest_commensurability(spin_ratio, body.max_order));
            body.offset = spin_ratio - this->p_reference_ratio(body);
            if (body.thermal) {
                body.temperature_slot = slot;
                body.temperature0.resize(body.num_layers);
                for (std::size_t i = 0; i < body.num_layers; ++i) {
                    const double temperature = body.world_ptr->get_layer(i)->get_temperature();
                    if (!(std::isfinite(temperature) && (temperature > 0.0))) {
                        throw std::invalid_argument(
                            "TidalPy: System.evolve with evolve_thermal needs every layer's temperature finite and "
                            "positive [K]; world '" + body.world_ptr->get_name() + "' layer "
                            + std::to_string(i) + " is " + std::to_string(temperature) + " K.");
                    }
                    body.temperature0[i] = temperature;
                }
                slot += body.num_layers;
            }
        }
        for (c_EvolveBody& body : this->p_bodies) {
            if (!body.rigid) { body.spin_slot = slot++; }
            if (body.thermal && !system_ptr->has_star()) {
                TIDALPY_LOG_WARN(
                    "TidalPy: System.evolve evolves the layer temperatures of '{}' in a system without a star, so no "
                    "heat leaves its surface. Add the star for its insolation temperature.",
                    body.world_ptr->get_name());
            }
        }
        this->p_num_y = slot;
    }

    // Integrates from `t_start` to `t_end` [s] and leaves the system at the final state. A failure stops the run and
    // returns what was integrated, with success false and the reason. Throws std::invalid_argument for a span that is
    // not finite and increasing.
    c_PairEvolutionRecord run(double t_start, double t_end) {
        if (!(std::isfinite(t_start) && std::isfinite(t_end) && (t_end > t_start))) {
            throw std::invalid_argument(
                "TidalPy: System.evolve needs a finite time span (t_start, t_end) with t_end > t_start [s].");
        }
        const auto started = std::chrono::steady_clock::now();
        this->p_next_progress = started;
        this->p_progress_time = -TidalPyConstants::d_INF;
        this->p_stopped = false;
        this->p_has_deadline = std::isfinite(this->p_settings.max_wall_time);
        if (this->p_has_deadline) {
            this->p_deadline = started + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                std::chrono::duration<double>(std::max(this->p_settings.max_wall_time, 0.0)));
        }
        this->p_record = c_PairEvolutionRecord();
        this->p_record.bodies.resize(2);
        for (std::size_t b = 0; b < 2; ++b) {
            c_BodyEvolutionRecord& body_record = this->p_record.bodies[b];
            body_record.world_index = this->p_bodies[b].index;
            body_record.rigid       = this->p_bodies[b].rigid;
            body_record.thermal     = this->p_bodies[b].thermal;
            body_record.num_layers  = this->p_bodies[b].thermal ? this->p_bodies[b].num_layers : 0;
        }

        // Each body's Love solves run at the evolution's radial tolerances, and its own settings come back afterwards.
        c_RestoreRadialOverrides restore[2];
        for (std::size_t b = 0; b < 2; ++b) {
            c_BaseWorld* world_ptr = this->p_bodies[b].world_ptr;
            c_RadialSolverOverrides overrides = world_ptr->get_radial_solver_overrides();
            restore[b].saved = overrides;
            restore[b].world_ptr = world_ptr;
            overrides.set("rtol", this->p_settings.radial_rtol);
            overrides.set("atol", this->p_settings.radial_atol);
            world_ptr->set_radial_solver_overrides(overrides);
        }
        std::string message;
        const bool success = this->p_integrate(t_start / d_EVOLVE_TIME_UNIT, t_end / d_EVOLVE_TIME_UNIT, message);
        this->p_record.success = success;
        this->p_record.message = message;

        // Leave the system at the final state: the last stored step, or the start when none was.
        if (!this->p_last_y.empty()) {
            this->p_has_deadline = false;
            this->p_stopped = false;
            try {
                std::vector<double> scratch(this->p_num_y + C_RECORD_SIZE);
                this->p_e_reference = this->p_last_e_reference;
                for (std::size_t b = 0; b < 2; ++b) {
                    this->p_set_reference(this->p_bodies[b], this->p_last_reference[b]);
                }
                this->p_evaluate(this->p_last_time, this->p_last_y.data(), scratch.data());
                for (c_EvolveBody& body : this->p_bodies) {
                    if (body.rigid) { continue; }
                    body.world_ptr->set_spin_frequency(this->p_spin_ratio(body) * this->p_orbital_frequency);
                }
            } catch (const std::exception& error) {
                this->p_record.success = false;
                this->p_record.message += std::string(" The system could not be set to the final state: ")
                    + error.what();
            }
        }
        for (std::size_t b = 0; b < 2; ++b) {
            this->p_record.bodies[b].num_tide_solves = this->p_bodies[b].num_tide_solves;
            this->p_record.bodies[b].num_eos_solves  = this->p_bodies[b].num_eos_solves;
        }
        if (!this->p_last_y.empty()) { this->p_report_progress(this->p_last_time * d_EVOLVE_TIME_UNIT, true); }
        this->p_record.elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        TIDALPY_LOG_INFO(
            "TidalPy: System.evolve: {} {} steps, {} segments, {} right-hand sides, {} Jacobians, {} and {} tide "
            "solves, {:.1f} s.",
            this->p_record.message, this->p_record.time.size(), this->p_record.segments.size(),
            this->p_record.num_rhs_calls, this->p_record.num_jacobians, this->p_bodies[0].num_tide_solves,
            this->p_bodies[1].num_tide_solves, this->p_record.elapsed);
        return std::move(this->p_record);
    }

private:
    // A body's rates (da/dt, de/dt, dn/dt, dspin/dt, heating, then each layer's dT/dt [SI]) at a set of inputs.
    struct c_TideEntry {
        bool                valid = false;
        double              key[3] = {0.0, 0.0, 0.0};
        int                 numerator = 0;
        int                 denominator = 1;
        std::size_t         eos_serial = 0;
        std::vector<double> rates;
    };

    struct c_EvolveBody {
        std::size_t         index      = 0;          // in the system
        c_BaseWorld*        world_ptr  = nullptr;
        bool                rigid      = false;      // no tide model
        bool                thermal    = false;      // its layer temperatures are evolved
        std::size_t         num_layers = 0;
        std::size_t         temperature_slot = 0;    // its first temperature in the state
        std::size_t         spin_slot  = 0;          // its spin offset in the state (dissipating bodies)
        std::vector<double> temperature0;            // [K]
        double              spin_frequency0 = 0.0;   // [rad s-1]; a rigid body's throughout
        int                 max_order   = 2;         // the largest m of its modes (its largest tidal degree)
        int                 numerator   = 0;         // the segment's commensurability j / m
        int                 denominator = 1;
        double              gap_below   = 1.0;       // the distances to its neighboring commensurabilities
        double              gap_above   = 1.0;
        double              offset      = 0.0;       // delta at the segment's start
        bool                armed       = false;

        // Its EOS state: the time [Myr], scaled temperatures, and surface temperature it was solved at.
        bool                eos_valid = false;
        double              eos_time  = TidalPyConstants::d_NAN;
        double              eos_surface = TidalPyConstants::d_NAN;
        std::vector<double> eos_temperatures;
        std::size_t         eos_serial = 0;

        // Its rates at the latest sets of the inputs they depend on (a / a0, |e|, the offset, the reference, and the
        // EOS state): the base and the two sides of a central difference, so a Jacobian column that moves only the
        // other body reuses this one's base rates.
        c_TideEntry         tides[3];
        std::size_t         newest = 0;

        std::size_t         num_tide_solves = 0;
        std::size_t         num_eos_solves  = 0;
    };

    // The step quantities each right-hand side call returns as CyRK extra outputs (stored with every step): per body
    // (spin ratio, heating, da/dt, de/dt, dspin/dt), then the summed da/dt, de/dt, dn/dt.
    static constexpr std::size_t C_RECORD_PER_BODY = 5;
    static constexpr std::size_t C_RECORD_SIZE = 2 * C_RECORD_PER_BODY + 3;

    static ODEMethod p_evolve_method(const std::string& name) {
        const std::string lower = c_to_lower(name);
        if (lower == "lsoda") { return ODEMethod::LSODA; }
        if (lower == "bdf")   { return ODEMethod::BDF; }
        if (lower == "radau") { return ODEMethod::RADAU; }
        throw std::invalid_argument(
            "TidalPy: System.evolve's method must be 'LSODA', 'BDF', or 'Radau' (an implicit method; the spins are "
            "stiff near a lock); got '" + name + "'.");
    }

    static std::string p_format_time(double time) {
        char buffer[64];
        std::snprintf(buffer, sizeof(buffer), "%.10g", time);
        return buffer;
    }

    void p_check_settings() const {
        const c_PairEvolveSettings& settings = this->p_settings;
        const std::pair<const char*, double> tolerances[7] = {
            {"semi_major_axis_rtol", settings.semi_major_axis_rtol}, {"eccentricity_rtol", settings.eccentricity_rtol},
            {"eccentricity_atol", settings.eccentricity_atol}, {"spin_rtol", settings.spin_rtol},
            {"thermal_rtol", settings.thermal_rtol}, {"radial_rtol", settings.radial_rtol},
            {"radial_atol", settings.radial_atol}};
        for (const auto& [name, value] : tolerances) {
            if (!(std::isfinite(value) && (value > 0.0))) {
                throw std::invalid_argument(std::string("TidalPy: System.evolve's ") + name + " must be positive.");
            }
        }
    }

    void p_check_deadline() const {
        if (this->p_stopped) { throw c_EvolveTimeout("TidalPy: System.evolve was interrupted."); }
        if (this->p_has_deadline && (std::chrono::steady_clock::now() > this->p_deadline)) {
            throw c_EvolveTimeout("TidalPy: System.evolve reached its max_wall_time.");
        }
    }

    // The spin offset at which a body's slowest resonant mode (m = 1) reaches its world's continuation frequency
    // (c_BaseWorld::calc_continuation_frequency), at the current state.
    double p_offset_scale(const c_EvolveBody& body) const {
        const double floor = body.world_ptr->calc_continuation_frequency();
        return std::max(floor, TidalPyConstants::d_EPS * this->p_orbital_frequency) / this->p_orbital_frequency;
    }

    double p_reference_ratio(const c_EvolveBody& body) const {
        return static_cast<double>(body.numerator) / static_cast<double>(body.denominator);
    }

    void p_set_reference(c_EvolveBody& body, const std::pair<int, int>& reference) const {
        body.numerator   = reference.first;
        body.denominator = reference.second;
        body.gap_below   = c_commensurability_gap(body.numerator, body.denominator, body.max_order, -1);
        body.gap_above   = c_commensurability_gap(body.numerator, body.denominator, body.max_order, 1);
    }

    double p_spin_ratio(const c_EvolveBody& body) const {
        if (body.rigid) { return body.spin_frequency0 / this->p_orbital_frequency; }
        return this->p_reference_ratio(body) + this->p_offset(body);
    }

    double p_offset(const c_EvolveBody& body) const { return this->p_current_offset[&body - this->p_bodies]; }

    // =================================================================================================================
    // Right-hand side
    // =================================================================================================================
    // Puts the system at state `y` at `time` [Myr]: the orbit (e extended oddly through zero, so the state crosses it
    // smoothly), each thermal body's layer temperatures and EOS (solved again only when its time, temperatures, or
    // surface temperature changed; without a star or for the star itself no flow leaves the surface), and each spin
    // offset. Throws std::runtime_error when an EOS solve fails.
    void p_set_state(double time, const double* y) {
        this->p_check_deadline();
        const double eccentricity = std::min(std::abs(this->p_e_reference + y[1]), d_EVOLVE_MAX_ECCENTRICITY);
        if (!((y[0] == this->p_orbit_key[0]) && (eccentricity == this->p_orbit_key[1]))) {
            this->p_system_ptr->set_semi_major_axis(this->p_orbit_index, (1.0 + y[0]) * this->p_semi_major_axis0);
            this->p_system_ptr->set_eccentricity(this->p_orbit_index, eccentricity);
            this->p_orbital_frequency = this->p_system_ptr->calc_orbital_frequency(this->p_orbit_index);
            this->p_orbit_key[0] = y[0];
            this->p_orbit_key[1] = eccentricity;
        }
        for (std::size_t b = 0; b < 2; ++b) {
            c_EvolveBody& body = this->p_bodies[b];
            this->p_current_offset[b] = body.rigid ? 0.0 : y[body.spin_slot];
            if (!body.thermal) { continue; }
            const double* temperatures = y + body.temperature_slot;
            const double surface = this->p_surface_temperature(body);
            const bool same_surface = (surface == body.eos_surface)
                || (std::isnan(surface) && std::isnan(body.eos_surface));
            if (body.eos_valid && (time == body.eos_time) && same_surface
                    && std::equal(temperatures, temperatures + body.num_layers, body.eos_temperatures.begin())) {
                continue;
            }
            body.eos_valid = false;   // Until its solve succeeds
            for (std::size_t i = 0; i < body.num_layers; ++i) {
                body.world_ptr->get_layer(i)->set_temperature(temperatures[i] * body.temperature0[i]);
            }
            body.world_ptr->clear_tidal_heating();
            c_WorldEOSSolveConfig config = body.world_ptr->make_eos_solve_config();
            config.solve_temperature   = true;
            config.surface_temperature = surface;
            config.time                = time * d_EVOLVE_TIME_UNIT;
            const c_WorldEOSReport report = body.world_ptr->solve_eos_report(config);
            ++body.num_eos_solves;
            if (!report.success) {
                throw std::runtime_error(
                    "TidalPy: the EOS solve of world '" + body.world_ptr->get_name() + "' failed: " + report.message);
            }
            body.eos_time    = time;
            body.eos_surface = surface;
            body.eos_temperatures.assign(temperatures, temperatures + body.num_layers);
            body.eos_valid   = true;
            ++body.eos_serial;
        }
    }

    double p_surface_temperature(const c_EvolveBody& body) const {
        if (!this->p_system_ptr->has_star() || (static_cast<int>(body.index) == this->p_system_ptr->get_star_index())) {
            return TidalPyConstants::d_NAN;
        }
        return this->p_system_ptr->calc_equilibrium_temperature(body.index);
    }

    // A body's rates at the current state (cached on its exact inputs): da/dt, de/dt, dn/dt, dspin/dt, heating, then
    // each layer's dT/dt, in SI. A rigid body contributes nothing but its layers' rates.
    const std::vector<double>& p_body_rates(std::size_t body_index) {
        c_EvolveBody& body = this->p_bodies[body_index];
        const double key[3] = {this->p_orbit_key[0], this->p_orbit_key[1], this->p_current_offset[body_index]};
        for (const c_TideEntry& entry : body.tides) {
            if (entry.valid && std::equal(key, key + 3, entry.key) && (entry.numerator == body.numerator)
                    && (entry.denominator == body.denominator) && (entry.eos_serial == body.eos_serial)) {
                return entry.rates;
            }
        }
        body.newest = (body.newest + 1) % 3;   // The oldest entry is replaced
        c_TideEntry& entry = body.tides[body.newest];
        entry.valid = false;
        std::vector<double>& rates = entry.rates;
        rates.assign(5 + (body.thermal ? body.num_layers : 0), 0.0);
        if (!body.rigid) {
            const double offset = this->p_current_offset[body_index];
            body.world_ptr->set_spin_frequency((this->p_reference_ratio(body) + offset) * this->p_orbital_frequency);
            c_WorldEvolution evolution;
            try {
                const c_TidalDissipation dissipation = this->p_system_ptr->p_dissipation(
                    body.index, this->p_bodies[1 - body_index].index, this->p_orbit_index,
                    body.numerator, body.denominator, offset);
                evolution = this->p_system_ptr->p_evolution(dissipation);
            } catch (const std::exception& error) {
                char state[160];
                std::snprintf(state, sizeof(state), " (at spin ratio %d/%d %+.6e, a/a0 - 1 %+.6e, e %.6e)",
                              body.numerator, body.denominator, offset, this->p_orbit_key[0], this->p_orbit_key[1]);
                throw std::runtime_error("TidalPy: the tides of '" + body.world_ptr->get_name() + "' failed" + state
                                         + ": " + error.what());
            }
            rates[0] = evolution.da_dt;
            rates[1] = evolution.de_dt;
            rates[2] = evolution.dn_dt;
            rates[3] = evolution.dspin_dt;
            rates[4] = evolution.tidal_heating;
            ++body.num_tide_solves;
        }
        if (body.thermal) {
            for (std::size_t i = 0; i < body.num_layers; ++i) {
                rates[5 + i] = body.world_ptr->calc_layer_temperature_rate(i);
            }
        }
        for (const double rate : rates) {
            if (!std::isfinite(rate)) {
                throw std::runtime_error(
                    "TidalPy: the rates of '" + body.world_ptr->get_name() + "' at spin ratio "
                    + std::to_string(this->p_spin_ratio(body)) + " are not finite.");
            }
        }
        std::copy(key, key + 3, entry.key);
        entry.numerator   = body.numerator;
        entry.denominator = body.denominator;
        entry.eos_serial  = body.eos_serial;
        entry.valid      = true;
        return rates;
    }

    // The scaled rates `dy` [per Myr] of state `y` at `time` [Myr], followed by the step quantities (C_RECORD_SIZE
    // extra outputs).
    void p_evaluate(double time, const double* y, double* dy) {
        ++this->p_record.num_rhs_calls;
        // The integrator stores a segment's steps only when the segment ends, so progress follows its evaluations.
        this->p_report_progress(time * d_EVOLVE_TIME_UNIT, false);
        this->p_set_state(time, y);
        const std::vector<double>* rates[2] = {&this->p_body_rates(0), &this->p_body_rates(1)};
        const double da_dt = (*rates[0])[0] + (*rates[1])[0];
        // The odd extension through e = 0: a negative e has the rate of |e| with its sign turned.
        const double de_dt = std::copysign(1.0, this->p_e_reference + y[1]) * ((*rates[0])[1] + (*rates[1])[1]);
        const double dn_dt = (*rates[0])[2] + (*rates[1])[2];
        dy[0] = da_dt / this->p_semi_major_axis0 * d_EVOLVE_TIME_UNIT;
        dy[1] = de_dt * d_EVOLVE_TIME_UNIT;
        double* record = dy + this->p_num_y;
        for (std::size_t b = 0; b < 2; ++b) {
            const c_EvolveBody& body = this->p_bodies[b];
            const std::vector<double>& body_rates = *rates[b];
            if (body.thermal) {
                for (std::size_t i = 0; i < body.num_layers; ++i) {
                    dy[body.temperature_slot + i] = body_rates[5 + i] / body.temperature0[i] * d_EVOLVE_TIME_UNIT;
                }
            }
            const double spin_ratio = this->p_spin_ratio(body);
            if (!body.rigid) {
                dy[body.spin_slot] = (body_rates[3] - spin_ratio * dn_dt) / this->p_orbital_frequency
                    * d_EVOLVE_TIME_UNIT;
            }
            double* entry = record + b * C_RECORD_PER_BODY;
            entry[0] = spin_ratio;
            entry[1] = body_rates[4];
            entry[2] = body_rates[0];
            entry[3] = body_rates[1];
            entry[4] = body_rates[3];
        }
        record[2 * C_RECORD_PER_BODY]     = da_dt;
        record[2 * C_RECORD_PER_BODY + 1] = de_dt;
        record[2 * C_RECORD_PER_BODY + 2] = dn_dt;
    }

    // The difference Jacobian (column-major) of the scaled rates at `y`, each column's step scaled to its
    // variable (see the file header). The spin columns come first, at the base EOS states, then the orbit, which can
    // move a thermal body's surface temperature, then the temperatures, each solving its body's EOS once.
    void p_jacobian(double time, const double* y, double* jacobian) {
        ++this->p_record.num_jacobians;
        const std::size_t num_y = this->p_num_y;
        std::vector<double> base(num_y + C_RECORD_SIZE);
        std::vector<double> moved(num_y + C_RECORD_SIZE);
        std::vector<double> y_step(y, y + num_y);
        this->p_evaluate(time, y, base.data());
        struct c_Column { std::size_t index; double step; bool central; };
        std::vector<c_Column> columns;
        for (const c_EvolveBody& body : this->p_bodies) {
            if (!body.rigid) {
                const std::size_t j = body.spin_slot;
                columns.push_back(c_Column{
                    j, d_EVOLVE_JACOBIAN_SPIN_STEP * std::max(std::abs(y[j]), this->p_offset_scale(body)), true});
            }
        }
        columns.push_back(c_Column{0, d_EVOLVE_JACOBIAN_ORBIT_STEP * std::max(std::abs(1.0 + y[0]), 1.0), false});
        columns.push_back(c_Column{
            1, d_EVOLVE_JACOBIAN_ORBIT_STEP * std::max(std::abs(this->p_e_reference + y[1]), d_EVOLVE_JACOBIAN_E_FLOOR),
            false});
        for (const c_EvolveBody& body : this->p_bodies) {
            for (std::size_t i = 0; body.thermal && (i < body.num_layers); ++i) {
                columns.push_back(c_Column{
                    body.temperature_slot + i,
                    d_EVOLVE_JACOBIAN_TEMPERATURE_STEP * std::abs(y[body.temperature_slot + i]), false});
            }
        }
        std::vector<double> lower(num_y + C_RECORD_SIZE);
        for (const c_Column& column : columns) {
            const std::size_t j = column.index;
            y_step[j] = y[j] + column.step;
            const double upper_step = y_step[j] - y[j];   // The steps as represented
            this->p_evaluate(time, y_step.data(), moved.data());
            double lower_step = 0.0;
            if (column.central) {
                y_step[j] = y[j] - column.step;
                lower_step = y[j] - y_step[j];
                this->p_evaluate(time, y_step.data(), lower.data());
            }
            const std::vector<double>& below = column.central ? lower : base;
            for (std::size_t i = 0; i < num_y; ++i) {
                jacobian[j * num_y + i] = (moved[i] - below[i]) / (upper_step + lower_step);
            }
            y_step[j] = y[j];
        }
    }

    // CyRK's view of the right-hand side and Jacobian. Nothing may escape the solver, so a failure is kept, stops the
    // solve, and returns zeros.
    static void p_diffeq(double* dy, double time, double* y, char* args, PreEvalFunc) {
        c_PairEvolver* self_ptr = nullptr;
        std::memcpy(&self_ptr, args, sizeof(self_ptr));
        self_ptr->p_guarded(dy, self_ptr->p_num_y + C_RECORD_SIZE, [&]() { self_ptr->p_evaluate(time, y, dy); });
    }

    static void p_jacobian_callback(double* jacobian, double time, double* y, char* args, PreEvalFunc) {
        c_PairEvolver* self_ptr = nullptr;
        std::memcpy(&self_ptr, args, sizeof(self_ptr));
        self_ptr->p_guarded(jacobian, self_ptr->p_num_y * self_ptr->p_num_y,
                            [&]() { self_ptr->p_jacobian(time, y, jacobian); });
    }

    template <class Call>
    void p_guarded(double* out, std::size_t size, Call&& call) noexcept {
        if (this->p_error.empty()) {
            try {
                call();
                return;
            } catch (const c_EvolveTimeout& error) {
                this->p_fail(error.what(), true);
            } catch (const std::exception& error) {
                this->p_fail(error.what(), false);
            } catch (...) {
                this->p_fail("TidalPy: System.evolve's right-hand side failed with an unknown error.", false);
            }
        }
        std::fill(out, out + size, 0.0);
        this->p_stop_solve();
    }

    void p_fail(const std::string& message, bool timeout) noexcept {
        if (this->p_error.empty()) {
            this->p_error   = message;
            this->p_timeout = timeout;
        }
    }

    void p_stop_solve() noexcept {
        if ((this->p_result_ptr != nullptr) && this->p_result_ptr->solver_uptr) {
            this->p_result_ptr->solver_uptr->set_external_error(CyrkErrorCodes::OTHER_ERROR);
        }
    }

    // =================================================================================================================
    // Integration
    // =================================================================================================================
    // The events of a segment (see the file header), each positive at its start and fired on a downward crossing.
    // Without `crossings`, a spin's crossing and re-arming events are left out (after an event fired at its segment's
    // start, see p_integrate).
    std::vector<Event> p_build_events(bool crossings) {
        std::vector<Event> events;
        this->p_event_is_crossing.clear();
        this->p_event_is_decay.clear();
        // Below its absolute tolerance an eccentricity's decay needs no relative accuracy.
        if (std::abs(this->p_e_reference) > 10.0 * this->p_settings.eccentricity_atol) {
            const double reference = this->p_e_reference;
            Event& decayed = events.emplace_back(nullptr, 1, -1);
            decayed.check = [reference](Event*, double, double* y, char*) -> double {
                return std::abs(reference + y[1]) - d_EVOLVE_ECCENTRICITY_REBASE * std::abs(reference);
            };
            this->p_event_is_crossing.push_back(false);
            this->p_event_is_decay.push_back(true);
        }
        for (c_EvolveBody& body : this->p_bodies) {
            if (body.rigid) { continue; }
            const double offset_scale = this->p_offset_scale(body);
            const std::size_t slot = body.spin_slot;
            const double above = d_EVOLVE_REBASE_SPACINGS * body.gap_above;
            const double below = d_EVOLVE_REBASE_SPACINGS * body.gap_below;
            Event& rebase = events.emplace_back(nullptr, 1, -1);
            rebase.check = [slot, above, below](Event*, double, double* y, char*) -> double {
                return (y[slot] >= 0.0) ? (above - y[slot]) : (below + y[slot]);
            };
            this->p_event_is_crossing.push_back(false);
            this->p_event_is_decay.push_back(false);
            const double arm_distance = d_EVOLVE_ARM_FACTOR * offset_scale;
            body.armed = crossings && (std::abs(body.offset) > arm_distance);
            if (!crossings) { continue; }
            Event& crossing = events.emplace_back(nullptr, 1, -1);
            this->p_event_is_crossing.push_back(body.armed);
            this->p_event_is_decay.push_back(false);
            if (body.armed) {
                const double side = std::copysign(1.0, body.offset);
                crossing.check = [slot, side](Event*, double, double* y, char*) -> double { return side * y[slot]; };
            } else {
                crossing.check = [slot, arm_distance](Event*, double, double* y, char*) -> double {
                    return 2.0 * arm_distance - std::abs(y[slot]);
                };
            }
        }
        return events;
    }

    bool p_integrate(double t_start, double t_end, std::string& message) {
        double time = t_start;
        std::vector<double> y_now(this->p_num_y, 1.0);
        y_now[0] = 0.0;
        y_now[1] = 0.0;
        this->p_e_reference = this->p_system_ptr->get_eccentricity(this->p_orbit_index);
        for (const c_EvolveBody& body : this->p_bodies) {
            if (!body.rigid) { y_now[body.spin_slot] = body.offset; }
        }
        this->p_keep_last(t_start, y_now.data());
        double first_step = d_EVOLVE_FIRST_STEP_FRACTION * (t_end - t_start);
        std::size_t failed = 0;
        std::size_t stalled = 0;
        bool first_segment = true;
        // A remainder of the span within rounding of its end is done (a segment needs room for its first step).
        const double end_rounding = 4.0 * std::numeric_limits<double>::epsilon() * std::abs(t_end);
        while (t_end - time > end_rounding) {
            if (stalled >= C_EVOLVE_MAX_STALLED_SEGMENTS) {
                message = "TidalPy: System.evolve stalled: " + std::to_string(stalled)
                    + " segments in a row ended where they began, at t = " + p_format_time(time) + " Myr.";
                return false;
            }
            if (failed >= C_EVOLVE_MAX_FAILED_SEGMENTS) {
                message = "TidalPy: System.evolve stopped at t = " + p_format_time(time) + " Myr: "
                    + std::to_string(failed) + " segments in a row ended on a failed evaluation (" + this->p_last_error
                    + ").";
                return false;
            }
            // The eccentricity's reference moves to its value now after it decayed to a tenth of the last one.
            if (this->p_rebase_eccentricity) {
                this->p_e_reference += y_now[1];
                y_now[1] = 0.0;
                this->p_rebase_eccentricity = false;
            }
            // Each spin's reference: the commensurability nearest it, and its offset from there.
            for (c_EvolveBody& body : this->p_bodies) {
                if (body.rigid) { continue; }
                const double spin_ratio = this->p_reference_ratio(body) + y_now[body.spin_slot];
                const std::pair<int, int> nearest = c_nearest_commensurability(spin_ratio, body.max_order);
                if ((nearest.first != body.numerator) || (nearest.second != body.denominator)) {
                    this->p_set_reference(body, nearest);
                    y_now[body.spin_slot] = spin_ratio - this->p_reference_ratio(body);
                }
                body.offset = y_now[body.spin_slot];
            }
            // The orbit at the segment's start sets the offset scale of its events and tolerances. A failure there (the
            // wall-clock cap, a solve failing at a state the run reached) ends the run.
            std::vector<Event> events;
            try {
                this->p_set_state(time, y_now.data());
                // An event that fired at its segment's start (a spin passing its arming distance within the time
                // resolution of the event search) would fire again at once: the next segment leaves the spins'
                // crossing events out until another event ends it.
                events = this->p_build_events(!this->p_event_at_start);
            } catch (const std::exception& error) {
                message = error.what();
                return false;
            }
            c_PairSegment entry;
            entry.time = time * d_EVOLVE_TIME_UNIT;
            for (const c_EvolveBody& body : this->p_bodies) {
                entry.bodies.push_back(c_BodySegment{
                    body.rigid ? TidalPyConstants::d_NAN : this->p_reference_ratio(body), body.rigid, body.armed});
            }
            std::vector<double> rtols(this->p_num_y, this->p_settings.thermal_rtol);
            std::vector<double> atols(this->p_num_y, d_EVOLVE_THERMAL_ATOL_FACTOR * this->p_settings.thermal_rtol);
            rtols[0] = this->p_settings.semi_major_axis_rtol;
            atols[0] = d_EVOLVE_SEMI_MAJOR_AXIS_ATOL;
            rtols[1] = this->p_settings.eccentricity_rtol;
            atols[1] = this->p_settings.eccentricity_atol;
            for (const c_EvolveBody& body : this->p_bodies) {
                if (body.rigid) { continue; }
                rtols[body.spin_slot] = this->p_settings.spin_rtol;
                atols[body.spin_slot] = d_EVOLVE_SPIN_ATOL_FACTOR * this->p_offset_scale(body);
            }

            first_step = std::min(first_step, d_EVOLVE_MAX_FIRST_STEP_FRACTION * (t_end - time));
            std::vector<double> t_eval;
            std::vector<char> args(sizeof(c_PairEvolver*));
            c_PairEvolver* self_ptr = this;
            std::memcpy(args.data(), &self_ptr, sizeof(self_ptr));
            auto result = std::make_unique<CySolverResult>(this->p_method);
            this->p_result_ptr = result.get();
            this->p_error.clear();
            this->p_timeout = false;
            try {
                baseline_cysolve_ivp_noreturn(
                    result.get(), &c_PairEvolver::p_diffeq, time, t_end, y_now, C_EVOLVE_EXPECTED_STEPS,
                    C_RECORD_SIZE,  // The step quantities, as extra outputs
                    args,
                    0,              // Max number of steps: from the RAM limit
                    C_EVOLVE_MAX_RAM_MB,
                    false,          // No dense output
                    t_eval, nullptr, events, rtols, atols, std::numeric_limits<double>::infinity(), first_step,
                    false,          // The solver is not needed after the segment
                    &c_PairEvolver::p_jacobian_callback);
            } catch (const std::exception& error) {
                this->p_fail(std::string("TidalPy: the integrator failed: ") + error.what(), false);
            }
            this->p_result_ptr = nullptr;
            this->p_record.segments.push_back(entry);
            if ((result->steps_taken == 0) || (result->size == 0)) {
                message = "TidalPy: System.evolve took no step at t = " + p_format_time(time) + " Myr ("
                    + (this->p_error.empty() ? result->message : this->p_error) + ").";
                return false;
            }
            const std::string segment_error = this->p_error;
            const bool segment_timeout = this->p_timeout;
            this->p_error.clear();

            // The stored steps, each with its step quantities as extra outputs.
            const std::size_t num_stored = this->p_num_y + C_RECORD_SIZE;
            const std::size_t size = result->size;
            for (std::size_t i = (first_segment ? 0 : 1); i < size; ++i) {
                const double* stored = result->solution.data() + i * num_stored;
                this->p_append_step(result->time_domain_vec[i], stored, stored + this->p_num_y);
            }
            first_segment = false;
            this->p_event_at_start = result->event_terminated && !(result->time_domain_vec[size - 1] > time);
            stalled = (result->time_domain_vec[size - 1] > time) ? 0 : (stalled + 1);
            time = result->time_domain_vec[size - 1];
            const double* y_end = result->solution.data() + (size - 1) * num_stored;
            y_now.assign(y_end, y_end + this->p_num_y);

            if (segment_timeout) {
                message = segment_error;
                return false;
            }
            if (!segment_error.empty()) {
                // A failed evaluation (a trial state a solver rejected): retry from the last stored step, shorter.
                this->p_record.segments.back().ended = segment_error;
                this->p_last_error = segment_error;
                TIDALPY_LOG_INFO("TidalPy: System.evolve segment ended at {:.6g} Myr: {}", time, segment_error);
                first_step *= d_EVOLVE_RETRY_STEP_FACTOR;
                ++failed;
                continue;
            }
            failed = 0;
            first_step = d_EVOLVE_FIRST_STEP_FRACTION * (t_end - t_start);
            const bool fired = result->event_terminated
                && (result->event_terminate_index < this->p_event_is_crossing.size());
            const bool crossed = fired && this->p_event_is_crossing[result->event_terminate_index];
            this->p_rebase_eccentricity = fired && this->p_event_is_decay[result->event_terminate_index];
            if (size > 1) {
                first_step = std::max(first_step,
                    (crossed ? d_EVOLVE_CROSSING_STEP_FACTOR : d_EVOLVE_RESTART_STEP_FACTOR)
                    * (result->time_domain_vec[size - 1] - result->time_domain_vec[size - 2]));
            }
            if (result->event_terminated) {
                this->p_record.segments.back().ended = "event";
                continue;
            }
            if (!result->success) {
                message = "TidalPy: System.evolve failed at t = " + p_format_time(time) + " Myr: " + result->message
                    + " (CyRK status " + std::to_string(static_cast<int>(result->status)) + ").";
                return false;
            }
        }
        message = "TidalPy: System.evolve reached the end of its span.";
        return true;
    }

    // Keeps a state and the references it is measured from, as the one the run leaves the system at.
    void p_keep_last(double time, const double* y) {
        this->p_last_time = time;
        this->p_last_y.assign(y, y + this->p_num_y);
        this->p_last_e_reference = this->p_e_reference;
        for (std::size_t b = 0; b < 2; ++b) {
            this->p_last_reference[b] = {this->p_bodies[b].numerator, this->p_bodies[b].denominator};
        }
    }

    // Reports the furthest time reached [s] to the settings' progress callback, at most every
    // d_EVOLVE_PROGRESS_INTERVAL of wall time unless `force`; a false return stops the run at the next right-hand side
    // (p_check_deadline).
    void p_report_progress(double time, bool force) {
        if (this->p_settings.progress == nullptr) { return; }
        if (!force) {
            if (!(time > this->p_progress_time)) { return; }
            this->p_progress_time = time;
        }
        const auto now = std::chrono::steady_clock::now();
        if (!force && (now < this->p_next_progress)) { return; }
        this->p_next_progress = now + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
            std::chrono::duration<double>(d_EVOLVE_PROGRESS_INTERVAL));
        if (!this->p_settings.progress(this->p_settings.progress_context, time)) { this->p_stopped = true; }
    }

    // Appends a stored step from its state `y` and its step quantities `record`.
    void p_append_step(double time, const double* y, const double* record) {
        c_PairEvolutionRecord& out = this->p_record;
        // The spin frequency follows the orbital motion of the step (Kepler's third law about the total mass).
        const double orbital_motion = this->p_orbital_frequency0 * std::pow(1.0 + y[0], -1.5);
        this->p_keep_last(time, y);
        out.time.push_back(time * d_EVOLVE_TIME_UNIT);
        out.semi_major_axis.push_back((1.0 + y[0]) * this->p_semi_major_axis0);
        out.eccentricity.push_back(std::abs(this->p_e_reference + y[1]));
        out.da_dt.push_back(record[2 * C_RECORD_PER_BODY]);
        out.de_dt.push_back(record[2 * C_RECORD_PER_BODY + 1] * std::copysign(1.0, this->p_e_reference + y[1]));
        out.dn_dt.push_back(record[2 * C_RECORD_PER_BODY + 2]);
        for (std::size_t b = 0; b < 2; ++b) {
            const c_EvolveBody& body = this->p_bodies[b];
            c_BodyEvolutionRecord& body_record = out.bodies[b];
            const double* entry = record + b * C_RECORD_PER_BODY;
            body_record.spin_ratio.push_back(entry[0]);
            body_record.spin_frequency.push_back(body.rigid ? body.spin_frequency0 : entry[0] * orbital_motion);
            body_record.tidal_heating.push_back(entry[1]);
            body_record.da_dt.push_back(entry[2]);
            body_record.de_dt.push_back(entry[3]);
            body_record.dspin_dt.push_back(entry[4]);
            if (body.thermal) {
                for (std::size_t i = 0; i < body.num_layers; ++i) {
                    body_record.temperature.push_back(y[body.temperature_slot + i] * body.temperature0[i]);
                }
            }
        }
    }

    c_System*              p_system_ptr = nullptr;
    c_PairEvolveSettings   p_settings;
    ODEMethod              p_method = ODEMethod::LSODA;
    std::size_t            p_orbit_index = 0;
    double                 p_semi_major_axis0   = TidalPyConstants::d_NAN;
    double                 p_orbital_frequency0 = TidalPyConstants::d_NAN;
    double                 p_orbital_frequency  = TidalPyConstants::d_NAN;   // at the current state
    double                 p_orbit_key[2] = {TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
    double                 p_e_reference = 0.0;        // the eccentricity the segment's state is measured from
    double                 p_last_e_reference = 0.0;   // that of the last stored step
    std::pair<int, int>    p_last_reference[2];        // each body's commensurability at the last stored step
    c_EvolveBody           p_bodies[2];
    double                 p_current_offset[2] = {0.0, 0.0};
    std::vector<bool>      p_event_is_crossing;   // per event of the segment: a spin crossing its commensurability
    std::vector<bool>      p_event_is_decay;      // per event of the segment: e falling to a tenth of its reference
    bool                   p_rebase_eccentricity = false;
    bool                   p_event_at_start = false;   // the last segment ended on an event at its own start
    std::size_t            p_num_y = 2;

    CySolverResult*        p_result_ptr = nullptr;
    std::string            p_error;
    std::string            p_last_error;
    bool                   p_timeout = false;

    bool                                  p_has_deadline = false;
    std::chrono::steady_clock::time_point p_deadline;
    std::chrono::steady_clock::time_point p_next_progress;
    double                                p_progress_time = -TidalPyConstants::d_INF;   // furthest time reported [s]
    bool                                  p_stopped = false;   // the progress report asked to stop
    c_PairEvolutionRecord  p_record;
    double                 p_last_time = TidalPyConstants::d_NAN;
    std::vector<double>    p_last_y;
};

// Evolves the pair of `world_index` and its tidal host over [t_start, t_end] [s] (System.evolve). The record is shared
// with the Python results that view it.
inline std::shared_ptr<c_PairEvolutionRecord> c_evolve_pair(
        c_System* system_ptr,
        std::size_t world_index,
        double t_start,
        double t_end,
        const c_PairEvolveSettings& settings) {
    c_PairEvolver evolver(system_ptr, world_index, settings);
    return std::make_shared<c_PairEvolutionRecord>(evolver.run(t_start, t_end));
}

// An empty record of a pair (two bodies), shared, for a result rebuilt from its arrays (unpickling).
inline std::shared_ptr<c_PairEvolutionRecord> c_new_pair_evolution_record() {
    std::shared_ptr<c_PairEvolutionRecord> record = std::make_shared<c_PairEvolutionRecord>();
    record->bodies.resize(2);
    return record;
}

} // namespace tidalpy
