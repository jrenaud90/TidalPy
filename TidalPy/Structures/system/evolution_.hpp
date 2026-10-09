#pragma once
/*
 * evolution_.hpp - the coupled orbital, spin, and thermal evolution of two worlds, one the tidal host of the other
 * (usually each other's), each spin tracking its spin-orbit equilibria (System.evolve).
 *
 * Near a stable spin-orbit equilibrium a world's spin relaxes onto it within years to kyr, while the orbit and the
 * interiors change over Myr to Gyr. Integrated directly, the spin forces an implicit integrator to steps of its
 * relaxation time, and an equilibrium narrower than the integrator's tolerance on the spin (about 1e-6 in spin ratio
 * in a cold mantle) is stepped over. Here each body's spin is *free* (a state variable) between equilibria and
 * *tracked* on one: while a stable root s*(x) of the body's spin balance b(s; x) = ds/dt exists for the slow state x
 * (the orbit and the layer temperatures), the spin sits on it and only x is integrated (a quasi-static spin;
 * Walterova and Behounkova 2020). Each tracked right-hand side finds s*(x) from the last root and returns the
 * Filippov combination of the rates on the root's two sides, so the spin and the orbit take one torque. A body with
 * no tide model is rigid: its spin frequency stays as it is.
 *
 * The two bodies are treated alike. A body's tide depends on its own spin, the shared orbit, and its interior, not on
 * the other body's spin, so a body's tide is solved again only when its own spin moves (probes are cached per slow
 * state). The two spin balances couple only through the mean-motion rate dn/dt; two tracked spins are found by
 * alternating root searches.
 *
 * Transitions are CyRK terminal events, each a function of the state alone and positive at its segment's start:
 *   - a free spin entering the band |s - k/2| < capture_band of the next commensurability, or coming within
 *     capture_margin of a stable root found ahead of it, runs a capture test: the first sign change of b in the
 *     direction the spin moves brackets the stable root it reaches;
 *   - a tracked spin's window (an interval about its root in which sample points spaced by factors of 4 show no other
 *     sign change) has an edge whose balance loses its restoring sign: the root is solved again there and the window
 *     rebuilt about it, or, with no root left nearby (the equilibrium vanished with its unstable neighbor), the spin
 *     is released.
 *
 * Assumptions:
 *   - A tracked spin sits exactly on its equilibrium; the neglected lag is its relaxation time over the evolution time.
 *   - The worlds have no permanent (triaxial) figure and no rotational flattening; the spins follow the tidal torques.
 *   - Layer boundaries are fixed; a two-phase layer melts and freezes inside its own radii (its zones move).
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "cysolve.hpp"     // CyRK: baseline_cysolve_ivp_noreturn, CySolverResult
#include "c_events.hpp"    // CyRK: Event

#include "system_.hpp"     // c_System, c_TidalDissipation, c_WorldEvolution
#include "constants_.hpp"  // TidalPyConstants, tidalpy_config_ptr
#include "../../Utilities/logging/logger_.hpp"

namespace tidalpy {

// The driver integrates in Myr, so the balances it compares are of order one.
inline constexpr double d_EVOLVE_TIME_UNIT = TidalPyConstants::d_SECONDS_PER_MYR;   // [s]
// Half the spacing of the half-integer commensurabilities k/2 [spin ratio].
inline constexpr double d_EVOLVE_HALF_SPACING = 0.25;
// The finest search offset keeps the resonant mode's frequency 2 n |s - k/2| this many times above [numerical]
// minimum_frequency, below which the mode is static.
inline constexpr double d_EVOLVE_STATIC_MARGIN = 30.0;
// Spin-ratio spacing of the four probes that sample the balance's noise at a root.
inline constexpr double d_EVOLVE_NOISE_STEP = 1.0e-13;
// Consecutive segments that end where they began, or that end on a failed evaluation, before the run stops.
inline constexpr std::size_t C_EVOLVE_MAX_STALLED_SEGMENTS = 5;
inline constexpr std::size_t C_EVOLVE_MAX_FAILED_SEGMENTS  = 10;
// Alternating root searches of two tracked spins (their balances couple only through dn/dt).
inline constexpr std::size_t C_EVOLVE_MAX_TRACKED_ROUNDS = 4;
// Relative slack of the band and capture-margin tests, so a spin an event stopped on an edge counts as inside.
inline constexpr double d_EVOLVE_EDGE_SLACK = 1.0e-6;
// The largest eccentricity a trial state is clipped to.
inline constexpr double d_EVOLVE_MAX_ECCENTRICITY = 1.0 - 1.0e-12;
// Growth per trial of a tracked root's warm-start offset, and the cap on that offset [spin ratio].
inline constexpr double d_EVOLVE_WARM_GROWTH   = 8.0;
inline constexpr double d_EVOLVE_MAX_WARM_STEP = 1.0e-2;
// Growth per trial of a window's sample distance from its root.
inline constexpr double d_EVOLVE_WINDOW_GROWTH = 4.0;
// A balance within this many times the spread of its noise samples has no reliable sign.
inline constexpr double d_EVOLVE_NOISE_FACTOR = 3.0;
// An Illinois trial keeps this fraction of the tolerance from either end of its bracket.
inline constexpr double d_EVOLVE_ILLINOIS_MARGIN = 0.25;
// After a failed evaluation a segment restarts with a first step of this fraction of the time reached, the time
// taken as at least the floor [Myr].
inline constexpr double d_EVOLVE_RESTART_STEP_FRACTION = 1.0e-6;
inline constexpr double d_EVOLVE_RESTART_TIME_FLOOR    = 1.0e-3;
// CyRK storage of a segment: the expected number of steps (storage grows past it) and the RAM cap [MB].
inline constexpr std::size_t C_EVOLVE_EXPECTED_STEPS = 256;
inline constexpr std::size_t C_EVOLVE_MAX_RAM_MB     = 2000;

// A tide solve at one spin ratio failed or returned non-finite rates.
class c_ProbeFailed : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// The run reached its wall-clock cap.
class c_EvolveTimeout : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// =====================================================================================================================
// Roots of a spin balance
// =====================================================================================================================
// The rates of a fixed slow state at one spin ratio: the balance ds/dt [Myr-1], the rates a root's two sides combine
// (Filippov), and the tidal heating [W]. `weakest` marks the weakest balance found where no root exists.
struct c_SpinProbe {
    double              spin_ratio = TidalPyConstants::d_NAN;
    double              balance    = TidalPyConstants::d_NAN;
    std::vector<double> rates;
    double              heating    = TidalPyConstants::d_NAN;
    bool                weakest    = false;
};

// An interval [lower, upper] in spin ratio about a stable root: the balance is positive at `lower`, negative at
// `upper`, and at its sample points between them it changes sign only at the root. `noise` is the balance's noise at
// the root [Myr-1].
struct c_SpinWindow {
    double lower = TidalPyConstants::d_NAN;
    double upper = TidalPyConstants::d_NAN;
    double noise = 0.0;
};

// The Filippov combination of a root's two sides: the spin ratio, the rates, and the heating.
struct c_CombinedSides {
    double              spin_ratio = TidalPyConstants::d_NAN;
    std::vector<double> rates;
    double              heating    = TidalPyConstants::d_NAN;
};

using c_SpinBracket = std::pair<c_SpinProbe, c_SpinProbe>;

// The combination of a final bracket's two sides whose weights zero the balance. A bracket of one probe (the balance
// exactly zero) is that probe.
inline c_CombinedSides c_combine_spin_sides(const c_SpinProbe& lower, const c_SpinProbe& upper) {
    c_CombinedSides out;
    if ((lower.spin_ratio == upper.spin_ratio) || (lower.balance == upper.balance)) {
        out.spin_ratio = lower.spin_ratio;
        out.rates      = lower.rates;
        out.heating    = lower.heating;
        return out;
    }
    const double weight = upper.balance / (upper.balance - lower.balance);
    out.spin_ratio = weight * lower.spin_ratio + (1.0 - weight) * upper.spin_ratio;
    out.heating    = weight * lower.heating + (1.0 - weight) * upper.heating;
    out.rates.resize(lower.rates.size());
    for (std::size_t i = 0; i < lower.rates.size(); ++i) {
        out.rates[i] = weight * lower.rates[i] + (1.0 - weight) * upper.rates[i];
    }
    return out;
}

// Shrinks a bracket below `tolerance` in spin ratio by the Illinois method, bisecting whenever two trials in a row fail
// to halve it; a trial whose tide solve fails is replaced by the midpoint once. The bracket has
// lower.spin_ratio < upper.spin_ratio and lower.balance > 0 > upper.balance. Returns the final (lower, upper), or one
// probe twice where the balance is exactly zero.
template <class Prober>
c_SpinBracket c_refine_spin_root(Prober& prober, c_SpinProbe lower, c_SpinProbe upper, double tolerance) {
    double lower_value = lower.balance;
    double upper_value = upper.balance;
    int retained = 0;
    int slow_trials = 0;
    while ((upper.spin_ratio - lower.spin_ratio) > tolerance) {
        const double width = upper.spin_ratio - lower.spin_ratio;
        double trial;
        if (slow_trials >= 2) {
            trial = lower.spin_ratio + 0.5 * width;
            slow_trials = 0;
        } else {
            trial = upper.spin_ratio - upper_value * width / (upper_value - lower_value);
            const double margin = d_EVOLVE_ILLINOIS_MARGIN * tolerance;
            trial = std::min(std::max(trial, lower.spin_ratio + margin), upper.spin_ratio - margin);
        }
        c_SpinProbe probe;
        try {
            probe = prober.probe(trial);
        } catch (const c_ProbeFailed&) {
            probe = prober.probe(lower.spin_ratio + 0.5 * width);
        }
        if (probe.balance > 0.0) {
            lower = probe;
            lower_value = probe.balance;
            if (retained == 1) { upper_value *= 0.5; }
            retained = 1;
        } else if (probe.balance < 0.0) {
            upper = probe;
            upper_value = probe.balance;
            if (retained == -1) { lower_value *= 0.5; }
            retained = -1;
        } else {
            return {probe, probe};
        }
        slow_trials = ((upper.spin_ratio - lower.spin_ratio) > 0.5 * width) ? (slow_trials + 1) : 0;
    }
    return {lower, upper};
}

// Points from `start` toward the far side of `center` (up to d_EVOLVE_HALF_SPACING past it), in the order a spin moving
// that way meets them, spaced by factors of 2 in their distance from `center` from `finest` out.
// Throws std::invalid_argument for a finest offset that is not finite and positive.
inline std::vector<double> c_spin_search_points(double start, double center, bool moving_up, double finest) {
    if (!(std::isfinite(finest) && (finest > 0.0))) {
        throw std::invalid_argument("TidalPy: the finest spin-ratio search offset must be finite and positive.");
    }
    // Truncated toward zero: an offset past the half spacing still takes one point each side.
    const int num_distances = static_cast<int>(std::trunc(std::log2(d_EVOLVE_HALF_SPACING / finest))) + 1;
    std::vector<double> points;
    points.reserve(2 * static_cast<std::size_t>(std::max(num_distances, 0)) + 1);
    for (int k = 0; k < num_distances; ++k) {
        const double distance = finest * std::ldexp(1.0, k);
        points.push_back(center - distance);
        points.push_back(center + distance);
    }
    points.push_back(center + (moving_up ? d_EVOLVE_HALF_SPACING : -d_EVOLVE_HALF_SPACING));
    std::vector<double> ahead;
    ahead.reserve(points.size());
    for (const double point : points) {
        if (moving_up ? (point > start) : (point < start)) { ahead.push_back(point); }
    }
    if (moving_up) {
        std::sort(ahead.begin(), ahead.end());
    } else {
        std::sort(ahead.begin(), ahead.end(), std::greater<double>());
    }
    return ahead;
}

// The nearest commensurability k/2 to a spin ratio (ties to even, as the search grid was laid out).
inline double c_nearest_commensurability(double spin_ratio) {
    return 0.5 * std::nearbyint(2.0 * spin_ratio);
}

// The stable root a spin at `spin_ratio` reaches at the current slow state: the first sign change of the balance in
// the direction it points, searched up to d_EVOLVE_HALF_SPACING past the nearest commensurability. Returns the refined
// bracket, or nothing. A search point whose tide solve fails is skipped.
// Assumes a pair of roots closer together than the search spacing near their location is not resolved.
template <class Prober>
std::optional<c_SpinBracket> c_locate_spin_root(Prober& prober, double spin_ratio, double tolerance, double finest) {
    const c_SpinProbe start = prober.probe(spin_ratio);
    if (start.balance == 0.0) {
        return c_SpinBracket{start, start};
    }
    const bool moving_up = start.balance > 0.0;
    c_SpinProbe previous = start;
    for (const double point : c_spin_search_points(
            spin_ratio, c_nearest_commensurability(spin_ratio), moving_up, prober.finest_offset(finest))) {
        c_SpinProbe probe;
        try {
            probe = prober.probe(point);
        } catch (const c_ProbeFailed&) {
            continue;
        }
        if ((probe.balance > 0.0) != moving_up) {
            return moving_up
                ? c_refine_spin_root(prober, previous, probe, tolerance)
                : c_refine_spin_root(prober, probe, previous, tolerance);
        }
        previous = probe;
    }
    return std::nullopt;
}

// The window about a stable root, or nothing when either side has no reliably restoring point.
//
// The net torque near a root is a small difference of large mode torques, so the balance carries the radial solve's
// error amplified by about the quality factor. Its noise is measured from four probes d_EVOLVE_NOISE_STEP apart at the
// root, and a balance within three times their spread has no reliable sign. Points step out from `resolution` past
// the root by factors of 4 (up to d_EVOLVE_HALF_SPACING) until one has a reliably wrong sign; each edge is the
// restoring point one step inside the last one, keeping a margin from the sign change beyond. A rebuilt window passes
// the edge that did not fire in `kept_lower` or `kept_upper` (NaN for none), reused while it still restores reliably.
// Points whose tide solve fails are skipped.
// Assumes a root pair between two sample points is not detected.
template <class Prober>
std::optional<c_SpinWindow> c_build_spin_window(
        Prober& prober,
        double root_ratio,
        double resolution,
        double kept_lower = TidalPyConstants::d_NAN,
        double kept_upper = TidalPyConstants::d_NAN) {
    std::vector<double> samples;
    for (int j = 0; j < 4; ++j) {
        try {
            samples.push_back(prober.probe(root_ratio + j * d_EVOLVE_NOISE_STEP).balance);
        } catch (const c_ProbeFailed&) {
            continue;
        }
    }
    double noise = 0.0;
    if (samples.size() > 1) {
        const auto [minimum, maximum] = std::minmax_element(samples.begin(), samples.end());
        noise = d_EVOLVE_NOISE_FACTOR * (*maximum - *minimum);
    }
    double edges[2] = {TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
    const double kept[2] = {kept_lower, kept_upper};
    for (int side = 0; side < 2; ++side) {
        const double direction = (side == 0) ? -1.0 : 1.0;
        if (std::isfinite(kept[side]) && (direction * (kept[side] - root_ratio) > resolution)) {
            try {
                const double balance = prober.probe(kept[side]).balance;
                if (((direction < 0.0) ? balance : -balance) > noise) {
                    edges[side] = kept[side];
                    continue;
                }
            } catch (const c_ProbeFailed&) {
            }
        }
        double distance = resolution;
        std::vector<double> restoring;
        while (distance <= d_EVOLVE_HALF_SPACING) {
            c_SpinProbe probe;
            try {
                probe = prober.probe(root_ratio + direction * distance);
            } catch (const c_ProbeFailed&) {
                distance *= d_EVOLVE_WINDOW_GROWTH;
                continue;
            }
            const double signed_balance = (direction < 0.0) ? probe.balance : -probe.balance;   // Positive: restoring
            if (signed_balance > noise) {
                restoring.push_back(probe.spin_ratio);
            } else if (signed_balance < -noise) {
                break;
            }
            distance *= d_EVOLVE_WINDOW_GROWTH;
        }
        if (restoring.empty()) {
            return std::nullopt;
        }
        edges[side] = (restoring.size() > 1) ? restoring[restoring.size() - 2] : restoring.back();
    }
    c_SpinWindow window;
    window.lower = edges[0];
    window.upper = edges[1];
    window.noise = noise;
    return window;
}

// A balance given by a caller's function (the synthetic balances of the tests): the rates are [s, balance] and the
// heating is s. A non-finite value is a failed probe.
struct c_CallbackSpinBalance {
    double (*function_ptr)(void*, double) = nullptr;
    void*       context_ptr = nullptr;
    std::size_t num_probes  = 0;

    c_SpinProbe probe(double spin_ratio) {
        const double value = this->function_ptr(this->context_ptr, spin_ratio);
        ++this->num_probes;
        if (!std::isfinite(value)) {
            throw c_ProbeFailed("TidalPy: the balance is not finite at spin ratio " + std::to_string(spin_ratio));
        }
        c_SpinProbe out;
        out.spin_ratio = spin_ratio;
        out.balance    = value;
        out.rates      = {spin_ratio, value};
        out.heating    = spin_ratio;
        return out;
    }

    double finest_offset(double finest) const noexcept { return finest; }
};

// =====================================================================================================================
// Settings and results
// =====================================================================================================================
// The settings of a pair evolution (System.evolve documents each).
struct c_PairEvolveSettings {
    bool   evolve_thermal = true;
    double orbit_rtol     = 1.0e-9;
    double thermal_rtol   = 1.0e-4;
    double spin_rtol      = 1.0e-4;
    double atol           = 1.0e-10;
    double capture_band   = 1.0e-2;
    double capture_margin = 1.0e-3;
    double root_tolerance = 1.0e-10;
    double resolution     = 1.0e-8;
    double max_wall_time  = TidalPyConstants::d_NAN;   // [s]; not finite: no cap
};

// How a body's spin is handled in a segment.
enum class c_SpinMode : int {
    Rigid    = 0,   // no tide model: the spin frequency stays as it is
    Free     = 1,   // a state variable
    Approach = 2,   // free, heading for a stable root found ahead of it
    Tracked  = 3,   // set to its equilibrium at every state
};

// One body's mode at the start of a segment, with its root and window (tracked), target (approach), or the reason it
// is free; NaN where not used.
struct c_BodySegment {
    int         mode         = static_cast<int>(c_SpinMode::Free);
    double      spin_ratio   = TidalPyConstants::d_NAN;
    double      root         = TidalPyConstants::d_NAN;
    double      window_lower = TidalPyConstants::d_NAN;
    double      window_upper = TidalPyConstants::d_NAN;
    double      target       = TidalPyConstants::d_NAN;
    std::string reason;
};

// One integration segment: its start time [s], each body's mode, and how it ended when it did not reach the end
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
    std::vector<unsigned char> tracked;          // 1 where the spin sat on an equilibrium
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
    bool                               success = false;
    std::string                        message;
    double                             elapsed = 0.0;     // wall clock [s]
};

// =====================================================================================================================
// The driver
// =====================================================================================================================
// Integrates the shared orbit, the two spins (each rigid, free, or tracked), and the layer temperatures of each body
// with layers (with evolve_thermal) over a span, in CyRK LSODA segments (see the file header). Construct, then run.
class c_PairEvolver {
public:
    // A body's balance at the current slow state, its partner's dn/dt held: the prober the root tools take.
    struct c_BodyProber {
        c_PairEvolver* evolver_ptr = nullptr;
        std::size_t    body        = 0;
        double         partner_dn  = 0.0;

        c_SpinProbe probe(double spin_ratio) {
            return this->evolver_ptr->p_probe(this->body, spin_ratio, this->partner_dn);
        }
        double finest_offset(double finest) const { return this->evolver_ptr->p_finest_offset(finest); }
    };

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
        const std::size_t partner_index = static_cast<std::size_t>(system_ptr->get_tidal_host_index(world_index));
        this->p_orbit_index = world_index;
        this->p_orbital_frequency0 = system_ptr->calc_orbital_frequency(world_index);
        if (!(std::isfinite(this->p_orbital_frequency0) && (this->p_orbital_frequency0 > 0.0))) {
            throw std::invalid_argument(
                "TidalPy: System.evolve needs a usable orbit (a positive semi-major axis) about the tidal host.");
        }
        this->p_semi_major_axis0 = system_ptr->get_semi_major_axis(world_index);
        this->p_orbital_frequency = this->p_orbital_frequency0;
        this->p_check_settings();

        const std::size_t indices[2] = {world_index, partner_index};
        std::size_t slot = 2;
        for (std::size_t b = 0; b < 2; ++b) {
            c_EvolveBody& body = this->p_bodies[b];
            body.index     = indices[b];
            body.world_ptr = system_ptr->p_worlds[indices[b]].get();
            body.rigid     = !body.world_ptr->get_tide_model_set();
            body.num_layers = body.world_ptr->get_num_layers();
            body.thermal   = settings.evolve_thermal && (body.num_layers > 0);
            const double spin_frequency = body.world_ptr->get_spin_frequency();
            body.rigid_spin_frequency = spin_frequency;
            if (!body.rigid && !(std::isfinite(spin_frequency) && (spin_frequency >= 0.0))) {
                throw std::invalid_argument(
                    "TidalPy: System.evolve needs each dissipating world spinning prograde (a finite spin_frequency of "
                    "at least zero); '" + body.world_ptr->get_name() + "' does not.");
            }
            body.spin_ratio = spin_frequency / this->p_orbital_frequency0;
            body.mode = body.rigid ? c_SpinMode::Rigid : c_SpinMode::Free;
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
        this->p_num_slow = slot;
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

        std::string message;
        const bool success = this->p_integrate(t_start / d_EVOLVE_TIME_UNIT, t_end / d_EVOLVE_TIME_UNIT, message);
        this->p_record.success = success;
        this->p_record.message = message;

        // Leave the system at the final state, each spin on its last ratio.
        if (!this->p_record.time.empty()) {
            this->p_has_deadline = false;
            try {
                this->p_set_slow_state(this->p_last_time, this->p_last_slow.data());
                for (std::size_t b = 0; b < 2; ++b) {
                    const c_EvolveBody& body = this->p_bodies[b];
                    if (!body.rigid) {
                        body.world_ptr->set_spin_frequency(
                            this->p_record.bodies[b].spin_ratio.back() * this->p_orbital_frequency);
                    }
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
        this->p_record.elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        TIDALPY_LOG_INFO(
            "TidalPy: System.evolve: {} {} steps, {} segments, {} and {} tide solves, {:.1f} s.",
            this->p_record.message, this->p_record.time.size(), this->p_record.segments.size(),
            this->p_bodies[0].num_tide_solves, this->p_bodies[1].num_tide_solves, this->p_record.elapsed);
        return std::move(this->p_record);
    }

private:
    // The probe of body `body` at `spin_ratio` and the current slow state, its partner's dn/dt `partner_dn` [rad s-2]:
    // one tide solve (cached per slow state). Throws c_ProbeFailed when the solve fails or its rates are not finite,
    // and c_EvolveTimeout past the wall-clock cap.
    c_SpinProbe p_probe(std::size_t body_index, double spin_ratio, double partner_dn) {
        this->p_check_deadline();
        c_EvolveBody& body = this->p_bodies[body_index];
        c_SpinProbe probe;
        const auto found = body.cache.find(spin_ratio);
        if (found != body.cache.end()) {
            probe = found->second;
        } else {
            probe.spin_ratio = spin_ratio;
            probe.heating    = 0.0;
            probe.rates.assign(4 + (body.thermal ? body.num_layers : 0), 0.0);
            try {
                if (!body.rigid) {
                    body.world_ptr->set_spin_frequency(spin_ratio * this->p_orbital_frequency);
                    const c_TidalDissipation dissipation = this->p_system_ptr->p_dissipation(
                        body.index, this->p_bodies[1 - body_index].index, this->p_orbit_index);
                    const c_WorldEvolution evolution = this->p_system_ptr->p_evolution(dissipation);
                    probe.rates[0] = evolution.da_dt;
                    probe.rates[1] = evolution.de_dt;
                    probe.rates[2] = evolution.dn_dt;
                    probe.rates[3] = evolution.dspin_dt;
                    probe.heating  = evolution.tidal_heating;
                    ++body.num_tide_solves;
                }
                if (body.thermal) {
                    for (std::size_t i = 0; i < body.num_layers; ++i) {
                        probe.rates[4 + i] = body.world_ptr->calc_layer_temperature_rate(i);
                    }
                }
            } catch (const std::exception& error) {
                std::string reason = error.what();
                const std::size_t last_line = reason.find_last_of('\n', reason.find_last_not_of("\n "));
                if (last_line != std::string::npos) { reason = reason.substr(last_line + 1); }
                throw c_ProbeFailed(
                    "TidalPy: the tide solve of '" + body.world_ptr->get_name() + "' failed at spin ratio "
                    + p_format(spin_ratio) + ": " + reason);
            }
            bool finite = std::isfinite(probe.heating);
            for (const double rate : probe.rates) { finite = finite && std::isfinite(rate); }
            if (!finite) {
                throw c_ProbeFailed(
                    "TidalPy: the rates of '" + body.world_ptr->get_name() + "' at spin ratio "
                    + p_format(spin_ratio) + " are not finite.");
            }
            body.cache.emplace(spin_ratio, probe);
        }
        probe.balance = body.rigid ? 0.0
            : (probe.rates[3] - spin_ratio * (probe.rates[2] + partner_dn)) / this->p_orbital_frequency
              * d_EVOLVE_TIME_UNIT;
        return probe;
    }

    // The finest spin-ratio offset searched about a commensurability: `finest`, or more where the resonant mode's
    // frequency 2 n |s - k/2| would come within d_EVOLVE_STATIC_MARGIN of [numerical] minimum_frequency.
    double p_finest_offset(double finest) const {
        const double minimum_frequency = (tidalpy_config_ptr == nullptr) ? 0.0 : tidalpy_config_ptr->d_MIN_FREQUENCY;
        return std::max(finest, d_EVOLVE_STATIC_MARGIN * minimum_frequency / (2.0 * this->p_orbital_frequency));
    }

    struct c_EvolveBody {
        std::size_t        index      = 0;          // in the system
        c_BaseWorld*       world_ptr  = nullptr;
        bool               rigid      = false;      // no tide model
        bool               thermal    = false;      // its layer temperatures are evolved
        std::size_t        num_layers = 0;
        std::size_t        temperature_slot = 0;    // its first temperature in the state
        std::vector<double> temperature0;           // [K]
        double             rigid_spin_frequency = 0.0;   // [rad s-1]

        c_SpinMode         mode       = c_SpinMode::Free;
        std::size_t        spin_slot  = 0;          // free or approach: its spin ratio in the segment state
        double             spin_ratio = TidalPyConstants::d_NAN;
        double             root       = TidalPyConstants::d_NAN;   // tracked: the last root (the warm start)
        double             step       = TidalPyConstants::d_NAN;   // tracked: the warm start's first offset
        c_SpinWindow       window;
        double             target     = TidalPyConstants::d_NAN;   // approach: the root ahead
        double             edge_spin  = TidalPyConstants::d_NAN;   // tracked: its root at the segment's start
        std::string        reason;                                 // free: why

        // After a segment: whether its mode is decided again, and a tracked window to rebuild about a found root.
        bool               redecide   = true;
        bool               recenter   = false;
        c_SpinBracket      recenter_bracket;
        double             kept_lower = TidalPyConstants::d_NAN;
        double             kept_upper = TidalPyConstants::d_NAN;

        std::unordered_map<double, c_SpinProbe> cache;   // its probes at the current slow state, by spin ratio
        // A thermal body's EOS key: the time [Myr] and the scaled temperatures it was solved at.
        bool               key_valid = false;
        double             key_time  = TidalPyConstants::d_NAN;
        std::vector<double> key_temperatures;
        c_CombinedSides    current;                      // its contribution at the last evaluation
        bool               current_weakest = false;
        std::size_t        num_tide_solves = 0;
        std::size_t        num_eos_solves  = 0;
    };

    // What a segment's event slot watches.
    enum class c_EventKind : int { BandBelow, BandAbove, Approach, EdgeLower, EdgeUpper };
    struct c_EventSlot {
        std::size_t body      = 0;
        c_EventKind kind      = c_EventKind::BandBelow;
        double      threshold = 0.0;   // band: the band edge; approach: the root
        double      side      = 1.0;   // approach: the side of the root the spin starts on
    };

    // The step quantities each right-hand side call returns as CyRK extra outputs (stored with every step, at the
    // accepted state and at an event's end): per body (spin ratio, heating, da/dt, de/dt, dspin/dt), then the summed
    // da/dt, de/dt, dn/dt.
    static constexpr std::size_t C_RECORD_PER_BODY = 5;
    static constexpr std::size_t C_RECORD_SIZE = 2 * C_RECORD_PER_BODY + 3;

    static std::string p_format(double value) {
        char buffer[64];
        std::snprintf(buffer, sizeof(buffer), "%.15f", value);
        return buffer;
    }

    static std::string p_format_time(double time) {
        char buffer[64];
        std::snprintf(buffer, sizeof(buffer), "%.10g", time);
        return buffer;
    }

    void p_check_settings() const {
        const c_PairEvolveSettings& settings = this->p_settings;
        const std::pair<const char*, double> checked[4] = {
            {"capture_band", settings.capture_band}, {"capture_margin", settings.capture_margin},
            {"root_tolerance", settings.root_tolerance}, {"resolution", settings.resolution}};
        for (const auto& [name, value] : checked) {
            if (!((value > 0.0) && (value < d_EVOLVE_HALF_SPACING))) {
                throw std::invalid_argument(
                    std::string("TidalPy: System.evolve's ") + name + " must be in (0, 0.25).");
            }
        }
        if (!((settings.root_tolerance < settings.resolution) && (settings.resolution < settings.capture_margin)
              && (settings.capture_margin <= settings.capture_band))) {
            throw std::invalid_argument(
                "TidalPy: System.evolve needs root_tolerance < resolution < capture_margin <= capture_band.");
        }
        const std::pair<const char*, double> tolerances[4] = {
            {"orbit_rtol", settings.orbit_rtol}, {"thermal_rtol", settings.thermal_rtol},
            {"spin_rtol", settings.spin_rtol}, {"atol", settings.atol}};
        for (const auto& [name, value] : tolerances) {
            if (!(std::isfinite(value) && (value > 0.0))) {
                throw std::invalid_argument(std::string("TidalPy: System.evolve's ") + name + " must be positive.");
            }
        }
    }

    void p_check_deadline() const {
        if (this->p_has_deadline && (std::chrono::steady_clock::now() > this->p_deadline)) {
            throw c_EvolveTimeout("TidalPy: System.evolve reached its max_wall_time.");
        }
    }

    // =================================================================================================================
    // The slow state
    // =================================================================================================================
    // Puts the system at slow state `slow_state` (scaled: a / a0, e, each thermal body's temperatures over their
    // initial values) at `time` [Myr]: the orbit, and for each thermal body its layer temperatures and its EOS, solved
    // with its temperature profile, its surface at its insolation temperature (no flow through it without a star or for
    // the star itself), and its heat sources at the time. Only what changed is solved again: a body's probes and EOS
    // depend on the orbit and, for a thermal body, on the time and its own temperatures, not on the other body's.
    // Throws std::runtime_error when an EOS solve fails.
    void p_set_slow_state(double time, const double* slow_state) {
        const bool orbit_changed = !(this->p_orbit_key_valid && (slow_state[0] == this->p_orbit_key[0])
                                     && (slow_state[1] == this->p_orbit_key[1]));
        if (orbit_changed) {
            this->p_orbit_key_valid = false;   // Until the orbit is set
            this->p_system_ptr->set_semi_major_axis(this->p_orbit_index, slow_state[0] * this->p_semi_major_axis0);
            this->p_system_ptr->set_eccentricity(
                this->p_orbit_index, std::min(std::max(slow_state[1], 0.0), d_EVOLVE_MAX_ECCENTRICITY));
            this->p_orbital_frequency = this->p_system_ptr->calc_orbital_frequency(this->p_orbit_index);
            this->p_orbit_key[0] = slow_state[0];
            this->p_orbit_key[1] = slow_state[1];
            this->p_orbit_key_valid = true;
        }
        for (c_EvolveBody& body : this->p_bodies) {
            if (!body.thermal) {
                if (orbit_changed) {
                    body.cache.clear();
                }
                continue;
            }
            const double* temperatures = slow_state + body.temperature_slot;
            if (!orbit_changed && body.key_valid && (time == body.key_time)
                    && std::equal(temperatures, temperatures + body.num_layers, body.key_temperatures.begin())) {
                continue;
            }
            body.key_valid = false;   // Until its solve succeeds
            body.cache.clear();
            for (std::size_t i = 0; i < body.num_layers; ++i) {
                body.world_ptr->get_layer(i)->set_temperature(temperatures[i] * body.temperature0[i]);
            }
            body.world_ptr->clear_tidal_heating();
            c_WorldEOSSolveConfig config = body.world_ptr->make_eos_solve_config();
            config.solve_temperature   = true;
            config.surface_temperature = this->p_surface_temperature(body);
            config.time                = time * d_EVOLVE_TIME_UNIT;
            const c_WorldEOSReport report = body.world_ptr->solve_eos_report(config);
            ++body.num_eos_solves;
            if (!report.success) {
                throw std::runtime_error(
                    "TidalPy: the EOS solve of world '" + body.world_ptr->get_name() + "' failed: " + report.message);
            }
            body.key_time = time;
            body.key_temperatures.assign(temperatures, temperatures + body.num_layers);
            body.key_valid = true;
        }
    }

    double p_surface_temperature(const c_EvolveBody& body) const {
        if (!this->p_system_ptr->has_star() || (static_cast<int>(body.index) == this->p_system_ptr->get_star_index())) {
            return TidalPyConstants::d_NAN;
        }
        return this->p_system_ptr->calc_equilibrium_temperature(body.index);
    }

    double p_partner_dn(std::size_t body_index) const {
        const c_EvolveBody& partner = this->p_bodies[1 - body_index];
        return partner.current.rates.empty() ? 0.0 : partner.current.rates[2];
    }

    c_BodyProber p_prober(std::size_t body_index) {
        c_BodyProber prober;
        prober.evolver_ptr = this;
        prober.body       = body_index;
        prober.partner_dn = this->p_partner_dn(body_index);
        return prober;
    }

    // Sets a body's contribution at its current spin (a rigid body's from its fixed frequency, a tracked body's at its
    // last root) at the current slow state.
    void p_refresh_current(std::size_t body_index) {
        c_EvolveBody& body = this->p_bodies[body_index];
        double spin_ratio = body.spin_ratio;
        if (body.mode == c_SpinMode::Rigid) {
            spin_ratio = body.rigid_spin_frequency / this->p_orbital_frequency;
        } else if (body.mode == c_SpinMode::Tracked) {
            spin_ratio = body.root;
        }
        const c_SpinProbe probe = this->p_probe(body_index, spin_ratio, 0.0);
        body.spin_ratio      = spin_ratio;
        body.current         = c_CombinedSides{probe.spin_ratio, probe.rates, probe.heating};
        body.current_weakest = false;
    }

    // =================================================================================================================
    // Tracked spins
    // =================================================================================================================
    // The tracked root's final bracket at the current slow state. The search steps out from the last root by factors of
    // 8, toward the side the balance points, as far as the window's edge, and refines the first sign change. A trial
    // state past an edge event (an implicit step can reach well beyond one) can have an edge with the wrong sign; the
    // search then follows the root past the window on a grid refined by factors of 2. With no root at all (past the
    // state where the equilibrium vanishes) it returns the weakest balance found, flagged, which keeps the right-hand
    // side defined; the edge event ends the step before that state.
    c_SpinBracket p_tracked_root(std::size_t body_index) {
        c_EvolveBody& body = this->p_bodies[body_index];
        c_BodyProber prober = this->p_prober(body_index);
        const double tolerance = this->p_settings.root_tolerance;
        const c_SpinProbe center = prober.probe(body.root);
        if (center.balance == 0.0) {
            return {center, center};
        }
        const bool moving_up = center.balance > 0.0;
        const double limit = moving_up ? body.window.upper : body.window.lower;
        double step = body.step;
        c_SpinProbe previous = center;
        while (true) {
            const double point = moving_up ? (body.root + step) : (body.root - step);
            const bool at_edge = moving_up ? (point >= limit) : (point <= limit);
            const c_SpinProbe probe = prober.probe(at_edge ? limit : point);
            if ((probe.balance > 0.0) != moving_up) {
                return moving_up
                    ? c_refine_spin_root(prober, previous, probe, tolerance)
                    : c_refine_spin_root(prober, probe, previous, tolerance);
            }
            if (at_edge) {
                break;
            }
            previous = probe;
            step *= d_EVOLVE_WARM_GROWTH;
        }
        previous = center;
        c_SpinProbe weakest = center;
        for (const double point : c_spin_search_points(
                body.root, body.root, moving_up, prober.finest_offset(tolerance))) {
            c_SpinProbe probe;
            try {
                probe = prober.probe(point);
            } catch (const c_ProbeFailed&) {
                continue;
            }
            if ((probe.balance > 0.0) != moving_up) {
                return moving_up
                    ? c_refine_spin_root(prober, previous, probe, tolerance)
                    : c_refine_spin_root(prober, probe, previous, tolerance);
            }
            if (std::abs(probe.balance) < std::abs(weakest.balance)) {
                weakest = probe;
            }
            previous = probe;
        }
        weakest.weakest = true;
        return {weakest, weakest};
    }

    // A tracked body's contribution at the current slow state: the Filippov combination of its root's sides. Only a
    // root inside the window becomes the next warm start, which keeps every search starting inside it.
    void p_solve_tracked(std::size_t body_index) {
        c_EvolveBody& body = this->p_bodies[body_index];
        const c_SpinBracket bracket = this->p_tracked_root(body_index);
        const c_CombinedSides combined = c_combine_spin_sides(bracket.first, bracket.second);
        if (!bracket.first.weakest && (body.window.lower < combined.spin_ratio)
                && (combined.spin_ratio < body.window.upper)) {
            // The next warm start brackets about as far as this root moved, within the tolerance and a hundredth.
            body.step = std::min(
                std::max(2.0 * std::abs(combined.spin_ratio - body.root), 4.0 * this->p_settings.root_tolerance),
                d_EVOLVE_MAX_WARM_STEP);
            body.root = combined.spin_ratio;
        }
        body.current         = combined;
        body.current_weakest = bracket.first.weakest;
        body.spin_ratio      = combined.spin_ratio;
    }

    // =================================================================================================================
    // Right-hand side
    // =================================================================================================================
    // The scaled rates `dy` of the segment state `y` at `time` [Myr]: [a / a0, e, temperatures, free spin ratios],
    // followed by the step quantities (C_RECORD_SIZE extra outputs). Each body's contribution is its probe at its spin
    // (free, approach, rigid) or its tracked combination; two tracked spins alternate root searches until neither
    // moves past root_tolerance.
    void p_evaluate(double time, const double* y, double* dy) {
        this->p_set_slow_state(time, y);
        std::vector<std::size_t> tracked;
        for (std::size_t b = 0; b < 2; ++b) {
            c_EvolveBody& body = this->p_bodies[b];
            if ((body.mode == c_SpinMode::Free) || (body.mode == c_SpinMode::Approach)) {
                body.spin_ratio = y[body.spin_slot];
            }
            this->p_refresh_current(b);
            if (body.mode == c_SpinMode::Tracked) {
                tracked.push_back(b);
            }
        }
        for (std::size_t round = 0; round < C_EVOLVE_MAX_TRACKED_ROUNDS; ++round) {
            double moved = 0.0;
            for (const std::size_t b : tracked) {
                const double before = this->p_bodies[b].spin_ratio;
                this->p_solve_tracked(b);
                moved = std::max(moved, std::abs(this->p_bodies[b].spin_ratio - before));
            }
            if ((tracked.size() < 2) || ((round > 0) && (moved <= this->p_settings.root_tolerance))) {
                break;
            }
            if ((round + 1 == C_EVOLVE_MAX_TRACKED_ROUNDS) && !this->p_round_cap_logged) {
                this->p_round_cap_logged = true;
                TIDALPY_LOG_INFO(
                    "TidalPy: System.evolve: the two tracked spins still moved {:.3g} in spin ratio after {} "
                    "alternating root searches at t = {:.10g} Myr; the last roots are used (shown once per run).",
                    moved, C_EVOLVE_MAX_TRACKED_ROUNDS, time);
            }
        }

        double da_dt = 0.0;
        double de_dt = 0.0;
        double dn_dt = 0.0;
        for (const c_EvolveBody& body : this->p_bodies) {
            da_dt += body.current.rates[0];
            de_dt += body.current.rates[1];
            dn_dt += body.current.rates[2];
        }
        dy[0] = da_dt / this->p_semi_major_axis0 * d_EVOLVE_TIME_UNIT;
        dy[1] = de_dt * d_EVOLVE_TIME_UNIT;
        for (const c_EvolveBody& body : this->p_bodies) {
            if (body.thermal) {
                for (std::size_t i = 0; i < body.num_layers; ++i) {
                    dy[body.temperature_slot + i] = body.current.rates[4 + i] / body.temperature0[i]
                        * d_EVOLVE_TIME_UNIT;
                }
            }
            if ((body.mode == c_SpinMode::Free) || (body.mode == c_SpinMode::Approach)) {
                dy[body.spin_slot] = (body.current.rates[3] - body.spin_ratio * dn_dt) / this->p_orbital_frequency
                    * d_EVOLVE_TIME_UNIT;
            }
        }

        double* record = dy + this->p_num_y;
        for (std::size_t b = 0; b < 2; ++b) {
            const c_EvolveBody& body = this->p_bodies[b];
            double* entry = record + b * C_RECORD_PER_BODY;
            entry[0] = body.spin_ratio;
            entry[1] = body.current.heating;
            entry[2] = body.current.rates[0];
            entry[3] = body.current.rates[1];
            entry[4] = body.current.rates[3];
        }
        record[2 * C_RECORD_PER_BODY]     = da_dt;
        record[2 * C_RECORD_PER_BODY + 1] = de_dt;
        record[2 * C_RECORD_PER_BODY + 2] = dn_dt;
    }

    // CyRK's view of the right-hand side. Nothing may escape the solver, so a failure is kept, stops the solve, and
    // returns zero rates.
    static void p_diffeq(double* dy, double time, double* y, char* args, PreEvalFunc) {
        c_PairEvolver* self_ptr = nullptr;
        std::memcpy(&self_ptr, args, sizeof(self_ptr));
        if (!self_ptr->p_error.empty()) {
            std::fill(dy, dy + self_ptr->p_num_y + C_RECORD_SIZE, 0.0);
            self_ptr->p_stop_solve();
            return;
        }
        try {
            self_ptr->p_evaluate(time, y, dy);
        } catch (const c_EvolveTimeout& error) {
            self_ptr->p_fail(error.what(), true);
            std::fill(dy, dy + self_ptr->p_num_y + C_RECORD_SIZE, 0.0);
        } catch (const std::exception& error) {
            self_ptr->p_fail(error.what(), false);
            std::fill(dy, dy + self_ptr->p_num_y + C_RECORD_SIZE, 0.0);
        } catch (...) {
            self_ptr->p_fail("TidalPy: System.evolve's right-hand side failed with an unknown error.", false);
            std::fill(dy, dy + self_ptr->p_num_y + C_RECORD_SIZE, 0.0);
        }
    }

    void p_fail(const std::string& message, bool timeout) {
        if (this->p_error.empty()) {
            this->p_error   = message;
            this->p_timeout = timeout;
        }
        this->p_stop_solve();
    }

    void p_stop_solve() {
        if ((this->p_result_ptr != nullptr) && this->p_result_ptr->solver_uptr) {
            this->p_result_ptr->solver_uptr->set_external_error(CyrkErrorCodes::OTHER_ERROR);
        }
    }

    // =================================================================================================================
    // Events
    // =================================================================================================================
    // The band events of a free spin at `spin_ratio`: entering the band of the next commensurability below or above, a
    // band it starts in (or on the edge of) excluded, so both are positive at the segment's start.
    void p_add_band_events(std::size_t body_index, double spin_ratio) {
        const double width = this->p_settings.capture_band;
        const double nearest = c_nearest_commensurability(spin_ratio);
        double below;
        double above;
        if (std::abs(spin_ratio - nearest) <= width * (1.0 + d_EVOLVE_EDGE_SLACK)) {
            below = nearest - 0.5;
            above = nearest + 0.5;
        } else {
            below = 0.5 * std::floor(2.0 * spin_ratio);
            above = 0.5 * std::ceil(2.0 * spin_ratio);
        }
        this->p_event_slots.push_back(c_EventSlot{body_index, c_EventKind::BandBelow, below + width, 1.0});
        this->p_event_slots.push_back(c_EventSlot{body_index, c_EventKind::BandAbove, above - width, 1.0});
    }

    // The value of one event slot at (time, y); every event is positive at its segment's start and fires crossing zero
    // downward.
    double p_event_value(const c_EventSlot& slot, double time, const double* y) {
        const c_EvolveBody& body = this->p_bodies[slot.body];
        switch (slot.kind) {
            case c_EventKind::BandBelow:
                return y[body.spin_slot] - slot.threshold;
            case c_EventKind::BandAbove:
                return slot.threshold - y[body.spin_slot];
            case c_EventKind::Approach:
                return slot.side * (y[body.spin_slot] - slot.threshold) - this->p_settings.capture_margin;
            case c_EventKind::EdgeLower:
                return this->p_edge_balance(slot.body, time, y, body.window.lower);
            case c_EventKind::EdgeUpper:
                return -this->p_edge_balance(slot.body, time, y, body.window.upper);
        }
        return 1.0;
    }

    // A tracked body's balance at a window edge at (time, y). An event must be a function of the state alone (the
    // integrator locates its root by evaluating it again), so the partner's dn/dt comes from its probe at a spin fixed
    // for the segment: its state value when free, its fixed frequency when rigid, and its root at the segment's start
    // when tracked (dn/dt barely depends on it).
    double p_edge_balance(std::size_t body_index, double time, const double* y, double edge) {
        this->p_set_slow_state(time, y);
        const c_EvolveBody& partner = this->p_bodies[1 - body_index];
        double partner_spin = partner.edge_spin;
        if ((partner.mode == c_SpinMode::Free) || (partner.mode == c_SpinMode::Approach)) {
            partner_spin = y[partner.spin_slot];
        }
        return this->p_edge_balance_at(body_index, edge, partner_spin);
    }

    // A body's balance at `edge` at the current slow state, its partner at spin ratio `partner_spin` (ignored for a
    // rigid partner, whose ratio follows from its fixed frequency).
    double p_edge_balance_at(std::size_t body_index, double edge, double partner_spin) {
        const std::size_t partner_index = 1 - body_index;
        const c_EvolveBody& partner = this->p_bodies[partner_index];
        if (partner.mode == c_SpinMode::Rigid) {
            partner_spin = partner.rigid_spin_frequency / this->p_orbital_frequency;
        }
        const double partner_dn = this->p_probe(partner_index, partner_spin, 0.0).rates[2];
        return this->p_probe(body_index, edge, partner_dn).balance;
    }

    // The partner spin ratio a body's window edges are evaluated with at the current state (see p_edge_balance).
    double p_edge_partner_spin(std::size_t body_index) const {
        const c_EvolveBody& partner = this->p_bodies[1 - body_index];
        return (partner.mode == c_SpinMode::Tracked) ? partner.root : partner.spin_ratio;
    }

    // Whether both edges of a tracked body's window restore at the current slow state (each is positive there, as an
    // event must be at its segment's start).
    bool p_window_holds(std::size_t body_index) {
        const c_EvolveBody& body = this->p_bodies[body_index];
        const double partner_spin = this->p_edge_partner_spin(body_index);
        const double lower_edge = body.window.lower;
        const double upper_edge = body.window.upper;
        return (this->p_edge_balance_at(body_index, lower_edge, partner_spin) > 0.0)
            && (-this->p_edge_balance_at(body_index, upper_edge, partner_spin) > 0.0);
    }

    std::vector<Event> p_build_events() {
        std::vector<Event> events;
        events.reserve(this->p_event_slots.size());
        for (std::size_t i = 0; i < this->p_event_slots.size(); ++i) {
            Event& event = events.emplace_back(nullptr, 1, -1);
            event.check = [this, i](Event*, double time, double* y, char*) -> double {
                if (!this->p_error.empty()) {
                    return 1.0;
                }
                try {
                    return this->p_event_value(this->p_event_slots[i], time, y);
                } catch (const c_EvolveTimeout& error) {
                    this->p_fail(error.what(), true);
                } catch (const std::exception& error) {
                    this->p_fail(error.what(), false);
                } catch (...) {
                    this->p_fail("TidalPy: a System.evolve event failed with an unknown error.", false);
                }
                return 1.0;
            };
        }
        return events;
    }

    // =================================================================================================================
    // Transitions
    // =================================================================================================================
    // The next mode of a body from the current slow state: tracked (with its root bracket in `bracket`), approach (with
    // its target), or free (with the reason).
    c_SpinMode p_decide(std::size_t body_index, c_SpinBracket& bracket) {
        c_EvolveBody& body = this->p_bodies[body_index];
        const double spin_ratio = body.spin_ratio;
        // A band event stops the spin on the band's edge, which counts as inside.
        if (std::abs(spin_ratio - c_nearest_commensurability(spin_ratio))
                > this->p_settings.capture_band * (1.0 + d_EVOLVE_EDGE_SLACK)) {
            body.reason = "outside the bands";
            return c_SpinMode::Free;
        }
        this->p_refresh_current(1 - body_index);
        c_BodyProber prober = this->p_prober(body_index);
        const std::optional<c_SpinBracket> found = c_locate_spin_root(
            prober, spin_ratio, this->p_settings.root_tolerance, this->p_settings.root_tolerance);
        if (!found) {
            body.reason = "no stable root ahead";
            return c_SpinMode::Free;
        }
        bracket = *found;
        const double root_ratio = c_combine_spin_sides(found->first, found->second).spin_ratio;
        if (std::abs(root_ratio - spin_ratio) <= this->p_settings.capture_margin * (1.0 + d_EVOLVE_EDGE_SLACK)) {
            return c_SpinMode::Tracked;   // An approach event ends on it
        }
        body.target = root_ratio;
        return c_SpinMode::Approach;
    }

    // Sets each body's mode for the next segment at slow state `slow_state` and `time` [Myr]. The bodies are decided in
    // turn, so a window built or kept for one may predate the other's decision: every window is then checked against
    // the final state of both, and a body whose window does not hold is decided again (twice at most).
    void p_decide_modes(double time, const double* slow_state) {
        this->p_set_slow_state(time, slow_state);
        for (std::size_t b = 0; b < 2; ++b) {
            c_EvolveBody& body = this->p_bodies[b];
            if (body.mode == c_SpinMode::Rigid) {
                continue;
            }
            if (body.recenter) {
                body.recenter = false;
                this->p_start_tracking(b, body.recenter_bracket, body.kept_lower, body.kept_upper);
                continue;
            }
            if (body.redecide) {
                this->p_redecide(b);
            }
        }
        for (std::size_t pass = 0; pass < 2; ++pass) {
            bool changed = false;
            for (std::size_t b = 0; b < 2; ++b) {
                if ((this->p_bodies[b].mode == c_SpinMode::Tracked) && !this->p_window_holds(b)) {
                    this->p_bodies[b].spin_ratio = this->p_bodies[b].root;
                    this->p_redecide(b);
                    changed = true;
                }
            }
            if (!changed) {
                break;
            }
        }
        for (std::size_t b = 0; b < 2; ++b) {
            c_EvolveBody& body = this->p_bodies[b];
            if ((body.mode == c_SpinMode::Tracked) && !this->p_window_holds(b)) {
                body.mode = c_SpinMode::Free;
                body.reason = "no window about the root";
            }
        }
    }

    // Decides a body's mode afresh from its current spin.
    void p_redecide(std::size_t body_index) {
        c_EvolveBody& body = this->p_bodies[body_index];
        body.reason.clear();
        body.target = TidalPyConstants::d_NAN;
        body.mode = c_SpinMode::Free;   // While deciding, the spin is the one given
        c_SpinBracket bracket;
        const c_SpinMode mode = this->p_decide(body_index, bracket);
        if (mode == c_SpinMode::Tracked) {
            this->p_start_tracking(body_index, bracket, TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
        } else {
            body.mode = mode;
        }
    }

    // Tracks a body on the root of `bracket`, with a window built about it (passing the kept edges), or leaves it free
    // when no window holds.
    void p_start_tracking(std::size_t body_index, const c_SpinBracket& bracket, double kept_lower, double kept_upper) {
        c_EvolveBody& body = this->p_bodies[body_index];
        const double root_ratio = c_combine_spin_sides(bracket.first, bracket.second).spin_ratio;
        this->p_refresh_current(1 - body_index);
        c_BodyProber prober = this->p_prober(body_index);
        const std::optional<c_SpinWindow> window = c_build_spin_window(
            prober, root_ratio, this->p_settings.resolution, kept_lower, kept_upper);
        if (!window) {
            body.mode = c_SpinMode::Free;
            body.reason = "no window about the root";
            return;
        }
        body.mode       = c_SpinMode::Tracked;
        body.root       = root_ratio;
        body.spin_ratio = root_ratio;
        body.window     = *window;
        body.step       = 4.0 * this->p_settings.root_tolerance;
    }

    // =================================================================================================================
    // Integration
    // =================================================================================================================
    bool p_integrate(double t_start, double t_end, std::string& message) {
        double time = t_start;
        std::vector<double> slow_state(this->p_num_slow, 1.0);
        slow_state[1] = this->p_system_ptr->get_eccentricity(this->p_orbit_index);
        double first_step = 0.0;   // LSODA's own choice, except after a failed evaluation
        std::size_t stalled = 0;
        std::size_t failed  = 0;
        bool first_segment = true;
        while (time < t_end) {
            if (stalled >= C_EVOLVE_MAX_STALLED_SEGMENTS) {
                message = "TidalPy: System.evolve stalled: " + std::to_string(stalled)
                    + " segments in a row ended at t = " + p_format_time(time) + " Myr.";
                return false;
            }
            if (failed >= C_EVOLVE_MAX_FAILED_SEGMENTS) {
                message = "TidalPy: System.evolve stopped at t = " + p_format_time(time) + " Myr: "
                    + std::to_string(failed) + " segments in a row ended on a failed evaluation.";
                return false;
            }
            // Modes, the state layout, and the events of the segment.
            try {
                this->p_decide_modes(time, slow_state.data());
            } catch (const std::exception& error) {   // The cap, or a solve failing at a reached state
                message = error.what();
                return false;
            }
            std::vector<double> y0 = slow_state;
            this->p_event_slots.clear();
            c_PairSegment entry;
            entry.time = time * d_EVOLVE_TIME_UNIT;
            entry.bodies.resize(2);
            std::vector<double> rtols(this->p_num_slow, this->p_settings.thermal_rtol);
            rtols[0] = this->p_settings.orbit_rtol;
            rtols[1] = this->p_settings.orbit_rtol;
            for (std::size_t b = 0; b < 2; ++b) {
                c_EvolveBody& body = this->p_bodies[b];
                c_BodySegment& body_entry = entry.bodies[b];
                body_entry.mode = static_cast<int>(body.mode);
                body_entry.spin_ratio = body.spin_ratio;
                body.redecide = true;
                if ((body.mode == c_SpinMode::Free) || (body.mode == c_SpinMode::Approach)) {
                    body.spin_slot = y0.size();
                    y0.push_back(body.spin_ratio);
                    rtols.push_back(this->p_settings.spin_rtol);
                    if (body.mode == c_SpinMode::Approach) {
                        body_entry.target = body.target;
                        this->p_event_slots.push_back(c_EventSlot{
                            b, c_EventKind::Approach, body.target,
                            std::copysign(1.0, body.spin_ratio - body.target)});
                    } else {
                        body_entry.reason = body.reason;
                    }
                    this->p_add_band_events(b, body.spin_ratio);
                } else if (body.mode == c_SpinMode::Tracked) {
                    body.edge_spin = body.root;
                    body_entry.root = body.root;
                    body_entry.window_lower = body.window.lower;
                    body_entry.window_upper = body.window.upper;
                    this->p_event_slots.push_back(c_EventSlot{b, c_EventKind::EdgeLower, 0.0, 1.0});
                    this->p_event_slots.push_back(c_EventSlot{b, c_EventKind::EdgeUpper, 0.0, 1.0});
                }
            }
            this->p_num_y = y0.size();

            // One LSODA segment. The wall-clock cap or a failed tide or EOS solve (a trial state the solvers reject)
            // ends it, its stored steps kept.
            std::vector<Event> events = this->p_build_events();
            std::vector<double> atols(1, this->p_settings.atol);
            std::vector<double> t_eval;
            std::vector<char> args(sizeof(c_PairEvolver*));
            c_PairEvolver* self_ptr = this;
            std::memcpy(args.data(), &self_ptr, sizeof(self_ptr));
            auto result = std::make_unique<CySolverResult>(ODEMethod::LSODA);
            this->p_result_ptr = result.get();
            this->p_error.clear();
            this->p_timeout = false;
            try {
                baseline_cysolve_ivp_noreturn(
                    result.get(),
                    &c_PairEvolver::p_diffeq,
                    time,
                    t_end,
                    y0,
                    C_EVOLVE_EXPECTED_STEPS,
                    C_RECORD_SIZE,  // The step quantities, as extra outputs
                    args,
                    0,              // Max number of steps: from the RAM limit
                    C_EVOLVE_MAX_RAM_MB,
                    false,          // No dense output
                    t_eval,
                    nullptr,
                    events,
                    rtols,
                    atols,
                    std::numeric_limits<double>::infinity(),
                    first_step,
                    false,          // The solver is not needed after the segment
                    nullptr);
            } catch (const std::exception& error) {
                // The integrator itself failed (its event location, for one): the segment ends on a failed
                // evaluation, its stored steps kept.
                if (this->p_error.empty()) {
                    this->p_error = std::string("TidalPy: the integrator failed: ") + error.what();
                }
            }
            this->p_result_ptr = nullptr;
            this->p_record.segments.push_back(entry);
            TIDALPY_LOG_DEBUG(
                "TidalPy: System.evolve segment at {:.6g} Myr: modes {} and {}, spin ratios {:.10f} and {:.10f}.",
                time, static_cast<int>(this->p_bodies[0].mode), static_cast<int>(this->p_bodies[1].mode),
                this->p_bodies[0].spin_ratio, this->p_bodies[1].spin_ratio);

            // A segment that took no step failed at its own start, a state the run already reached: the run ends.
            if ((result->steps_taken == 0) || (result->size == 0)) {
                message = "TidalPy: System.evolve took no step at t = " + p_format_time(time) + " Myr ("
                    + (this->p_error.empty() ? result->message : this->p_error) + ").";
                return false;
            }
            const std::string segment_error = this->p_error;
            const bool segment_timeout = this->p_timeout;
            this->p_error.clear();

            // The stored steps, each with its step quantities as extra outputs. The last step first settles the
            // tracked roots at the end state and is evaluated there again.
            const std::size_t num_y = this->p_num_y;
            const std::size_t num_stored = num_y + C_RECORD_SIZE;
            const std::size_t size = result->size;
            std::vector<double> scratch(num_stored);
            int fired_slot = -1;
            if (segment_error.empty() && result->event_terminated
                    && (result->event_terminate_index < this->p_event_slots.size())) {
                fired_slot = static_cast<int>(result->event_terminate_index);
            }
            const double segment_start = time;
            if (segment_timeout) {
                this->p_has_deadline = false;   // Its stored steps are still recorded
            }
            for (std::size_t i = 0; i < size; ++i) {
                const double step_time = result->time_domain_vec[i];
                const double* y = result->solution.data() + i * num_stored;
                const double* record = y + num_y;
                if (i + 1 == size) {
                    try {
                        this->p_settle_end_state(step_time, y, fired_slot);
                        this->p_evaluate(step_time, y, scratch.data());
                    } catch (const std::exception& error) {
                        message = error.what();
                        return false;
                    }
                    // A tracked spin ends on its last real root (kept when the end state has none).
                    for (std::size_t b = 0; b < 2; ++b) {
                        c_EvolveBody& body = this->p_bodies[b];
                        if (body.mode == c_SpinMode::Tracked) {
                            body.spin_ratio = body.root;
                            scratch[num_y + b * C_RECORD_PER_BODY] = body.root;
                        }
                    }
                    record = scratch.data() + num_y;
                }
                if (first_segment || (i > 0)) {
                    this->p_append_step(step_time, y, record);
                }
            }
            first_segment = false;
            time = result->time_domain_vec[size - 1];
            const double* y_end = result->solution.data() + (size - 1) * num_stored;
            slow_state.assign(y_end, y_end + this->p_num_slow);
            stalled = (time <= segment_start) ? (stalled + 1) : 0;

            if (segment_timeout) {
                message = segment_error;
                return false;
            }
            // A failed evaluation: retest from the last stored step, LSODA restarting with a short first step.
            if (!segment_error.empty()) {
                this->p_record.segments.back().ended = segment_error;
                TIDALPY_LOG_INFO("TidalPy: System.evolve segment ended at {:.6g} Myr: {}", time, segment_error);
                for (c_EvolveBody& body : this->p_bodies) {
                    body.redecide = true;
                    body.recenter = false;
                }
                first_step = d_EVOLVE_RESTART_STEP_FRACTION * std::max(time, d_EVOLVE_RESTART_TIME_FLOOR);
                ++failed;
                continue;
            }
            first_step = 0.0;
            failed = 0;
            if (result->event_terminated) {
                this->p_record.segments.back().ended = "event";
                continue;
            }
            if (!result->success) {
                message = "TidalPy: System.evolve failed at t = " + p_format_time(time) + " Myr: " + result->message;
                return false;
            }
        }
        message = "TidalPy: System.evolve reached the end of its span.";
        return true;
    }

    // At a segment's end state: free spins take their final value; each tracked root is solved again (the last
    // right-hand side call can be a trial state past the end), from the last real root (kept when the end state has
    // none). A tracked body whose own window edge fired, with its root still within one window width of the window,
    // has its window rebuilt about it; one whose edge did not fire keeps its window while its root is inside it. Any
    // other (a root that moved farther or vanished) is left to the capture test.
    void p_settle_end_state(double time, const double* y, int fired_slot) {
        this->p_set_slow_state(time, y);
        for (std::size_t b = 0; b < 2; ++b) {
            c_EvolveBody& body = this->p_bodies[b];
            if ((body.mode == c_SpinMode::Free) || (body.mode == c_SpinMode::Approach)) {
                body.spin_ratio = y[body.spin_slot];
            }
        }
        for (std::size_t b = 0; b < 2; ++b) {
            c_EvolveBody& body = this->p_bodies[b];
            if (body.mode != c_SpinMode::Tracked) {
                body.redecide = (body.mode != c_SpinMode::Rigid);
                continue;
            }
            const bool fired = (fired_slot >= 0) && (this->p_event_slots[fired_slot].body == b);
            body.redecide = true;
            this->p_refresh_current(1 - b);
            const c_SpinBracket bracket = this->p_tracked_root(b);
            if (bracket.first.weakest) {
                continue;
            }
            const double root_ratio = c_combine_spin_sides(bracket.first, bracket.second).spin_ratio;
            if (!fired) {
                // The window is kept only while its root is inside it.
                if ((body.window.lower < root_ratio) && (root_ratio < body.window.upper)) {
                    body.root = root_ratio;
                    body.redecide = false;
                }
                continue;
            }
            const double width = body.window.upper - body.window.lower;
            if (((body.window.lower - width) < root_ratio) && (root_ratio < (body.window.upper + width))) {
                const bool fired_lower = this->p_event_slots[fired_slot].kind == c_EventKind::EdgeLower;
                body.root = root_ratio;
                body.recenter = true;
                body.redecide = false;
                body.recenter_bracket = bracket;
                body.kept_lower = fired_lower ? TidalPyConstants::d_NAN : body.window.lower;
                body.kept_upper = fired_lower ? body.window.upper : TidalPyConstants::d_NAN;
            }
        }
    }

    // Appends a stored step from its state `y` and its step quantities `record`: the time, the orbit, each body's spin,
    // heating, contributions, tracked flag, and temperatures, and the summed rates.
    void p_append_step(double time, const double* y, const double* record) {
        c_PairEvolutionRecord& out = this->p_record;
        const double semi_major_axis = y[0] * this->p_semi_major_axis0;
        // The spin frequency follows the orbital motion of the step (Kepler's third law about the total mass).
        const double orbital_motion = this->p_orbital_frequency0 * std::pow(y[0], -1.5);
        this->p_last_time = time;
        this->p_last_slow.assign(y, y + this->p_num_slow);
        out.time.push_back(time * d_EVOLVE_TIME_UNIT);
        out.semi_major_axis.push_back(semi_major_axis);
        out.eccentricity.push_back(y[1]);
        out.da_dt.push_back(record[2 * C_RECORD_PER_BODY]);
        out.de_dt.push_back(record[2 * C_RECORD_PER_BODY + 1]);
        out.dn_dt.push_back(record[2 * C_RECORD_PER_BODY + 2]);
        for (std::size_t b = 0; b < 2; ++b) {
            const c_EvolveBody& body = this->p_bodies[b];
            c_BodyEvolutionRecord& body_record = out.bodies[b];
            const double* entry = record + b * C_RECORD_PER_BODY;
            body_record.spin_ratio.push_back(entry[0]);
            body_record.spin_frequency.push_back(
                body.rigid ? body.rigid_spin_frequency : entry[0] * orbital_motion);
            body_record.tidal_heating.push_back(entry[1]);
            body_record.da_dt.push_back(entry[2]);
            body_record.de_dt.push_back(entry[3]);
            body_record.dspin_dt.push_back(entry[4]);
            body_record.tracked.push_back(body.mode == c_SpinMode::Tracked ? 1 : 0);
            if (body.thermal) {
                for (std::size_t i = 0; i < body.num_layers; ++i) {
                    body_record.temperature.push_back(y[body.temperature_slot + i] * body.temperature0[i]);
                }
            }
        }
    }

    c_System*              p_system_ptr = nullptr;
    c_PairEvolveSettings   p_settings;
    std::size_t            p_orbit_index = 0;
    double                 p_semi_major_axis0   = TidalPyConstants::d_NAN;
    double                 p_orbital_frequency0 = TidalPyConstants::d_NAN;
    double                 p_orbital_frequency  = TidalPyConstants::d_NAN;   // at the current slow state
    c_EvolveBody           p_bodies[2];
    std::size_t            p_num_slow = 2;
    std::size_t            p_num_y    = 2;

    bool                   p_orbit_key_valid = false;
    double                 p_orbit_key[2]    = {TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
    bool                   p_round_cap_logged = false;

    std::vector<c_EventSlot>                p_event_slots;
    CySolverResult*        p_result_ptr = nullptr;
    std::string            p_error;
    bool                   p_timeout = false;

    bool                                  p_has_deadline = false;
    std::chrono::steady_clock::time_point p_deadline;
    c_PairEvolutionRecord  p_record;
    double                 p_last_time = TidalPyConstants::d_NAN;
    std::vector<double>    p_last_slow;
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

// The root tools on a balance given by a caller's function (System.evolve's tests use synthetic balances). Each
// returns its probes count in `num_probes`.
struct c_SpinRootResult {
    bool        found       = false;
    double      lower       = TidalPyConstants::d_NAN;
    double      upper       = TidalPyConstants::d_NAN;
    double      spin_ratio  = TidalPyConstants::d_NAN;   // the combination of the two sides
    double      balance     = TidalPyConstants::d_NAN;   // the combined balance (zero at a root)
    double      noise       = TidalPyConstants::d_NAN;   // a window's noise
    std::size_t num_probes  = 0;
};

inline c_SpinRootResult c_spin_root_result(const c_SpinBracket& bracket, std::size_t num_probes) {
    c_SpinRootResult out;
    const c_CombinedSides combined = c_combine_spin_sides(bracket.first, bracket.second);
    out.found      = true;
    out.lower      = bracket.first.spin_ratio;
    out.upper      = bracket.second.spin_ratio;
    out.spin_ratio = combined.spin_ratio;
    out.balance    = combined.rates.size() > 1 ? combined.rates[1] : TidalPyConstants::d_NAN;
    out.num_probes = num_probes;
    return out;
}

inline c_SpinRootResult c_locate_callback_root(
        double (*function_ptr)(void*, double),
        void* context_ptr,
        double spin_ratio,
        double tolerance,
        double finest) {
    c_CallbackSpinBalance balance{function_ptr, context_ptr, 0};
    const std::optional<c_SpinBracket> found = c_locate_spin_root(balance, spin_ratio, tolerance, finest);
    if (!found) {
        c_SpinRootResult out;
        out.num_probes = balance.num_probes;
        return out;
    }
    return c_spin_root_result(*found, balance.num_probes);
}

inline c_SpinRootResult c_refine_callback_root(
        double (*function_ptr)(void*, double),
        void* context_ptr,
        double lower,
        double upper,
        double tolerance) {
    c_CallbackSpinBalance balance{function_ptr, context_ptr, 0};
    const c_SpinProbe lower_probe = balance.probe(lower);
    const c_SpinProbe upper_probe = balance.probe(upper);
    return c_spin_root_result(c_refine_spin_root(balance, lower_probe, upper_probe, tolerance), balance.num_probes);
}

inline c_SpinRootResult c_build_callback_window(
        double (*function_ptr)(void*, double),
        void* context_ptr,
        double root_ratio,
        double resolution) {
    c_CallbackSpinBalance balance{function_ptr, context_ptr, 0};
    const std::optional<c_SpinWindow> window = c_build_spin_window(balance, root_ratio, resolution);
    c_SpinRootResult out;
    out.num_probes = balance.num_probes;
    if (window) {
        out.found = true;
        out.lower = window->lower;
        out.upper = window->upper;
        out.noise = window->noise;
    }
    return out;
}

} // namespace tidalpy
