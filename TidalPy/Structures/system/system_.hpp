#pragma once
/*
 * system_.hpp - c_System: a gravitationally bound set of worlds.
 *
 * A system links two or more worlds (stars, planets, moons) so that tides, orbital evolution, and
 * stellar insolation can be computed between them. Two roles are tracked independently:
 *
 *   - a tidal host per world: the body that raises that world's tides. Each world names its own host (or none)
 *     and carries a two-body orbit about it (semi-major axis + eccentricity). The Moon's host is the Earth, and
 *     the Earth's can be the Moon: two worlds that host each other share one orbit, so one of them may leave its
 *     elements unset and take its partner's, and when both carry them they have to agree.
 *   - the star: the body that supplies insolation. Each world also carries a separate two-body orbit
 *     about the star, because the star need not be its tidal host. For the Earth-Moon system the Moon's tidal
 *     host is the Earth but the star is the Sun; for an exoplanet orbiting its star the star is also the
 *     tidal host and the two orbits coincide.
 *
 * A world interacts only with its tidal host and the star. The system is also the tide-state provider of its
 * worlds (c_TideStateProvider): a world it holds can ask it for the orbital state its tides are raised in.
 *
 * The system owns its worlds through shared_ptr so the Python world wrappers and the system can co-own
 * the same underlying C++ world. It also holds the orbital rate engine (c_OrbitSolver) that turns each
 * dissipating world's tidal-potential derivatives into orbital rates.
 */

#include <cmath>
#include <cstdint>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../worlds/base_.hpp"                 // c_BaseWorld
#include "../worlds/stellar_.hpp"              // c_StarWorld (for the star's luminosity in insolation)
#include "../worlds/factory_.hpp"              // c_world_from_binary (world binary-dispatch factory)
#include "../../Dynamics/orbit_solver_.hpp"    // c_OrbitSolver / c_OrbitState (orbital rate engine)
#include "../../Utilities/math/numerics_.hpp"  // c_isclose
#include "../../Utilities/conversions/conversions_.hpp"  // c_semi_a2orbital_motion, c_orbital_motion2semi_a
#include "constants_.hpp"                      // TidalPyConstants::d_EPS / d_NAN / d_PI, tidalpy_config_ptr->d_G

namespace tidalpy {

// The two members of a mutual pair describe one orbit, so a disagreement past rounding is a contradiction
// in the input.
inline constexpr double d_SHARED_ORBIT_RTOL = 1.0e-12;

// The two-body orbital elements of one world about another, used both for its orbit about the tidal host
// and for its orbit about the star. The reference body's own entry is unused. NaN leaves an element unset, so two
// element sets that describe one orbit (a mutual pair, a world hosted by the star) can each fill what the other
// leaves out (c_merge_orbit_elements); an eccentricity still unset reads as 0 (c_resolve_orbit_elements).
struct c_OrbitElements {
    double semi_major_axis = TidalPyConstants::d_NAN;   // a [m]
    double eccentricity    = TidalPyConstants::d_NAN;   // e [dimensionless]
};

// A semi-major axis is positive and an eccentricity lies in [0, 1), NaN leaving either unset: anything else is no
// bound orbit, and would reach the rate equations as NaN or a negative square root.
inline void c_check_orbit(double semi_major_axis, double eccentricity, const std::string& world_name) {
    if (!std::isnan(semi_major_axis) && !(std::isfinite(semi_major_axis) && (semi_major_axis > 0.0))) {
        throw std::invalid_argument(
            "TidalPy: world '" + world_name + "' was given a semi-major axis of " + std::to_string(semi_major_axis)
            + " m; it must be positive.");
    }
    if (!std::isnan(eccentricity) && !((eccentricity >= 0.0) && (eccentricity < 1.0))) {
        throw std::invalid_argument(
            "TidalPy: world '" + world_name + "' was given an eccentricity of " + std::to_string(eccentricity)
            + "; a bound orbit has 0 <= e < 1.");
    }
}

// The elements as the rate equations and insolation read them: an unset eccentricity is a circular orbit.
inline c_OrbitElements c_resolve_orbit_elements(const c_OrbitElements& elements) noexcept {
    c_OrbitElements resolved = elements;
    if (std::isnan(resolved.eccentricity)) { resolved.eccentricity = 0.0; }
    return resolved;
}

// One orbit from two element sets that both describe it, element by element: an element given by either set is
// kept, and one given by both must agree (to d_SHARED_ORBIT_RTOL), or `conflict` is set and `first` is returned.
inline c_OrbitElements c_merge_orbit_elements(
        const c_OrbitElements& first,
        const c_OrbitElements& second,
        bool& conflict) noexcept {
    conflict = false;
    c_OrbitElements merged = first;
    if (std::isnan(first.semi_major_axis)) {
        merged.semi_major_axis = second.semi_major_axis;
    } else if (!std::isnan(second.semi_major_axis)
               && !c_isclose(first.semi_major_axis, second.semi_major_axis, d_SHARED_ORBIT_RTOL, 0.0)) {
        conflict = true;
    }
    if (std::isnan(first.eccentricity)) {
        merged.eccentricity = second.eccentricity;
    } else if (!std::isnan(second.eccentricity)
               && !c_isclose(first.eccentricity, second.eccentricity, d_SHARED_ORBIT_RTOL, d_SHARED_ORBIT_RTOL)) {
        conflict = true;
    }
    return conflict ? first : merged;
}

// The largest spin-lock tolerance, in spin ratio: half the spacing of the half-integer spin-orbit commensurabilities,
// so a probe never passes the midpoint to the next one.
inline constexpr double d_SPIN_LOCK_TOLERANCE_MAX = 0.25;

// Spin-orbit resonance locks for the evolution calls (calc_world_evolution, calc_pair_evolution,
// calc_system_evolution). With use_locks, a dissipating body whose spin ratio s = spin / n sits within lock_tolerance
// of a stable zero of its spin balance is held there (c_System::p_hold_spin). NaN leaves lock_tolerance unset, which
// only use_locks off allows. The Python wrappers take both from the [dynamics] table of the TidalPy configuration.
struct c_SpinLockConfig {
    bool use_locks = false;
    double lock_tolerance = TidalPyConstants::d_NAN;  // [dimensionless spin ratio]
};

// Throws std::invalid_argument for a lock_tolerance that is not finite and in (0, d_SPIN_LOCK_TOLERANCE_MAX), unless
// locks are off and the tolerance is unset.
inline void c_validate_spin_lock_config(const c_SpinLockConfig& lock_config) {
    const double tolerance = lock_config.lock_tolerance;
    if (!lock_config.use_locks && std::isnan(tolerance)) { return; }
    if (!(std::isfinite(tolerance) && (tolerance > 0.0) && (tolerance < d_SPIN_LOCK_TOLERANCE_MAX))) {
        std::ostringstream message;
        message << "TidalPy: lock_tolerance must be finite and in (0, " << d_SPIN_LOCK_TOLERANCE_MAX
                << "), in units of the spin ratio (spin / mean motion); got " << tolerance << ".";
        throw std::invalid_argument(message.str());
    }
}

// The tidal, orbital, and spin rates of one orbiting world for a single tidal solve, with the state used
// and the raw tidal outputs so the energy balance can be checked. evolved is false when the world has no
// tidal host or no usable orbit about it, and the numeric fields are then unset. has_tide_model is false for a
// rigid world (no tide model attached): it raises no tide, so its rates and energy terms are zero while evolved
// stays true. A spin held by a spin-orbit resonance lock (c_SpinLockConfig) reports the combined rates of its two sides
// (c_System::p_hold_spin).
struct c_WorldEvolution {
    std::size_t world_index    = 0;
    bool        evolved        = false;
    bool        has_tide_model = false;

    // The state the solve used.
    double orbital_frequency = TidalPyConstants::d_NAN;  // mean motion n            [rad s-1]
    double semi_major_axis   = TidalPyConstants::d_NAN;  // a about the host         [m]
    double eccentricity      = 0.0;                      // e about the host         [dimensionless]
    double spin_frequency    = 0.0;                      // world spin rate          [rad s-1]
    double host_mass         = 0.0;                      // tidal host mass          [kg]
    double target_mass       = 0.0;                      // dissipating world mass   [kg]

    double tidal_heating = 0.0;  // total tidal heating                     [W]
    double dU_dM         = 0.0;  // potential derivative wrt mean anomaly   [J kg-1 rad-1]
    double dU_dw         = 0.0;  // potential derivative wrt arg pericenter [J kg-1 rad-1]
    double dU_dO         = 0.0;  // potential derivative wrt node longitude [J kg-1 rad-1]
    // dU_dM - dU_dw summed mode by mode [J kg-1 rad-1] (c_GlobalTideResult::dU_dM_minus_dw).
    double dU_dM_minus_dw = 0.0;

    // Rates.
    double da_dt    = 0.0;  // semi-major-axis rate                 [m s-1]
    double de_dt    = 0.0;  // eccentricity rate                    [s-1]
    double dn_dt    = 0.0;  // mean-motion rate                     [rad s-2]
    double dspin_dt = 0.0;  // spin rate (0 for a rigid world)      [rad s-2]

    double moment_of_inertia = TidalPyConstants::d_NAN;  // world MoI [kg m2] (NaN for a rigid world)
    bool   has_spin          = false;                    // set when the world's spin was evolved

    // The tidal heating is drawn from the orbit and the spin, so under conservation
    // energy_residual = tidal_heating + dE_orbit_dt + dE_spin_dt is about zero.
    double dE_orbit_dt     = 0.0;
    double dE_spin_dt      = 0.0;
    double energy_residual = 0.0;

    // Set when a spin-orbit resonance lock held the spin. The rates are then w (current spin) + (1 - w) (probe), w
    // being spin_lock_weight, and spin_lock_bracket is the spin ratio (spin / n) of the probe. Both are NaN for a free
    // spin.
    bool spin_locked = false;
    double spin_lock_weight = TidalPyConstants::d_NAN;   // [dimensionless]
    double spin_lock_bracket = TidalPyConstants::d_NAN;  // [dimensionless]
};

// The dual-body tidal evolution of an orbiting world and its tidal host. Both raise a tide on the shared
// orbit, so the top-level rates and energy balance are the sum of each body's single-body contribution,
// held in `world` and `host`. has_tide_model is true when at least one of the two bodies carries a tide model;
// false means both are rigid and every rate is zero.
struct c_PairEvolution {
    std::size_t world_index    = 0;      // the orbiting world
    std::size_t host_index     = 0;      // its tidal host
    bool        evolved        = false;
    bool        has_tide_model = false;  // either body can dissipate

    // Shared two-body orbit state.
    double orbital_frequency = TidalPyConstants::d_NAN;  // mean motion n [rad s-1]
    double semi_major_axis   = TidalPyConstants::d_NAN;  // a             [m]
    double eccentricity      = 0.0;                      // e             [dimensionless]

    double da_dt = 0.0;  // [m s-1]
    double de_dt = 0.0;  // [s-1]
    double dn_dt = 0.0;  // [rad s-2]

    // Each body's single-body dissipation contribution to the shared orbit.
    c_WorldEvolution world;   // the orbiting world (tide raised by the host)
    c_WorldEvolution host;    // the host (tide raised by the orbiting world)

    // Total heating against the shared-orbit energy loss and both spins.
    double tidal_heating_total = 0.0;  // world + host heating
    double dE_orbit_dt         = 0.0;  // from the combined da/dt
    double dE_spin_dt_total    = 0.0;  // both spins
    double energy_residual     = 0.0;  // heating_total + dE_orbit_dt + dE_spin_dt_total (~0 conserved)
};

class c_System : public c_TidalPyBaseClass, public c_TideStateProvider {
public:
    c_System() = default;
    explicit c_System(const std::string& name) : p_name(name) {}

    // The worlds can outlive the system, the Python wrappers co-owning them, so they stop pointing at it.
    ~c_System() override { this->p_release_worlds(); }

    // Worlds hold a pointer to the system they belong to, so it is not copied.
    c_System(const c_System&) = delete;
    c_System& operator=(const c_System&) = delete;

    const std::string& get_name() const noexcept { return this->p_name; }
    void set_name(const std::string& name) { this->p_name = name; }

    // The semi_major_axis and eccentricity here describe the world's orbit about its tidal host, named
    // afterwards with set_tidal_host, since a host may be added after the worlds it hosts; its orbit about
    // the star is set separately. NaN leaves an element unset (see c_OrbitElements). The last world added with
    // is_star is the star.
    std::size_t add_world(
            std::shared_ptr<c_BaseWorld> world,
            bool is_star = false,
            double semi_major_axis = TidalPyConstants::d_NAN,
            double eccentricity = TidalPyConstants::d_NAN) {
        if (world == nullptr) {
            throw std::invalid_argument("TidalPy: c_System::add_world - world is null");
        }
        for (const std::shared_ptr<c_BaseWorld>& member : this->p_worlds) {
            if (member == world) {
                throw std::invalid_argument(
                    "TidalPy: world '" + world->get_name() + "' is already a member of system '" + this->p_name
                    + "'.");
            }
        }
        // Worlds are found by name, so two of one name would leave the second unreachable.
        if (this->find_world_index(world->get_name()) >= 0) {
            throw std::invalid_argument(
                "TidalPy: system '" + this->p_name + "' already has a world named '" + world->get_name()
                + "'; give each world its own name.");
        }
        c_check_orbit(semi_major_axis, eccentricity, world->get_name());
        const std::size_t index = this->p_worlds.size();
        world->set_tide_state_provider(this, index);
        this->p_worlds.push_back(std::move(world));
        this->p_orbits.push_back(c_OrbitElements{semi_major_axis, eccentricity});
        this->p_stellar_orbits.push_back(c_OrbitElements{});
        this->p_host_index_byworld.push_back(-1);
        {
            const std::lock_guard<std::mutex> lock(this->p_warning_mutex);
            this->p_no_tide_model_warned.push_back(0);
        }
        if (is_star) {
            this->p_star_index = static_cast<int>(index);
        }
        return index;
    }

    std::size_t get_num_worlds() const noexcept { return this->p_worlds.size(); }

    const std::shared_ptr<c_BaseWorld>& get_world(std::size_t index) const {
        this->check_index(index);
        return this->p_worlds[index];
    }

    // Case-sensitive; -1 when no world matches.
    int find_world_index(const std::string& name) const noexcept {
        for (std::size_t i = 0; i < this->p_worlds.size(); ++i) {
            if (this->p_worlds[i]->get_name() == name) {
                return static_cast<int>(i);
            }
        }
        return -1;
    }

    bool has_tidal_host(std::size_t index) const {
        this->check_index(index);
        return this->p_host_index_byworld[index] >= 0;
    }

    // -1 for a world with none.
    int get_tidal_host_index(std::size_t index) const {
        this->check_index(index);
        return this->p_host_index_byworld[index];
    }

    void set_tidal_host(std::size_t index, std::size_t host_index) {
        this->check_index(index);
        this->check_index(host_index);
        if (index == host_index) {
            throw std::invalid_argument("TidalPy: c_System::set_tidal_host - a world cannot be its own tidal host");
        }
        const int previous_host = this->p_host_index_byworld[index];
        this->p_host_index_byworld[index] = static_cast<int>(host_index);
        try {
            this->p_reconcile_star_hosted_orbit(index);
        } catch (...) {
            this->p_host_index_byworld[index] = previous_host;
            throw;
        }
    }

    void clear_tidal_host(std::size_t index) {
        this->check_index(index);
        this->p_host_index_byworld[index] = -1;
    }

    const std::shared_ptr<c_BaseWorld>& get_tidal_host(std::size_t index) const {
        if (!this->has_tidal_host(index)) {
            throw std::runtime_error(
                "TidalPy: c_System - world '" + this->p_worlds[index]->get_name()
                + "' has no tidal host (call set_tidal_host)");
        }
        return this->p_worlds[static_cast<std::size_t>(this->p_host_index_byworld[index])];
    }

    double get_tidal_host_mass(std::size_t index) const {
        return this->get_tidal_host(index)->get_mass();
    }

    // The world and its tidal host host each other, so the two share one orbit.
    bool is_mutual_pair(std::size_t index) const {
        this->check_index(index);
        const int host_index = this->p_host_index_byworld[index];
        return (host_index >= 0)
            && (this->p_host_index_byworld[static_cast<std::size_t>(host_index)] == static_cast<int>(index));
    }

    // The star is the insolation source, and may or may not also be the tidal host.
    bool has_star() const noexcept {
        return this->p_star_index >= 0 && static_cast<std::size_t>(this->p_star_index) < this->p_worlds.size();
    }
    int get_star_index() const noexcept { return this->p_star_index; }

    void set_star(std::size_t index) {
        this->check_index(index);
        const int previous_star = this->p_star_index;
        this->p_star_index = static_cast<int>(index);
        try {
            for (std::size_t world_i = 0; world_i < this->p_worlds.size(); ++world_i) {
                this->p_reconcile_star_hosted_orbit(world_i);
            }
        } catch (...) {
            this->p_star_index = previous_star;
            throw;
        }
    }

    const std::shared_ptr<c_BaseWorld>& get_star() const {
        if (!this->has_star()) {
            throw std::runtime_error("TidalPy: c_System has no star world (call set_star / add_world with is_star=true)");
        }
        return this->p_worlds[static_cast<std::size_t>(this->p_star_index)];
    }

    double get_star_mass() const {
        return this->get_star()->get_mass();
    }

    // NaN when the star world is not a c_StarWorld.
    double get_star_luminosity() const {
        const c_StarWorld* star = dynamic_cast<const c_StarWorld*>(this->get_star().get());
        if (star == nullptr) {
            return TidalPyConstants::d_NAN;
        }
        return star->get_luminosity();
    }

    // For a world whose tidal host is the star the two element sets are one orbit, so these also set the stellar
    // elements (as the stellar setters set the tidal ones); a stale stellar copy would otherwise conflict with the
    // tidal orbit once the host changes away from the star and back. The two members of a mutual pair share one orbit
    // too, so these set the partner's elements as well; otherwise the old value the partner still carries would
    // conflict with the new one at the next read (get_host_orbit).
    void set_semi_major_axis(std::size_t index, double semi_major_axis) {
        this->check_index(index);
        c_check_orbit(semi_major_axis, 0.0, this->p_worlds[index]->get_name());
        this->p_orbits[index].semi_major_axis = semi_major_axis;
        if (this->is_hosted_by_star(index)) { this->p_stellar_orbits[index].semi_major_axis = semi_major_axis; }
        if (this->is_mutual_pair(index)) { this->p_partner_orbit(index).semi_major_axis = semi_major_axis; }
    }
    void set_eccentricity(std::size_t index, double eccentricity) {
        this->check_index(index);
        c_check_orbit(TidalPyConstants::d_NAN, eccentricity, this->p_worlds[index]->get_name());
        this->p_orbits[index].eccentricity = eccentricity;
        if (this->is_hosted_by_star(index)) { this->p_stellar_orbits[index].eccentricity = eccentricity; }
        if (this->is_mutual_pair(index)) { this->p_partner_orbit(index).eccentricity = eccentricity; }
    }

protected:
    // The tidal elements of a mutual pair's other member (is_mutual_pair).
    c_OrbitElements& p_partner_orbit(std::size_t index) {
        return this->p_orbits[static_cast<std::size_t>(this->p_host_index_byworld[index])];
    }

    // Mean motion n = sqrt(mu / a^3) [rad s-1]; NaN for a non-finite mu or a non-positive semi-major axis.
    static double p_mean_motion(double mu, double semi_major_axis) {
        if (!std::isfinite(mu) || !std::isfinite(semi_major_axis) || semi_major_axis <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return c_semi_a2orbital_motion(semi_major_axis, mu);
    }

    // A world whose tidal host is the star has one orbit, stored in both element sets. When it becomes star-hosted
    // after its elements were set (a builder sets them before the roles), the two sets merge element by element, so
    // an element given on either is kept whichever set gave the semi-major axis; an element the two give different
    // values for throws std::invalid_argument, and nothing changes, rather than one value being dropped.
    void p_reconcile_star_hosted_orbit(std::size_t index) {
        if (!this->is_hosted_by_star(index)) { return; }
        bool conflict = false;
        const c_OrbitElements merged =
            c_merge_orbit_elements(this->p_orbits[index], this->p_stellar_orbits[index], conflict);
        if (conflict) {
            throw std::invalid_argument(
                "TidalPy: c_System - world '" + this->p_worlds[index]->get_name() + "' has the star as its "
                "tidal host, so its orbit about the star is its tidal orbit, but the two sets of elements "
                "give different values for the same element; give one orbit.");
        }
        this->p_orbits[index]         = merged;
        this->p_stellar_orbits[index] = merged;
    }

public:
    // True when the world's tidal host is the star: the two orbits are then one, and the stellar elements are its
    // tidal elements.
    bool is_hosted_by_star(std::size_t index) const {
        this->check_index(index);
        return this->has_star() && (this->p_host_index_byworld[index] == this->p_star_index);
    }

    // The world's orbit about the star, the source of its insolation; an unset eccentricity reads as 0.
    c_OrbitElements get_stellar_orbit(std::size_t index) const {
        this->check_index(index);
        return this->is_hosted_by_star(index)
            ? this->get_host_orbit(index) : c_resolve_orbit_elements(this->p_stellar_orbits[index]);
    }

    // The elements of the world's orbit about its tidal host; an unset eccentricity reads as 0. The two members of
    // a mutual pair describe the same orbit, so their elements merge element by element: an element either one
    // carries is kept, and one they both carry with different values throws std::invalid_argument.
    c_OrbitElements get_host_orbit(std::size_t index) const {
        this->check_index(index);
        const c_OrbitElements& own = this->p_orbits[index];
        if (!this->is_mutual_pair(index)) {
            return c_resolve_orbit_elements(own);
        }
        const c_OrbitElements& partner =
            this->p_orbits[static_cast<std::size_t>(this->p_host_index_byworld[index])];
        bool conflict = false;
        const c_OrbitElements merged = c_merge_orbit_elements(own, partner, conflict);
        if (conflict) {
            throw std::invalid_argument(
                "TidalPy: c_System - worlds '" + this->p_worlds[index]->get_name() + "' and '"
                + this->get_tidal_host(index)->get_name()
                + "' host each other, so they share one orbit, but state different orbital elements for it. "
                  "Set the semi-major axis and eccentricity on one of them, or the same values on both.");
        }
        return c_resolve_orbit_elements(merged);
    }
    double get_semi_major_axis(std::size_t index) const { return this->get_host_orbit(index).semi_major_axis; }
    double get_eccentricity(std::size_t index) const { return this->get_host_orbit(index).eccentricity; }

    // Standard gravitational parameter mu = G (M_host + M_world) [m^3 s-2] for the world's orbit about its
    // tidal host. NaN for a world with no tidal host.
    double calc_gravitational_parameter(std::size_t index) const {
        this->check_index(index);
        if (!this->has_tidal_host(index) || tidalpy_config_ptr == nullptr) {
            return TidalPyConstants::d_NAN;
        }
        const double total_mass = this->get_tidal_host_mass(index) + this->p_worlds[index]->get_mass();
        return tidalpy_config_ptr->d_G * total_mass;
    }

    // Mean motion n = sqrt(mu / a^3) [rad s-1] for the world's two-body orbit about its tidal host.
    // Returns NaN for a non-positive/degenerate semi-major axis or a world with no tidal host.
    double calc_orbital_frequency(std::size_t index) const {
        return p_mean_motion(this->calc_gravitational_parameter(index), this->get_host_orbit(index).semi_major_axis);
    }

    // Semi-major axis a = (mu / n^2)^(1/3) [m] from a mean motion (the inverse of calc_orbital_frequency).
    // Returns NaN for a non-positive frequency or a world with no tidal host.
    double calc_semi_major_axis_from_frequency(std::size_t index, double orbital_frequency) const {
        const double mu = this->calc_gravitational_parameter(index);
        if (!std::isfinite(mu) || orbital_frequency <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return c_orbital_motion2semi_a(orbital_frequency, mu);
    }

    // Orbital elements about the star (per world, by index; the source of insolation)
    //
    // A world's orbit about the star can differ from its orbit about the tidal host. For a moon these
    // are the moon-about-planet (tidal) and the moon-about-star (roughly the planet's heliocentric
    // orbit) ellipses; for a planet whose tidal host is the star they coincide and can be set to the
    // same values. For a world whose tidal host is the star they are the same orbit: the stellar elements read
    // the tidal ones, and setting either sets both, so an evolution loop that moves the orbit moves the
    // insolation with it. The star's own entry is unused.
    void set_stellar_semi_major_axis(std::size_t index, double semi_major_axis) {
        this->check_index(index);
        c_check_orbit(semi_major_axis, 0.0, this->p_worlds[index]->get_name());
        this->p_stellar_orbits[index].semi_major_axis = semi_major_axis;
        if (this->is_hosted_by_star(index)) { this->p_orbits[index].semi_major_axis = semi_major_axis; }
    }
    void set_stellar_eccentricity(std::size_t index, double eccentricity) {
        this->check_index(index);
        c_check_orbit(TidalPyConstants::d_NAN, eccentricity, this->p_worlds[index]->get_name());
        this->p_stellar_orbits[index].eccentricity = eccentricity;
        if (this->is_hosted_by_star(index)) { this->p_orbits[index].eccentricity = eccentricity; }
    }
    double get_stellar_semi_major_axis(std::size_t index) const {
        return this->get_stellar_orbit(index).semi_major_axis;
    }
    double get_stellar_eccentricity(std::size_t index) const {
        return this->get_stellar_orbit(index).eccentricity;
    }

    // mu = G (M_star + M_world) for the world's orbit about the star; NaN for the star's own index.
    double calc_stellar_gravitational_parameter(std::size_t index) const {
        this->check_index(index);
        if (!this->has_star()) {
            throw std::runtime_error("TidalPy: c_System::calc_stellar_gravitational_parameter - no star set");
        }
        if (static_cast<int>(index) == this->p_star_index) {
            return TidalPyConstants::d_NAN;
        }
        if (tidalpy_config_ptr == nullptr) {
            return TidalPyConstants::d_NAN;
        }
        const double total_mass = this->get_star_mass() + this->p_worlds[index]->get_mass();
        return tidalpy_config_ptr->d_G * total_mass;
    }

    // n = sqrt(mu / a^3) for the world's orbit about the star.
    double calc_stellar_orbital_frequency(std::size_t index) const {
        return p_mean_motion(
            this->calc_stellar_gravitational_parameter(index), this->get_stellar_orbit(index).semi_major_axis);
    }

    // Orbit-averaged incident stellar flux [W m-2] (c_orbit_averaged_flux), with a and e the world's orbital elements
    // about the star; the incident flux, before the world's own albedo and emissivity are applied.
    double calc_insolation_flux(std::size_t index) const {
        this->check_index(index);
        if (!this->has_star()) {
            throw std::runtime_error("TidalPy: c_System::calc_insolation_flux - no star world set");
        }
        if (static_cast<int>(index) == this->p_star_index) {
            return TidalPyConstants::d_NAN;
        }
        const c_OrbitElements stellar_orbit = this->get_stellar_orbit(index);
        return c_orbit_averaged_flux(
            this->get_star_luminosity(), stellar_orbit.semi_major_axis, stellar_orbit.eccentricity);
    }

    // From stellar insolation alone, by gray-body radiative balance on the world's albedo and emissivity:
    // T = ((1-A) F / (4 eps sigma))^(1/4).
    double calc_equilibrium_temperature(std::size_t index) const {
        const double flux = this->calc_insolation_flux(index);
        if (!std::isfinite(flux)) {
            return TidalPyConstants::d_NAN;
        }
        return this->p_worlds[index]->calc_equilibrium_temperature(flux);
    }

    // What a world of this system is told about its own tides.
    bool get_tide_state(std::size_t world_index, c_TideSolveConfig& state_out) const override {
        if (world_index >= this->p_worlds.size() || !this->has_tidal_host(world_index)) {
            return false;
        }
        const double orbital_frequency = this->calc_orbital_frequency(world_index);
        if (!std::isfinite(orbital_frequency)) {
            return false;
        }
        const c_OrbitElements orbit = this->get_host_orbit(world_index);
        const c_BaseWorld* world_ptr = this->p_worlds[world_index].get();
        state_out.orbital_frequency = orbital_frequency;
        state_out.spin_frequency    = world_ptr->get_spin_frequency();
        state_out.eccentricity      = orbit.eccentricity;
        state_out.obliquity         = world_ptr->get_obliquity();
        state_out.semi_major_axis   = orbit.semi_major_axis;
        state_out.host_mass         = this->get_tidal_host_mass(world_index);
        return true;
    }

    double get_equilibrium_temperature(std::size_t world_index) const override {
        if (world_index >= this->p_worlds.size() || !this->has_star()) {
            return TidalPyConstants::d_NAN;
        }
        return this->calc_equilibrium_temperature(world_index);
    }

    // Single-body tidal dissipation: run one world's global tidal solve in the current system state, then
    // turn the tidal-potential derivatives into the orbital rates and the spin rate. Only this world raises
    // tides, its host being a point mass; calc_pair_evolution adds the host's own tide. A rigid world (no tide
    // model) comes back evolved with zero rates and has_tide_model false, and is warned about once. With locks on
    // (c_SpinLockConfig), a spin at a stable zero of its spin balance is held there (p_hold_spin).
    c_WorldEvolution calc_world_evolution(
            std::size_t index,
            const c_SpinLockConfig& lock_config = c_SpinLockConfig()) {
        this->check_index(index);
        c_validate_spin_lock_config(lock_config);
        c_WorldEvolution out;
        out.world_index    = index;
        out.has_tide_model = this->p_worlds[index]->get_tide_model_set();

        // A world with no tidal host, or no usable orbit about it, has no tide to evolve under.
        if (!this->has_tidal_host(index)) {
            return out;
        }
        const double orbital_frequency = this->calc_orbital_frequency(index);
        if (!std::isfinite(orbital_frequency)) {
            return out;
        }
        const c_OrbitElements orbit = this->get_host_orbit(index);
        out = this->calc_dissipation(
            index,
            this->get_tidal_host_mass(index),
            orbital_frequency,
            orbit.semi_major_axis,
            orbit.eccentricity);
        if (!out.has_tide_model) {
            this->p_warn_no_tide_model(index);
        }
        return this->p_hold_spin(out, lock_config, 0.0);
    }

    // Single-body dissipation for every world, in index order. The two members of a mutual pair each get a
    // row: their contributions to the orbit they share add.
    std::vector<c_WorldEvolution> calc_system_evolution(const c_SpinLockConfig& lock_config = c_SpinLockConfig()) {
        c_validate_spin_lock_config(lock_config);
        std::vector<c_WorldEvolution> results;
        results.reserve(this->p_worlds.size());
        for (std::size_t i = 0; i < this->p_worlds.size(); ++i) {
            results.push_back(this->calc_world_evolution(i, lock_config));
        }
        return results;
    }

    // Dual-body tidal evolution. Each body dissipates as a self-consistent single-body problem with the
    // other as the tide raiser, masses swapped, so the shared-orbit rates are the sum of the two and each
    // body evolves its own spin. The energy balance is the sum of the two single-body balances:
    //   heating_world + heating_host = -(dE_orbit/dt + dE_spin_world/dt + dE_spin_host/dt).
    // A body with no tide model is rigid and contributes nothing. A rigid host beside a dissipating world is a
    // normal setup (a star treated as a point mass), so only a pair in which neither body can dissipate is warned
    // about, once per world, and reports has_tide_model false. With locks on (c_SpinLockConfig), either body's spin
    // can be held (p_hold_spin); each body's hold test takes the orbit's total dn/dt with the other body's
    // contribution as solved at its own spin, a coupling of second order in the two bodies' rates.
    c_PairEvolution calc_pair_evolution(
            std::size_t index,
            const c_SpinLockConfig& lock_config = c_SpinLockConfig()) {
        this->check_index(index);
        c_validate_spin_lock_config(lock_config);
        c_PairEvolution out;
        out.world_index = index;
        if (!this->has_tidal_host(index)) {
            return out;
        }
        out.host_index = static_cast<std::size_t>(this->p_host_index_byworld[index]);
        const double orbital_frequency = this->calc_orbital_frequency(index);
        if (!std::isfinite(orbital_frequency)) {
            return out;
        }
        const c_OrbitElements orbit = this->get_host_orbit(index);
        const double a = orbit.semi_major_axis;
        const double e = orbit.eccentricity;
        const double world_mass = this->p_worlds[index]->get_mass();
        const double host_mass  = this->get_tidal_host_mass(index);

        out.orbital_frequency = orbital_frequency;
        out.semi_major_axis   = a;
        out.eccentricity      = e;

        // Each body dissipates on the shared orbit with the other body as the tide raiser.
        out.world = this->calc_dissipation(index, host_mass, orbital_frequency, a, e);
        out.host  = this->calc_dissipation(out.host_index, world_mass, orbital_frequency, a, e);
        const double world_dn_dt_free = out.world.dn_dt;
        out.world = this->p_hold_spin(out.world, lock_config, out.host.dn_dt);
        out.host  = this->p_hold_spin(out.host, lock_config, world_dn_dt_free);

        // Both are linear in each body's da/dt, so they add.
        out.da_dt = out.world.da_dt + out.host.da_dt;
        out.de_dt = out.world.de_dt + out.host.de_dt;
        out.dn_dt = out.world.dn_dt + out.host.dn_dt;
        out.tidal_heating_total = out.world.tidal_heating + out.host.tidal_heating;
        out.dE_orbit_dt         = out.world.dE_orbit_dt + out.host.dE_orbit_dt;
        out.dE_spin_dt_total    = out.world.dE_spin_dt + out.host.dE_spin_dt;
        out.energy_residual     = out.tidal_heating_total + out.dE_orbit_dt + out.dE_spin_dt_total;
        out.has_tide_model      = out.world.has_tide_model || out.host.has_tide_model;
        out.evolved             = true;
        if (!out.has_tide_model) {
            this->p_warn_no_tide_model(index);
            this->p_warn_no_tide_model(out.host_index);
        }
        return out;
    }

    // E_orbit = -G M_host M_world / (2 a), so dE_orbit/dt = G M_host M_world / (2 a^2) da/dt.
    double calc_orbital_energy_derivative(const c_WorldEvolution& evolution) const {
        if (tidalpy_config_ptr == nullptr) {
            return TidalPyConstants::d_NAN;
        }
        const double a = evolution.semi_major_axis;
        if (!std::isfinite(a) || std::abs(a) <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return tidalpy_config_ptr->d_G * evolution.host_mass * evolution.target_mass
             / (2.0 * a * a) * evolution.da_dt;
    }

    // E_spin = (1/2) I spin^2, so dE_spin/dt = I spin dspin/dt. Zero for a rigid world.
    double calc_spin_energy_derivative(const c_WorldEvolution& evolution) const noexcept {
        if (!evolution.has_spin || !std::isfinite(evolution.moment_of_inertia)) {
            return 0.0;
        }
        return evolution.moment_of_inertia * evolution.spin_frequency * evolution.dspin_dt;
    }

    // tidal_heating + dE_orbit/dt + dE_spin/dt, which conservation makes about zero: every watt dissipated
    // is drawn from the orbit and the spin.
    double calc_energy_residual(const c_WorldEvolution& evolution) const {
        return evolution.tidal_heating
             + this->calc_orbital_energy_derivative(evolution)
             + this->calc_spin_energy_derivative(evolution);
    }

    // One body's tidal-dissipation contribution to a two-body orbit: the shared primitive behind
    // calc_world_evolution, where the companion is the host, and calc_pair_evolution, which runs each body
    // in turn. A body with no tide model is rigid and contributes nothing: its result is evolved with zero rates
    // and has_tide_model false. This primitive does not warn; its callers decide when a rigid body is a problem.
    c_WorldEvolution calc_dissipation(
            std::size_t dissipator_index,
            double companion_mass,
            double orbital_frequency,
            double semi_major_axis,
            double eccentricity) {
        this->check_index(dissipator_index);
        c_WorldEvolution out;
        out.world_index = dissipator_index;

        c_BaseWorld* world_ptr      = this->p_worlds[dissipator_index].get();
        // The world's call lock is held from its tide solve through the reads of the result, the spin model, and the
        // moment of inertia, so a concurrent evolution call or setter on the same world takes its turn instead of
        // replacing the result between the solve and the read. The locked world calls below take it again, which
        // the recursive lock allows; one world is locked at a time, so a pair never waits on itself.
        const c_WorldCallLock call_lock(world_ptr->get_call_mutex());
        const double target_mass    = world_ptr->get_mass();
        const double spin_frequency = world_ptr->get_spin_frequency();

        out.orbital_frequency = orbital_frequency;
        out.semi_major_axis   = semi_major_axis;
        out.eccentricity      = eccentricity;
        out.spin_frequency    = spin_frequency;
        out.host_mass         = companion_mass;
        out.target_mass       = target_mass;

        // Rigid: no tide raised, nothing contributed.
        out.has_tide_model = world_ptr->get_tide_model_set();
        if (!out.has_tide_model) {
            out.evolved = true;
            return out;
        }

        // With the companion as the tide raiser.
        world_ptr->calc_tides(this->p_tide_state(out, world_ptr->get_obliquity(), spin_frequency));

        // calc_tides throws on failure, so the tide result is populated here, and the call lock keeps it this solve's.
        // Every world carries a spin model: it uses the EOS moment of inertia after solve_eos and its
        // moment_of_inertia_factor before (the [worlds] default of its world type unless set).
        out.moment_of_inertia = world_ptr->get_moment_of_inertia();
        this->p_fill_rates(out, world_ptr->get_tide_result(), world_ptr->get_spin_model());
        return out;
    }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::System); }

    // What load_binary reads a file into first, so a bad file never reaches this system or its worlds
    // (c_TidalPyBaseClass::make_binary_scratch).
    std::unique_ptr<c_TidalPyBaseClass> make_binary_scratch() const override {
        return std::make_unique<c_System>();
    }

protected:
    // The container state, then every world's complete binary record; p_read_payload rebuilds the
    // heterogeneous world list through c_world_from_binary, each world with its models and settings; only solved
    // state (the EOS profiles) is recomputed after load. c_OrbitSolver is stateless, so it needs no serialized
    // state.
    void p_write_payload(std::ostream& out) const override {
        const auto num_worlds = static_cast<uint64_t>(this->p_worlds.size());
        write_binary_string(out, this->p_name);
        const int32_t star_index = this->p_star_index;
        out.write(reinterpret_cast<const char*>(&star_index), sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&num_worlds), sizeof(uint64_t));
        for (const int host_index_value : this->p_host_index_byworld) {
            const int32_t host_index = host_index_value;
            out.write(reinterpret_cast<const char*>(&host_index), sizeof(int32_t));
        }
        for (const c_OrbitElements& orbit : this->p_orbits) {
            out.write(reinterpret_cast<const char*>(&orbit.semi_major_axis), sizeof(double));
            out.write(reinterpret_cast<const char*>(&orbit.eccentricity),    sizeof(double));
        }
        for (const c_OrbitElements& orbit : this->p_stellar_orbits) {
            out.write(reinterpret_cast<const char*>(&orbit.semi_major_axis), sizeof(double));
            out.write(reinterpret_cast<const char*>(&orbit.eccentricity),    sizeof(double));
        }
        for (const std::shared_ptr<c_BaseWorld>& world : this->p_worlds) {
            world->write_binary(out);
        }
        if (!out) {
            throw std::runtime_error("TidalPy: failed to write System binary data");
        }
    }

    // The whole record is read into locals and committed only once it is complete, so a corrupt or truncated file
    // throws with the system unchanged rather than half replaced. The saved roles and orbits went through the same
    // checks as add_world, set_tidal_host, and the orbit setters, so a tidal host index that names no other world,
    // an out-of-range star index, a duplicate world name, or an unbound orbit can only come from a corrupt file and
    // throws std::runtime_error rather than being dropped.
    void p_read_payload(std::istream& in, bool force) override {
        std::string name = read_binary_string(in);

        int32_t star_index = -1;
        in.read(reinterpret_cast<char*>(&star_index), sizeof(int32_t));
        uint64_t num_worlds = 0;
        in.read(reinterpret_cast<char*>(&num_worlds), sizeof(uint64_t));
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read System binary data");
        }
        check_binary_count(in, num_worlds, sizeof(int32_t) + 4 * sizeof(double), "world");

        // Per-world tidal host, then the orbital elements about it and about the star.
        std::vector<int> host_index_byworld(num_worlds, -1);
        for (uint64_t i = 0; i < num_worlds; ++i) {
            int32_t host_index = -1;
            in.read(reinterpret_cast<char*>(&host_index), sizeof(int32_t));
            host_index_byworld[i] = host_index;
        }
        std::vector<c_OrbitElements> orbits(num_worlds, c_OrbitElements{});
        for (uint64_t i = 0; i < num_worlds; ++i) {
            in.read(reinterpret_cast<char*>(&orbits[i].semi_major_axis), sizeof(double));
            in.read(reinterpret_cast<char*>(&orbits[i].eccentricity),    sizeof(double));
        }
        std::vector<c_OrbitElements> stellar_orbits(num_worlds, c_OrbitElements{});
        for (uint64_t i = 0; i < num_worlds; ++i) {
            in.read(reinterpret_cast<char*>(&stellar_orbits[i].semi_major_axis), sizeof(double));
            in.read(reinterpret_cast<char*>(&stellar_orbits[i].eccentricity),    sizeof(double));
        }
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read System binary data");
        }
        if ((star_index < -1) || ((star_index >= 0) && (static_cast<uint64_t>(star_index) >= num_worlds))) {
            throw std::runtime_error("TidalPy: corrupt System binary data: the star index is out of range");
        }
        // -1 marks a world with no tidal host; any other value must name a different world of this system.
        for (uint64_t i = 0; i < num_worlds; ++i) {
            const int host_index = host_index_byworld[i];
            const bool names_other_world = (host_index >= 0) && (static_cast<uint64_t>(host_index) < num_worlds)
                                           && (static_cast<uint64_t>(host_index) != i);
            if ((host_index != -1) && !names_other_world) {
                throw std::runtime_error(
                    "TidalPy: corrupt System binary data: world " + std::to_string(i) + " has the tidal host index "
                    + std::to_string(host_index) + ", which names no other world of this "
                    + std::to_string(num_worlds) + "-world system.");
            }
        }

        // Each world's concrete type is recovered from its own record.
        std::vector<std::shared_ptr<c_BaseWorld>> worlds;
        worlds.reserve(num_worlds);
        for (uint64_t i = 0; i < num_worlds; ++i) {
            worlds.push_back(c_world_from_binary(in, force));
        }
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read System binary data");
        }

        // Checked once the worlds are read so the errors can name them.
        for (uint64_t i = 0; i < num_worlds; ++i) {
            const std::string& world_name = worlds[i]->get_name();
            for (uint64_t j = 0; j < i; ++j) {
                if (worlds[j]->get_name() == world_name) {
                    throw std::runtime_error(
                        "TidalPy: corrupt System binary data: two worlds are named '" + world_name + "'.");
                }
            }
            try {
                c_check_orbit(orbits[i].semi_major_axis, orbits[i].eccentricity, world_name);
                c_check_orbit(stellar_orbits[i].semi_major_axis, stellar_orbits[i].eccentricity, world_name);
            } catch (const std::invalid_argument& error) {
                throw std::runtime_error(std::string("TidalPy: corrupt System binary data: ") + error.what());
            }
        }

        // Commit.
        this->p_release_worlds();
        this->p_name               = std::move(name);
        this->p_worlds             = std::move(worlds);
        this->p_orbits             = std::move(orbits);
        this->p_stellar_orbits     = std::move(stellar_orbits);
        this->p_host_index_byworld = std::move(host_index_byworld);
        this->p_star_index         = star_index;
        for (std::size_t i = 0; i < this->p_worlds.size(); ++i) {
            this->p_worlds[i]->set_tide_state_provider(this, i);
        }
        // The loaded worlds are new instances, so each can be warned about once more.
        const std::lock_guard<std::mutex> lock(this->p_warning_mutex);
        this->p_no_tide_model_warned.assign(this->p_worlds.size(), 0);
    }

    // Bounds check shared by the index-based accessors.
    void check_index(std::size_t index) const {
        if (index >= this->p_worlds.size()) {
            throw std::out_of_range("TidalPy: c_System world index out of range");
        }
    }

    // The tidal state of one dissipating body: the orbit and companion of `evolution`, the body's obliquity [rad], and
    // a spin frequency [rad s-1].
    static c_TideSolveConfig p_tide_state(
            const c_WorldEvolution& evolution,
            double obliquity,
            double spin_frequency) noexcept {
        c_TideSolveConfig state;
        state.orbital_frequency = evolution.orbital_frequency;
        state.spin_frequency    = spin_frequency;
        state.eccentricity      = evolution.eccentricity;
        state.obliquity         = obliquity;
        state.semi_major_axis   = evolution.semi_major_axis;
        state.host_mass         = evolution.host_mass;
        return state;
    }

    // The rates and energy terms of one dissipating body from one tidal solve's result. `out` carries the state the
    // solve used and the body's moment of inertia; `spin_model` is the body's.
    void p_fill_rates(c_WorldEvolution& out, const c_GlobalTideResult& tide, const c_Spin& spin_model) const {
        out.tidal_heating  = tide.tidal_heating;
        out.dU_dM          = tide.dU_dM;
        out.dU_dw          = tide.dU_dw;
        out.dU_dO          = tide.dU_dO;
        out.dU_dM_minus_dw = tide.dU_dM_minus_dw;

        // From the tidal-potential derivatives, this body being the dissipator.
        c_OrbitState orbit_state;
        orbit_state.orbital_frequency = out.orbital_frequency;
        orbit_state.semi_major_axis   = out.semi_major_axis;
        orbit_state.eccentricity      = out.eccentricity;
        orbit_state.target_mass       = out.target_mass;
        orbit_state.host_mass         = out.host_mass;
        const c_OrbitDerivatives rates =
            this->p_orbit_solver.calc_derivatives(orbit_state, out.dU_dM, out.dU_dw, tide.dU_dM_minus_dw);
        out.da_dt = rates.da_dt;
        out.de_dt = rates.de_dt;
        out.dn_dt = rates.dn_dt;

        // From this body's own spin model, under the torque from the companion.
        out.dspin_dt = spin_model.calc_dspin_dt(out.host_mass, out.dU_dO, out.moment_of_inertia);
        out.has_spin = true;

        out.dE_orbit_dt     = this->calc_orbital_energy_derivative(out);
        out.dE_spin_dt      = this->calc_spin_energy_derivative(out);
        out.energy_residual = out.tidal_heating + out.dE_orbit_dt + out.dE_spin_dt;
        out.evolved         = true;
    }

    // The spin-orbit resonance lock of one dissipating body. `free_evolution` is its calc_dissipation result at its own
    // spin, and `dn_dt_other` the other body's mean-motion rate [rad s-2] in a pair (0 for a single body). With the
    // spin ratio s = spin / n and the spin balance b = dspin/dt - s dn/dt, dn/dt being the orbit's total, one more
    // tidal solve probes the side the spin moves toward: s + lock_tolerance when b0 = b(s) > 0, s - lock_tolerance
    // when b0 < 0, giving b1. The spin is held when the two straddle zero the restoring way, b(lower) > 0 > b(upper).
    // Its balance is then B = b0^2 / (b0 - b1), the free balance b0 times the probe side's Filippov weight
    // b0 / (b0 - b1): B lies between b1 and b0, equals b0 at the outer edge of the hold (b1 -> 0), so the spin rate is
    // continuous there, and is 0 with a vanishing slope at the equilibrium (b0 -> 0), so the held spin relaxes onto the
    // equilibrium with no stiffness. Every rate and the heating are the combination w (current) + (1 - w) (probe) with
    // w b0 + (1 - w) b1 = B, and dspin/dt is the spin model's from the combined dU/dO, so the spin and the orbit take
    // one torque and dspin/dt = s dn/dt + B. Both balances are taken at the state's own s (both vector fields at one
    // point). A body at an exact zero (b0 = 0) is held on its own solve. Otherwise, with locks off, or for a rigid or
    // unevolved body, `free_evolution` comes back unchanged. The probe commits nothing to the world: its tide result,
    // layer heating, and tidal heat source stay those of its own spin.
    c_WorldEvolution p_hold_spin(
            const c_WorldEvolution& free_evolution,
            const c_SpinLockConfig& lock_config,
            double dn_dt_other) {
        const double orbital_frequency = free_evolution.orbital_frequency;
        if (!lock_config.use_locks || !free_evolution.evolved || !free_evolution.has_tide_model
                || !(orbital_frequency > 0.0)) {
            return free_evolution;
        }
        const double spin_ratio = free_evolution.spin_frequency / orbital_frequency;
        const double balance_free = free_evolution.dspin_dt - spin_ratio * (free_evolution.dn_dt + dn_dt_other);
        if (!std::isfinite(balance_free)) { return free_evolution; }
        c_WorldEvolution held = free_evolution;
        if (balance_free == 0.0) {
            held.spin_locked       = true;
            held.spin_lock_weight  = 1.0;
            held.spin_lock_bracket = spin_ratio;
            return held;
        }

        const bool probe_above = balance_free > 0.0;
        const double probe_offset = probe_above ? lock_config.lock_tolerance : -lock_config.lock_tolerance;
        const double probe_ratio = spin_ratio + probe_offset;
        c_BaseWorld* world_ptr = this->p_worlds[free_evolution.world_index].get();
        const c_WorldCallLock call_lock(world_ptr->get_call_mutex());
        c_WorldEvolution probe = free_evolution;
        const c_TideSolveOutcome outcome = world_ptr->calc_tides_probe(
            this->p_tide_state(free_evolution, world_ptr->get_obliquity(), probe_ratio * orbital_frequency));
        this->p_fill_rates(probe, outcome.tide_result, world_ptr->get_spin_model());
        const double balance_probe = probe.dspin_dt - spin_ratio * (probe.dn_dt + dn_dt_other);

        const double balance_lower = probe_above ? balance_free : balance_probe;
        const double balance_upper = probe_above ? balance_probe : balance_free;
        if (!((balance_lower > 0.0) && (balance_upper < 0.0))) { return free_evolution; }
        // The straddle makes b0 - b1 nonzero, and the held balance lies between b1 and b0.
        const double balance_spread = balance_free - balance_probe;
        const double balance_held = balance_free * balance_free / balance_spread;
        const double weight = (balance_held - balance_probe) / balance_spread;
        const auto combine = [weight](double current_value, double probe_value) {
            return weight * current_value + (1.0 - weight) * probe_value;
        };
        c_GlobalTideResult combined;
        combined.tidal_heating  = combine(free_evolution.tidal_heating, probe.tidal_heating);
        combined.dU_dM          = combine(free_evolution.dU_dM, probe.dU_dM);
        combined.dU_dw          = combine(free_evolution.dU_dw, probe.dU_dw);
        combined.dU_dO          = combine(free_evolution.dU_dO, probe.dU_dO);
        combined.dU_dM_minus_dw = combine(free_evolution.dU_dM_minus_dw, probe.dU_dM_minus_dw);
        this->p_fill_rates(held, combined, world_ptr->get_spin_model());
        held.spin_locked       = true;
        held.spin_lock_weight  = weight;
        held.spin_lock_bracket = probe_ratio;
        return held;
    }

    // Warns, once per world of this system, that a world with no tide model is rigid, so the evolution it enters
    // has zero rates. The flags are shared by concurrent evolution calls, hence the lock.
    void p_warn_no_tide_model(std::size_t index) {
        {
            const std::lock_guard<std::mutex> lock(this->p_warning_mutex);
            if (this->p_no_tide_model_warned[index] != 0) {
                return;
            }
            this->p_no_tide_model_warned[index] = 1;
        }
        TIDALPY_LOG_WARN(
            "TidalPy: world '{}' of system '{}' has no tide model, so it is rigid: it raises no tide and the "
            "evolution rates it enters are zero (has_tide_model is false in the result). Attach one with "
            "set_tide_model. Shown once per world.",
            this->p_worlds[index]->get_name(), this->p_name);
    }

    // The worlds stop asking this system for their tide state (those still pointing at it, that is: a world
    // added to a second system points at that one).
    void p_release_worlds() noexcept {
        for (const std::shared_ptr<c_BaseWorld>& world : this->p_worlds) {
            if (world && (world->get_tide_state_provider() == static_cast<const c_TideStateProvider*>(this))) {
                world->set_tide_state_provider(nullptr, 0);
            }
        }
    }

    std::string p_name;
    std::vector<std::shared_ptr<c_BaseWorld>> p_worlds;         // owned worlds (shared with the Python wrappers)
    std::vector<c_OrbitElements>              p_orbits;         // orbit about each world's tidal host
    std::vector<c_OrbitElements>              p_stellar_orbits; // orbit about the star; star entry unused
    std::vector<int> p_host_index_byworld;                      // each world's tidal host in p_worlds, or -1
    int p_star_index = -1;                                      // index into p_worlds, or -1 if unset
    c_OrbitSolver p_orbit_solver;                              // stateless engine turning dU/dX into orbital rates
    std::vector<uint8_t> p_no_tide_model_warned;                // per world: the rigid-world warning was shown
    std::mutex p_warning_mutex;                                 // guards p_no_tide_model_warned
};

} // namespace tidalpy
