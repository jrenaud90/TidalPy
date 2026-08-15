#pragma once
/*
 * system_.hpp - c_System: a gravitationally bound set of worlds.
 *
 * A system links two or more worlds (stars, planets, moons) so that tides, orbital evolution, and
 * stellar insolation can be computed between them. Two roles are tracked independently:
 *
 *   - the tidal host: the body a world raises tides on / orbits tidally (e.g. the Moon orbits the
 *     Earth). Each world carries a two-body orbit about the host (semi-major axis + eccentricity).
 *   - the star: the body that supplies insolation. Each world also carries a separate two-body orbit
 *     about the star, because the star need not be the tidal host. For the Earth-Moon system the tidal
 *     host is the Earth but the star is the Sun; for an exoplanet orbiting its star the star is also the
 *     tidal host and the two orbits coincide.
 *
 * Worlds do not interact with one another, only with the host / the star.
 *
 * The system owns its worlds through shared_ptr so the Python world wrappers and the system can co-own
 * the same underlying C++ world. It also holds the orbital rate engine (c_OrbitSolver) that turns each
 * dissipating world's tidal-potential derivatives into orbital rates; that wiring is added on top of
 * this container.
 *
 * All quantities MKS: masses [kg], radii [m], semi-major axes [m], eccentricities dimensionless,
 * orbital frequencies [rad s-1], gravitational parameters [m^3 s-2], insolation flux [W m-2],
 * temperatures [K].
 */

#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../worlds/base_.hpp"                // c_BaseWorld
#include "../worlds/stellar_.hpp"             // c_StarWorld (for the star's luminosity in insolation)
#include "../worlds/layered_.hpp"             // c_LayeredWorld (rheology tidal solve + spin model)
#include "../worlds/factory_.hpp"             // c_world_from_binary (world binary-dispatch factory)
#include "../../dynamics_x/orbit_solver_.hpp" // c_OrbitSolver / c_OrbitState (orbital rate engine)
#include "constants_.hpp"                     // TidalPyConstants::d_EPS / d_NAN / d_PI, tidalpy_config_ptr->d_G

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_OrbitElements - the two-body orbital elements (semi-major axis + eccentricity) of one world about
// another body. Used both for a world's orbit about the tidal host and for its orbit about the star.
// The reference body's own entry is unused.
// -------------------------------------------------------------------------------
struct c_OrbitElements {
    double semi_major_axis = TidalPyConstants::d_NAN;   // a [m]
    double eccentricity    = 0.0;                       // e [dimensionless]
};

// -------------------------------------------------------------------------------
// c_WorldEvolution - the tidal + orbital + spin rates of one orbiting world for a single tidal solve,
// together with the orbital/spin state used and the raw tidal outputs so the energy balance can be
// checked. Produced by c_System::calc_world_evolution. evolved is false when the world does not orbit
// the host (the host's own entry, or a world with no usable orbit); the numeric fields are then unset.
// -------------------------------------------------------------------------------
struct c_WorldEvolution {
    std::size_t world_index = 0;
    bool        evolved     = false;

    // Orbital + spin state used for the solve.
    double orbital_frequency = TidalPyConstants::d_NAN;  // mean motion n            [rad s-1]
    double semi_major_axis   = TidalPyConstants::d_NAN;  // a about the host         [m]
    double eccentricity      = 0.0;                      // e about the host         [dimensionless]
    double spin_frequency    = 0.0;                      // world spin rate          [rad s-1]
    double host_mass         = 0.0;                      // tidal host mass          [kg]
    double target_mass       = 0.0;                      // dissipating world mass   [kg]

    // Raw global-tidal-solve outputs for the world.
    double tidal_heating = 0.0;  // total tidal heating                     [W]
    double dU_dM         = 0.0;  // potential derivative wrt mean anomaly   [J kg-1 rad-1]
    double dU_dw         = 0.0;  // potential derivative wrt arg pericenter [J kg-1 rad-1]
    double dU_dO         = 0.0;  // potential derivative wrt node longitude [J kg-1 rad-1]

    // Rates.
    double da_dt    = 0.0;  // semi-major-axis rate                 [m s-1]
    double de_dt    = 0.0;  // eccentricity rate                    [s-1]
    double dn_dt    = 0.0;  // mean-motion rate                     [rad s-2]
    double dspin_dt = 0.0;  // spin rate (0 without a spin model)   [rad s-2]

    double moment_of_inertia = TidalPyConstants::d_NAN;  // world MoI [kg m2] (NaN without a spin model)
    bool   has_spin          = false;                    // set when the world carries a spin model

    // Energy-balance terms [W]: the tidal heating is drawn from the orbit + the spin, so
    // energy_residual = tidal_heating + dE_orbit_dt + dE_spin_dt is ~0 under conservation.
    double dE_orbit_dt     = 0.0;
    double dE_spin_dt      = 0.0;
    double energy_residual = 0.0;
};

// -------------------------------------------------------------------------------
// c_PairEvolution - the two-body (dual-body) tidal evolution of an orbiting world together with its
// tidal host. Both bodies raise a tide on the shared orbit, so the orbital rates are the sum of the two
// bodies' individual dissipation contributions, and each body evolves its own spin. Produced by
// c_System::calc_pair_evolution. `world` and `host` hold each body's full single-body contribution (its
// own solve, spin, heating, and its share of the orbital rates + energy); the combined orbital rates and
// energy balance are the top-level fields. evolved is false when there is no host / no usable orbit.
// -------------------------------------------------------------------------------
struct c_PairEvolution {
    std::size_t world_index = 0;   // the orbiting world (its p_orbits entry holds the shared orbit)
    std::size_t host_index  = 0;   // the tidal host
    bool        evolved     = false;

    // Shared two-body orbit state.
    double orbital_frequency = TidalPyConstants::d_NAN;  // mean motion n [rad s-1]
    double semi_major_axis   = TidalPyConstants::d_NAN;  // a             [m]
    double eccentricity      = 0.0;                      // e             [dimensionless]

    // Combined orbital rates (sum of both bodies' contributions).
    double da_dt = 0.0;  // [m s-1]
    double de_dt = 0.0;  // [s-1]
    double dn_dt = 0.0;  // [rad s-2]

    // Each body's single-body dissipation contribution to the shared orbit.
    c_WorldEvolution world;   // the orbiting world (tide raised by the host)
    c_WorldEvolution host;    // the host (tide raised by the orbiting world)

    // Combined energy balance [W]: total heating vs the shared-orbit energy loss + both spins.
    double tidal_heating_total = 0.0;  // world + host heating
    double dE_orbit_dt         = 0.0;  // from the combined da/dt
    double dE_spin_dt_total    = 0.0;  // both spins
    double energy_residual     = 0.0;  // heating_total + dE_orbit_dt + dE_spin_dt_total (~0 conserved)
};

// -------------------------------------------------------------------------------
// c_System
// -------------------------------------------------------------------------------
class c_System : public c_TidalPyBaseClass {
public:
    c_System() = default;
    explicit c_System(const std::string& name) : p_name(name) {}
    ~c_System() override = default;

    // -----------------------------------------------------------------------
    // Identity
    // -----------------------------------------------------------------------
    const std::string& get_name() const noexcept { return this->p_name; }
    void set_name(const std::string& name) { this->p_name = name; }

    // -----------------------------------------------------------------------
    // World membership
    //
    // Add a world to the system. is_host / is_star are independent roles (a world can be both, e.g. an
    // exoplanet's star that is also the tidal host); the last world added with each flag wins. The
    // semi_major_axis / eccentricity here describe the world's orbit about the tidal host; its orbit
    // about the star is set separately (set_stellar_semi_major_axis / set_stellar_eccentricity).
    // -----------------------------------------------------------------------
    std::size_t add_world(
            std::shared_ptr<c_BaseWorld> world,
            bool is_host = false,
            bool is_star = false,
            double semi_major_axis = TidalPyConstants::d_NAN,
            double eccentricity = 0.0) {
        if (world == nullptr) {
            throw std::invalid_argument("TidalPy: c_System::add_world - world is null");
        }
        const std::size_t index = this->p_worlds.size();
        this->p_worlds.push_back(std::move(world));
        this->p_orbits.push_back(c_OrbitElements{semi_major_axis, eccentricity});
        this->p_stellar_orbits.push_back(c_OrbitElements{});
        if (is_host) {
            this->p_host_index = static_cast<int>(index);
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

    // Index of the world whose name matches (case-sensitive), or -1 if none.
    int find_world_index(const std::string& name) const noexcept {
        for (std::size_t i = 0; i < this->p_worlds.size(); ++i) {
            if (this->p_worlds[i]->get_name() == name) {
                return static_cast<int>(i);
            }
        }
        return -1;
    }

    // -----------------------------------------------------------------------
    // Host
    // -----------------------------------------------------------------------
    bool has_host() const noexcept {
        return this->p_host_index >= 0 && static_cast<std::size_t>(this->p_host_index) < this->p_worlds.size();
    }
    int get_host_index() const noexcept { return this->p_host_index; }

    void set_host(std::size_t index) {
        this->check_index(index);
        this->p_host_index = static_cast<int>(index);
    }

    const std::shared_ptr<c_BaseWorld>& get_host() const {
        if (!this->has_host()) {
            throw std::runtime_error("TidalPy: c_System has no host world (call set_host / add_world with is_host=true)");
        }
        return this->p_worlds[static_cast<std::size_t>(this->p_host_index)];
    }

    double get_host_mass() const {
        return this->get_host()->get_mass();
    }

    // -----------------------------------------------------------------------
    // Star (the insolation source; may or may not be the tidal host)
    // -----------------------------------------------------------------------
    bool has_star() const noexcept {
        return this->p_star_index >= 0 && static_cast<std::size_t>(this->p_star_index) < this->p_worlds.size();
    }
    int get_star_index() const noexcept { return this->p_star_index; }

    void set_star(std::size_t index) {
        this->check_index(index);
        this->p_star_index = static_cast<int>(index);
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

    // The star's luminosity [W]. Returns NaN if the star world is not a c_StarWorld.
    double get_star_luminosity() const {
        const c_StarWorld* star = dynamic_cast<const c_StarWorld*>(this->get_star().get());
        if (star == nullptr) {
            return TidalPyConstants::d_NAN;
        }
        return star->get_luminosity();
    }

    // -----------------------------------------------------------------------
    // Orbital elements about the tidal host (per orbiting world, by index)
    // -----------------------------------------------------------------------
    void set_semi_major_axis(std::size_t index, double semi_major_axis) {
        this->check_index(index);
        this->p_orbits[index].semi_major_axis = semi_major_axis;
    }
    void set_eccentricity(std::size_t index, double eccentricity) {
        this->check_index(index);
        this->p_orbits[index].eccentricity = eccentricity;
    }
    double get_semi_major_axis(std::size_t index) const {
        this->check_index(index);
        return this->p_orbits[index].semi_major_axis;
    }
    double get_eccentricity(std::size_t index) const {
        this->check_index(index);
        return this->p_orbits[index].eccentricity;
    }

    // Standard gravitational parameter mu = G (M_host + M_world) [m^3 s-2] for the world's orbit about
    // the host. Throws if no host is set; returns NaN for the host's own index.
    double calc_gravitational_parameter(std::size_t index) const {
        this->check_index(index);
        if (!this->has_host()) {
            throw std::runtime_error("TidalPy: c_System::calc_gravitational_parameter - no tidal host world set");
        }
        if (static_cast<int>(index) == this->p_host_index) {
            return TidalPyConstants::d_NAN;
        }
        if (tidalpy_config_ptr == nullptr) {
            return TidalPyConstants::d_NAN;
        }
        const double total_mass = this->get_host_mass() + this->p_worlds[index]->get_mass();
        return tidalpy_config_ptr->d_G * total_mass;
    }

    // Mean motion n = sqrt(mu / a^3) [rad s-1] for the world's two-body orbit about the host.
    // Returns NaN for a non-positive/degenerate semi-major axis or the host's own index.
    double calc_orbital_frequency(std::size_t index) const {
        const double mu = this->calc_gravitational_parameter(index);
        const double semi_major_axis = this->p_orbits[index].semi_major_axis;
        if (!std::isfinite(mu) || !std::isfinite(semi_major_axis) || semi_major_axis <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return std::sqrt(mu / (semi_major_axis * semi_major_axis * semi_major_axis));
    }

    // Semi-major axis a = (mu / n^2)^(1/3) [m] from a mean motion (the inverse of calc_orbital_frequency).
    // Returns NaN for a non-positive frequency or the host's own index.
    double calc_semi_major_axis_from_frequency(std::size_t index, double orbital_frequency) const {
        const double mu = this->calc_gravitational_parameter(index);
        if (!std::isfinite(mu) || orbital_frequency <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return std::cbrt(mu / (orbital_frequency * orbital_frequency));
    }

    // -----------------------------------------------------------------------
    // Orbital elements about the star (per world, by index; the source of insolation)
    //
    // A world's orbit about the star can differ from its orbit about the tidal host. For a moon these
    // are the moon-about-planet (tidal) and the moon-about-star (roughly the planet's heliocentric
    // orbit) ellipses; for a planet whose tidal host is the star they coincide and can be set to the
    // same values. The star's own entry is unused.
    // -----------------------------------------------------------------------
    void set_stellar_semi_major_axis(std::size_t index, double semi_major_axis) {
        this->check_index(index);
        this->p_stellar_orbits[index].semi_major_axis = semi_major_axis;
    }
    void set_stellar_eccentricity(std::size_t index, double eccentricity) {
        this->check_index(index);
        this->p_stellar_orbits[index].eccentricity = eccentricity;
    }
    double get_stellar_semi_major_axis(std::size_t index) const {
        this->check_index(index);
        return this->p_stellar_orbits[index].semi_major_axis;
    }
    double get_stellar_eccentricity(std::size_t index) const {
        this->check_index(index);
        return this->p_stellar_orbits[index].eccentricity;
    }

    // Standard gravitational parameter mu = G (M_star + M_world) [m^3 s-2] for the world's orbit about
    // the star. Throws if no star is set; returns NaN for the star's own index.
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

    // Mean motion n = sqrt(mu / a^3) [rad s-1] for the world's orbit about the star.
    // Returns NaN for a non-positive/degenerate stellar semi-major axis or the star's own index.
    double calc_stellar_orbital_frequency(std::size_t index) const {
        const double mu = this->calc_stellar_gravitational_parameter(index);
        const double semi_major_axis = this->p_stellar_orbits[index].semi_major_axis;
        if (!std::isfinite(mu) || !std::isfinite(semi_major_axis) || semi_major_axis <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return std::sqrt(mu / (semi_major_axis * semi_major_axis * semi_major_axis));
    }

    // -----------------------------------------------------------------------
    // Insolation (stellar irradiation of a world)
    //
    // The orbit-averaged incident stellar flux [W m-2] at a world, F = L_star / (4 pi a^2 sqrt(1-e^2)),
    // where a, e are the world's orbital elements about the star (the sqrt(1-e^2) is the time-average of
    // 1/r^2 over an eccentric orbit; Mendez & Rivera-Valentin 2017). This is the incident flux before
    // the world's own albedo/emissivity are applied. Returns NaN if there is no star, the index is the
    // star's own, the stellar semi-major axis is unset, or the star has no luminosity.
    // -----------------------------------------------------------------------
    double calc_insolation_flux(std::size_t index) const {
        this->check_index(index);
        if (!this->has_star()) {
            throw std::runtime_error("TidalPy: c_System::calc_insolation_flux - no star world set");
        }
        if (static_cast<int>(index) == this->p_star_index) {
            return TidalPyConstants::d_NAN;
        }
        const double luminosity = this->get_star_luminosity();
        const double semi_major_axis = this->p_stellar_orbits[index].semi_major_axis;
        const double eccentricity = this->p_stellar_orbits[index].eccentricity;
        if (!std::isfinite(luminosity) || !std::isfinite(semi_major_axis)
                || semi_major_axis <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        const double ecc_factor = std::sqrt(1.0 - eccentricity * eccentricity);
        const double denom = 4.0 * TidalPyConstants::d_PI * semi_major_axis * semi_major_axis * ecc_factor;
        if (std::abs(denom) <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return luminosity / denom;
    }

    // Surface equilibrium temperature [K] of a world from stellar insolation alone (grey-body radiative
    // balance using the world's albedo + emissivity: T = ((1-A) F / (4 eps sigma))^(1/4)). Delegates the
    // radiative balance to the world. Returns NaN if the insolation flux is unavailable.
    double calc_equilibrium_temperature(std::size_t index) const {
        const double flux = this->calc_insolation_flux(index);
        if (!std::isfinite(flux)) {
            return TidalPyConstants::d_NAN;
        }
        return this->p_worlds[index]->calc_equilibrium_temperature(flux);
    }

    // -----------------------------------------------------------------------
    // Orbital + spin evolution (single-body tidal dissipation)
    //
    // For one orbiting world, run its global tidal solve in the current system state (mean motion from
    // Kepler's third law, spin + obliquity from the world, eccentricity + semi-major axis from its orbit
    // about the host, host mass from the host world), then turn the resulting tidal-potential derivatives
    // into the orbital rates (da/dt, de/dt, dn/dt) via the orbital rate engine and the spin rate
    // (dspin/dt) via the world's own spin model. Only this world raises tides; the host is treated as a
    // point mass (its own tide is added by calc_pair_evolution). The returned struct carries the state and
    // raw tidal outputs alongside the rates so the orbit + spin energy balance can be checked.
    // -----------------------------------------------------------------------
    c_WorldEvolution calc_world_evolution(std::size_t index) {
        this->check_index(index);
        c_WorldEvolution out;
        out.world_index = index;

        // The host does not orbit itself; a world without a usable orbit about the host cannot evolve.
        if (!this->has_host() || static_cast<int>(index) == this->p_host_index) {
            return out;
        }
        const double orbital_frequency = this->calc_orbital_frequency(index);
        if (!std::isfinite(orbital_frequency)) {
            return out;
        }
        return this->calc_dissipation(
            index, this->get_host_mass(), orbital_frequency,
            this->p_orbits[index].semi_major_axis, this->p_orbits[index].eccentricity);
    }

    // Evolve every world in the system (single-body dissipation), returning one c_WorldEvolution per
    // world in index order. The host's own entry and any world without a usable orbit come back with
    // evolved = false.
    std::vector<c_WorldEvolution> calc_system_evolution() {
        std::vector<c_WorldEvolution> results;
        results.reserve(this->p_worlds.size());
        for (std::size_t i = 0; i < this->p_worlds.size(); ++i) {
            results.push_back(this->calc_world_evolution(i));
        }
        return results;
    }

    // -----------------------------------------------------------------------
    // Dual-body tidal evolution
    //
    // Both the orbiting world and its tidal host raise a tide on their shared orbit. Each body dissipates
    // as a self-consistent single-body problem with the other body as the tide raiser (masses swapped),
    // so the shared-orbit rates are the sum of the two contributions and each body evolves its own spin.
    // The energy balance is the sum of the two single-body balances:
    //   heating_world + heating_host = -(dE_orbit/dt + dE_spin_world/dt + dE_spin_host/dt).
    //
    // A body with no tide model attached is rigid and contributes nothing (the single-body limit is
    // recovered when the host has no tide model). Returns evolved = false with zero rates for the host's
    // own entry, a hostless system, or a world with no usable orbit.
    // -----------------------------------------------------------------------
    c_PairEvolution calc_pair_evolution(std::size_t index) {
        this->check_index(index);
        c_PairEvolution out;
        out.world_index = index;
        if (!this->has_host() || static_cast<int>(index) == this->p_host_index) {
            return out;
        }
        out.host_index = static_cast<std::size_t>(this->p_host_index);
        const double orbital_frequency = this->calc_orbital_frequency(index);
        if (!std::isfinite(orbital_frequency)) {
            return out;
        }
        const double a = this->p_orbits[index].semi_major_axis;
        const double e = this->p_orbits[index].eccentricity;
        const double world_mass = this->p_worlds[index]->get_mass();
        const double host_mass  = this->get_host_mass();

        out.orbital_frequency = orbital_frequency;
        out.semi_major_axis   = a;
        out.eccentricity      = e;

        // Each body dissipates on the shared orbit with the other body as the tide raiser.
        out.world = this->calc_dissipation(index, host_mass, orbital_frequency, a, e);
        out.host  = this->calc_dissipation(out.host_index, world_mass, orbital_frequency, a, e);

        // Combined orbital rates + energy (both linear in each body's da/dt, so they add).
        out.da_dt = out.world.da_dt + out.host.da_dt;
        out.de_dt = out.world.de_dt + out.host.de_dt;
        out.dn_dt = out.world.dn_dt + out.host.dn_dt;
        out.tidal_heating_total = out.world.tidal_heating + out.host.tidal_heating;
        out.dE_orbit_dt         = out.world.dE_orbit_dt + out.host.dE_orbit_dt;
        out.dE_spin_dt_total    = out.world.dE_spin_dt + out.host.dE_spin_dt;
        out.energy_residual     = out.tidal_heating_total + out.dE_orbit_dt + out.dE_spin_dt_total;
        out.evolved             = true;
        return out;
    }

    // Rate of change of a world's Keplerian two-body orbital energy [W] from its semi-major-axis rate:
    //   E_orbit = -G M_host M_world / (2 a)   ->   dE_orbit/dt = G M_host M_world / (2 a^2) da/dt.
    // Returns NaN for an unusable semi-major axis or a null config pointer.
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

    // Rate of change of a world's rotational (spin) energy [W]:
    //   E_spin = (1/2) I spin^2   ->   dE_spin/dt = I spin dspin/dt.
    // Returns 0 for a world with no spin model.
    double calc_spin_energy_derivative(const c_WorldEvolution& evolution) const noexcept {
        if (!evolution.has_spin || !std::isfinite(evolution.moment_of_inertia)) {
            return 0.0;
        }
        return evolution.moment_of_inertia * evolution.spin_frequency * evolution.dspin_dt;
    }

    // Energy-balance residual [W] for a single evolved world: tidal_heating + dE_orbit/dt + dE_spin/dt.
    // Conservation (all dissipated tidal power is drawn from the orbit + spin) makes this ~0.
    double calc_energy_residual(const c_WorldEvolution& evolution) const {
        return evolution.tidal_heating
             + this->calc_orbital_energy_derivative(evolution)
             + this->calc_spin_energy_derivative(evolution);
    }

    // Compute one body's tidal-dissipation contribution to a two-body orbit: solve its global tides for
    // the shared orbital state (its own spin + obliquity, the companion mass as the tide raiser), then its
    // orbital rates (as the dissipating body) via the rate engine and its spin rate via its own spin
    // model. A body with no tide model attached is rigid: it raises no tide and contributes zero rates /
    // heating / spin. This is the shared primitive behind calc_world_evolution (companion = the host) and
    // calc_pair_evolution (each body in turn, companion = the other body).
    //
    // c_LayeredWorld hides the base analytic calc_tides with the rheology + layer-distribution path and
    // owns the spin model, so the concrete type is resolved here to run the right solve and reach the
    // spin rate; a layerless world (e.g. a star) uses the base analytic solve and contributes no spin.
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
        const double target_mass    = world_ptr->get_mass();
        const double spin_frequency = world_ptr->get_spin_frequency();

        out.orbital_frequency = orbital_frequency;
        out.semi_major_axis   = semi_major_axis;
        out.eccentricity      = eccentricity;
        out.spin_frequency    = spin_frequency;
        out.host_mass         = companion_mass;
        out.target_mass       = target_mass;

        // A body without a tide model is rigid: it raises no tide and contributes nothing.
        if (!world_ptr->get_tide_model_set()) {
            out.evolved = true;
            return out;
        }

        // Global tidal solve for this body with the companion as the tide raiser.
        c_TideSolveConfig state;
        state.orbital_frequency = orbital_frequency;
        state.spin_frequency    = spin_frequency;
        state.eccentricity      = eccentricity;
        state.obliquity         = world_ptr->get_obliquity();
        state.semi_major_axis   = semi_major_axis;
        state.host_mass         = companion_mass;

        c_LayeredWorld* layered = dynamic_cast<c_LayeredWorld*>(world_ptr);
        if (layered != nullptr) {
            layered->calc_tides(state);
        } else {
            world_ptr->calc_tides(state);
        }

        // calc_tides succeeded above (it throws on failure), so the tide result is populated.
        const c_GlobalTideResult& tide = world_ptr->get_tide_result();
        out.tidal_heating = tide.tidal_heating;
        out.dU_dM         = tide.dU_dM;
        out.dU_dw         = tide.dU_dw;
        out.dU_dO         = tide.dU_dO;

        // Orbital rates from the tidal-potential derivatives (this body is the dissipator).
        c_OrbitState orbit_state;
        orbit_state.orbital_frequency = orbital_frequency;
        orbit_state.semi_major_axis   = semi_major_axis;
        orbit_state.eccentricity      = eccentricity;
        orbit_state.target_mass       = target_mass;
        orbit_state.host_mass         = companion_mass;
        const c_OrbitDerivatives rates =
            this->p_orbit_solver.calc_derivatives(orbit_state, out.dU_dM, out.dU_dw);
        out.da_dt = rates.da_dt;
        out.de_dt = rates.de_dt;
        out.dn_dt = rates.dn_dt;

        // Spin rate from this body's own spin model (torque from the companion).
        if (layered != nullptr) {
            out.moment_of_inertia = layered->get_moment_of_inertia();
            out.dspin_dt          = layered->calc_spin_derivative(companion_mass);
            out.has_spin          = true;
        }

        // Energy-balance bookkeeping for this body.
        out.dE_orbit_dt     = this->calc_orbital_energy_derivative(out);
        out.dE_spin_dt      = this->calc_spin_energy_derivative(out);
        out.energy_residual = out.tidal_heating + out.dE_orbit_dt + out.dE_spin_dt;
        out.evolved         = true;
        return out;
    }

    // -----------------------------------------------------------------------
    // Binary I/O
    //
    // write_binary records the system's own container state (name, host/star indices, and each world's
    // orbital elements about both the host and the star) followed by every world's complete binary
    // record. read_binary reverses it, rebuilding the heterogeneous world list from the stream via the
    // world binary-dispatch factory (c_world_from_binary). Physics sub-models a world does not serialize
    // (the star's luminosity model, the layer EOS profile data, the orbit-solver/spin/tide models) are
    // reattached after load, exactly as for a directly-loaded world.
    // -----------------------------------------------------------------------
    void write_binary(std::ostream& out) const override {
        const auto num_worlds = static_cast<uint64_t>(this->p_worlds.size());
        uint64_t payload =
            binary_string_bytes(this->p_name)
            + sizeof(int32_t) * 2                   // host index, star index
            + sizeof(uint64_t)                      // world count
            + num_worlds * (4 * sizeof(double));    // per-world host-orbit + star-orbit elements
        write_binary_header(out, static_cast<uint32_t>(BinaryClassID::System), payload);
        write_binary_string(out, this->p_name);
        const int32_t host_index = this->p_host_index;
        const int32_t star_index = this->p_star_index;
        out.write(reinterpret_cast<const char*>(&host_index), sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&star_index), sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&num_worlds), sizeof(uint64_t));
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

    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);
        this->p_name = read_binary_string(in);

        int32_t host_index = -1;
        int32_t star_index = -1;
        in.read(reinterpret_cast<char*>(&host_index), sizeof(int32_t));
        in.read(reinterpret_cast<char*>(&star_index), sizeof(int32_t));
        uint64_t num_worlds = 0;
        in.read(reinterpret_cast<char*>(&num_worlds), sizeof(uint64_t));
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read System binary data");
        }

        // Per-world orbital elements about the host, then about the star (same order write_binary used).
        this->p_orbits.assign(num_worlds, c_OrbitElements{});
        for (uint64_t i = 0; i < num_worlds; ++i) {
            in.read(reinterpret_cast<char*>(&this->p_orbits[i].semi_major_axis), sizeof(double));
            in.read(reinterpret_cast<char*>(&this->p_orbits[i].eccentricity),    sizeof(double));
        }
        this->p_stellar_orbits.assign(num_worlds, c_OrbitElements{});
        for (uint64_t i = 0; i < num_worlds; ++i) {
            in.read(reinterpret_cast<char*>(&this->p_stellar_orbits[i].semi_major_axis), sizeof(double));
            in.read(reinterpret_cast<char*>(&this->p_stellar_orbits[i].eccentricity),    sizeof(double));
        }

        // Rebuild the heterogeneous world list; each world's concrete type is recovered from its record.
        this->p_worlds.clear();
        this->p_worlds.reserve(num_worlds);
        for (uint64_t i = 0; i < num_worlds; ++i) {
            this->p_worlds.push_back(c_world_from_binary(in, force));
        }

        this->p_host_index = host_index;
        this->p_star_index = star_index;
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read System binary data");
        }
    }

protected:
    // Bounds check shared by the index-based accessors.
    void check_index(std::size_t index) const {
        if (index >= this->p_worlds.size()) {
            throw std::out_of_range("TidalPy: c_System world index out of range");
        }
    }

    std::string p_name;
    std::vector<std::shared_ptr<c_BaseWorld>> p_worlds;         // owned worlds (shared with the Python wrappers)
    std::vector<c_OrbitElements>              p_orbits;         // orbit about the tidal host; host entry unused
    std::vector<c_OrbitElements>              p_stellar_orbits; // orbit about the star; star entry unused
    int p_host_index = -1;                                      // index into p_worlds, or -1 if unset
    int p_star_index = -1;                                      // index into p_worlds, or -1 if unset
    c_OrbitSolver p_orbit_solver;                              // stateless engine turning dU/dX into orbital rates
};

} // namespace tidalpy
