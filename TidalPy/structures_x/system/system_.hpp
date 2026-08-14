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

#include "../worlds/base_.hpp"      // c_BaseWorld
#include "../worlds/stellar_.hpp"   // c_StarWorld (for the star's luminosity in insolation)
#include "constants_.hpp"           // TidalPyConstants::d_EPS / d_NAN / d_PI, tidalpy_config_ptr->d_G

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
    // Binary I/O
    //
    // write_binary records the system's own container state (name, host/star indices, and each world's
    // orbital elements about both the host and the star) followed by every world's complete binary
    // record. Loading a system back (reconstructing the heterogeneous world list from the stream) is
    // handled separately once the world binary-dispatch factory is available.
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
};

} // namespace tidalpy
