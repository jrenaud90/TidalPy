#pragma once
/*
 * structure_base_.hpp — c_StructureBase: spherical geometry base class.
 *
 * Stores radius [m] and mass [kg].
 * All calc_* methods are const and take explicit arguments — they are pure
 * functions that do not depend on the object's stored radius/mass.  This
 * matches the functional-style design principle (no hidden state in calcs).
 *
 * Binary format (20-byte header + 16 bytes payload):
 *   header: class_id = BinaryClassID::StructureBase (2)
 *   payload: p_radius (double, 8 bytes) | p_mass (double, 8 bytes)
 */

#include <cmath>
#include <stdexcept>

#include "constants_.hpp"
#include "tidalpy_base_.hpp"

namespace tidalpy {

class c_StructureBase : public c_TidalPyBaseClass {
public:
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------
    c_StructureBase() = default;

    c_StructureBase(double radius_m, double mass_kg)
        : p_radius(radius_m), p_mass(mass_kg) {}

    ~c_StructureBase() override = default;

    // -----------------------------------------------------------------------
    // Getters
    // -----------------------------------------------------------------------
    double get_radius() const noexcept { return p_radius; }
    double get_mass()   const noexcept { return p_mass; }

    // -----------------------------------------------------------------------
    // Geometry calculations (all const; all MKS inputs and outputs)
    // -----------------------------------------------------------------------
    // Surface area of a sphere [m^2]
    double calc_surface_area(double radius) const noexcept {
        return 4.0 * TidalPyConstants::d_PI * radius * radius;
    }

    // Volume of a solid sphere [m^3]
    double calc_volume_sphere(double radius) const noexcept {
        return (4.0 / 3.0) * TidalPyConstants::d_PI * radius * radius * radius;
    }

    // Volume of a spherical shell [m^3]
    double calc_volume_shell(double radius_outer, double radius_inner) const noexcept {
        return calc_volume_sphere(radius_outer) - calc_volume_sphere(radius_inner);
    }

    // Surface gravitational acceleration [m/s^2]
    double calc_surface_gravity(double mass, double radius) const noexcept {
        if (radius == 0.0) return 0.0;
        return tidalpy_config_ptr->d_G * mass / (radius * radius);
    }

    // Mean density [kg/m^3]
    double calc_mean_density(double mass, double volume) const noexcept {
        if (volume == 0.0) return 0.0;
        return mass / volume;
    }

    // Escape velocity [m/s]
    double calc_escape_velocity(double mass, double radius) const noexcept {
        if (radius == 0.0) return 0.0;
        return std::sqrt(2.0 * tidalpy_config_ptr->d_G * mass / radius);
    }

    // -----------------------------------------------------------------------
    // Binary I/O
    // -----------------------------------------------------------------------
    void write_binary(std::ostream& out) const override {
        constexpr uint64_t payload = 2 * sizeof(double);
        write_binary_header(
            out,
            static_cast<uint32_t>(BinaryClassID::StructureBase),
            payload);
        out.write(reinterpret_cast<const char*>(&p_radius), sizeof(double));
        out.write(reinterpret_cast<const char*>(&p_mass),   sizeof(double));
        if (!out) {
            throw std::runtime_error(
                "TidalPy: failed to write StructureBase binary data");
        }
    }

    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);
        in.read(reinterpret_cast<char*>(&p_radius), sizeof(double));
        in.read(reinterpret_cast<char*>(&p_mass),   sizeof(double));
        if (!in) {
            throw std::runtime_error(
                "TidalPy: failed to read StructureBase binary data");
        }
    }

protected:
    double p_radius = 0.0;  // [m]
    double p_mass   = 0.0;  // [kg]
};

} // namespace tidalpy
