#pragma once
/* Spherical geometry base class.
 *
 * The calc_* methods take explicit arguments rather than reading the stored radius and mass, so they stay
 * pure functions with no hidden state.
 *
 * Binary payload: radius then mass.
 */

#include <cmath>
#include <stdexcept>

#include "constants_.hpp"
#include "tidalpy_base_.hpp"

namespace tidalpy {

class c_StructureBase : public c_TidalPyBaseClass {
public:
    c_StructureBase() = default;

    c_StructureBase(double radius, double mass)
        : p_radius(radius), p_mass(mass) {}

    ~c_StructureBase() override = default;

    double get_radius() const noexcept { return p_radius; }
    double get_mass()   const noexcept { return p_mass; }

    // Geometry; all MKS in and out.
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
        return c_get_G() * mass / (radius * radius);
    }

    // Mean density [kg/m^3]
    double calc_mean_density(double mass, double volume) const noexcept {
        if (volume == 0.0) return 0.0;
        return mass / volume;
    }

    // Escape velocity [m/s]
    double calc_escape_velocity(double mass, double radius) const noexcept {
        if (radius == 0.0) return 0.0;
        return std::sqrt(2.0 * c_get_G() * mass / radius);
    }

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
