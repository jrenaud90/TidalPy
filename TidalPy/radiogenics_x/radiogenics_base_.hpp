#pragma once
/*
 * radiogenics_base_.hpp - c_RadiogenicsBase: abstract base for all TidalPy radiogenics models.
 *
 * Inherits c_PhysicsBase (Utilities_x/classes_x/physics_base_.hpp).
 *
 * Defines the abstract interface that every radiogenics model must satisfy:
 *   calc_heating(time, mass) -> double [W]
 *
 * The three implemented models (Off, Isotope, Fixed) live in radiogenics_.hpp.
 *
 * All calc_* methods are const and operate in MKS units (time [s], mass [kg],
 * heating [W]).
 */

#include <stdexcept>
#include <string>
#include <vector>

#include "physics_base_.hpp"

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_RadiogenicsBase
// -------------------------------------------------------------------------------
class c_RadiogenicsBase : public c_PhysicsBase {
public:
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------
    c_RadiogenicsBase() = default;

    explicit c_RadiogenicsBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_RadiogenicsBase() override = default;

    // -----------------------------------------------------------------------
    // Radiogenic heating (pure virtual) [W]
    //
    // Each model overrides this to implement its decay law and return the total
    // radiogenic heating produced by the given mass at the given time.
    // c_SolidLiquidLayer calls this through calc_radiogenic_heating.
    //
    // Parameters
    // ----------
    // time_s : elapsed time [s] (must share its zero point with the reference
    //          time stored on the model)
    // mass_kg : mass of the radiogenic material [kg]
    //
    // Returns
    // -------
    // double : radiogenic heating [W]
    //
    // Assumptions
    // -----------
    // - Exponential radioactive decay from a reference time.
    // - All inputs and outputs are MKS.
    // -----------------------------------------------------------------------
    virtual double calc_heating(double time_s, double mass_kg) const = 0;

    // -----------------------------------------------------------------------
    // Vectorized heating — vary time at constant mass.
    //
    // Evaluates calc_heating over each time at the single constant mass.
    // out_heating is resized to the time vector length.
    // -----------------------------------------------------------------------
    void calc_heating_vectorize_time(
            const std::vector<double>& time_s,
            double mass_kg,
            std::vector<double>& out_heating) const {
        const std::size_t n = time_s.size();
        out_heating.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_heating[i] = this->calc_heating(time_s[i], mass_kg);
        }
    }

    // -----------------------------------------------------------------------
    // Vectorized heating — vary mass at constant time.
    //
    // Evaluates calc_heating over each mass at the single constant time.
    // out_heating is resized to the mass vector length.
    // -----------------------------------------------------------------------
    void calc_heating_vectorize_mass(
            double time_s,
            const std::vector<double>& mass_kg,
            std::vector<double>& out_heating) const {
        const std::size_t n = mass_kg.size();
        out_heating.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_heating[i] = this->calc_heating(time_s, mass_kg[i]);
        }
    }

    // -----------------------------------------------------------------------
    // Vectorized heating — vary both time and mass element-wise.
    //
    // The two input vectors must have the same length N; out_heating is resized
    // to N.  Throws std::invalid_argument if the input vectors differ in length.
    // -----------------------------------------------------------------------
    void calc_heating_vectorize_all(
            const std::vector<double>& time_s,
            const std::vector<double>& mass_kg,
            std::vector<double>& out_heating) const {
        if (time_s.size() != mass_kg.size()) {
            throw std::invalid_argument(
                "TidalPy: calc_heating_vectorize_all — time and mass vectors must "
                "have the same length");
        }
        const std::size_t n = time_s.size();
        out_heating.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_heating[i] = this->calc_heating(time_s[i], mass_kg[i]);
        }
    }
};

} // namespace tidalpy
