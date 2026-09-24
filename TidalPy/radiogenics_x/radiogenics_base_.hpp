#pragma once
/* Abstract base for TidalPy radiogenics models. Concrete models live in radiogenics_.hpp. All MKS. */

#include <stdexcept>
#include <string>
#include <vector>

#include "physics_base_.hpp"

namespace tidalpy {

class c_RadiogenicsBase : public c_PhysicsBase {
public:
    c_RadiogenicsBase() = default;

    explicit c_RadiogenicsBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_RadiogenicsBase() override = default;

    // Heating [W] from `mass` [kg] at `time` [s]. Time shares its zero point with the model's
    // reference time. Assumes exponential decay from that reference time.
    virtual double calc_heating(double time, double mass) const = 0;

    // Time [s] the model's quoted abundances or rate apply at; zero for a model with no decay.
    virtual double get_ref_time() const noexcept { return 0.0; }

    // Over time at constant mass. The vectorized calls are virtual so a model can hoist per-call constants.
    virtual void calc_heating_vectorize_time(
            const std::vector<double>& time,
            double mass,
            std::vector<double>& out_heating) const {
        const std::size_t n = time.size();
        out_heating.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_heating[i] = this->calc_heating(time[i], mass);
        }
    }

    // Over mass at constant time.
    virtual void calc_heating_vectorize_mass(
            double time,
            const std::vector<double>& mass,
            std::vector<double>& out_heating) const {
        const std::size_t n = mass.size();
        out_heating.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_heating[i] = this->calc_heating(time, mass[i]);
        }
    }

    // Element-wise over time and mass.
    virtual void calc_heating_vectorize_all(
            const std::vector<double>& time,
            const std::vector<double>& mass,
            std::vector<double>& out_heating) const {
        if (time.size() != mass.size()) {
            throw std::invalid_argument(
                "TidalPy::calc_heating_vectorize_all: time and mass vectors must "
                "have the same length");
        }
        const std::size_t n = time.size();
        out_heating.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_heating[i] = this->calc_heating(time[i], mass[i]);
        }
    }
};

} // namespace tidalpy
