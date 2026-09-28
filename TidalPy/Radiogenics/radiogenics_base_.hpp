#pragma once
/* Abstract base for TidalPy radiogenics models. Concrete models live in radiogenics_.hpp. All MKS. */

#include <string>
#include <vector>

#include "broadcast_.hpp"
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

    // Element-wise over time and mass. Each holds one value per point or a single value used at every point.
    virtual void calc_heating_vectorize(
            const std::vector<double>& time,
            const std::vector<double>& mass,
            std::vector<double>& out_heating) const {
        const std::size_t num_points = c_broadcast_length({time.size(), mass.size()}, "calc_heating_vectorize");
        const std::size_t time_stride = c_broadcast_stride(time.size());
        const std::size_t mass_stride = c_broadcast_stride(mass.size());
        out_heating.resize(num_points);
        for (std::size_t i = 0; i < num_points; ++i) {
            out_heating[i] = this->calc_heating(time[i * time_stride], mass[i * mass_stride]);
        }
    }
};

} // namespace tidalpy
