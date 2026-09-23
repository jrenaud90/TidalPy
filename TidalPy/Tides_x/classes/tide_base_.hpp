#pragma once
/*
 * tide_base_.hpp - c_TideBase: abstract base for TidalPy's global (1D) tidal dissipation models.
 *
 * Inherits c_PhysicsBase. A tide model converts a per-mode Love number into the dissipation
 * multiplier used by the global mode collapse (tide_collapse_.hpp). The concrete models live in
 * tide_.hpp:
 *
 *   c_RheologyTide  (alias "rheology")           - k_l from the radial solver (frequency dependent).
 *   c_FixedQTide    (alias "cpl"/"fixed_q")      - constant phase lag, k_l*(1 - i/Q_l).
 *   c_FixedLagTide  (alias "ctl"/"fixed_dt")     - constant time lag,  k_l*(1 - i*omega*dt_l).
 *   c_CTLQTide      (alias "ctl_q"/"fixed_dt_q") - k_l*(1 - i*omega*dt_l/Q_l).
 *
 * The full c_LoveNumbers suite (k, h, l) is the transport type even though the collapse needs only
 * -Im[k_l(omega)], so the radial solver's displacement Love numbers are never thrown away. All
 * quantities MKS; frequencies in rad s-1.
 *
 * References
 * ----------
 * - Renaud et al. (2021, PSJ): global dual-body tidal dissipation (collapse form).
 * - Efroimsky and Makarov (2013): CPL and CTL frequency dependence of the dissipation.
 */

#include <complex>
#include <limits>
#include <string>

// Explicit relative path, not a bare "love_.hpp": the layered and world extensions also carry
// RadialSolver_x on their include path, which holds a different love_.hpp in the global namespace, so a
// bare include can resolve to the wrong file depending on include-dir order.
#include "../love/love_.hpp"   // tidalpy::c_LoveNumbers
#include "physics_base_.hpp"

namespace tidalpy {

class c_TideBase : public c_PhysicsBase {
public:
    c_TideBase() = default;

    explicit c_TideBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_TideBase() override = default;

    // Complex Love numbers at the forcing frequency magnitude |omega_lmpq| [rad s-1]. The analytic models
    // build k_l from their fixed per-degree parameters, set h and l to NaN, and ignore solver_love; the
    // rheology model returns solver_love unchanged.
    virtual c_LoveNumbers calc_love_numbers(
            int degree_l, double frequency, const c_LoveNumbers& solver_love) const = 0;

    // -Im[k_l]: the dissipation multiplier used in the mode collapse.
    double calc_neg_imk(int degree_l, double frequency, const c_LoveNumbers& solver_love) const {
        return -std::imag(this->calc_love_numbers(degree_l, frequency, solver_love).k);
    }

    // True when the world must run the radial solver to supply the Love numbers.
    virtual bool needs_radial_solve() const = 0;

    // Fixed per-degree quality factor and time lag [s] when the model carries them; NaN otherwise. The
    // world's cpl and ctl Love methods fall back on these when no explicit value is configured.
    virtual double get_fixed_q(int /*degree_l*/) const { return std::numeric_limits<double>::quiet_NaN(); }
    virtual double get_fixed_dt(int /*degree_l*/) const { return std::numeric_limits<double>::quiet_NaN(); }
};

} // namespace tidalpy
