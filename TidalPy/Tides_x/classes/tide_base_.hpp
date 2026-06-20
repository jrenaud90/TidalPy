#pragma once
/*
 * tide_base_.hpp — c_TideBase: abstract base for TidalPy's global (1D) tidal
 * dissipation models.
 *
 * Inherits c_PhysicsBase (Utilities_x/classes_x/physics_base_.hpp). A tide model is a
 * world-level physics model: it converts a per-mode Love number into the dissipation
 * multiplier used by the global mode collapse (tide_collapse_.hpp). The four concrete
 * models live in tide_.hpp:
 *
 *   c_RheologyTide  (alias "rheology")           — k_l from the radial solver (frequency dependent).
 *   c_FixedQTide    (alias "cpl"/"fixed_q")      — constant phase lag, k_l*(1 - i/Q_l).
 *   c_FixedLagTide  (alias "ctl"/"fixed_dt")     — constant time lag,  k_l*(1 - i*omega*dt_l).
 *   c_CTLQTide      (alias "ctl_q"/"fixed_dt_q") — k_l*(1 - i*omega*dt_l/Q_l).
 *
 * The collapse only needs the dissipation multiplier -Im[k_l(omega)], but the full
 * c_LoveNumbers suite (k, h, l) is always the transport type so the displacement Love
 * numbers from the radial solver are never thrown away. The potential terms (heating,
 * dU/dM, dU/dw, dU/dO) are produced once, model-independently, by c_global_potential
 * (potential/global_.hpp). All quantities MKS; frequencies in rad s-1. calc_* methods
 * are const.
 *
 * References
 * ----------
 * - Renaud et al. (2021, PSJ) — global dual-body tidal dissipation (collapse form).
 * - Efroimsky & Makarov (2013) — CPL/CTL frequency dependence of the dissipation.
 */

#include <complex>
#include <string>

#include "love_.hpp"          // tidalpy::c_LoveNumbers
#include "physics_base_.hpp"

namespace tidalpy {

class c_TideBase : public c_PhysicsBase {
public:
    c_TideBase() = default;

    explicit c_TideBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_TideBase() override = default;

    // -----------------------------------------------------------------------
    // Complex Love numbers (k, h, l) at tidal frequency [rad s-1] (pure virtual).
    //
    // The full c_LoveNumbers suite is always the transport type, even though only k is
    // needed for heating/orbital dynamics. Analytic models (FixedQ/FixedLag/CTLQ) build
    // k_l from their fixed per-degree parameters and frequency law and set h, l to NaN
    // (they cannot produce displacement Love numbers without a radial solution); they
    // ignore `solver_love`. The rheology model returns `solver_love` unchanged (the full
    // k/h/l from the radial solver, supplied by the world at this frequency).
    //
    // Assumptions
    // -----------
    // - frequency is the tidal forcing frequency magnitude |omega_lmpq| (>= 0).
    // -----------------------------------------------------------------------
    virtual c_LoveNumbers calc_love_numbers(
            int degree_l, double frequency, const c_LoveNumbers& solver_love) const = 0;

    // -Im[k_l] — the dissipation multiplier used in the mode collapse. Derived
    // from the k component of the full Love-number suite.
    double calc_neg_imk(int degree_l, double frequency, const c_LoveNumbers& solver_love) const {
        return -std::imag(this->calc_love_numbers(degree_l, frequency, solver_love).k);
    }

    // True if the model requires the radial solver to supply the Love numbers (rheology);
    // false for the analytic models. The world uses this to decide whether to run a radial
    // solve.
    virtual bool needs_radial_solve() const = 0;
};

} // namespace tidalpy
