#pragma once
/*
 * tide_collapse_.hpp — global (1D) tidal mode collapse.
 *
 * c_global_potential (potential/global_.hpp) produces, model-independently, the per-mode
 * potential terms (dU/dM, dU/dw, dU/dO, E_dot) carrying the common coefficient
 * G_lpq^2 * F_lmp^2 * (l-m)!/(l+m)! * (R/a)^(2l+1) * G*M_host/a, plus the unique-frequency
 * maps. This header collapses those terms with a tide model's per-mode dissipation
 * multiplier -Im[k_l(omega)] to give the world's global tidal heating and the
 * three orbital-potential derivatives.
 *
 * Collapse (per active mode lmpq, omega = |omega_lmpq|):
 *   neg_imk = -Im[ tide_model.calc_love_numbers(l, omega, solver_love_lmpq).k ]
 *   tidal_heating += E_dot_term * neg_imk      [W]
 *   dU/dX         += dU_dX_term * neg_imk       [J kg-1 rad-1]   (X = M, w, O)
 *
 * The full c_LoveNumbers suite (k, h, l) is the transport type even though only k drives
 * the collapse. Layer scaling (tidal_scale) is applied by the world afterward, not here;
 * the whole-body collapse uses the unscaled -Im[k]. The rheology model needs the
 * radial-solver Love numbers per mode (solver_love_by_lmpq); analytic models pass nullptr.
 *
 * Reference: Renaud et al. (2021, PSJ), Eq. 7 and surrounding (GlobalTidalDissipation).
 */

#include <complex>

#include "global_.hpp"        // c_GlobalPotentialStorage, c_GlobalPotentialResultAtMode
#include "intmap_.hpp"        // c_IntMap
#include "keys_.hpp"          // c_Key4
#include "tide_base_.hpp"     // tidalpy::c_TideBase
#include "tide_result_.hpp"   // c_GlobalTideResult

// Collapse the per-mode global potential terms with the tide model's dissipation
// multiplier. `solver_love_by_lmpq` supplies the radial-solver Love numbers (k, h, l) per
// mode for the rheology model; pass nullptr for the analytic (CPL/CTL/CTL_Q) models.
inline c_GlobalTideResult c_collapse_global_tides(
        const c_GlobalPotentialStorage& potential,
        const tidalpy::c_TideBase& tide_model,
        const c_IntMap<c_Key4, tidalpy::c_LoveNumbers>* solver_love_by_lmpq = nullptr)
{
    c_GlobalTideResult result;
    result.error_code = potential.error_code;
    if (potential.error_code != 0) {
        return result;
    }

    bool found = false;
    for (const auto& mode_entry : potential.potential_map) {
        const c_Key4& lmpq_key                       = mode_entry.first;
        const c_GlobalPotentialResultAtMode& terms   = mode_entry.second;
        const int degree_l                           = static_cast<int>(lmpq_key.a);

        // Resolve this mode's unique tidal frequency.
        found = false;
        const size_t freq_index = potential.unique_freq_index_map.get(found, lmpq_key);
        if (!found) {
            // Mode had no recorded (nonzero) frequency; it contributes nothing.
            continue;
        }
        const double frequency = potential.unique_freq_map[freq_index].frequency;

        // Radial-solver Love numbers (k, h, l) for this mode (rheology only).
        tidalpy::c_LoveNumbers solver_love;
        if (solver_love_by_lmpq != nullptr) {
            bool love_found = false;
            solver_love = solver_love_by_lmpq->get(love_found, lmpq_key);
        }

        const double neg_imk = tide_model.calc_neg_imk(degree_l, frequency, solver_love);

        result.tidal_heating += terms.E_dot * neg_imk;
        result.dU_dM         += terms.dU_dM * neg_imk;
        result.dU_dw         += terms.dU_dw * neg_imk;
        result.dU_dO         += terms.dU_dO * neg_imk;
        result.num_modes     += 1;
    }

    return result;
}
