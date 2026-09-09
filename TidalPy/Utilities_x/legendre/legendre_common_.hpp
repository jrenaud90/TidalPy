#pragma once
/*
 * legendre_common_.hpp - Shared type for the precomputed associated-Legendre tables.
 */

#include "constants_.hpp"

namespace tidalpy {

// One associated Legendre function value and its first/second colatitude (theta) derivatives.
struct c_LegendreValue {
    double p;             // P_lm(cos theta)
    double dp_dtheta;     // dP_lm/dtheta
    double d2p_dtheta2;   // d2P_lm/dtheta2
};

} // namespace tidalpy
