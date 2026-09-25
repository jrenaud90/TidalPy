#pragma once
/* Gauss-Legendre nodes and weights.
 *
 * An n-node rule integrates a polynomial of degree 2n - 1 exactly, which is why it serves both the
 * colatitude integrals of the 3D tides (polynomials in cos theta) and the radial integrals over a layer.
 */

#include <cmath>
#include <vector>

#include "../../constants_.hpp"

namespace tidalpy {

/// Nodes and weights on [-1, 1], by Newton-Raphson on the Legendre polynomial. Map to [a, b] with
/// x = 0.5 (b - a) node + 0.5 (a + b) and weight 0.5 (b - a).
inline void c_gauss_legendre_nodes(
        int num_nodes,
        std::vector<double>& nodes,
        std::vector<double>& weights) {
    nodes.resize(num_nodes);
    weights.resize(num_nodes);
    const double pi = TidalPyConstants::d_PI;
    const int half = (num_nodes + 1) / 2;
    for (int i = 0; i < half; ++i) {
        double x = std::cos(pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(num_nodes) + 0.5));
        double dp = 1.0;
        for (int iter = 0; iter < 100; ++iter) {
            double p0 = 1.0;   // P_0
            double p1 = x;     // P_1
            for (int k = 2; k <= num_nodes; ++k) {
                const double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / static_cast<double>(k);
                p0 = p1;
                p1 = p2;
            }
            dp = static_cast<double>(num_nodes) * (x * p1 - p0) / (x * x - 1.0);
            const double dx = -p1 / dp;
            x += dx;
            if (std::abs(dx) <= 1.0e-15) { break; }
        }
        nodes[i]                 = -x;
        nodes[num_nodes - 1 - i] = x;
        const double w = 2.0 / ((1.0 - x * x) * dp * dp);
        weights[i]                 = w;
        weights[num_nodes - 1 - i] = w;
    }
}

} // namespace tidalpy
