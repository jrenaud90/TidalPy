#pragma once
/*
 * eccentricity_exact_.hpp - the eccentricity functions from the exact Kepler orbit (the "exact" truncation).
 *
 * G_lpq(e) = X^{-(l+1), m}_k(e), m = l - 2p, k = m + q, is the k-th Fourier coefficient in mean anomaly M of
 * (r/a)^{-(l+1)} e^{i m f}, f the true anomaly (Kaula 1964). Sampling that function at N points in M (Kepler's equation
 * solved at each) and taking one FFT gives every G_lpq of one (l, p) at once, exactly up to aliasing. Nothing is
 * truncated in e, so the result holds at any e < 1 given enough modes, and a product of two functions is the plain
 * product (no cancellation between modes, unlike the tabulated series at large e).
 *
 * Mode range. The coefficients fall off as exp(-rho |k|) with rho = ln[(1 + sqrt(1 - e^2)) / e] - sqrt(1 - e^2), set by
 * the orbit's singularity in complex mean anomaly. The modes kept are the |q| <= Q for which the q^2-weighted tail sum
 * over |q| > Q of G_lpq^2 is below `tolerance` of the whole q^2-weighted sum for every p. That weight is the
 * synchronous constant-time-lag heating's, the most demanding of the tide models measured, so the tolerance bounds
 * its relative error; the q = 0 modes, which carry no tail, do not enter it. The FFT length leaves the aliased
 * coefficients at |k| >= N - k_max below double precision.
 */

#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "eccentricity_accuracy_.hpp"

// The largest mode range the exact functions build; an eccentricity that needs more (about e > 0.99) is refused.
inline constexpr int C_ECCENTRICITY_EXACT_MAX_Q = 20000;

// G_lpq for one degree: value(p, q) for 0 <= p <= l and |q| <= max_q.
struct c_ExactEccentricityModes {
    int degree_l = 0;
    int max_q = 0;
    std::vector<double> values;   // [p * (2 max_q + 1) + (q + max_q)]

    double value(int p, int q) const noexcept {
        return this->values[static_cast<size_t>(p * (2 * this->max_q + 1) + (q + this->max_q))];
    }
};

namespace eccentricity_detail {

// In-place iterative radix-2 FFT, forward sign: X_k = sum_j x_j exp(-2 pi i j k / N). data.size() is a power of 2.
inline void c_fft_radix2(std::vector<std::complex<double>>& data) {
    const size_t n = data.size();
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) { j ^= bit; }
        j ^= bit;
        if (i < j) { std::swap(data[i], data[j]); }
    }
    const double pi = std::acos(-1.0);
    for (size_t length = 2; length <= n; length <<= 1) {
        const double angle = -2.0 * pi / static_cast<double>(length);
        const std::complex<double> step(std::cos(angle), std::sin(angle));
        for (size_t start = 0; start < n; start += length) {
            std::complex<double> twiddle(1.0, 0.0);
            for (size_t k = 0; k < length / 2; ++k) {
                const std::complex<double> even = data[start + k];
                const std::complex<double> odd = data[start + k + length / 2] * twiddle;
                data[start + k] = even + odd;
                data[start + k + length / 2] = even - odd;
                twiddle *= step;
            }
        }
    }
}

// Eccentric anomaly from mean anomaly (Newton on Kepler's equation from E = M + e sin M).
inline double c_kepler_eccentric_anomaly(double mean_anomaly, double eccentricity) noexcept {
    double eccentric = mean_anomaly + eccentricity * std::sin(mean_anomaly);
    for (int iteration = 0; iteration < 100; ++iteration) {
        const double residual = eccentric - eccentricity * std::sin(eccentric) - mean_anomaly;
        const double step = residual / (1.0 - eccentricity * std::cos(eccentric));
        eccentric -= step;
        if (std::abs(step) < 1.0e-15) { break; }
    }
    return eccentric;
}

}  // namespace eccentricity_detail

// Every G_lpq of one degree at one eccentricity, with the mode range set by `tolerance`. Throws std::invalid_argument
// for an eccentricity outside [0, 1) or a tolerance outside (0, 1), and std::runtime_error when the eccentricity needs
// more than C_ECCENTRICITY_EXACT_MAX_Q modes.
inline c_ExactEccentricityModes c_exact_eccentricity_modes(double eccentricity, int degree_l, double tolerance) {
    if (!((eccentricity >= 0.0) && (eccentricity < 1.0))) {
        throw std::invalid_argument(
            "TidalPy: the exact eccentricity functions need 0 <= e < 1; got " + std::to_string(eccentricity) + ".");
    }
    if (!((tolerance > 0.0) && (tolerance < 1.0))) {
        throw std::invalid_argument(
            "TidalPy: the exact eccentricity tolerance must be in (0, 1); got " + std::to_string(tolerance) + ".");
    }
    c_ExactEccentricityModes modes;
    modes.degree_l = degree_l;
    const int num_p = degree_l + 1;

    // A circular orbit: G_lpq = 1 at q = 0 and 0 otherwise (X^{n,m}_k(0) = delta_km).
    if (eccentricity == 0.0) {
        modes.max_q = 0;
        modes.values.assign(static_cast<size_t>(num_p), 1.0);
        return modes;
    }

    const double root = std::sqrt(1.0 - eccentricity * eccentricity);
    const double decay_rate = std::log((1.0 + root) / eccentricity) - root;
    // The tolerance's tail sits near the Q where exp(-2 rho Q) Q^2 ~ tolerance; search well past it.
    const double q_search_d = (std::log(1.0 / tolerance) + 40.0) / (2.0 * decay_rate) + 2.0 * degree_l + 8.0;
    if (!(q_search_d < static_cast<double>(C_ECCENTRICITY_EXACT_MAX_Q))) {
        throw std::runtime_error(
            "TidalPy: eccentricity " + std::to_string(eccentricity) + " needs more than "
            + std::to_string(C_ECCENTRICITY_EXACT_MAX_Q) + " modes for the exact eccentricity functions.");
    }
    const int q_search = static_cast<int>(q_search_d);
    const int k_max = q_search + degree_l;
    // Aliased coefficients |k| >= N - k_max fall off as exp(-rho (N - k_max)): 40 / rho leaves them below 1e-17.
    const double length_needed = std::fmax(2.0 * k_max + 2.0, k_max + 40.0 / decay_rate);
    size_t num_samples = 64;
    while (static_cast<double>(num_samples) < length_needed) { num_samples <<= 1; }

    // (r/a) and e^{i f} at the samples.
    const double two_pi = 2.0 * std::acos(-1.0);
    std::vector<double> radius_ratio(num_samples);
    std::vector<std::complex<double>> exp_true(num_samples);
    for (size_t j = 0; j < num_samples; ++j) {
        const double mean_anomaly = two_pi * static_cast<double>(j) / static_cast<double>(num_samples);
        const double eccentric = eccentricity_detail::c_kepler_eccentric_anomaly(mean_anomaly, eccentricity);
        const double ratio = 1.0 - eccentricity * std::cos(eccentric);
        radius_ratio[j] = ratio;
        exp_true[j] = std::complex<double>(std::cos(eccentric) - eccentricity, root * std::sin(eccentric)) / ratio;
    }

    // One FFT per p; keep |q| <= q_search.
    const int width_search = 2 * q_search + 1;
    std::vector<double> search(static_cast<size_t>(num_p * width_search), 0.0);
    std::vector<std::complex<double>> samples(num_samples);
    const double n_power = -static_cast<double>(degree_l + 1);
    for (int p = 0; p < num_p; ++p) {
        const int order_m = degree_l - 2 * p;
        for (size_t j = 0; j < num_samples; ++j) {
            samples[j] = std::pow(radius_ratio[j], n_power) * std::pow(exp_true[j], order_m);
        }
        eccentricity_detail::c_fft_radix2(samples);
        for (int q = -q_search; q <= q_search; ++q) {
            const long k = static_cast<long>(order_m + q);
            const size_t index = static_cast<size_t>((k % static_cast<long>(num_samples) + static_cast<long>(num_samples))
                                                     % static_cast<long>(num_samples));
            search[static_cast<size_t>(p * width_search + (q + q_search))] =
                samples[index].real() / static_cast<double>(num_samples);
        }
    }

    // The smallest Q whose q^2-weighted tail is below the tolerance for every p.
    int max_q = 0;
    for (int p = 0; p < num_p; ++p) {
        double total = 0.0;
        for (int q = -q_search; q <= q_search; ++q) {
            const double g = search[static_cast<size_t>(p * width_search + (q + q_search))];
            total += static_cast<double>(q) * q * g * g;
        }
        if (!(total > 0.0)) { continue; }
        double tail = 0.0;
        int q_needed = 0;
        for (int q = q_search; q >= 1; --q) {
            const double g_plus = search[static_cast<size_t>(p * width_search + (q + q_search))];
            const double g_minus = search[static_cast<size_t>(p * width_search + (-q + q_search))];
            tail += static_cast<double>(q) * q * (g_plus * g_plus + g_minus * g_minus);
            if (tail >= tolerance * total) {
                q_needed = q;
                break;
            }
        }
        if (q_needed > max_q) { max_q = q_needed; }
    }

    modes.max_q = max_q;
    const int width = 2 * max_q + 1;
    modes.values.assign(static_cast<size_t>(num_p * width), 0.0);
    for (int p = 0; p < num_p; ++p) {
        for (int q = -max_q; q <= max_q; ++q) {
            modes.values[static_cast<size_t>(p * width + (q + max_q))] =
                search[static_cast<size_t>(p * width_search + (q + q_search))];
        }
    }
    return modes;
}
