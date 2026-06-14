#pragma once
/*
 * interp_.hpp — 1-D linear interpolation utilities (numpy.interp-style).
 *
 * The search routine is a binary search seeded with a guess (adapted from
 * NumPy's `compiled_interp`); the interpolation matches `numpy.interp`:
 * out-of-range queries clamp to the endpoint values, and a NaN slope falls back
 * to the value from the other bracket endpoint.
 *
 * The domain `x_domain` must be sorted ascending and have length >= 1.
 */

#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>

namespace tidalpy {

// ---------------------------------------------------------------------------
// c_binary_search_with_guess
// ---------------------------------------------------------------------------
// Find the index j such that array[j] <= key < array[j+1] for a sorted-ascending
// array, seeded with `guess` to accelerate near-sequential queries. Returns
// `length` when key is past the right end; sets `code = -1` (and returns 0) when
// key is left of the array. Requires length >= 3 (callers handle shorter arrays).
inline std::size_t c_binary_search_with_guess(
        double key,
        const double* array,
        std::size_t length,
        std::size_t guess,
        int& code) {
    constexpr std::size_t LIKELY_IN_CACHE_SIZE = 8;
    code = 0;

    if (key > array[length - 1]) { return length; }
    if (key < array[0])          { code = -1; return 0; }

    if (guess > (length - 3)) { guess = length - 3; }
    if (guess < 1)            { guess = 1; }

    std::size_t imin = 0;
    std::size_t imax = length;

    // Check the most likely values: guess - 1, guess, guess + 1.
    if (key < array[guess]) {
        if (key < array[guess - 1]) {
            imax = guess - 1;
            if ((guess > LIKELY_IN_CACHE_SIZE) && (key >= array[guess - LIKELY_IN_CACHE_SIZE])) {
                imin = guess - LIKELY_IN_CACHE_SIZE;
            }
        } else {
            return guess - 1;
        }
    } else {
        if (key < array[guess + 1]) {
            return guess;
        } else if (key < array[guess + 2]) {
            return guess + 1;
        } else {
            imin = guess + 2;
            if ((guess < (length - LIKELY_IN_CACHE_SIZE - 1)) &&
                (key < array[guess + LIKELY_IN_CACHE_SIZE])) {
                imax = guess + LIKELY_IN_CACHE_SIZE;
            }
        }
    }

    // Bisection over the restricted range.
    while (imin < imax) {
        const std::size_t imid = imin + ((imax - imin) >> 1);
        if (key >= array[imid]) { imin = imid + 1; } else { imax = imid; }
    }

    if (imin == 0) { code = -1; }
    return imin - 1;
}

// ---------------------------------------------------------------------------
// c_interp — real linear interpolation (numpy.interp-style)
// ---------------------------------------------------------------------------
// Interpolate `dependent_values` (sampled on the sorted-ascending `x_domain`) at
// `desired_x`. Out-of-range queries clamp to the endpoint values. `len_x` is the
// length of both arrays.
//
// `guess` seeds the binary search (pass the previous result for near-sequential
// queries; pass 0 otherwise). Returns NaN only for an empty domain.
inline double c_interp(
        double desired_x,
        const double* x_domain,
        const double* dependent_values,
        std::size_t len_x,
        std::size_t guess = 0) {
    if (len_x == 0) { return std::numeric_limits<double>::quiet_NaN(); }
    if (len_x == 1) { return dependent_values[0]; }

    // Endpoint clamping (matches numpy.interp's default left/right behaviour).
    if (desired_x <= x_domain[0])         { return dependent_values[0]; }
    if (desired_x >= x_domain[len_x - 1]) { return dependent_values[len_x - 1]; }

    std::size_t j;
    if (len_x == 2) {
        j = 0;  // only one interval; the search routine needs len >= 3.
    } else {
        int code = 0;
        j = c_binary_search_with_guess(desired_x, x_domain, len_x, guess, code);
        if (code == -1)     { return dependent_values[0]; }
        if (j >= len_x)     { return dependent_values[len_x - 1]; }
        if (j == len_x - 1) { return dependent_values[j]; }
    }

    const double xp_j   = x_domain[j];
    const double fp_j   = dependent_values[j];
    if (xp_j == desired_x) { return fp_j; }
    const double xp_jp1 = x_domain[j + 1];
    const double fp_jp1 = dependent_values[j + 1];
    const double slope  = (fp_jp1 - fp_j) / (xp_jp1 - xp_j);

    double result = slope * (desired_x - xp_j) + fp_j;
    // If we get NaN in one direction, try the other (numpy's robustness trick).
    if (std::isnan(result)) {
        result = slope * (desired_x - xp_jp1) + fp_jp1;
        if (std::isnan(result) && (fp_jp1 == fp_j)) { result = fp_j; }
    }
    return result;
}

// ---------------------------------------------------------------------------
// c_interp_complex — complex linear interpolation
// ---------------------------------------------------------------------------
// As c_interp, but the dependent values are complex. Real and imaginary parts are
// interpolated independently.
inline std::complex<double> c_interp_complex(
        double desired_x,
        const double* x_domain,
        const std::complex<double>* dependent_values,
        std::size_t len_x,
        std::size_t guess = 0) {
    if (len_x == 0) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return std::complex<double>(nan, nan);
    }
    if (len_x == 1) { return dependent_values[0]; }

    if (desired_x <= x_domain[0])         { return dependent_values[0]; }
    if (desired_x >= x_domain[len_x - 1]) { return dependent_values[len_x - 1]; }

    std::size_t j;
    if (len_x == 2) {
        j = 0;
    } else {
        int code = 0;
        j = c_binary_search_with_guess(desired_x, x_domain, len_x, guess, code);
        if (code == -1)     { return dependent_values[0]; }
        if (j >= len_x)     { return dependent_values[len_x - 1]; }
        if (j == len_x - 1) { return dependent_values[j]; }
    }

    const double xp_j = x_domain[j];
    if (xp_j == desired_x) { return dependent_values[j]; }
    const double xp_jp1   = x_domain[j + 1];
    const double inv_dx   = 1.0 / (xp_jp1 - xp_j);
    const std::complex<double> fp_j   = dependent_values[j];
    const std::complex<double> fp_jp1 = dependent_values[j + 1];
    const std::complex<double> slope  = (fp_jp1 - fp_j) * inv_dx;
    return slope * (desired_x - xp_j) + fp_j;
}

} // namespace tidalpy
