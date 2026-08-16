#pragma once
/*
 * interp_.hpp (utilities/arrays) — pointer-style front end over the shared interpolation
 * utilities in Utilities_x/arrays/interp_.hpp (the single implementation).
 *
 * The wrappers below keep the pointer/output-argument call style used by the EOS solution,
 * the radial-solver storage, and the world Love-number solve. `provided_j_ptr` seeds the
 * binary search (fast for near-sequential queries); the search always verifies the interval,
 * so a stale index can no longer select the wrong segment.
 */

#include <limits>
#include <cmath>
#include <complex>
#include <cstddef>

#include "../../Utilities_x/arrays/interp_.hpp"

// Find index j such that array[j] <= key < array[j+1]; returns `length` past the right end
// and sets code[0] = -1 (returning 0) left of the array.
inline std::size_t c_binary_search_with_guess(
        double key,
        double* array,
        std::size_t length,
        std::size_t guess,
        int* code)
{
    int status = 0;
    const std::size_t index = tidalpy::c_binary_search_with_guess(key, array, length, guess, status);
    code[0] = status;
    return index;
}

// Interpolate `dependent_values` at desired_x_ptr[0]; result written to result_ptr[0].
inline void c_interp(
        double* desired_x_ptr,
        double* x_domain_ptr,
        double* dependent_values_ptr,
        std::size_t len_x,
        std::size_t* provided_j_ptr,
        double* result_ptr)
{
    const std::size_t guess = (provided_j_ptr[0] < len_x) ? provided_j_ptr[0] : 0;
    result_ptr[0] = tidalpy::c_interp(desired_x_ptr[0], x_domain_ptr, dependent_values_ptr, len_x, guess);
}

// Complex variant: `dependent_values_ptr` is interleaved (re, im) pairs; the result's real and
// imaginary parts are written to result_ptr[0] and result_ptr[1].
inline void c_interp_complex(
        double desired_x,
        double* x_domain_ptr,
        double* dependent_values_ptr,
        std::size_t len_x,
        std::size_t* provided_j_ptr,
        double* result_ptr)
{
    const std::size_t guess = (provided_j_ptr[0] < len_x) ? provided_j_ptr[0] : 0;
    const std::complex<double>* values_ptr =
        reinterpret_cast<const std::complex<double>*>(dependent_values_ptr);
    const std::complex<double> result =
        tidalpy::c_interp_complex(desired_x, x_domain_ptr, values_ptr, len_x, guess);
    result_ptr[0] = result.real();
    result_ptr[1] = result.imag();
}
