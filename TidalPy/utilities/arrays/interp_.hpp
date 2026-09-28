#pragma once
/* 1-D linear interpolation matching numpy.interp, using a binary search seeded with a guess (adapted from NumPy's
 * compiled_interp). Out-of-range queries clamp to the endpoint values. In the real-valued `cf_interp`, a NaN result
 * from one bracket endpoint falls back to the other. The x domain must be sorted ascending.
 *
 * `provided_j_ptr` seeds the binary search (fast for near-sequential queries); it is read, not updated. The search
 * always verifies the interval, so a stale index cannot select the wrong segment.
 */

#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>

// Find index j such that array[j] <= key < array[j+1], seeded with `guess` to speed up near-sequential queries.
// Returns `length` past the right end; sets code[0] = -1 (returning 0) left of the array.
inline std::size_t cf_binary_search_with_guess(
        double key,
        double* array,
        std::size_t length,
        std::size_t guess,
        int* code)
{
    const std::size_t LIKELY_IN_CACHE_SIZE = 8;
    code[0] = 0;

    if (length == 0)
    {
        code[0] = -1;
        return 0;
    }
    if (key > array[length - 1])
    {
        return length;
    }
    if (key < array[0])
    {
        code[0] = -1;
        return 0;
    }

    // Too short for the guess fast paths, which read array[guess - 1] .. array[guess + 2]: scan the few intervals
    // directly. This also keeps the size_t subtraction below safe.
    if (length <= 4)
    {
        std::size_t interval = 0;
        while ((interval + 2 < length) && (key >= array[interval + 1]))
        {
            ++interval;
        }
        return interval;
    }

    if (guess > (length - 3))
    {
        guess = length - 3;
    }
    if (guess < 1)
    {
        guess = 1;
    }

    std::size_t imin = 0;
    std::size_t imax = length;

    // Check the most likely values: guess - 1, guess, guess + 1.
    if (key < array[guess])
    {
        if (key < array[guess - 1])
        {
            imax = guess - 1;
            // Last attempt to restrict the search to items in cache.
            if ((guess > LIKELY_IN_CACHE_SIZE) && (key >= array[guess - LIKELY_IN_CACHE_SIZE]))
            {
                imin = guess - LIKELY_IN_CACHE_SIZE;
            }
        }
        else
        {
            return guess - 1;
        }
    }
    else
    {
        if (key < array[guess + 1])
        {
            return guess;
        }
        else if (key < array[guess + 2])
        {
            return guess + 1;
        }
        else
        {
            imin = guess + 2;
            // Written as an addition rather than length - window - 1: the subtraction form underflows in unsigned
            // size_t for short arrays and would then read out of bounds.
            if ((guess + LIKELY_IN_CACHE_SIZE + 1 < length) && (key < array[guess + LIKELY_IN_CACHE_SIZE]))
            {
                imax = guess + LIKELY_IN_CACHE_SIZE;
            }
        }
    }

    // Find the index by bisection.
    while (imin < imax)
    {
        const std::size_t imid = imin + ((imax - imin) >> 1);
        if (key >= array[imid])
        {
            imin = imid + 1;
        }
        else
        {
            imax = imid;
        }
    }

    if (imin == 0)
    {
        code[0] = -1;
    }
    return imin - 1;
}

/* Interpolate `dependent_values_ptr` at desired_x_ptr[0]; the result is written to result_ptr[0].
 *
 * provided_j_ptr[0] seeds the search: pass an index near the expected interval for near-sequential queries, else 0.
 * A NaN query or an empty domain returns NaN, as numpy.interp does.
 */
inline void cf_interp(
        double* desired_x_ptr,
        double* x_domain_ptr,
        double* dependent_values_ptr,
        std::size_t len_x,
        std::size_t* provided_j_ptr,
        double* result_ptr)
{
    const double desired_x = desired_x_ptr[0];
    if ((len_x == 0) || std::isnan(desired_x))
    {
        result_ptr[0] = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    if (len_x == 1)
    {
        result_ptr[0] = dependent_values_ptr[0];
        return;
    }

    // Matches numpy.interp's default left/right behavior.
    if (desired_x <= x_domain_ptr[0])
    {
        result_ptr[0] = dependent_values_ptr[0];
        return;
    }
    if (desired_x >= x_domain_ptr[len_x - 1])
    {
        result_ptr[0] = dependent_values_ptr[len_x - 1];
        return;
    }

    std::size_t j;
    if (len_x == 2)
    {
        // Only one interval; the search routine needs len >= 3.
        j = 0;
    }
    else
    {
        const std::size_t guess = (provided_j_ptr[0] < len_x) ? provided_j_ptr[0] : 0;
        int code = 0;
        j = cf_binary_search_with_guess(desired_x, x_domain_ptr, len_x, guess, &code);
        if (code == -1)
        {
            result_ptr[0] = dependent_values_ptr[0];
            return;
        }
        if (j >= len_x - 1)
        {
            result_ptr[0] = dependent_values_ptr[len_x - 1];
            return;
        }
    }

    const double xp_j = x_domain_ptr[j];
    const double fp_j = dependent_values_ptr[j];
    if (xp_j == desired_x)
    {
        result_ptr[0] = fp_j;
        return;
    }
    const double xp_jp1 = x_domain_ptr[j + 1];
    const double fp_jp1 = dependent_values_ptr[j + 1];
    const double slope  = (fp_jp1 - fp_j) / (xp_jp1 - xp_j);

    double result = slope * (desired_x - xp_j) + fp_j;
    // NaN from one direction: try the other, as numpy does.
    if (std::isnan(result))
    {
        result = slope * (desired_x - xp_jp1) + fp_jp1;
        if (std::isnan(result) && (fp_jp1 == fp_j))
        {
            result = fp_j;
        }
    }
    result_ptr[0] = result;
}

/* Complex variant of `cf_interp`: `dependent_values_ptr` holds interleaved (real, imag) pairs, and the result's real
 * and imaginary parts are written to result_ptr[0] and result_ptr[1].
 */
inline void cf_interp_complex(
        double desired_x,
        double* x_domain_ptr,
        double* dependent_values_ptr,
        std::size_t len_x,
        std::size_t* provided_j_ptr,
        double* result_ptr)
{
    const std::complex<double>* values_ptr = reinterpret_cast<const std::complex<double>*>(dependent_values_ptr);
    std::complex<double> result;

    if ((len_x == 0) || std::isnan(desired_x))
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        result = std::complex<double>(nan, nan);
    }
    else if ((len_x == 1) || (desired_x <= x_domain_ptr[0]))
    {
        result = values_ptr[0];
    }
    else if (desired_x >= x_domain_ptr[len_x - 1])
    {
        result = values_ptr[len_x - 1];
    }
    else
    {
        std::size_t j    = 0;
        int code         = 0;
        bool at_endpoint = false;
        if (len_x > 2)
        {
            const std::size_t guess = (provided_j_ptr[0] < len_x) ? provided_j_ptr[0] : 0;
            j = cf_binary_search_with_guess(desired_x, x_domain_ptr, len_x, guess, &code);
            if (code == -1)
            {
                result      = values_ptr[0];
                at_endpoint = true;
            }
            else if (j >= len_x - 1)
            {
                result      = values_ptr[len_x - 1];
                at_endpoint = true;
            }
        }

        if (!at_endpoint)
        {
            const double xp_j = x_domain_ptr[j];
            if (xp_j == desired_x)
            {
                result = values_ptr[j];
            }
            else
            {
                // Real and imaginary parts are interpolated independently.
                const double inv_dx              = 1.0 / (x_domain_ptr[j + 1] - xp_j);
                const std::complex<double> slope = (values_ptr[j + 1] - values_ptr[j]) * inv_dx;
                result = slope * (desired_x - xp_j) + values_ptr[j];
            }
        }
    }
    result_ptr[0] = result.real();
    result_ptr[1] = result.imag();
}
