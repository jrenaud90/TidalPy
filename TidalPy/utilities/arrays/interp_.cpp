#include <cmath>
#include <complex>

size_t cf_binary_search_with_guess(
        double key,
        double* array,
        size_t length,
        size_t guess,
        int* code)
{
    const size_t LIKELY_IN_CACHE_SIZE = 8;
    
    // Set code to default
    code[0] = 0;

    if (key > array[length - 1]){
        return length;
    }
    else if (key < array[0]){
        // Tell output that the key was to the left of the array.
        code[0] = -1;
        return 0;
    }

    if (guess > (length - 3)){
        guess = length - 3;
    }
    if (guess < 1) {
        guess = 1;
    }

    size_t imin = 0;
    size_t imax = length;
    size_t imid = 0;

    /* Check most likely values: guess - 1, guess, guess + 1 */
    if (key < array[guess]){
        if (key < array[guess - 1]){
            imax = guess - 1;
            /* last attempt to restrict search to items in cache */
            if ((guess > LIKELY_IN_CACHE_SIZE) && (key >= array[guess - LIKELY_IN_CACHE_SIZE])){
                imin = guess - LIKELY_IN_CACHE_SIZE;
            }
        }
        else {
            return guess - 1;
        }
    }
    else {
        if (key < array[guess + 1]){
            return guess;
        }
        else {
            if (key < array[guess + 2]){
                return guess + 1;
            }
            else {
                imin = guess + 2;
                /* last attempt to restrict search to items in cache */
                if ((guess < (length - LIKELY_IN_CACHE_SIZE - 1)) && (key < array[guess + LIKELY_IN_CACHE_SIZE])){
                    imax = guess + LIKELY_IN_CACHE_SIZE;
                }
            }
        }
    }
    /* Finally, find index by bisection */
    while (imin < imax){
        imid = imin + ((imax - imin) >> 1);
        if (key >= array[imid]){
            imin = imid + 1;
        }
        else {
            imax = imid;
        }
    }

    if (imin == 0)
    {
        code[0] = -1;
    }
    return imin - 1;
}

void cf_interp(
        double* desired_x_ptr,
        double* x_domain_ptr,
        double* dependent_values_ptr,
        size_t len_x,
        size_t* provided_j_ptr,
        double* result_ptr
        )
{
    /*
    Interpolation function for floats.

    Provided a domain, `x_domain` and a dependent array `dependent_values` search domain for value closest to 
    `desired_x` and return the value of `dependent_values` at that location if it is defined. Otherwise, use local 
    slopes of `x_domain` and `dependent_values` to interpolate a value of `dependent_values` at `desired_x`.

    Based on `numpy`'s `interp` function.

    Parameters
    ----------
    desired_x_ptr : double*
        Location where `dependent_variables` is desired.
    x_domain : double*
        Domain to search for the correct location.
    dependent_values : double*
        Dependent values that are to be returned after search and interpolation.
    provided_j : int*, input & output
        Give a j index from a previous interpolation to improve performance.
    rseult_ptr : double*, output
        Desired value of `dependent_values`.
    */
    
    // TODO: Needs to be at least 3 item long array. Add exception here?

    const double left_value  = dependent_values_ptr[0];
    const double right_value = dependent_values_ptr[len_x - 1];
    const double left_x      = x_domain_ptr[0];
    const double right_x     = x_domain_ptr[len_x - 1];

    // Binary Search with Guess
    int b_search_code = 0;
    size_t j = 0;
        

    // Perform binary search with guess
    if (provided_j_ptr[0] >= 0)
    {
        // j was provided, use it for later calcs.
        j = provided_j_ptr[0];
    }
    else
    {
        // j was not provided (or it was provided as an unsupported negative) so we need to find it
        // Get a guess for where in the x-array we might be.
        const size_t j_guess = len_x * (size_t)std::fabs(std::floor(desired_x_ptr[0] / (right_x - left_x)));
        j = cf_binary_search_with_guess(desired_x_ptr[0], x_domain_ptr, len_x, j_guess, &b_search_code);
    }
    
    double slope;
    double fp_at_j;
    double xp_at_j;
    double fp_at_jp1;
    double xp_at_jp1;

    if (b_search_code == -1)
    {
        result_ptr[0] = left_value;
    }
    else if (j >= len_x)
    {
        result_ptr[0] = right_value;
    }
    else
    {
        fp_at_j = dependent_values_ptr[j];
        xp_at_j = x_domain_ptr[j];
        if (j == (len_x - 1))
        {
            result_ptr[0] = fp_at_j;
        }
        else if (xp_at_j == desired_x_ptr[0])
        {
            result_ptr[0] = fp_at_j;
        }
        else
        {
            fp_at_jp1 = dependent_values_ptr[j + 1];
            xp_at_jp1 = x_domain_ptr[j + 1];
            slope = (fp_at_jp1 - fp_at_j) / (xp_at_jp1 - xp_at_j);
            result_ptr[0] = slope * (desired_x_ptr[0] - xp_at_j) + fp_at_j;

            // If we get nan in one direction, try the other
            if (std::isnan<double>(result_ptr[0]))
            {
                result_ptr[0] = slope * (desired_x_ptr[0] - xp_at_jp1) + fp_at_jp1;
                if (std::isnan<double>(result_ptr[0]) && (fp_at_jp1 == fp_at_j))
                {
                    result_ptr[0] = fp_at_j;
                }
            }
        }
    }
}

void cf_interp_complex(
        double desired_x,
        double* x_domain_ptr,
        double* dependent_values_ptr,
        size_t len_x,
        size_t* provided_j_ptr,
        double* result_ptr
        )
{
    /*
    Interpolation function for floats.

    Provided a domain, `x_domain` and a dependent array `dependent_values` search domain for value closest to 
    `desired_x` and return the value of `dependent_values` at that location if it is defined. Otherwise, use local 
    slopes of `x_domain` and `dependent_values` to interpolate a value of `dependent_values` at `desired_x`.

    Based on `numpy`'s `interp` function.

    Parameters
    ----------
    desired_x : double
        Location where `dependent_variables` is desired.
    x_domain_ptr : double*
        Domain to search for the correct location.
    dependent_values_ptr : double complex*
        Dependent values that are to be returned after search and interpolation.
    len_x : size_t 
        Size of `x_domain_ptr`.
    provided_j_ptr : size_t*, input & output
        Give a j index from a previous interpolation to improve performance.
    rseult_ptr : double complex*, output
        Desired value of `dependent_values`.
    */
    
    // TODO: Needs to be at least 3 item long array. Add exception here?

    const std::complex<double> left_value  = 
        std::complex<double>(dependent_values_ptr[0], dependent_values_ptr[1]);
    const std::complex<double> right_value = 
        std::complex<double>(dependent_values_ptr[2 * len_x - 2], dependent_values_ptr[2 * len_x - 1]);
    const double left_x  = x_domain_ptr[0];
    const double right_x = x_domain_ptr[len_x - 1];

    // Binary Search with Guess
    int b_search_code = 0;
    size_t j = 0;

    // Perform binary search with guess
    if (provided_j_ptr[0] >= 0)
    {
        // j was provided, use it for later calcs.
        j = provided_j_ptr[0];
    }
    else
    {
        // j was not provided (or it was provided as an unsupported negative) so we need to find it
        // Get a guess for where in the x-array we might be.
        size_t j_guess = len_x * (size_t)std::fabs(std::floor(desired_x / (right_x - left_x)));
        j_guess = std::max<size_t>(std::min<size_t>(j_guess, len_x), 0);
        j = cf_binary_search_with_guess(desired_x, x_domain_ptr, len_x, j_guess, &b_search_code);
    }
    
    double slope_real;
    double slope_imag;
    double x_slope_inverse;

    double result_real;
    double result_imag;
    double fp_at_j_real;
    double fp_at_j_imag;
    double xp_at_j;
    double fp_at_jp1_real;
    double fp_at_jp1_imag;
    double xp_at_jp1;

    if (b_search_code == -1)
    {
        result_ptr[0] = left_value.real();
        result_ptr[1] = left_value.imag();
    }
    else if (j >= len_x)
    {
        result_ptr[0] = right_value.real();
        result_ptr[1] = right_value.imag();
    }
    else
    {
        fp_at_j_real = dependent_values_ptr[2 * j];
        fp_at_j_imag = dependent_values_ptr[2 * j + 1];
        xp_at_j      = x_domain_ptr[j];
        if (j == (len_x - 1))
        {
            result_ptr[0] = fp_at_j_real;
            result_ptr[0] = fp_at_j_imag;
        }
        else if (xp_at_j == desired_x)
        {
            result_ptr[0] = fp_at_j_real;
            result_ptr[0] = fp_at_j_imag;
        }
        else
        {
            fp_at_jp1_real  = dependent_values_ptr[2 * (j + 1)];
            fp_at_jp1_imag  = dependent_values_ptr[2 * (j + 1) + 1];
            xp_at_jp1       = x_domain_ptr[j + 1];
            x_slope_inverse = 1.0 / (xp_at_jp1 - xp_at_j);
            slope_real = (fp_at_jp1_real - fp_at_j_real) * x_slope_inverse;
            slope_imag = (fp_at_jp1_imag - fp_at_j_imag) * x_slope_inverse;

            // If we get nan in one direction try the other
            // Real Part
            result_ptr[0] = slope_real * (desired_x - xp_at_j) + fp_at_j_real;
            if (std::isnan<double>(result_ptr[0]))
            {
                result_ptr[0] = slope_real * (desired_x - xp_at_jp1) + fp_at_jp1_real;
                if (std::isnan<double>(result_ptr[0]) && (fp_at_jp1_real == fp_at_j_real))
                {
                    result_ptr[0] = fp_at_j_real;
                }
            }

            // Imaginary Part
            result_ptr[1] = slope_imag * (desired_x - xp_at_j) + fp_at_j_imag;
            if (std::isnan<double>(result_ptr[1]))
            {
                result_ptr[1] = slope_imag * (desired_x - xp_at_jp1) + fp_at_jp1_imag;
                if (std::isnan<double>(result_ptr[1]) && (fp_at_jp1_imag == fp_at_j_imag))
                {
                    result_ptr[1] = fp_at_j_imag;
                }
            }
        }
    }
}
