#pragma once
/* Power-series arithmetic shared by the tabulated eccentricity and obliquity functions.
 *
 * Dependency-free, so any extension can include it.
 */


// A run of power-series coefficients c_0 .. c_(count - 1).
struct c_SeriesSpan
{
    const double* coefficients;
    int count;
};


// sum_j c_j x^j over the first num_terms coefficients (Horner).
inline double c_series_horner(const double* coefficients, int num_terms, double x) noexcept
{
    double total = 0.0;
    for (int j = num_terms - 1; j >= 0; --j)
    {
        total = total * x + coefficients[j];
    }
    return total;
}


// x^power for a small non-negative integer power.
inline double c_series_power(double x, int power) noexcept
{
    double result = 1.0;
    for (int i = 0; i < power; ++i) { result *= x; }
    return result;
}


// The Cauchy product of two series cut after the x^max_s term, evaluated by Horner in x.
inline double c_series_cut_product(c_SeriesSpan series_a, c_SeriesSpan series_b, int max_s, double x) noexcept
{
    double total = 0.0;
    for (int s = max_s; s >= 0; --s)
    {
        double coefficient = 0.0;
        const int i_min = (s - (series_b.count - 1) > 0) ? (s - (series_b.count - 1)) : 0;
        const int i_max = (s < series_a.count - 1) ? s : (series_a.count - 1);
        for (int i = i_min; i <= i_max; ++i)
        {
            coefficient += series_a.coefficients[i] * series_b.coefficients[s - i];
        }
        total = total * x + coefficient;
    }
    return total;
}
