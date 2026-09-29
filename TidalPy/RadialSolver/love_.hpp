#pragma once

#include <cmath>
#include <complex>

#include "../constants_.hpp"    // TidalPyConstants (d_INF, d_PI)


struct c_LoveNumbers
{
    std::complex<double> k;
    std::complex<double> h;
    std::complex<double> l;

    c_LoveNumbers() :
        k(0.0, 0.0),
        h(0.0, 0.0),
        l(0.0, 0.0)
    {
    }

    c_LoveNumbers(const std::complex<double>& k_, const std::complex<double>& h_, const std::complex<double>& l_) :
        k(k_),
        h(h_),
        l(l_)
    {
    }

    c_LoveNumbers(const double k_, const double h_, const double l_) :
        k(k_, 0.0),
        h(h_, 0.0),
        l(l_, 0.0)
    {
    }

    double get_Q_k() const { return c_LoveNumbers::p_quality_factor(this->k); }
    double get_Q_h() const { return c_LoveNumbers::p_quality_factor(this->h); }
    double get_Q_l() const { return c_LoveNumbers::p_quality_factor(this->l); }

    double get_lag_k() const { return c_LoveNumbers::p_lag(this->k); }
    double get_lag_h() const { return c_LoveNumbers::p_lag(this->h); }
    double get_lag_l() const { return c_LoveNumbers::p_lag(this->l); }

private:
    // Sign of the elastic part of a Love number: -1 for a negative real part (loading k' and h'), else +1. A
    // response lags its forcing when its imaginary part has the opposite sign to its real part (Im k < 0 for a
    // tidal k), so the quality factor and lag below take this sign to be positive for a dissipative response of
    // either sign.
    static double p_real_sign(const std::complex<double>& love_number)
    {
        return (std::real(love_number) < 0.0) ? -1.0 : 1.0;
    }

    // Quality factor -s |n| / Im(n) of a complex Love number n, with s the sign of Re(n); infinite for a purely
    // elastic response.
    static double p_quality_factor(const std::complex<double>& love_number)
    {
        const double love_abs  = std::abs(love_number);
        const double love_imag = std::imag(love_number);

        if (love_imag == 0.0) [[unlikely]]
        {
            return TidalPyConstants::d_INF;
        }
        else
        {
            return -c_LoveNumbers::p_real_sign(love_number) * love_abs / love_imag;
        }
    }

    // Phase lag arctan2(-s Im(n), |Re(n)|) of a complex Love number n, with s the sign of Re(n); in
    // [-pi/2, pi/2], pi/2 for a purely imaginary n with Im(n) < 0.
    static double p_lag(const std::complex<double>& love_number)
    {
        const double love_real = std::real(love_number);
        const double love_imag = std::imag(love_number);

        if (love_imag == 0.0)
        {
            return 0.0;
        }
        else
        {
            return std::atan2(-c_LoveNumbers::p_real_sign(love_number) * love_imag, std::abs(love_real));
        }
    }
};


/// Love and Shida numbers from the y-values at the planet surface, in the Tobie et al. (2005) sign
/// convention for y5.
///
/// References
/// ----------
/// HH14: Henning & Hurford (2014), Eq. A9
/// RN08: Roberts & Nimmo (2008), Eq. A8
/// T05:  Tobie et al. (2005), Eqs. 9 & 36
///
/// Parameters
/// ----------
/// surface_solutions_ptr : std::complex<double>*
///     y-values at the planet surface: [y1, y2, y3, y4, y5, y6].
/// surface_gravity : double
///     Gravitational acceleration at the surface [m s-2].
///
/// Returns
/// -------
/// c_LoveNumbers : the complex k, h, and l.
inline c_LoveNumbers c_find_love(
        std::complex<double>* surface_solutions_ptr,
        double surface_gravity
        ) noexcept
{
    const std::complex<double> k = surface_solutions_ptr[4] - 1.0;              // k = y5 - 1
    const std::complex<double> h = surface_solutions_ptr[0] * surface_gravity;  // h = y1 * g
    const std::complex<double> l = surface_solutions_ptr[2] * surface_gravity;  // l = y3 * g

    return c_LoveNumbers(k, h, l);
}
