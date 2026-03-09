#pragma once

#include <complex>
#include "constants_.hpp"       // RadialSolver_x: C_MAX_NUM_Y, etc.
#include "../constants_.hpp"    // TidalPy: TidalPyConstants (d_INF, d_PI)


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

    double get_Q_k() const
    {
        const double k_abs  = std::abs(this->k);
        const double k_imag = std::imag(this->k);

        if (k_imag == 0.0) [[unlikely]]
        {
            return TidalPyConstants::d_INF;
        }
        else
        {
            return -k_abs / k_imag;
        }
    }

    double get_Q_h() const
    {
        const double h_abs  = std::abs(this->h);
        const double h_imag = std::imag(this->h);

        if (h_imag == 0.0) [[unlikely]]
        {
            return TidalPyConstants::d_INF;
        }
        else
        {
            return -h_abs / h_imag;
        }
    }

    double get_Q_l() const
    {
        const double l_abs  = std::abs(this->l);
        const double l_imag = std::imag(this->l);

        if (l_imag == 0.0) [[unlikely]]
        {
            return TidalPyConstants::d_INF;
        }
        else
        {
            return -l_abs / l_imag;
        }
    }

    double get_lag_k() const
    {
        const double k_real = std::real(this->k);
        const double k_imag = std::imag(this->k);
        
        if (k_imag == 0.0)
        {
            return 0.0;
        }
        else if (k_real == 0.0) [[unlikely]]
        {
            // Limit of arctan(inf)
            return TidalPyConstants::d_PI / 2.0;
        }
        else
        {
            return std::atan(-k_imag / k_real);
        }
    }

    double get_lag_h() const
    {
        const double h_real = std::real(this->h);
        const double h_imag = std::imag(this->h);
        
        if (h_imag == 0.0)
        {
            return 0.0;
        }
        else if (h_real == 0.0) [[unlikely]]
        {
            // Limit of arctan(inf)
            return TidalPyConstants::d_PI / 2.0;
        }
        else
        {
            return std::atan(-h_imag / h_real);
        }
    }

    double get_lag_l() const
    {
        const double l_real = std::real(this->l);
        const double l_imag = std::imag(this->l);
        
        if (l_imag == 0.0)
        {
            return 0.0;
        }
        else if (l_real == 0.0) [[unlikely]]
        {
            // Limit of arctan(inf)
            return TidalPyConstants::d_PI / 2.0;
        }
        else
        {
            return std::atan(-l_imag / l_real);
        }
    }
};


/// Compute Love and Shida numbers from the radial solution at the planet surface.
///
/// Uses the convention of Tobie et al. (2005) for y5 sign:
///   k = y5 - 1
///   h = y1 * surface_gravity
///   l = y3 * surface_gravity
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
/// Return
/// ------
/// c_LoveNumbers complex_love_numbers
///    C++ class containing complex Love numbers.
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
