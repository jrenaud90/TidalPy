#include <complex>
#include "constants_.hpp"


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

    c_LoveNumbers(std::complex<double>& k_, std::complex<double>& h_, std::complex<double>& l_) :
        k(k_),
        h(h_),
        l(l_)
    {
    }

    c_LoveNumbers(double k_, double h_, double l_) :
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
