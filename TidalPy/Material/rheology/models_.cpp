#include "models_.hpp"

/*
=======================================================================================================================
Elastic Rheology
=======================================================================================================================
*/
void Elastic::call(double* output_ptr, const double frequency, const double modulus, const double viscosity) const
{
    // Elastic rheology just returns the static modulus as the real part and 0 as the imaginary.
    output_ptr[0] = modulus;
    output_ptr[1] = 0.0;
}

/*
=======================================================================================================================
Newton Rheology
=======================================================================================================================
*/
void Newton::call(double* output_ptr, const double frequency, const double modulus, const double viscosity) const
{
    // TODO: Should frequency be abs here? Assuming so as the others are.
    const double frequency_abs = std::abs(frequency);

    // Check for extreme values of frequency. If found, use pre-calculated limits.
    if (frequency_abs < MIN_FREQUENCY) [[unlikely]]
    {
        output_ptr[0] = 0.0;
        output_ptr[1] = 0.0;
    }
    else if (frequency_abs > MAX_FREQUENCY) [[unlikely]]
    {
        output_ptr[0] = 0.0;
        output_ptr[1] = INF_DBL;
    }
    else [[likely]]
    {
        // Newtonian fluid
        output_ptr[0] = 0.0;
        output_ptr[1] = frequency_abs * viscosity;
    }
}

/*
=======================================================================================================================
Maxwell Rheology
=======================================================================================================================
*/
void Maxwell::call(double* output_ptr, const double frequency, const double modulus, const double viscosity) const
{
    const double frequency_abs = std::abs(frequency);

    // Check for extreme values of frequency. If found, use pre-calculated limits.
    if (frequency_abs < MIN_FREQUENCY) [[unlikely]]
    {
        output_ptr[0] = 0.0;
        output_ptr[1] = 0.0;
    }
    else if ((frequency_abs > MAX_FREQUENCY) || std::isinf<double>(frequency_abs)) [[unlikely]]
    {
        output_ptr[0] = modulus;
        output_ptr[1] = 0.0;
    }
    else if (modulus < MIN_MODULUS) [[unlikely]]
    {
        output_ptr[0] = 0.0;
        output_ptr[1] = 0.0;
    }
    else [[likely]]
    {
        // Maxwell rheology
        const double maxwell_time = viscosity / modulus;
        const std::complex<double> denom  = std::complex<double>(frequency_abs * maxwell_time, -1.0);
        const std::complex<double> result = (viscosity * frequency_abs) / denom;
        output_ptr[0] = result.real();
        output_ptr[1] = result.imag();
    }
}

/*
=======================================================================================================================
Voigt-Kelvin Rheology
=======================================================================================================================
*/
void Voigt::call(double* output_ptr, const double frequency, const double modulus, const double viscosity) const
{

    const double voigt_modulus   = this->parameter_storage_ptr[0] * modulus;
    const double voigt_viscosity = this->parameter_storage_ptr[1] * viscosity;

    const double frequency_abs = std::abs(frequency);

    // Check for extreme values of frequency. If found, use pre-calculated limits.
    if (frequency_abs < MIN_FREQUENCY) [[unlikely]]
    {
        output_ptr[0] = voigt_modulus;
        output_ptr[1] = 0.0;
    }
    else if ((frequency_abs > MAX_FREQUENCY) || std::isinf<double>(frequency_abs)) [[unlikely]]
    {
        output_ptr[0] = 0.0;
        output_ptr[1] = INF_DBL;
    }
    else if (modulus < MIN_MODULUS) [[unlikely]]
    {
        output_ptr[0] = 0.0;
        output_ptr[1] = voigt_viscosity * frequency_abs;
    }
    else [[likely]]
    {
        // Voigt-Kelvin rheology
        const std::complex<double> result = std::complex<double>(voigt_modulus, frequency_abs * voigt_viscosity);
        output_ptr[0] = result.real();
        output_ptr[1] = result.imag();
    }
}

/*
=======================================================================================================================
Burgers Rheology
=======================================================================================================================
*/
void Burgers::call(double* output_ptr, const double frequency, const double modulus, const double viscosity) const
{

    const double voigt_modulus   = this->parameter_storage_ptr[0] * modulus;
    const double voigt_viscosity = this->parameter_storage_ptr[1] * viscosity;

    const double frequency_abs = std::abs(frequency);

    // Check for extreme values of frequency. If found, use pre-calculated limits.
    if (frequency_abs < MIN_FREQUENCY) [[unlikely]]
    {
        output_ptr[0] = 0.0;
        output_ptr[1] = 0.0;
    }
    else if ((frequency_abs > MAX_FREQUENCY) || std::isinf<double>(frequency_abs)) [[unlikely]]
    {
        output_ptr[0] = modulus;
        output_ptr[1] = 0.0;
    }
    else if (modulus < MIN_MODULUS) [[unlikely]]
    {
        output_ptr[0] = 0.0;
        output_ptr[1] = 0.0;
    }
    else [[likely]]
    {
        // Burgers rheology (combination of Voigt-Kelvin and Maxwell)
        const double voigt_time   = voigt_viscosity / voigt_modulus;
        const double maxwell_time = viscosity / modulus;
        const double maxwell_parm = maxwell_time * frequency_abs;
        
        const std::complex<double> voigt_param = std::complex<double>(frequency_abs * voigt_time, -1.0);
        const std::complex<double> denom       = maxwell_parm * voigt_param + std::complex<double>(-1.0, -frequency_abs * (viscosity + voigt_viscosity) / voigt_modulus);
        const std::complex<double> result      = (viscosity * frequency_abs * voigt_param) / denom;
        output_ptr[0] = result.real();
        output_ptr[1] = result.imag();
    }
}

/*
=======================================================================================================================
Andrade Rheology
=======================================================================================================================
*/

void Andrade::update_parameters(const double* parameters_ptr)
{
    // Andrade has a few other constants that only change when parameters are updated. We will precalculate those here
    // for performance.
    
    // First, call base class to load in parameters.
    RheologyModelBaseCC::update_parameters(parameters_ptr);

    // Parameters are:
    // 0 : Andrade Alpha
    // 1 : Andrade Zeta

    // Pull out Andrade alpha and zeta
    const double alpha = this->parameter_storage_ptr[0];
    const double zeta  = this->parameter_storage_ptr[1];

    // Calculate other constants and put them back into storage.
    this->parameter_storage_ptr[2] = std::tgamma(alpha + 1.0);              // alpha_factorial
    this->parameter_storage_ptr[3] = std::cos(PI_DBL * alpha / 2.);         // sine_term_real
    this->parameter_storage_ptr[4] = -1.0 * std::sin(PI_DBL * alpha / 2.);  // sine_term_imag
}

void Andrade::call(double* output_ptr, const double frequency, const double modulus, const double viscosity) const
{

    const double alpha = this->parameter_storage_ptr[0];
    const double zeta  = this->parameter_storage_ptr[1];
    const double alpha_factorial = this->parameter_storage_ptr[2];
    const std::complex<double> sine_term = std::complex<double>(this->parameter_storage_ptr[3], this->parameter_storage_ptr[4]);

    const double frequency_abs = std::abs(frequency);

    // Check for extreme values of frequency. If found, use pre-calculated limits.
    if (frequency_abs < MIN_FREQUENCY) [[unlikely]]
    {
        output_ptr[0] = 0.0;
        output_ptr[1] = 0.0;
    }
    else if ((frequency_abs > MAX_FREQUENCY) || std::isinf<double>(frequency_abs)) [[unlikely]]
    {
        output_ptr[0] = modulus;
        output_ptr[1] = 0.0;
    }
    else if (modulus < MIN_MODULUS) [[unlikely]]
    {
        output_ptr[0] = 0.0;
        output_ptr[1] = 0.0;
    }
    else [[likely]]
    {
        // Andrade rheology (including Maxwell component)
        const double maxwell_time = viscosity / modulus;
        const double maxwell_parm = maxwell_time * frequency_abs;
        const double andrade_term = std::pow(maxwell_parm * zeta, alpha);
        
        const std::complex<double> denom  = maxwell_parm * alpha_factorial * sine_term + std::complex<double>(maxwell_parm * andrade_term, -andrade_term);
        const std::complex<double> result = (viscosity * frequency_abs * andrade_term) / denom;
        output_ptr[0] = result.real();
        output_ptr[1] = result.imag();
    }
}

/*
=======================================================================================================================
Sundberg-Cooper Rheology
=======================================================================================================================
*/

void SundbergCooper::update_parameters(const double* parameters_ptr)
{
    // The andrade portion of S-C has a few other constants that only change when parameters are updated.
    // We will precalculate those here for performance.
    
    // First, call base class to load in parameters.
    RheologyModelBaseCC::update_parameters(parameters_ptr);

    // Parameters are:
    // 0 : Voigt Modulus Scale
    // 1 : Voigt Viscosity Scale
    // 2 : Andrade Alpha
    // 3 : Andrade Zeta

    // Pull out Andrade alpha and zeta
    const double alpha = this->parameter_storage_ptr[2];
    const double zeta  = this->parameter_storage_ptr[3];

    // Calculate other constants and put them back into storage.
    this->parameter_storage_ptr[4] = std::tgamma(alpha + 1.0);              // alpha_factorial
    this->parameter_storage_ptr[5] = std::cos(PI_DBL * alpha / 2.);         // sine_term_real
    this->parameter_storage_ptr[6] = -1.0 * std::sin(PI_DBL * alpha / 2.);  // sine_term_imag
}

void Andrade::call(double* output_ptr, const double frequency, const double modulus, const double viscosity) const
{

    const double voigt_modulus   = this->parameter_storage_ptr[0] * modulus;
    const double voigt_viscosity = this->parameter_storage_ptr[1] * viscosity;
    const double alpha = this->parameter_storage_ptr[2];
    const double zeta  = this->parameter_storage_ptr[3];
    const double alpha_factorial = this->parameter_storage_ptr[4];
    const std::complex<double> sine_term = std::complex<double>(this->parameter_storage_ptr[5], this->parameter_storage_ptr[6]);

    const double frequency_abs = std::abs(frequency);

    // Check for extreme values of frequency. If found, use pre-calculated limits.
    if (frequency_abs < MIN_FREQUENCY) [[unlikely]]
    {
        output_ptr[0] = 0.0;
        output_ptr[1] = 0.0;
    }
    else if ((frequency_abs > MAX_FREQUENCY) || std::isinf<double>(frequency_abs)) [[unlikely]]
    {
        output_ptr[0] = modulus;
        output_ptr[1] = 0.0;
    }
    else if (modulus < MIN_MODULUS) [[unlikely]]
    {
        output_ptr[0] = 0.0;
        output_ptr[1] = 0.0;
    }
    else [[likely]]
    {
        // Andrade rheology (including Maxwell component)
        const double voigt_time   = voigt_viscosity / voigt_modulus;
        const double maxwell_time = viscosity / modulus;
        const double maxwell_parm = maxwell_time * frequency_abs;
        const double andrade_term = std::pow(maxwell_parm * zeta, alpha);
        
        const std::complex<double> voigt_param = std::complex<double>(frequency_abs * voigt_time, -1.0);
        const std::complex<double> sundberg_term = std::complex<double>(0.0, -1.0) * (andrade_term * voigt_param + andrade_term * maxwell_parm / this->parameter_storage_ptr[0]);
        const std::complex<double> denom = 
            maxwell_parm * alpha_factorial * sine_term * voigt_param + 
            maxwell_parm * andrade_term * voigt_param +
            sundberg_term;

        const std::complex<double> result = (viscosity * frequency_abs * andrade_term * voigt_param) / denom;
        output_ptr[0] = result.real();
        output_ptr[1] = result.imag();
    }
}
