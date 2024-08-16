#include <cstring>
#include <stdio.h>
#include <stdlib.h>

#include "base_.hpp"

RheologyModelBaseCC::RheologyModelBaseCC() :
        rheology_model_int(-1),
        num_parameters(0)

{
    std::strcpy(this->rheology_name_ptr, "RheologyModelBaseCC");
}

RheologyModelBaseCC::RheologyModelBaseCC(
    const double* parameters_ptr,
    int num_parameters,
    char rheology_model_int,
    char* rheology_name) :
        rheology_model_int(rheology_model_int),
        num_parameters(num_parameters)
{
    // Initialize base class
    std::strcpy(this->rheology_name_ptr, rheology_name);

    if (num_parameters > MAX_NUM_PARAMETERS)
    {
        printf("Unsupported number of parameters in RheologyModelBaseCC: %d of a max of %d.", num_parameters, MAX_NUM_PARAMETERS);
        exit(EXIT_FAILURE);
    }
    
    // Setup rheological model parameters
    if (parameters_ptr != NULL)
    {
        this->update_parameters(parameters_ptr);
    }
}

void RheologyModelBaseCC::update_parameters(
    const double* parameters_ptr)
{
    // Make a copy of parameters and load them into the base class.
    for (int i = 0; i < this->num_parameters; i++)
    {
        this->parameter_storage_ptr[i] = parameters_ptr[i];
    }
}

void RheologyModelBaseCC::call(
    double* output_ptr,
    double frequency,
    double modulus,
    double viscosity) const
{
    // Overriden by subclasses. The base class simply returns nans.
    // `output_ptr` should point to a double complex number so it has two double members.
    output_ptr[0] = NAN;
    output_ptr[1] = NAN;
}

void RheologyModelBaseCC::vectorize_frequency(
        double* output_ptr,
        size_t output_size,
        double* frequency_ptr,
        double modulus, 
        double viscosity) const
{
    /* Vectorize over frequency */
    double result[2] = { };
    double* result_ptr = &result[0];

    // OPT: Make parallelize if `output_size` is large?
    for (size_t i = 0; i < output_size; i++)
    {
        this->call(result_ptr, frequency_ptr[i], modulus, viscosity);
        // `output_ptr` is a double pointer but it is actually pointing to double complex data
        // So we need to store two results for each call.
        output_ptr[2 * i]     = result_ptr[0];
        output_ptr[2 * i + 1] = result_ptr[1];
    }
    
}

void RheologyModelBaseCC::vectorize_modulus_viscosity(
        double* output_ptr,
        size_t output_size,
        double frequency,
        double* modulus_ptr, 
        double* viscosity_ptr) const
{
    /* Vectorize over modulus and viscosity (they must have the same lengths) */
    double result[2] = { };
    double* result_ptr = &result[0];

    // OPT: Make parallelize if `output_size` is large?
    for (size_t i = 0; i < output_size; i++)
    {
        this->call(result_ptr, frequency, modulus_ptr[i], viscosity_ptr[i]);
        // `output_ptr` is a double pointer but it is actually pointing to double complex data
        // So we need to store two results for each call.
        output_ptr[2 * i]     = result_ptr[0];
        output_ptr[2 * i + 1] = result_ptr[1];
    }
}

void RheologyModelBaseCC::vectorize_all(
        double* output_ptr,
        size_t output_size,
        double* frequency_ptr,
        double* modulus_ptr, 
        double* viscosity_ptr) const
{
    /* Vectorize over modulus and viscosity (they must have the same lengths) */
    double result[2] = { };
    double* result_ptr = &result[0];

    // OPT: Make parallelize if `output_size` is large?
    for (size_t i = 0; i < output_size; i++)
    {
        this->call(result_ptr, frequency_ptr[i], modulus_ptr[i], viscosity_ptr[i]);
        // `output_ptr` is a double pointer but it is actually pointing to double complex data
        // So we need to store two results for each call.
        output_ptr[2 * i]     = result_ptr[0];
        output_ptr[2 * i + 1] = result_ptr[1];
    }
}
