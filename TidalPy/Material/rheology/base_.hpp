#pragma once

#include <cmath>

static const int MAX_NUM_PARAMETERS = 16;

class RheologyModelBaseCC
{
// Public Attributes
public:
    char rheology_model_int = -1;
    char num_parameters = 0;
    char* rheology_name_ptr = &rheology_name[0];

// Protected Attributes
protected:
    double* parameter_storage_ptr = &parameter_storage[0];
    double parameter_storage[MAX_NUM_PARAMETERS] = { };
    char rheology_name[32] = { };

// Public Methods
public:
    virtual ~RheologyModelBaseCC() {};
    RheologyModelBaseCC();
    RheologyModelBaseCC(
        const double* parameters,
        int num_parameters, 
        char rheology_model_int,
        char* rheology_name);

    virtual void update_parameters(
        const double* parameters);

    virtual void call(
        double* output_ptr,
        double frequency,
        double modulus,
        double viscosity) const;
    
    void vectorize_frequency(
        double* output_ptr,
        size_t output_size,
        double* frequency_ptr,
        double modulus, 
        double viscosity) const;

    void vectorize_modulus_viscosity(
        double* output_ptr,
        size_t output_size,
        double frequency,
        double* modulus_ptr, 
        double* viscosity_ptr) const;
    
    void vectorize_all(
        double* output_ptr,
        size_t output_size,
        double* frequency_ptr,
        double* modulus_ptr, 
        double* viscosity_ptr) const;

// Protected Methods
protected:
};