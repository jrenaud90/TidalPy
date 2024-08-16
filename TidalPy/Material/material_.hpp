
#pragma once

#include <stdio.h>
#include <memory>
#include "rheology\base_.hpp"
#include "rheology\models_.hpp"

class MaterialCC
{
public:
    // State properties
    double pressure    = 0.0;
    double temperature = 0.0;

    // Physical properties
    double gravity = 0.0;
    double density = 0.0;

    // Viscoelastic properties (set by EOS)
    double shear_viscosity      = 0.0;
    double bulk_viscosity       = 0.0;
    double static_shear_modulus = 0.0;
    double static_bulk_modulus  = 0.0;

    // Rheological properties (set by rheology)
    double complex_shear_modulus[2]   = {0.0, 0.0};
    double* complex_shear_modulus_ptr = &complex_shear_modulus[0];
    int shear_rheology_model_int = -1;
    std::shared_ptr<RheologyModelBaseCC> shear_rheology_ptr = nullptr;
    double complex_bulk_modulus[2]    = {0.0, 0.0};
    double* complex_bulk_modulus_ptr  = &complex_shear_modulus[0];
    int bulk_rheology_model_int = -1;
    std::shared_ptr<RheologyModelBaseCC> bulk_rheology_ptr = nullptr;

    virtual ~MaterialCC();
    MaterialCC() {};

    void build_shear_rheology(int rheology_model_int, const double* rheology_parameters);
    void build_bulk_rheology(int rheology_model_int, const double* rheology_parameters);
};