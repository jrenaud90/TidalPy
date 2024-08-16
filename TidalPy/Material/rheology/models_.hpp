#pragma once

#include <cmath>
#include <algorithm>
#include <complex>
#include <numbers>

#include <tpy_constants.hpp>
#include "base_.hpp"

/*
=======================================================================================================================
Elastic Rheology
=======================================================================================================================
*/
static const char ELASTIC_RHEOLOGY_INT   = 0;
static const char ELASTIC_NUM_PARAMETERS = 0;
static const double ELASTIC_DEFAULTS[1]  = {
    NAN_DBL // No real parameters, set the pointer equal to something though so it is pointing to real memory.
};
static const double* ELASTIC_DEFAULTS_PTR = &ELASTIC_DEFAULTS[0];

class Elastic : public RheologyModelBaseCC
{
// Public Methods
public:
    virtual ~Elastic() {};
    Elastic() : RheologyModelBaseCC(ELASTIC_DEFAULTS_PTR, ELASTIC_NUM_PARAMETERS, ELASTIC_RHEOLOGY_INT, "elastic") {};
    Elastic(const double* parameters_ptr) : RheologyModelBaseCC(parameters_ptr, ELASTIC_NUM_PARAMETERS, ELASTIC_RHEOLOGY_INT, "elastic") {};

    virtual void call(
        double* output_ptr,
        double frequency,
        double modulus,
        double viscosity) const override;
};

/*
=======================================================================================================================
Newton Rheology
=======================================================================================================================
*/
static const char NEWTON_RHEOLOGY_INT   = 1;
static const char NEWTON_NUM_PARAMETERS = 0;
static const double NEWTON_DEFAULTS[1]  = {
    NAN_DBL // No real parameters, set the pointer equal to something though so it is pointing to real memory.
};
static const double* NEWTON_DEFAULTS_PTR = &NEWTON_DEFAULTS[0];

class Newton : public RheologyModelBaseCC
{
// Public Methods
public:
    virtual ~Newton() {};
    Newton() : RheologyModelBaseCC(NEWTON_DEFAULTS_PTR, NEWTON_NUM_PARAMETERS, NEWTON_RHEOLOGY_INT, "newton") {};
    Newton(const double* parameters_ptr) : RheologyModelBaseCC(parameters_ptr, NEWTON_NUM_PARAMETERS, NEWTON_RHEOLOGY_INT, "newton") {};

    virtual void call(
        double* output_ptr,
        double frequency,
        double modulus,
        double viscosity) const override;
};

/*
=======================================================================================================================
Maxwell Rheology
=======================================================================================================================
*/
static const char MAXWELL_RHEOLOGY_INT   = 2;
static const char MAXWELL_NUM_PARAMETERS = 0;
static const double MAXWELL_DEFAULTS[1]  = {
    NAN_DBL // No real parameters, set the pointer equal to something though so it is pointing to real memory.
};
static const double* MAXWELL_DEFAULTS_PTR = &MAXWELL_DEFAULTS[0];

class Maxwell : public RheologyModelBaseCC
{
// Public Methods
public:
    virtual ~Maxwell() {};
    Maxwell() : RheologyModelBaseCC(MAXWELL_DEFAULTS_PTR, MAXWELL_NUM_PARAMETERS, MAXWELL_RHEOLOGY_INT, "maxwell") {};
    Maxwell(const double* parameters_ptr) : RheologyModelBaseCC(parameters_ptr, MAXWELL_NUM_PARAMETERS, MAXWELL_RHEOLOGY_INT, "maxwell") {};

    virtual void call(
        double* output_ptr,
        const double frequency,
        const double modulus,
        const double viscosity) const override;
};

/*
=======================================================================================================================
Voigt-Kelvin Rheology
=======================================================================================================================
*/
static const char VOIGT_RHEOLOGY_INT   = 3;
static const char VOIGT_NUM_PARAMETERS = 2;
static const double VOIGT_DEFAULTS[2]  = {
     5.0,  // Voigt Modulus Scaling Factor
     0.02  // Voigt Viscosity Scaling Factor
};
static const double* VOIGT_DEFAULTS_PTR = &VOIGT_DEFAULTS[0];

class Voigt : public RheologyModelBaseCC
{
// Public Methods
public:
    virtual ~Voigt() {};
    Voigt() : RheologyModelBaseCC(VOIGT_DEFAULTS_PTR, VOIGT_NUM_PARAMETERS, VOIGT_RHEOLOGY_INT, "voigt-kelvin") {};
    Voigt(const double* parameters_ptr) : RheologyModelBaseCC(parameters_ptr, VOIGT_NUM_PARAMETERS, VOIGT_RHEOLOGY_INT, "voigt-kelvin") {};

    virtual void call(
        double* output_ptr,
        const double frequency,
        const double modulus,
        const double viscosity) const override;
};

/*
=======================================================================================================================
Burgers Rheology
=======================================================================================================================
*/
static const char BURGERS_RHEOLOGY_INT   = 4;
static const char BURGERS_NUM_PARAMETERS = 2;
static const double BURGERS_DEFAULTS[2]  = {
     5.0,  // Voigt Modulus Scaling Factor
     0.02  // Voigt Viscosity Scaling Factor
};
static const double* BURGERS_DEFAULTS_PTR = &BURGERS_DEFAULTS[0];

class Burgers : public RheologyModelBaseCC
{
// Public Methods
public:
    virtual ~Burgers() {};
    Burgers() : RheologyModelBaseCC(BURGERS_DEFAULTS_PTR, BURGERS_NUM_PARAMETERS, BURGERS_RHEOLOGY_INT, "burgers") {};
    Burgers(const double* parameters_ptr) : RheologyModelBaseCC(parameters_ptr, BURGERS_NUM_PARAMETERS, BURGERS_RHEOLOGY_INT, "burgers") {};

    virtual void call(
        double* output_ptr,
        const double frequency,
        const double modulus,
        const double viscosity) const override;
};

/*
=======================================================================================================================
Andrade Rheology
=======================================================================================================================
*/
static const char ANDRADE_RHEOLOGY_INT   = 5;
static const char ANDRADE_NUM_PARAMETERS = 2;
static const double ANDRADE_DEFAULTS[2]  = {
     0.3,  // Andrade Alpha (exponential parameter)
     1.0   // Andrade Zeta (Andrade to Maxwell timescale ratio)
};
static const double* ANDRADE_DEFAULTS_PTR = &ANDRADE_DEFAULTS[0];

class Andrade : public RheologyModelBaseCC
{
// Public Methods
public:
    virtual ~Andrade() {};
    Andrade() : RheologyModelBaseCC(ANDRADE_DEFAULTS_PTR, ANDRADE_NUM_PARAMETERS, ANDRADE_RHEOLOGY_INT, "andrade") {};
    Andrade(const double* parameters_ptr) : RheologyModelBaseCC(parameters_ptr, ANDRADE_NUM_PARAMETERS, ANDRADE_RHEOLOGY_INT, "andrade") {};

    virtual void update_parameters(
        const double* parameters) override;

    virtual void call(
        double* output_ptr,
        const double frequency,
        const double modulus,
        const double viscosity) const override;
};

/*
=======================================================================================================================
Sundberg-Cooper Rheology
=======================================================================================================================
*/
static const char SUNDBERG_RHEOLOGY_INT   = 6;
static const char SUNDBERG_NUM_PARAMETERS = 4;
static const double SUNDBERG_DEFAULTS[4]  = {
     5.0,   // Voigt Modulus Scaling Factor
     0.02,  // Voigt Viscosity Scaling Factor
     0.3,   // Andrade Alpha (exponential parameter)
     1.0    // Andrade Zeta (Andrade to Maxwell timescale ratio)
};
static const double* SUNDBERG_DEFAULTS_PTR = &SUNDBERG_DEFAULTS[0];

class SundbergCooper : public RheologyModelBaseCC
{
// Public Methods
public:
    virtual ~SundbergCooper() {};
    SundbergCooper() : RheologyModelBaseCC(SUNDBERG_DEFAULTS_PTR, SUNDBERG_NUM_PARAMETERS, SUNDBERG_RHEOLOGY_INT, "sundberg-cooper") {};
    SundbergCooper(const double* parameters_ptr) : RheologyModelBaseCC(parameters_ptr, SUNDBERG_NUM_PARAMETERS, SUNDBERG_RHEOLOGY_INT, "sundberg-cooper") {};

    virtual void update_parameters(
        const double* parameters) override;

    virtual void call(
        double* output_ptr,
        const double frequency,
        const double modulus,
        const double viscosity) const override;
};