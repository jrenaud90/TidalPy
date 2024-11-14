#pragma once

#include <cstdio>
#include <vector>
#include <memory>
#include <string>

#include "cysolution.hpp"  // Part of the CyRK python package. Header should be included in setup using CyRK.get_include()

static const size_t EOS_Y_VALUES = 4;
static const size_t EOS_EXTRA_VALUES = 5;
static const size_t EOS_DY_VALUES = 9;
// There are 4 dependent y values in the EOS ODE:
// 0 : gravity
// 1 : pressure
// 2 : total mass
// 3 : moment of inertia
// Additionally there are 5 extra parameters that are tracked.
// 4 : density
// 5 : Real part of complex shear modulus
// 6 : Imag part of complex shear modulus
// 7 : Real part of complex bulk modulus
// 8 : Imag part of complex bulk modulus

class EOSSolutionCC
{

// Attributes
protected:
    char message[256] = { };

public:
    int iterations        = -1;
    bool success          = false;
    bool max_iters_hit    = false;
    bool radius_array_set = false;
    bool other_vecs_set   = false;

    char* message_ptr           = &message[0];
    size_t current_layers_saved = 0;
    size_t num_layers           = 0;
    size_t radius_array_size    = 0;

    double pressure_error   = NULL;
    double surface_gravity  = NULL;
    double surface_pressure = NULL;
    double mass             = NULL;
    double moi              = NULL;
    
    // Store results from CyRK's cysolve_ivp.
    std::vector<double> upper_radius_bylayer_vec = std::vector<double>();
    std::vector<std::shared_ptr<CySolverResult>> cysolver_results_bylayer_vec = std::vector<std::shared_ptr<CySolverResult>>();

    // Copy of user-provided radius array
    std::vector<double> radius_array_vec = std::vector<double>();

    // Store interpolated arrays based on provided radius array
    std::vector<double> gravity_array_vec  = std::vector<double>();
    std::vector<double> pressure_array_vec = std::vector<double>();
    std::vector<double> mass_array_vec     = std::vector<double>();
    std::vector<double> moi_array_vec      = std::vector<double>();
    std::vector<double> density_array_vec  = std::vector<double>();

    // These complex arrays are stored as double arrays with twice the length (Cython and C++ don't play nicely with complex across all systems)
    std::vector<double> complex_shear_array_vec = std::vector<double>();
    std::vector<double> complex_bulk_array_vec  = std::vector<double>();

// Methods
protected:

public:

    virtual ~EOSSolutionCC( ) { };
    EOSSolutionCC( ) { };
    EOSSolutionCC(
        double* upper_radius_bylayer_ptr,
        const size_t num_layers,
        double* radius_array_ptr,
        const size_t radius_array_size
        );

    // Save result of cysolve_ivp integration
    void save_cyresult(
        std::shared_ptr<CySolverResult> new_cysolver_result_sptr);
    
    // Call a specific layer's interpolate function.
    void call(
        const size_t layer_index,
        const double radius,
        double* y_interp_ptr);
    void call_vectorize(
        const size_t layer_index,
        const double* radius_array_ptr,
        size_t len_radius_array,
        double* y_interp_ptr);
    
    // Prepare storage vectors for interpolation
    void rest_radius_array(
        double* radius_array_ptr,
        const size_t radius_array_size);
    
    // Run full planet interpolation through each layer.
    void interpolate_full_planet();
}