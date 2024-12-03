#pragma once

#include <cstdio>
#include <cstring>
#include <vector>
#include <memory>

#include "cysolution.hpp"  // Part of the CyRK python package. Header should be included in setup using CyRK.get_include()

#include "nondimensional_.hpp" // Part of the TidalPy.utilities module

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
    int iterations             = -1;
    int error_code             = -100;
    int nondim_status          = 0;
    int solution_nondim_status = 0;
    bool success          = false;
    bool max_iters_hit    = false;
    bool radius_array_set = false;
    bool other_vecs_set   = false;

    char* message_ptr           = &message[0];
    size_t current_layers_saved = 0;
    size_t num_layers           = 0;
    size_t radius_array_size    = 0;
    size_t num_cyolver_calls    = 0;

    double pressure_error   = NULL;
    double surface_gravity  = NULL;
    double surface_pressure = NULL;
    double central_pressure = NULL;
    double radius           = NULL;
    double mass             = NULL;
    double moi              = NULL;

    double redim_length_scale  = NULL;
    double redim_gravity_scale = NULL;
    double redim_mass_scale    = NULL;
    double redim_density_scale = NULL;
    double redim_moi_scale     = NULL;
    double redim_pascal_scale  = NULL;
    
    // Store results from CyRK's cysolve_ivp.
    std::vector<double> upper_radius_bylayer_vec = std::vector<double>();
    std::vector<std::shared_ptr<CySolverResult>> cysolver_results_sptr_bylayer_vec = std::vector<std::shared_ptr<CySolverResult>>();
    std::vector<size_t> steps_taken_vec = std::vector<size_t>();

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

    virtual ~EOSSolutionCC( );
    EOSSolutionCC( );
    EOSSolutionCC(
            double* upper_radius_bylayer_ptr,
            size_t num_layers,
            double* radius_array_ptr,
            size_t radius_array_size
        );

    // Prepare storage vectors for interpolation
    void change_radius_array(
        double* new_radius_ptr,
        size_t new_radius_size);

    // Save result of cysolve_ivp integration
    void save_cyresult(std::shared_ptr<CySolverResult> new_cysolver_result_sptr);
    void save_steps_taken(size_t steps_taken);
    
    // Call a specific layer's interpolate function.
    void call(
        const size_t layer_index,
        const double radius,
        double* y_interp_ptr) const;
    void call_vectorize(
        const size_t layer_index,
        const double* radius_array_ptr,
        size_t len_radius_array,
        double* y_interp_ptr) const;
    
    // Run full planet interpolation through each layer.
    void interpolate_full_planet();
    void dimensionalize_data(
        NonDimensionalScalesCC* nondim_scales,
        bool redimensionalize);
};