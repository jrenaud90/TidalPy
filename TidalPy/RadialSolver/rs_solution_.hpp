#pragma once

#include <cstring>
#include <complex>
#include <vector>
#include <memory>

#include "eos_solution_.hpp"
#include "love_.hpp"
#include "nondimensional_.hpp"

const size_t MAX_NUM_Y      = 6;
const size_t MAX_NUM_Y_REAL = 12; // Maximum number of y values counting both the real and imaginary portions.


// Error Codes:
// -1 : Equation of State storage (EOSSolutionCC) could not be initialized.
// -2 : (set by python wrapper) Unknown / Unsupported boundary condition provided.
// -5 : There was a problem with the inputs to radial solver

// -1X : Error in shooting method
// -10 : Error in finding starting conditions
// -11 : Numerical integration failed
// -12 : Error using ZGESV solver with boundary condition
// 
// -2X : Error in propagation matrix method
// -20 : Unknown core starting conditions
// -21 : Error using ZGESV solver with boundary condition

class RadialSolutionStorageCC
{
// Attributes
protected:
    char message[256] = { };

public:
    bool success      = false;
    int error_code    = -100;
    int degree_l      = 0;
    char* message_ptr = &message[0];
    size_t num_ytypes = 0;
    size_t num_slices = 0;
    size_t num_layers = 0;
    size_t total_size = 0;
    
    // Equation of state solution
    std::unique_ptr<EOSSolutionCC> eos_solution_uptr = nullptr;

    // Radial solution results
    std::vector<double> full_solution_vec = std::vector<double>();

    // Love number attributes
    std::vector<double> complex_love_vec = std::vector<double>();

    // Diagnostic data
    std::vector<size_t> shooting_method_steps_taken_vec = std::vector<size_t>();

    // Constructors and methods
    virtual ~RadialSolutionStorageCC();
    RadialSolutionStorageCC();
    RadialSolutionStorageCC(
        size_t num_ytypes,
        double* upper_radius_bylayer_ptr,
        size_t num_layers,
        double* radius_array_ptr,
        size_t size_radius_array,
        int degree_l);
    
    EOSSolutionCC* get_eos_solution_ptr();
    void change_radius_array(
        double* new_radius_array_ptr,
        size_t new_size_radius_array,
        bool array_changed);
    void set_message(const char* new_message);
    void find_love();
    void dimensionalize_data(
        NonDimensionalScalesCC* nondim_scales,
        bool redimensionalize);
};
