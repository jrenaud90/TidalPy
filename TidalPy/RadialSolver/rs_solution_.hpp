#pragma once

#include <cstring>
#include <complex>
#include <vector>
#include <memory>

#include "eos_solution_.hpp"
#include "love_.hpp"

const int MAX_NUM_Y = 6;
const int MAX_NUM_Y_REAL = 12; // Maximum number of y values counting both the real and imaginary portions.


// Error Codes:
// -1 : Equation of State storage (EOSSolutionCC) could not be initialized.
// -2 : (set by python wrapper) Unknown / Unsupported boundary condition provided.

class RadialSolutionStorageCC
{
// Attributes
protected:
    char message[256] = { };

public:
    bool success      = false;
    int error_code    = 0;
    char num_ytypes   = 0;
    char* message_ptr = &message[0];
    size_t num_slices = 0;
    size_t num_layers = 0;
    size_t total_size = 0;
    
    // Equation of state solution
    std::shared_ptr<EOSSolutionCC> eos_solution_sptr = std::make_shared<EOSSolutionCC>();

    // Radial solution results
    std::vector<double> full_solution_vec = std::vector<double>();

    // Love number attributes
    std::vector<double> complex_love_vec = std::vector<double>();

    // Raw pointers to radial solver storage
    double* full_solution_ptr = nullptr;
    double* complex_love_ptr  = nullptr;

    // Raw pointers to equation of state storage
    double* radius_array_ptr   = nullptr;
    double* gravity_array_ptr  = nullptr;
    double* pressure_array_ptr = nullptr;
    double* mass_array_ptr     = nullptr;
    double* moi_array_ptr      = nullptr;
    double* density_array_ptr  = nullptr;
    // These complex arrays are stored as double arrays with twice the length (Cython and C++ don't play nicely with complex across all systems)
    double* complex_shear_array_ptr = nullptr;
    double* complex_bulk_array_ptr  = nullptr;

    // Constructors and methods
    virtual ~RadialSolutionStorageCC() { };
    RadialSolutionStorageCC() { };
    RadialSolutionStorageCC(
        char num_ytypes,
        double* upper_radius_bylayer_ptr,
        const size_t num_layers,
        double* radius_array_ptr,
        const size_t radius_array_size);
    
    void change_radius_array(
        double* radius_array_ptr,
        const size_t radius_array_size,
        cpp_bool array_changed)
    void set_message(const char* new_message);
    void find_love();
};
