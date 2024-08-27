#pragma once

#include <cstring>
#include <complex>
#include <vector>

#include "love_.hpp"

const int MAX_NUM_Y = 6;
const int MAX_NUM_Y_REAL = 12; // Maximum number of y values counting both the real and imaginary portions.


class RadialSolutionStorageCC
{
    // Public Attributes
public:
    size_t num_slices;
    size_t total_size;
    char message[256] = { };
    char* message_ptr = &message[0];
    bool success      = false;
    char num_ytypes   = 0;

    std::vector<double> complex_love_vec = std::vector<double>(0);
    double* complex_love_ptr = NULL;
    std::vector<double> full_solution_vec = std::vector<double>(0);
    double* full_solution_ptr = NULL;

    RadialSolutionStorageCC() { };
    virtual ~RadialSolutionStorageCC() { };
    RadialSolutionStorageCC(size_t num_slices, char num_ytypes);
    
    void set_message(const char* new_message);
    void find_love(double surface_gravity);
};
