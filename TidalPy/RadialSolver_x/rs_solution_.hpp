// rs_solution_.hpp - Radial solver solution storage class
// Ported from TidalPy/RadialSolver/rs_solution_.hpp + rs_solution_.cpp
#pragma once

#include <cstring>
#include <complex>
#include <vector>
#include <memory>
#include <string>

#include "love_.hpp"
#include "rs_constants_.hpp"
#include "../Material_x/eos/eos_solution_.hpp"
#include "../../constants_.hpp"
#include "../utilities/dimensions/nondimensional_.hpp"


// Error Codes:
// -1 : Equation of State storage (c_EOSSolution) could not be initialized.
// -2 : (set by python wrapper) Unknown / Unsupported boundary condition provided.
// -5 : There was a problem with the inputs to radial solver
//
// -1X : Error in shooting method
// -10 : Error in finding starting conditions
// -11 : Numerical integration failed
// -12 : Error using ZGESV solver with boundary condition
//
// -2X : Error in propagation matrix method
// -20 : Unknown core starting conditions
// -21 : Error using ZGESV solver with boundary condition

class c_RadialSolutionStorage
{
public:
    bool success        = false;
    int error_code      = -100;
    int degree_l        = 0;
    std::string message = "No Message Set.";
    size_t num_ytypes   = 0;
    size_t num_slices   = 0;
    size_t num_layers   = 0;
    size_t total_size   = 0;

    // Equation of state solution
    std::unique_ptr<c_EOSSolution> eos_solution_uptr = nullptr;

    // Radial solution results (stores double-pairs for complex values)
    std::vector<double> full_solution_vec = std::vector<double>();

    // Love number attributes (stores double-pairs for complex values)
    std::vector<c_LoveNumbers> complex_love_vec = std::vector<c_LoveNumbers>();

    // Diagnostic data
    std::vector<size_t> shooting_method_steps_taken_vec = std::vector<size_t>();

    // Default constructor
    c_RadialSolutionStorage() = default;

    // Main constructor
    c_RadialSolutionStorage(
        size_t num_ytypes,
        double* upper_radius_bylayer_ptr,
        size_t num_layers,
        double* radius_array_ptr,
        size_t size_radius_array,
        int degree_l) :
            success(false),
            error_code(0),
            degree_l(degree_l),
            num_ytypes(num_ytypes),
            num_slices(size_radius_array),
            num_layers(num_layers)
    {
        // Create equation of state class instance
        this->eos_solution_uptr = std::make_unique<c_EOSSolution>(
            upper_radius_bylayer_ptr,
            num_layers,
            radius_array_ptr,
            this->num_slices
            );

        // Setup diagnostic array. Size = max possible number of solutions per layer (3) * num layers
        this->shooting_method_steps_taken_vec.resize(3 * this->num_layers);
        for (size_t layer_i = 0; layer_i < this->num_layers; ++layer_i)
        {
            this->shooting_method_steps_taken_vec[3 * layer_i]     = 0;
            this->shooting_method_steps_taken_vec[3 * layer_i + 1] = 0;
            this->shooting_method_steps_taken_vec[3 * layer_i + 2] = 0;
        }

        // Setup radius array based vectors
        if (this->eos_solution_uptr.get())
        {
            this->change_radius_array(
                radius_array_ptr,
                size_radius_array,
                false  // We are in initialization, this is not an array change.
                );

            this->message = "Radial solution storage initialized successfully.";
        }
        else
        {
            this->error_code = -1;
            this->message = "c_RadialSolutionStorage:: Could not initialize equation of state storage.";
        }
    }

    virtual ~c_RadialSolutionStorage()
    {
        this->eos_solution_uptr.reset();
    }

    c_EOSSolution* get_eos_solution_ptr()
    {
        return this->eos_solution_uptr.get();
    }

    void change_radius_array(
        double* new_radius_array_ptr,
        size_t new_size_radius_array,
        bool array_changed)
    {
        if (this->error_code == 0)
        {
            if (array_changed)
            {
                if (this->eos_solution_uptr.get())
                {
                    this->eos_solution_uptr->change_radius_array(new_radius_array_ptr, new_size_radius_array);
                }

                this->message = "Radius array changed. Radial solution reset.";
                this->success = false;
            }

            this->num_slices = new_size_radius_array;
            this->total_size = static_cast<size_t>(C_MAX_NUM_Y_REAL) * this->num_slices * this->num_ytypes;

            this->full_solution_vec.resize(this->total_size);

            // Love number structure for each ytype
            this->complex_love_vec.resize(this->num_ytypes);
        }
    }

    void find_love()
    {
        if (this->success && this->eos_solution_uptr->success && this->error_code == 0) [[likely]]
        {
            const size_t top_slice_i   = this->num_slices - 1;
            const size_t num_output_ys = C_MAX_NUM_Y_REAL * this->num_ytypes;

            // For each y-type, use the solution at the surface to find each Love number.
            // The full_solution_vec stores values as double pairs (real, imag alternating).
            // We need to convert to std::complex<double> for c_find_love.
            std::complex<double> surface_solutions[C_MAX_NUM_Y];
            std::complex<double> love_numbers[3];

            for (size_t ytype_i = 0; ytype_i < this->num_ytypes; ++ytype_i)
            {
                // Pull out surface solutions for this y-type
                for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i)
                {
                    const size_t lhs_y_index = ytype_i * C_MAX_NUM_Y_REAL + y_i * 2;
                    double real_part = this->full_solution_vec[top_slice_i * num_output_ys + lhs_y_index];
                    double imag_part = this->full_solution_vec[top_slice_i * num_output_ys + lhs_y_index + 1];
                    surface_solutions[y_i] = std::complex<double>(real_part, imag_part);
                }

                this->complex_love_vec[ytype_i] = c_find_love(
                    surface_solutions,
                    this->eos_solution_uptr->surface_gravity
                );
            }
        }
    }

    void dimensionalize_data(
        c_NonDimensionalScales* nondim_scales,
        bool redimensionalize)
    {
        // Perform dimensionalization on the EOS solution first.
        double* full_solution_ptr       = this->full_solution_vec.data();
        c_EOSSolution* eos_solution_ptr = this->get_eos_solution_ptr();
        eos_solution_ptr->dimensionalize_data(nondim_scales, redimensionalize);

        const double displacement_scale = (nondim_scales->second2_conversion / nondim_scales->length_conversion);
        const double stress_scale       = (nondim_scales->mass_conversion / nondim_scales->length3_conversion);
        const double potential_scale    = (1.0 / nondim_scales->length_conversion);

        if (this->success)
        {
            for (size_t solver_i = 0; solver_i < this->num_ytypes; ++solver_i)
            {
                const size_t bc_stride = solver_i * C_MAX_NUM_Y_REAL;
                for (size_t slice_i = 0; slice_i < this->num_slices; ++slice_i)
                {
                    const size_t slice_stride = bc_stride + slice_i * C_MAX_NUM_Y_REAL * this->num_ytypes;
                    // y1 (real and imag)
                    full_solution_ptr[slice_stride + 0] *= displacement_scale;
                    full_solution_ptr[slice_stride + 1] *= displacement_scale;

                    // y3 (real and imag)
                    full_solution_ptr[slice_stride + 4] *= displacement_scale;
                    full_solution_ptr[slice_stride + 5] *= displacement_scale;

                    // y2 (real and imag)
                    full_solution_ptr[slice_stride + 2] *= stress_scale;
                    full_solution_ptr[slice_stride + 3] *= stress_scale;

                    // y4 (real and imag)
                    full_solution_ptr[slice_stride + 6] *= stress_scale;
                    full_solution_ptr[slice_stride + 7] *= stress_scale;

                    // y5 is unitless - no conversion needed

                    // y6 (real and imag)
                    full_solution_ptr[slice_stride + 10] *= potential_scale;
                    full_solution_ptr[slice_stride + 11] *= potential_scale;
                }
            }
        }
    }
};
