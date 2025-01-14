
#include "rs_solution_.hpp"

RadialSolutionStorageCC::RadialSolutionStorageCC( )
{
}

RadialSolutionStorageCC::RadialSolutionStorageCC(
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
    this->eos_solution_uptr = std::make_unique<EOSSolutionCC>(
        upper_radius_bylayer_ptr,
        num_layers,
        radius_array_ptr,
        this->num_slices
        );
    
    // Setup diagnostic array. Size = max possible number of solutions per layer (3) * num layers
    this->shooting_method_steps_taken_vec.resize(3 * this->num_layers);
    // Set steps taken to 0
    for (size_t layer_i = 0; layer_i < this->num_layers; layer_i++)
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
        
        this->set_message("Radial solution storage initialized successfully.");
    }
    else
    {
        // Equation of state solution storage could not be initialized.
        this->error_code = -1;
        this->set_message("RadialSolutionStorageCC:: Could not initialize equation of state storage.");
    }
}

RadialSolutionStorageCC::~RadialSolutionStorageCC( )
{
    this->eos_solution_uptr.reset();
}

EOSSolutionCC* RadialSolutionStorageCC::get_eos_solution_ptr()
{
    return this->eos_solution_uptr.get();
}

void RadialSolutionStorageCC::change_radius_array(
        double* new_radius_array_ptr,
        size_t new_size_radius_array,
        bool array_changed)
{
    if (this->error_code == 0)
    {
        if (array_changed)
        {
            // If the array has changed since class initialization then we need to pass the new array to the equation 
            // of state solver.
            if (this->eos_solution_uptr.get())
            {
                this->eos_solution_uptr->change_radius_array(new_radius_array_ptr, new_size_radius_array);
            }

            this->set_message("Radius array changed. Radial solution reset.");
            this->success = false;
        }

        // Update properties held by this class instance
        this->num_slices = new_size_radius_array;
        this->total_size = (size_t)MAX_NUM_Y_REAL * this->num_slices * (size_t)this->num_ytypes;

        // Reserve space in the storage vectors.
        // These are double vectors but store double complex values so we need to double the amount of storage.
        // This is done already by using MAX_NUM_Y_REAL over MAX_NUM_Y
        this->full_solution_vec.resize(this->total_size);

        // Three Love numbers are stored for each requested y-type.
        // These are also double vectors storing double complex values, double the storage.
        this->complex_love_vec.resize(2 * 3 * this->num_ytypes);
    }
}

void RadialSolutionStorageCC::set_message(const char* new_message_ptr)
{
    std::strcpy(this->message_ptr, new_message_ptr);
}


void RadialSolutionStorageCC::find_love()
{   
    // Uses the equation of state results to calculate the Love numbers.
    if (this->success && this->eos_solution_uptr->success && this->error_code==0) [[likely]]
    {
        const size_t top_slice_i   = this->num_slices - 1;
        const size_t num_output_ys = MAX_NUM_Y_REAL * this->num_ytypes;

        // For each y-type, use the solution at the surface of the planet to find each Love number.
        double surface_solutions[MAX_NUM_Y_REAL] = { };
        double* surface_solutions_ptr            = &surface_solutions[0];

        for (size_t ytype_i = 0; ytype_i < this->num_ytypes; ytype_i++)
        {
            // Pull out surface solutions for this y-type
            for (size_t y_i = 0; y_i < MAX_NUM_Y_REAL; y_i++)
            {
                const size_t lhs_y_index = ytype_i * MAX_NUM_Y_REAL + y_i;
                surface_solutions_ptr[y_i] = this->full_solution_vec[top_slice_i * num_output_ys + lhs_y_index];
            }
            double* love_ptr = &this->complex_love_vec[2 * 3 * ytype_i];
            find_love_cf(
                love_ptr,
                surface_solutions_ptr,
                this->eos_solution_uptr->surface_gravity);
        }
    }
    else
    {
        // Could pass a new message to update the state but it will overwrite any error message that is already there.
        // TODO: Think about doing an append or logging system in the future.
        // this->set_message("Can not update Love number values when solution is not complete or is unsuccessful.")
    }
}

void RadialSolutionStorageCC::dimensionalize_data(
            NonDimensionalScalesCC* nondim_scales,
            bool redimensionalize)
{

    // Perform dimensionalization on the EOS solution first.
    double* full_solution_ptr       = this->full_solution_vec.data();
    EOSSolutionCC* eos_solution_ptr = this->get_eos_solution_ptr();
    eos_solution_ptr->dimensionalize_data(nondim_scales, redimensionalize);

    const double displacement_scale = (nondim_scales->second2_conversion / nondim_scales->length_conversion);
    const double stress_scale       = (nondim_scales->mass_conversion / nondim_scales->length3_conversion);
    const double potential_scale    = (1. / nondim_scales->length_conversion);

    // Redimensionalize the radial solutions
    if (this->success)
    {
        for (size_t solver_i = 0; solver_i < this->num_ytypes; solver_i++)
        {
            const size_t bc_stride = solver_i * MAX_NUM_Y_REAL;
            for (size_t slice_i = 0; slice_i < this->num_slices; slice_i++)
            {
                const size_t slice_stride = bc_stride + slice_i * MAX_NUM_Y_REAL * this->num_ytypes;
                //    y1, y3 are the radial and tangential displacements with units of [s2 m-1]
                //    y2, y4 are the radial and tangential stresses with units of [kg m-3]
                //    y5 is the tidal potential which is unitless and thus needs no conversion.
                //    y6 is a "potential stress" with units of [m-1]
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

                // y5 (real and imag)
                // full_solution_ptr[slice_stride + 8] *= 1.0;
                // full_solution_ptr[slice_stride + 9] *= 1.0;

                // y6 (real and imag)
                full_solution_ptr[slice_stride + 10] *= potential_scale;
                full_solution_ptr[slice_stride + 11] *= potential_scale;
            }
        }
    }
}