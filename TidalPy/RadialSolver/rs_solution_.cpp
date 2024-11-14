
#include "rs_solution_.hpp"


RadialSolutionStorageCC::RadialSolutionStorageCC(
        char num_ytypes,
        double* upper_radius_bylayer_ptr,
        const size_t num_layers,
        double* radius_array_ptr,
        const size_t radius_array_size) :
            success(false),
            num_ytypes(num_ytypes),
            num_layers(num_layers),
            num_slices(radius_array_size)
{
    // Find total data size
    this->total_size = (size_t)MAX_NUM_Y_REAL * this->num_slices * (size_t)this->num_ytypes;

    // Create equation of state class instance
    this->eos_solution_sptr = std::make_shared<EOSSolutionCC>(
        upper_radius_bylayer_ptr,
        num_layers,
        radius_array_ptr,
        radius_array_size
        );
    
    if (this->eos_solution_sptr.get())
    {
        // Setup equation of state pointers
        this->radius_array_ptr   = &this->eos_solution_sptr->radius_array_vec[0];
        this->gravity_array_ptr  = &this->eos_solution_sptr->gravity_array_vec[0];
        this->pressure_array_ptr = &this->eos_solution_sptr->pressure_array_vec[0];
        this->mass_array_ptr     = &this->eos_solution_sptr->mass_array_vec[0];
        this->moi_array_ptr      = &this->eos_solution_sptr->moi_array_vec[0];
        this->density_array_ptr  = &this->eos_solution_sptr->density_array_vec[0];
        // These complex arrays are stored as double arrays with twice the length (Cython and C++ don't play nicely with complex across all systems)
        this->complex_shear_array_ptr = &this->eos_solution_sptr->complex_shear_array_vec[0];
        this->complex_bulk_array_ptr  = &this->eos_solution_sptr->complex_bulk_array_vec[0];
    }

    // Reserve space in the storage vectors.
    // These are double vectors but store double complex values so we need to double the amount of storage.
    // This is done already by using MAX_NUM_Y_REAL over MAX_NUM_Y
    this->full_solution_vec.resize(this->total_size);
    this->full_solution_ptr = &this->full_solution_vec[0];

    // Three Love numbers are stored for each requested y-type.
    // These are also double vectors storing double complex values, double the storage.
    this->complex_love_vec.resize(2 * 3 * this->num_ytypes);
    this->complex_love_ptr = &this->complex_love_vec[0];
}


void RadialSolutionStorageCC::set_message(const char* new_message_ptr)
{
    std::strcpy(this->message_ptr, new_message_ptr);
}


void RadialSolutionStorageCC::find_love()
{   
    // Uses the equation of state results to calculate the Love numbers.
    if (this->success && this->eos_solution_sptr->success) [[likely]]
    {
        const size_t top_slice_i   = this->num_slices - 1;
        const size_t num_output_ys = MAX_NUM_Y_REAL * this->num_ytypes;

        // For each y-type, use the solution at the surface of the planet to find each Love number.
        double surface_solutions[MAX_NUM_Y_REAL] = { };
        double* surface_solutions_ptr            = &surface_solutions[0];

        for (char ytype_i = 0; ytype_i < this->num_ytypes; ytype_i++)
        {
            // Pull out surface solutions for this y-type
            for (int y_i = 0; y_i < MAX_NUM_Y_REAL; y_i++)
            {
                const size_t lhs_y_index = ytype_i * MAX_NUM_Y_REAL + y_i;
                surface_solutions_ptr[y_i] = this->full_solution_ptr[top_slice_i * num_output_ys + lhs_y_index];
            }
            find_love_cf(
                &this->complex_love_ptr[2 * 3 * ytype_i],
                surface_solutions_ptr,
                this->eos_solution_sptr->surface_gravity);
        }
    }
    else
    {
        // Could pass a new message to update the state but it will overwrite any error message that is already there.
        // TODO: Think about doing an append or logging system in the future.
        // this->set_message("Can not update Love number values when solution is not complete or is unsuccessful.")
    }
}
