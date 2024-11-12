
#include "solutions_.hpp"

RadialSolutionStorageCC::RadialSolutionStorageCC(
        size_t num_slices,
        char num_ytypes) :
            num_slices(num_slices),
            num_ytypes(num_ytypes)
{
    // Initialize status
    this->success = false;
    this->set_message("RadialSolutionStorageCC has not been updated with data yet.");

    // Find total data size
    this->total_size = (size_t)MAX_NUM_Y_REAL * this->num_slices * (size_t)this->num_ytypes;

    // Reserve space in the storage vectors.
    // These are double vectors but store double complex values so we need to double the amount of storage.
    // This is done already by using MAX_NUM_Y_REAL over MAX_NUM_Y
    this->full_solution_vec.reserve(this->total_size);
    this->full_solution_ptr = &this->full_solution_vec[0];

    // Three Love numbers are stored for each requested y-type.
    // These are also double vectors storing double complex values, double the storage.
    this->complex_love_vec.reserve(2 * 3 * this->num_ytypes);
    this->complex_love_ptr = &this->complex_love_vec[0];

    // The radial solution class will also store the output of the equation of state solver. 
    // Setup the eos solution vector and associated pointers.
    // The vector will be 7 x total slices
    // Array 0: Gravity
    // Array 1: Pressure
    // Array 2: Density
    // Double Array 3a + 3b: Shear Mod (real, imag)
    // Double Array 4a + 4b: Bulk Mod (real, imag)
    this->eos_properties_vec.reserve(7 * this->num_slices);
    this->eos_properties_ptr = &this->full_solution_vec[0];
    double* gravity_ptr   = &this->eos_properties_ptr[0];
    double* pressure_ptr  = &this->eos_properties_ptr[this->num_slices];
    double* density_ptr   = &this->eos_properties_ptr[2 * this->num_slices];
    double* shear_mod_ptr = &this->eos_properties_ptr[3 * this->num_slices];
    double* bulk_mod_ptr  = &this->eos_properties_ptr[5 * this->num_slices]; // The jump from 3 to 5 here is because shear_mod_ptr (and this bulk_mod_ptr) is twice as long as the others.
}

void RadialSolutionStorageCC::find_love(double surface_gravity)
{   
    if (this->success) [[likely]]
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
            find_love_cf(&this->complex_love_ptr[2 * 3 * ytype_i], surface_solutions_ptr, surface_gravity);
        }
    }
    else
    {
        // Could pass a new message to update the state but it will overwrite any error message that is already there.
        // TODO: Think about doing an append or logging system in the future.
        // this->set_message("Can not update Love number values when solution is not complete or successful.")
    }
}

void RadialSolutionStorageCC::set_message(const char* new_message_ptr)
{
    std::strcpy(this->message_ptr, new_message_ptr);
}
