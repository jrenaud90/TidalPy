#include "eos_solution_.hpp"


EOSSolutionCC::EOSSolutionCC(
        double* upper_radius_bylayer_ptr,
        const size_t num_layers,
        double* radius_array_ptr,
        const size_t radius_array_size
        ) :
            error_code(0),
            num_layers(num_layers)
{
    this->cysolver_results_bylayer_vec.reserve(num_layers);
    this->upper_radius_bylayer_vec.reserve(num_layers);

    // Store the upper radius of each layer to make it easier to call interpolators later
    this->upper_radius_bylayer_vec.insert(
        this->upper_radius_bylayer_vec.begin(),
        upper_radius_bylayer_ptr,
        upper_radius_bylayer_ptr + num_layers
        );

    // Use provided radius array to setup various storage vectors
    this->rest_radius_array(radius_array_ptr, radius_array_size);
}


void EOSSolutionCC::save_cyresult(
        std::shared_ptr<CySolverResult> new_cysolver_result_sptr)
{
    this->cysolver_results_bylayer_vec.push_back(new_cysolver_result_sptr);
    this->current_layers_saved++;
}


void EOSSolutionCC::call(
        const size_t layer_index,
        const double radius,
        double* y_interp_ptr)
{
    if (layer_index < this->current_layers_saved)
    {
        this->cysolver_results_bylayer_vec[layer_index]->call(radius, y_interp_ptr);
    }
    else
    {
        // TODO: Better error handling
        printf("Error! EOSSolution::call was asked to interpolate a layer that it has not saved.");
    }
}


void EOSSolutionCC::call_vectorize(
        const size_t layer_index,
        const double* radius_array_ptr,
        size_t len_radius_array,
        double* y_interp_ptr)
{
    if (layer_index < this->current_layers_saved)
    {
        this->cysolver_results_bylayer_vec[layer_index]->call_vectorize(radius_array_ptr, len_radius_array, y_interp_ptr);
    }
    else
    {
        // TODO: Better error handling
        printf("Error! EOSSolution::call was asked to interpolate a layer that it has not saved.");
    }
}


void EOSSolutionCC::change_radius_array(
        double* radius_array_ptr,
        const size_t radius_array_size)
{
    this->radius_array_size = radius_array_size;
    if (this->radius_array_set)
    {
        // Radius array was already set before. Reset the storage array vectors.
        this->radius_array_vec.clear();
        this->gravity_array_vec.clear();
        this->pressure_array_vec.clear();
        this->mass_array_vec.clear();
        this->moi_array_vec.clear();
        this->density_array_vec.clear();
        this->complex_shear_array_vec.clear();
        this->complex_bulk_array_vec.clear();

        // Indicate that all vectors are cleared.
        this->other_vecs_set = false;
    }
    this->radius_array_set = true;

    // Reserve capacity in vectors for new radius array size
    this->radius_array_vec.reserve(radius_array_size);
    this->gravity_array_vec.reserve(radius_array_size);
    this->pressure_array_vec.reserve(radius_array_size);
    this->mass_array_vec.reserve(radius_array_size);
    this->moi_array_vec.reserve(radius_array_size);
    this->density_array_vec.reserve(radius_array_size);
    
    // These complex arrays are stored as double arrays with twice the length (Cython and C++ don't play nicely with complex across all systems)
    this->complex_shear_array_vec.reserve(2 * radius_array_size);
    this->complex_bulk_array_vec.reserve(2 * radius_array_size);

    // Copy over the radius array values
    this->radius_array_vec.insert(
        this->radius_array_vec.begin(),
        radius_array_ptr,
        radius_array_ptr + radius_array_size
        );
}


void EOSSolutionCC::interpolate_full_planet()
{
    size_t current_layer_index        = 0;
    double current_layer_upper_radius = this->upper_radius_bylayer_vec[0];

    double y_interp_arr[EOS_DY_VALUES];
    double* y_interp_ptr = &y_interp_arr[0];

    for (size_t radius_i = 0; radius_i < this->radius_array_size; radius_i++)
    {
        const double radius = this->radius_array_vec[radius_i];

        if (radius > current_layer_upper_radius)
        {
            // No longer in current layer
            current_layer_index++;
            if (current_layer_index > (this->num_layers - 1))
            {
                // Outside of the planet
                break;
            }
            current_layer_upper_radius = this->upper_radius_bylayer_vec[current_layer_index];
        }

        // Call interpolate using temp array as holder.
        this->cysolver_results_bylayer_vec[current_layer_index]->call(radius, y_interp_ptr);
        
        // Store results
        this->gravity_array_vec.push_back(y_interp_ptr[0]);
        this->pressure_array_vec.push_back(y_interp_ptr[1]);
        this->mass_array_vec.push_back(y_interp_ptr[2]);
        this->moi_array_vec.push_back(y_interp_ptr[3]);
        this->density_array_vec.push_back(y_interp_ptr[4]);

        // These complex arrays are stored as double arrays with twice the length (Cython and C++ don't play nicely with complex across all systems)
        this->complex_shear_array_vec.push_back(y_interp_ptr[5]);
        this->complex_shear_array_vec.push_back(y_interp_ptr[6]);
        this->complex_bulk_array_vec.push_back(y_interp_ptr[7]);
        this->complex_bulk_array_vec.push_back(y_interp_ptr[8]);
    }

    // At the end we should be at the planet's surface. Use the last interpolated values to set surface constants
    this->surface_gravity  = y_interp_ptr[0];
    this->surface_pressure = y_interp_ptr[1];
    this->mass             = y_interp_ptr[2];
    this->moi              = y_interp_ptr[3];

    // Finished
    this->other_vecs_set = true;
}