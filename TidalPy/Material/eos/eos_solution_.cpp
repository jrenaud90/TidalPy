#include "eos_solution_.hpp"
#include <exception>
#include <cmath>


bool isclose(double a, double b)
{
    // TODO: Make the numerics.pyx utility a c file that we can import here instead of reimplementing it here.
    const double rtol = 1.0e-9;
    const double atol = 0.0;

    if (isnan(a))
    {
        return false;
    }

    if (isnan(b))
    {
        return false;
    }

    if (a == b)
    {
        return true;
    }
    
    const double lhs = std::fabs(a - b);
    const double rhs = std::fmax(rtol * std::fmax(std::fabs(a), std::fabs(b)), atol);

    return lhs <= rhs;
}


EOSSolutionCC::EOSSolutionCC( )
{
    
}

EOSSolutionCC::EOSSolutionCC(
        double* upper_radius_bylayer_ptr,
        size_t num_layers,
        double* radius_array_ptr,
        size_t radius_array_size
        ) :
            error_code(0),
            current_layers_saved(0),
            num_layers(num_layers)
{
    this->cysolver_results_sptr_bylayer_vec.reserve(num_layers);
    this->upper_radius_bylayer_vec.resize(num_layers);

    // Store the upper radius of each layer to make it easier to call interpolators later
    std::memcpy(this->upper_radius_bylayer_vec.data(), upper_radius_bylayer_ptr, this->num_layers * sizeof(double));

    // Use provided radius array to setup various storage vectors
    this->change_radius_array(radius_array_ptr, radius_array_size);
}

EOSSolutionCC::~EOSSolutionCC( )
{
    // Reset each shared pointer in the cysolver vector
    for (size_t i = 0; i < this->cysolver_results_sptr_bylayer_vec.size(); i++)
    {
        this->cysolver_results_sptr_bylayer_vec[i]->dense_vec.clear();
        // Reset this class reference to the shared pointer.
        this->cysolver_results_sptr_bylayer_vec[i].reset();
    }
    // Clear all vectors
    this->cysolver_results_sptr_bylayer_vec.clear();
    this->upper_radius_bylayer_vec.clear();
    this->radius_array_vec.clear();
    this->gravity_array_vec.clear();
    this->pressure_array_vec.clear();
    this->mass_array_vec.clear();
    this->moi_array_vec.clear();
    this->density_array_vec.clear();
    this->complex_shear_array_vec.clear();
    this->complex_bulk_array_vec.clear();
}

void EOSSolutionCC::save_cyresult(std::shared_ptr<CySolverResult> new_cysolver_result_sptr)
{
    // We will save a copy of the shared pointer to ensure that the underlying object does not get deconstructed as long
    // as this object is alive. We will also save the raw pointer for ease of access and performance. 
    this->cysolver_results_sptr_bylayer_vec.push_back(new_cysolver_result_sptr);
    this->current_layers_saved++;
}

void EOSSolutionCC::save_steps_taken(size_t steps_taken)
{
    this->steps_taken_vec.push_back(steps_taken);
    this->num_cyolver_calls++;
}

void EOSSolutionCC::call(
        const size_t layer_index,
        const double radius,
        double* y_interp_ptr) const
{
    if (layer_index < this->current_layers_saved) [[likely]]
    {
        this->cysolver_results_sptr_bylayer_vec[layer_index]->call(radius, y_interp_ptr);

        if (this->nondim_status == 1)
        {
            // Acceleration due to Gravity
            y_interp_ptr[0] *= this->redim_gravity_scale;

            // Pressure
            y_interp_ptr[1] *= this->redim_pascal_scale;

            // These assume sphereical symmetry
            // Total mass
            y_interp_ptr[2] *= this->redim_mass_scale;

            // Moment of inertia (r^2 multipled by dm which was found above)
            y_interp_ptr[3] *= this->redim_moi_scale;

            // Density
            y_interp_ptr[4] *= this->redim_density_scale;
            
            // Shear modulus (real and complex)
            y_interp_ptr[5] *= this->redim_pascal_scale;
            y_interp_ptr[6] *= this->redim_pascal_scale;

            // Bulk modulus (real and complex)
            y_interp_ptr[7] *= this->redim_pascal_scale;
            y_interp_ptr[8] *= this->redim_pascal_scale;
        }
        else if (this->nondim_status == -1)
        {
            // Acceleration due to Gravity
            y_interp_ptr[0] /= this->redim_gravity_scale;

            // Pressure
            y_interp_ptr[1] /= this->redim_pascal_scale;

            // These assume sphereical symmetry
            // Total mass
            y_interp_ptr[2] /= this->redim_mass_scale;

            // Moment of inertia (r^2 multipled by dm which was found above)
            y_interp_ptr[3] /= this->redim_moi_scale;

            // Density
            y_interp_ptr[4] /= this->redim_density_scale;
            
            // Shear modulus (real and complex)
            y_interp_ptr[5] /= this->redim_pascal_scale;
            y_interp_ptr[6] /= this->redim_pascal_scale;

            // Bulk modulus (real and complex)
            y_interp_ptr[7] /= this->redim_pascal_scale;
            y_interp_ptr[8] /= this->redim_pascal_scale;
        }
    }
    else
    {
        // TODO: Better error handling
        std::exception();
    }
}


void EOSSolutionCC::change_radius_array(
        double* new_radius_ptr,
        size_t new_radius_size)
{
    this->radius_array_size = new_radius_size;
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
        for (size_t i = 0; i < this->cysolver_results_sptr_bylayer_vec.size(); i++)
        {
            this->cysolver_results_sptr_bylayer_vec[i]->dense_vec.clear();
            this->cysolver_results_sptr_bylayer_vec[i].reset();
        }
        this->cysolver_results_sptr_bylayer_vec.clear();
        this->current_layers_saved = 0;

        // Indicate that all vectors are cleared.
        this->other_vecs_set = false;
    }
    this->radius_array_set = true;

    // Reserve capacity in vectors for new radius array size
    this->gravity_array_vec.reserve(this->radius_array_size);
    this->pressure_array_vec.reserve(this->radius_array_size);
    this->mass_array_vec.reserve(this->radius_array_size);
    this->moi_array_vec.reserve(this->radius_array_size);
    this->density_array_vec.reserve(this->radius_array_size);
    
    // These complex arrays are stored as double arrays with twice the length (Cython and C++ don't play nicely with complex across all systems)
    this->complex_shear_array_vec.reserve(2 * this->radius_array_size);
    this->complex_bulk_array_vec.reserve(2 * this->radius_array_size);

    // Copy over the radius array values
    this->radius_array_vec.resize(this->radius_array_size);
    std::memcpy(this->radius_array_vec.data(), new_radius_ptr, new_radius_size * sizeof(double));

    // Get constants
    this->radius = this->radius_array_vec.back();
}


void EOSSolutionCC::interpolate_full_planet()
{
    this->solution_nondim_status = this->nondim_status;

    if (this->current_layers_saved == 0)
    {
        // No layers have been saved, we can't perform the interpolation.
        throw std::exception();
    }
    
    size_t current_layer_index        = 0;
    double current_layer_upper_radius = this->upper_radius_bylayer_vec[0];

    double y_interp_arr[EOS_DY_VALUES];
    double* y_interp_ptr = &y_interp_arr[0];

    bool ready_for_next_layer = false;

    size_t interface_check = 0;

    for (size_t radius_i = 0; radius_i < this->radius_array_size; radius_i++)
    {
        const double radius = this->radius_array_vec[radius_i];

        if (isclose(radius, current_layer_upper_radius))
        {
            // At the layer's radius. We want to capture it once at interfaces (there will be two of the same radii for interface layers)
            if (interface_check == 1)
            {
                // No longer in current layer
                ready_for_next_layer = true;
            }
            interface_check++;
        }
        else if (radius > current_layer_upper_radius)
        {
            // No longer in current layer
            ready_for_next_layer = true;
        }

        if (ready_for_next_layer)
        {
            current_layer_index++;
            if (current_layer_index > (this->num_layers - 1))
            {
                // Outside of the planet
                break;
            }
            else
            {
                interface_check = 0;
                current_layer_upper_radius = this->upper_radius_bylayer_vec[current_layer_index];
                ready_for_next_layer = false;
            }
        }

        // Call interpolate using temp array as holder.
        this->cysolver_results_sptr_bylayer_vec[current_layer_index]->call(radius, y_interp_ptr);
        
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

        // Record central pressure
        if (current_layer_index == 0 && radius_i == 0)
        {
            this->central_pressure = y_interp_ptr[1];
        }
        
    }

    // At the end we should be at the planet's surface. Use the last interpolated values to set surface constants
    this->surface_gravity  = y_interp_ptr[0];
    this->surface_pressure = y_interp_ptr[1];
    this->mass             = y_interp_ptr[2];
    this->moi              = y_interp_ptr[3];

    // Finished
    this->other_vecs_set = true;
}

void EOSSolutionCC::dimensionalize_data(NonDimensionalScalesCC* nondim_scales, bool redimensionalize)
{   
    // Save scalers
    this->redim_length_scale  = nondim_scales->length_conversion;
    this->redim_gravity_scale = nondim_scales->length_conversion / nondim_scales->second2_conversion;
    this->redim_mass_scale    = nondim_scales->mass_conversion;
    this->redim_density_scale = nondim_scales->density_conversion;
    this->redim_moi_scale     = nondim_scales->mass_conversion * nondim_scales->length_conversion * nondim_scales->length_conversion;
    this->redim_pascal_scale  = nondim_scales->pascal_conversion;

    // Figure out how to set flags used during eos solution calls
    if (this->solution_nondim_status == 0)
    {
        // EOS was not explicitly non-dim'd or re-dim'd when solution was found. 
        if (redimensionalize)
        {
            // User is now "redimensionalizing." Thus we assume that the solution was non-dim'd
            this->nondim_status = 1;
        }
        else
        {
            // User is now "dimensionalizing." Thus we assume that the solution was dim'd
            this->nondim_status = -1;
        }
    }
    else
    {
        // TODO: Deal with this case! For now push the problem to the user when they try to call...
    }

    // Update other constants
    if (redimensionalize)
    {
        this->pressure_error *= this->redim_pascal_scale;
    }
    else
    {
        this->pressure_error /= this->redim_pascal_scale;
    }

    // Update layer data
    for (size_t layer_i = 0; layer_i < this->num_layers; layer_i++)
    {
        if (redimensionalize)
        {
            this->upper_radius_bylayer_vec[layer_i] *= this->redim_length_scale;
        }
        else
        {
            this->upper_radius_bylayer_vec[layer_i] /= this->redim_length_scale;
        }
    }
    
    if (this->other_vecs_set)
    {
        // Work through arrays
        for (size_t slice_i = 0; slice_i < this->radius_array_size; slice_i++)
        {
            if (redimensionalize)
            {
                this->radius_array_vec[slice_i]   *= this->redim_length_scale;
                this->gravity_array_vec[slice_i]  *= this->redim_gravity_scale;
                this->pressure_array_vec[slice_i] *= this->redim_pascal_scale;
                this->mass_array_vec[slice_i]     *= this->redim_mass_scale;
                this->moi_array_vec[slice_i]      *= this->redim_moi_scale;
                this->density_array_vec[slice_i]  *= this->redim_density_scale;

                // These complex arrays are stored as double arrays with twice the length (Cython and C++ don't play nicely with complex across all systems)
                this->complex_shear_array_vec[2 * slice_i]     *= this->redim_pascal_scale;
                this->complex_shear_array_vec[2 * slice_i + 1] *= this->redim_pascal_scale;
                this->complex_bulk_array_vec[2 * slice_i]      *= this->redim_pascal_scale;
                this->complex_bulk_array_vec[2 * slice_i + 1]  *= this->redim_pascal_scale;
            }
            else
            {
                this->radius_array_vec[slice_i]   /= this->redim_length_scale;
                this->gravity_array_vec[slice_i]  /= this->redim_gravity_scale;
                this->pressure_array_vec[slice_i] /= this->redim_pascal_scale;
                this->mass_array_vec[slice_i]     /= this->redim_mass_scale;
                this->moi_array_vec[slice_i]      /= this->redim_moi_scale;
                this->density_array_vec[slice_i]  /= this->redim_density_scale;

                // These complex arrays are stored as double arrays with twice the length (Cython and C++ don't play nicely with complex across all systems)
                this->complex_shear_array_vec[2 * slice_i]     /= this->redim_pascal_scale;
                this->complex_shear_array_vec[2 * slice_i + 1] /= this->redim_pascal_scale;
                this->complex_bulk_array_vec[2 * slice_i]      /= this->redim_pascal_scale;
                this->complex_bulk_array_vec[2 * slice_i + 1]  /= this->redim_pascal_scale;
            }
        }

        // Update global constants
        this->radius           = this->radius_array_vec[this->radius_array_size - 1];
        this->surface_gravity  = this->gravity_array_vec[this->radius_array_size - 1];
        this->surface_pressure = this->pressure_array_vec[this->radius_array_size - 1];
        this->mass             = this->mass_array_vec[this->radius_array_size - 1];
        this->moi              = this->moi_array_vec[this->radius_array_size - 1];
        this->central_pressure = this->pressure_array_vec[0];
    }
}