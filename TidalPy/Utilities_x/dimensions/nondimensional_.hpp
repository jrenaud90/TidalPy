#pragma once
/*
 * nondimensional_.hpp - c_NonDimensionalScales: unit-conversion scales for non-dimensionalized solves.
 *
 * The scheme follows Martens (2016, PhD thesis, CalTech, DOI: 10.7907/Z9N29TX7, ~p. 99): the time
 * scale is set by 1/(pi G rho_bulk) (frequency independent), the length scale by the mean radius,
 * and the density scale by the bulk density; the mass and pascal scales follow from those. Used by
 * the radial solver and the EOS solution to convert between SI and solve units.
 */

#include <cmath>
#include "../../constants_.hpp"

class c_NonDimensionalScales
{
public:
    double second2_conversion;
    double second_conversion;
    double length_conversion;
    double length3_conversion;
    double density_conversion;
    double mass_conversion;
    double pascal_conversion;

    c_NonDimensionalScales() = default;
    ~c_NonDimensionalScales() = default;
    c_NonDimensionalScales(const c_NonDimensionalScales&) = default;
    c_NonDimensionalScales(c_NonDimensionalScales&&) = default;
    c_NonDimensionalScales(
        double frequency,
        double mean_radius,
        double bulk_density
    )
    {
        if (tidalpy_config_ptr != nullptr)
        {
            this->second2_conversion = 1. / (TidalPyConstants::d_PI * tidalpy_config_ptr->d_G * bulk_density);
        }
        else
        {
            // Config is not initialized. Default to standard value.
            const double d_G = 6.674300000000e-11;
            this->second2_conversion = 1. / (TidalPyConstants::d_PI * d_G * bulk_density);
        }
        this->second_conversion  = std::sqrt(this->second2_conversion);
        this->length_conversion  = mean_radius;
        this->length3_conversion = mean_radius * mean_radius * mean_radius;
        this->density_conversion = bulk_density;
        this->mass_conversion    = bulk_density * this->length3_conversion;
        this->pascal_conversion  = \
            this->mass_conversion / (this->length_conversion * this->second2_conversion);
    }
};
